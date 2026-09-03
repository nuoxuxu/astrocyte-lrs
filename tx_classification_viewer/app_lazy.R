suppressPackageStartupMessages({
  library(ggtranscript)
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(readr)
  library(glue)
  library(shiny)
  library(memoise)
})

cat("Initializing app (lazy loading enabled)...\n")

# Data directory
data_dir <- "data"
if (!dir.exists(data_dir)) {
  stop("ERROR: 'data/' directory not found!\n",
       "Run 'Rscript setup_data_optimized.R' first.")
}

# Load data lazily - only load when first accessed
# This reduces memory footprint at startup
.data_cache <- new.env()

load_data <- function(name) {
  if (is.null(.data_cache[[name]])) {
    cat("Loading", name, "...\n")
    .data_cache[[name]] <- readRDS(file.path(data_dir, paste0(name, ".rds")))
  }
  .data_cache[[name]]
}

# Lazy getters
get_gencode <- function() load_data("gencode")
get_ribotie <- function() load_data("ribotie")
get_orfanage <- function() load_data("orfanage")
get_pbid_to_pr_transcripts <- function() load_data("pbid_to_pr_transcripts")
get_pbid_to_orfanage_template <- function() load_data("pbid_to_orfanage_template")
get_orf_type_lookup <- function() load_data("orf_type_lookup")
get_pclass_lookup <- function() load_data("pclass_lookup")

# Cache combined annotation (computed once after first data load)
.annotation_cache <- NULL
get_annotation_gtf <- function() {
  if (is.null(.annotation_cache)) {
    cat("Building annotation GTF...\n")
    .annotation_cache <<- bind_rows(
      get_gencode(),
      get_ribotie(),
      get_orfanage()
    )
    .annotation_cache$source <<- factor(.annotation_cache$source,
                                        levels = c("RiboTIE", "ORFanage", "GENCODE"))
  }
  .annotation_cache
}

cat("App initialized (data will load on first access)\n\n")

# Extracted and memoized difference-finding function
get_differences <- memoise(function(tx_of_interest, reference_to_highlight) {
  annotation_gtf <- get_annotation_gtf()
  pbid_to_orfanage_template <- get_pbid_to_orfanage_template()
  pbid_to_pr_transcripts <- get_pbid_to_pr_transcripts()

  orfanage_template <- pbid_to_orfanage_template %>%
    filter(transcript_id == sub("_.*$", "", tx_of_interest)) %>%
    pull(orfanage_template)
  pr_transcripts <- pbid_to_pr_transcripts %>%
    filter(isoform_id == tx_of_interest) %>%
    pull(pr_transcripts)

  annotation_per_tx <- annotation_gtf %>%
    dplyr::filter(
      case_when(
        source == "GENCODE" ~ transcript_id %in% c(pr_transcripts, orfanage_template),
        source == "RiboTIE" ~ transcript_id == tx_of_interest,
        source == "ORFanage" ~ transcript_id == sub("_.*$", "", tx_of_interest)
      )
    )

  ribotie_tx <- annotation_per_tx %>%
    filter(source == "RiboTIE", type == "CDS") %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  pick_diff <- function(a, b) {
    d <- GenomicRanges::setdiff(a, b)
    if (length(d) != 0) {
      d_ranges <- a[subjectHits(findOverlaps(d, a))]
    } else {
      d <- GenomicRanges::setdiff(b, a)
      d_ranges <- b[subjectHits(findOverlaps(d, b))]
    }
    d_ranges
  }

  if (reference_to_highlight == "ORFanage prediction") {
    orfanage_prediction_tx <- annotation_per_tx %>%
      filter(source == "ORFanage", type == "CDS") %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    pick_diff(ribotie_tx, orfanage_prediction_tx)
  } else if (reference_to_highlight == "ORFanage template (GENCODE)") {
    orfanage_template_tx <- annotation_per_tx %>%
      filter(transcript_id == orfanage_template, type == "CDS") %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    pick_diff(ribotie_tx, orfanage_template_tx)
  } else {
    pr_transcripts_tx <- annotation_per_tx %>%
      filter(transcript_id %in% pr_transcripts, type == "CDS") %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    pick_diff(ribotie_tx, pr_transcripts_tx)
  }
})

get_orfanage_template <- function(tx_of_interest) {
  get_pbid_to_orfanage_template() %>%
    filter(transcript_id == sub("_.*$", "", tx_of_interest)) %>%
    pull(orfanage_template)
}

get_pr_transcripts <- function(tx_of_interest) {
  get_pbid_to_pr_transcripts() %>%
    filter(isoform_id == tx_of_interest) %>%
    pull(pr_transcripts)
}

# Plotting function
set_theme(theme_minimal())
plot_tx <- function(tx_of_interest, reference_to_highlight, focused = FALSE, diff_index = 1) {
  annotation_gtf <- get_annotation_gtf()
  orf_type_lookup <- get_orf_type_lookup()
  pclass_lookup <- get_pclass_lookup()

  orfanage_template <- get_orfanage_template(tx_of_interest)
  pr_transcripts <- get_pr_transcripts(tx_of_interest)

  annotation_per_tx <- annotation_gtf %>%
    dplyr::filter(
      case_when(
        source == "GENCODE" ~ transcript_id %in% c(pr_transcripts, orfanage_template),
        source == "RiboTIE" ~ transcript_id == tx_of_interest,
        source == "ORFanage" ~ transcript_id == sub("_.*$", "", tx_of_interest)
      )
    )

  exons <- annotation_per_tx %>% dplyr::filter(type == "exon")
  CDS <- annotation_per_tx %>% dplyr::filter(type == "CDS")
  transcript_order <- unique(exons$transcript_id)

  pclass_val <- pclass_lookup[tx_of_interest]
  if (is.na(pclass_val)) pclass_val <- NA_character_
  ORF_biotype <- orf_type_lookup[tx_of_interest]
  if (is.na(ORF_biotype)) ORF_biotype <- NA_character_

  d_ranges <- get_differences(tx_of_interest, reference_to_highlight)

  if (orfanage_template == pr_transcripts) {
    title_a <- glue("The RiboTIE ORF has the biotype of {ORF_biotype} in respect to ORFanage-predicted ORF \nand {pclass_val} in respect to GENCODE ORF, matched by SQANTI3 algorithm. \nThe ORFanage template and SQANTI3 protein-matched GENCODE transcript are the same.")
  } else if (pr_transcripts == "novel") {
    title_a <- glue("The RiboTIE ORF has the biotype of {ORF_biotype} in respect to ORFanage-predicted ORF \nand {pclass_val} in respect to GENCODE ORF, matched by SQANTI3 algorithm. \n No corresponding GENCODE transcript matched by SQANTI3.")
  } else {
    title_a <- glue("The RiboTIE ORF has the biotype of {ORF_biotype} in respect to ORFanage-predicted ORF \nand {pclass_val} in respect to GENCODE ORF, matched by SQANTI3 algorithm.")
  }

  exons <- exons %>%
    mutate(transcript_id = factor(transcript_id, levels = transcript_order))

  if (reference_to_highlight == "ORFanage prediction") {
    highlight_df <- exons %>%
      filter(source %in% c("ORFanage", "RiboTIE"))
  } else if (reference_to_highlight == "ORFanage template (GENCODE)") {
    highlight_df <- exons %>%
      filter((transcript_id == orfanage_template) | (source == "RiboTIE"))
  } else if (reference_to_highlight == "SQANTI3 protein-matched GENCODE transcript") {
    highlight_df <- exons %>%
      filter((transcript_id == pr_transcripts) | (source == "RiboTIE"))
  }

  highlight_df <- highlight_df %>%
    group_by(transcript_id) %>%
    summarize(xmin = min(start), xmax = max(end), .groups = "drop") %>%
    mutate(
      y_num = as.numeric(factor(transcript_id, levels = transcript_order)),
      ymin = y_num - 0.4,
      ymax = y_num + 0.4
    )

  p <- exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_rect(
      data = highlight_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "yellow", alpha = 0.2, inherit.aes = FALSE
    ) +
    geom_range(
      height = 0.25
    ) +
    geom_range(
      data = CDS,
      aes(fill = source)
    ) +
    geom_intron(
      data = to_intron(exons, "transcript_id"),
      aes(strand = strand),
      arrow.min.intron.length = 500
    ) +
    labs(title = title_a)

  if (length(d_ranges) != 0) {
    if (focused && length(d_ranges) > 1) {
      diff_index <- min(diff_index, length(d_ranges))
      diff_index <- max(diff_index, 1)
      current_diff <- d_ranges[diff_index]
      p <- p +
        geom_vline(xintercept = start(current_diff), linetype = "dashed", color = "red") +
        geom_vline(xintercept = end(current_diff), linetype = "dashed", color = "red") +
        ggplot2::coord_cartesian(xlim = c(start(current_diff) - 100, end(current_diff) + 100))
    } else {
      p <- p +
        geom_vline(xintercept = start(d_ranges), linetype = "dashed", color = "red") +
        geom_vline(xintercept = end(d_ranges), linetype = "dashed", color = "red")
      if (focused) {
        p <- p + ggplot2::coord_cartesian(xlim = c(min(start(d_ranges)) - 100, max(end(d_ranges)) + 100))
      }
    }
  }
  p
}

# Shiny UI
ui <- fluidPage(
  titlePanel("Transcript Classification Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput(
        "tx_of_interest",
        "Select Transcript:",
        choices = NULL,
        selected = "PB.11906.164_253",
        options = list(maxOptions = 5000)
      ),
      radioButtons(
        "reference_to_highlight",
        "Reference to Highlight:",
        choices = c("ORFanage prediction", "SQANTI3 protein-matched GENCODE transcript", "ORFanage template (GENCODE)"),
        selected = "ORFanage prediction"
      ),
      checkboxInput(
        "focused",
        "Focused View (zoom on differences)",
        value = FALSE
      ),
      br(),
      div(
        style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px;",
        textOutput("diff_counter"),
        div(
          style = "display: flex; gap: 10px; margin-top: 10px; justify-content: center;",
          actionButton("prev_diff", "← Previous", class = "btn-sm"),
          actionButton("next_diff", "Next →", class = "btn-sm")
        )
      )
    ),
    mainPanel(
      plotOutput("plot", height = "600px")
    )
  )
)

# Shiny Server
server <- function(input, output, session) {
  # Server-side selectize for 27k transcripts (massive performance improvement)
  shiny::updateSelectizeInput(
    session,
    "tx_of_interest",
    choices = sort(unique(get_ribotie()$transcript_id)),
    server = TRUE
  )

  diff_index <- reactiveVal(1)

  get_n_diffs <- reactive({
    n_diffs <- length(get_differences(input$tx_of_interest, input$reference_to_highlight))
    n_diffs
  })

  observeEvent(input$tx_of_interest, {
    diff_index(1)
  })

  observeEvent(input$reference_to_highlight, {
    diff_index(1)
  })

  observeEvent(input$next_diff, {
    n_diffs <- get_n_diffs()
    if (n_diffs > 1) {
      new_index <- diff_index() + 1
      if (new_index <= n_diffs) {
        diff_index(new_index)
      }
    }
  })

  observeEvent(input$prev_diff, {
    n_diffs <- get_n_diffs()
    if (n_diffs > 1) {
      new_index <- diff_index() - 1
      if (new_index >= 1) {
        diff_index(new_index)
      }
    }
  })

  output$diff_counter <- renderText({
    n_diffs <- get_n_diffs()
    if (n_diffs == 0) {
      "No differences"
    } else if (n_diffs == 1) {
      "1 difference found"
    } else {
      paste("Difference", diff_index(), "of", n_diffs)
    }
  })

  output$plot <- renderPlot({
    req(input$tx_of_interest)

    # Force garbage collection before rendering to free memory
    gc(verbose = FALSE)

    plot_tx(
      tx_of_interest = input$tx_of_interest,
      reference_to_highlight = input$reference_to_highlight,
      focused = input$focused,
      diff_index = diff_index()
    )
  })
}

# Run the app
shinyApp(ui, server)
