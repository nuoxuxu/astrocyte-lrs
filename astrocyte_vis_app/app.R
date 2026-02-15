library(ggplot2)
library(ggtranscript)
library(dplyr)
library(shiny)
library(bslib)
library(patchwork)
library(dplyr)
library(ggnewscale)
library(IsoformSwitchAnalyzeR)

my_theme <- theme_bw() + theme(
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_blank(),
  strip.text.y = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.ticks.y = element_blank(),
  legend.position = "top",
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 14)
)

theme_set(my_theme)

colorVector <- c(
    "FSM" = "#009E73",
    "ISM" = "#0072B2",
    "NIC" = "#D55E00",
    "NNC" = "#E69F00",
    "Other" = "#000000"
)

pastelColorVector <- c(
    "FSM" = "#7FCEB9",
    "ISM" = "#7FB8D8",
    "NIC" = "#EAAE7F",
    "NNC" = "#F2CF7F",
    "Other" = "#7F7F7F"
)

gtf <- read.csv("data/gtf.csv")
lr_log2_cpm <- read.csv("data/lr_log2_cpm.csv") %>%
  as_tibble() %>%
  mutate(
    structural_category = factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC", "Other")),
    condition = factor(condition, levels = c("mean_Unstim", "mean_Stim"))
  )

switchlist_low <- readRDS("data/switchlist_low.rds")
switchlist_high <- readRDS("data/switchlist_high.rds")
switch_genes <- sort(unique(na.omit(switchlist_low$isoformFeatures$gene_name)))

plot_ggtranscript <- function(gene_of_interest) {
  exons <- gtf %>%
    filter(
      type == "exon",
      associated_gene == {{ gene_of_interest }}
    )
  CDS <- gtf %>%
    filter(
      type == "CDS",
      associated_gene == {{ gene_of_interest }}
    )
  exons %>%
    ggplot(
      aes(xstart = start, xend = end, y = transcript_id, strand = strand)
    ) +
    geom_intron(
      data = to_intron(exons, "transcript_id")
    ) +
    # hide exon legend
    geom_range(
      aes(fill = structural_category),
      height = 0.25,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = pastelColorVector, guide = "none") +
    ggnewscale::new_scale_fill() +
    # keep only CDS legend
    geom_range(
      data = CDS,
      aes(fill = structural_category),
      color = "black"
    ) +
    scale_fill_manual("Structural Category", values = colorVector) +
    facet_grid(structural_category ~ ., scales = "free", space = "free")
}

plot_abundance <- function(gene_of_interest) {
  lr_log2_cpm %>%
    filter(associated_gene == gene_of_interest) %>%
    ggplot(aes(x = abundance, y = isoform, fill = condition)) +
    geom_col(position = "dodge") +
    xlab("mean log2(CPM + 1)") +
    ylab(gene_of_interest) +
    theme(
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 14),
    ) +
    scale_fill_manual(breaks = c("mean_Unstim", "mean_Stim"), values = c("#e0a19c", "#e14bd5")) +
    facet_grid(structural_category ~ ., scales = "free", space = "free")
}

plot_combined <- function(gene_of_interest) {
  ggtranscript_plot <- plot_ggtranscript(gene_of_interest)
  abundance <- plot_abundance(gene_of_interest)
  ggtranscript_plot | abundance + plot_layout(widths = c(2, 1))
}

get_n_transcripts <- function(gene_of_interest) {
  gtf %>%
    filter(type == "transcript") %>%
    filter(associated_gene == gene_of_interest) %>%
    nrow()
}

ui <- navbarPage(
  "Astrocyte Full-length Transcriptome",
  tabPanel(
    "Transcript View",
    helpText("Explore transcript structures and abundances for genes of interest. Transcript structural categories are based on comparison to GENCODE v47."),
    sidebarLayout(
      sidebarPanel(
        selectizeInput(
          "select_genes",
          "Select gene of interest:",
          choices = unique(gtf$associated_gene),
          selected = "AGO1"
        )
      ),
      mainPanel(
        plotOutput("combined_plot")
      )
    )
  ),
  tabPanel(
    "Isoform Switching",
    helpText("Visualize isoform switching events using IsoformSwitchAnalyzeR switchPlot."),
    sidebarLayout(
      sidebarPanel(
        selectizeInput(
          "switch_gene",
          "Select gene of interest:",
          choices = switch_genes,
          selected = switch_genes[1]
        ),
        radioButtons(
          "stringency",
          "Stringency:",
          choices = c("Low" = "low", "High" = "high"),
          selected = "low"
        )
      ),
      mainPanel(
        plotOutput("switch_plot", height = "600px")
      )
    )
  )
)

server <- function(input, output, session) {
  output$combined_plot <- renderPlot(
    {
      validate(need(input$select_genes, 'Choose a gene!'))
      plot_combined(input$select_genes)
    },
    height = reactive({
      35 * sum(get_n_transcripts(input$select_genes)) + 100
    })
  )

  selected_switchlist <- reactive({
    if (input$stringency == "high") switchlist_high else switchlist_low
  })

  output$switch_plot <- renderPlot({
    validate(need(input$switch_gene, 'Choose a gene!'))
    switchPlot(
      selected_switchlist(),
      gene = input$switch_gene,
      condition1 = "Unstim",
      condition2 = "Stim"
    )
  })
}

shinyApp(ui, server)
