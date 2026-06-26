library(shiny)
library(bslib)
library(ggplot2)
library(grid)
library(GenomicRanges)
library(IRanges)

source("classes.R")
source("tools.R")
source("isoform_plots.R")
source("switch_plot_from_tables.R")

app_env <- new.env()
load("data/minlen/switchPlotFromTables.RData", envir = app_env)
app_genes <- sort(unique(na.omit(app_env$isoformFeatures$gene_name)))

ui <- navbarPage(
  "Astrocyte Full-length Transcriptome",

  tabPanel(
    "Isoform Switching",
    helpText("Visualize isoform switching events using IsoformSwitchAnalyzeR switchPlot."),
    sidebarLayout(
      sidebarPanel(
        uiOutput("switch_gene_ui"),
        hr(),
        numericInput(
          "IF_cutoff",
          "IFcutoff — minimum isoform fraction to display:",
          value = 0.05, min = 0, max = 1, step = 0.01
        ),
        numericInput(
          "dIF_cutoff",
          "dIFcutoff — minimum |dIF| to call a switch:",
          value = 0.1, min = 0, max = 1, step = 0.01
        ),
        hr(),
        helpText("To view the genomic region of the selected gene in the UCSC Genome Browser, click the link below:"),
        uiOutput("ucsc_link")
      ),
      mainPanel(
        uiOutput("switch_plot_ui")
      )
    )
  ),

)

server <- function(input, output, session) {

  output$switch_gene_ui <- renderUI({
    current <- isolate(input$switch_gene)
    selected <- if (!is.null(current) && current %in% app_genes) current else app_genes[1]
    selectizeInput(
      "switch_gene",
      "Select gene of interest:",
      choices  = app_genes,
      selected = selected
    )
  })

  n_isoforms <- reactive({
    req(input$switch_gene)
    feats <- app_env$isoformFeatures[app_env$isoformFeatures$gene_name == input$switch_gene, ]
    if (!all(is.na(feats$IF1))) {
      feats <- feats[pmax(feats$IF1, feats$IF2, na.rm = TRUE) >= input$IF_cutoff, ]
    }
    length(unique(feats$isoform_id))
  })

  output$switch_plot_ui <- renderUI({
    n <- max(n_isoforms(), 1L)
    input$dIF_cutoff  # ensure container re-renders when dIF threshold changes
    height <- max(500, 50 + n * 50)
    plotOutput("switch_plot", height = paste0(height, "px"))
  })

  output$ucsc_link <- renderUI({
    req(input$switch_gene)
    gene_exons <- app_env$exons[app_env$exons$gene_name == input$switch_gene]
    validate(need(length(gene_exons) > 0, "No exon data for this gene."))
    chrom <- as.character(seqnames(gene_exons))[1]
    pos_start <- min(start(gene_exons))
    pos_end   <- max(end(gene_exons))
    url <- paste0("https://genome.ucsc.edu/s/nuoxuxu/asotryctes_ribo_seq?position=", chrom, ":", pos_start, "-", pos_end)
    tags$a(href = url, "View in UCSC Genome Browser", target = "_blank")
  })

  output$switch_plot <- renderPlot({
    validate(need(input$switch_gene, "Choose a gene!"))
    switchPlotFromTables(
      isoformFeatures = app_env$isoformFeatures,
      exons           = app_env$exons,
      conditions      = app_env$conditions,
      gene            = input$switch_gene,
      condition1      = "Unstim",
      condition2      = "Stim",
      orfAnalysis     = app_env$orfAnalysis,
      domainAnalysis  = app_env$domainAnalysis,
      IFcutoff        = input$IF_cutoff,
      dIFcutoff       = input$dIF_cutoff
    )
  })

}

shinyApp(ui, server)
