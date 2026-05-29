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

# Pre-load both datasets into isolated environments
orf_env <- new.env()
load("data/ORF/switchPlotFromTables.RData", envir = orf_env)

all_env <- new.env()
load("data/ALL/switchPlotFromTables.RData", envir = all_env)

orf_genes <- sort(unique(na.omit(orf_env$isoformFeatures$gene_name)))
all_genes  <- sort(unique(na.omit(all_env$isoformFeatures$gene_name)))

version_descriptions <- list(
  ORF = paste(
    "Contains only transcripts predicted to be protein-coding by",
    "ORFanage + RiboTIE + Ribo-seq reads filtering."
  ),
  ALL = "Contains all transcripts from the Iso-seq pipeline."
)

# Global version selector shown in the navbar header
version_bar <- div(
  style = paste(
    "background-color: #f1f3f5; padding: 8px 20px;",
    "border-bottom: 1px solid #dee2e6; display: flex;",
    "align-items: center; gap: 30px;"
  ),
  div(
    style = "font-weight: 600; white-space: nowrap;",
    "Analysis Version:"
  ),
  radioButtons(
    "version", label = NULL,
    choices = c("ORF" = "ORF", "ALL" = "ALL"),
    selected = "ORF", inline = TRUE
  ),
  uiOutput("version_note")
)

ui <- navbarPage(
  "Astrocyte Full-length Transcriptome",
  header = version_bar,

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
        )
      ),
      mainPanel(
        uiOutput("switch_plot_ui")
      )
    )
  ),

)

server <- function(input, output, session) {

  dataset <- reactive({
    if (input$version == "ORF") orf_env else all_env
  })

  gene_list <- reactive({
    if (input$version == "ORF") orf_genes else all_genes
  })

  output$version_note <- renderUI({
    div(
      style = "color: #495057; font-size: 0.92em;",
      tags$b(input$version), "—", version_descriptions[[input$version]]
    )
  })

  output$switch_gene_ui <- renderUI({
    genes <- gene_list()
    current <- isolate(input$switch_gene)
    selected <- if (!is.null(current) && current %in% genes) current else genes[1]
    selectizeInput(
      "switch_gene",
      "Select gene of interest:",
      choices  = genes,
      selected = selected
    )
  })

  n_isoforms <- reactive({
    req(input$switch_gene)
    env <- dataset()
    feats <- env$isoformFeatures[env$isoformFeatures$gene_name == input$switch_gene, ]
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

  output$switch_plot <- renderPlot({
    validate(need(input$switch_gene, "Choose a gene!"))
    env <- dataset()
    switchPlotFromTables(
      isoformFeatures = env$isoformFeatures,
      exons           = env$exons,
      conditions      = env$conditions,
      gene            = input$switch_gene,
      condition1      = "Unstim",
      condition2      = "Stim",
      orfAnalysis     = env$orfAnalysis,
      domainAnalysis  = env$domainAnalysis,
      IFcutoff        = input$IF_cutoff,
      dIFcutoff       = input$dIF_cutoff
    )
  })

}

shinyApp(ui, server)
