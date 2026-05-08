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

load("data/switchPlotFromTables.RData")
switch_genes <- sort(unique(na.omit(isoformFeatures$gene_name)))

ui <- navbarPage(
  "Astrocyte Full-length Transcriptome",
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
        )
      ),
      mainPanel(
        plotOutput("switch_plot", height = "600px")
      )
    )
  )
)

server <- function(input, output, session) {
  output$switch_plot <- renderPlot({
    validate(need(input$switch_gene, 'Choose a gene!'))
    switchPlotFromTables(
        isoformFeatures = isoformFeatures,
        exons           = exons,
        conditions      = conditions,
        gene            = input$switch_gene,
        condition1      = "Unstim",
        condition2      = "Stim",
        orfAnalysis     = orfAnalysis,
        domainAnalysis  = domainAnalysis
    )
  })
}

shinyApp(ui, server)