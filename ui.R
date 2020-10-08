# example from http://shiny.rstudio.com/gallery/kmeans-example.html

library(shiny)

ui <- fluidPage(
  titlePanel("Genomic Exploration in HAPI-FHIR"),

  sidebarLayout(
      sidebarPanel(
          style = "position:fixed;width:inherit;",
          width = 3,
          textInput("gene", label = "Search for mutations in Gene", value = "EGFR"),
          textInput("region", label = "Or search in genomic region (GRCh38)"),
          helpText("For example, 7:55019020-55211628"),
          actionButton("search", "Search"),
          br(),
          icon("arrow-down"),
          hr(),
          actionButton("annot", "Variant Annotation"),
          helpText("To predict amino acid changes for each gene mutation. Note: this step will take some time."),
          icon("arrow-down"),
          hr(),
          actionButton("lplot", "Hotspot Plot"),
          helpText("Genetic mutation data in protein domains using a lollipop-diagram"),
          icon("arrow-down"),
          hr(),
          actionButton("dgi", "Search for Drug Interaction"),
          helpText("Search in DGIdb, the Drug Gene Interaction Database")
    ),

    mainPanel(
        h4("Variants in Patients"),
        withSpinner(DT::dataTableOutput("vars")),
        hr(),
        h4("Variants Annotation"),
        withSpinner(DT::dataTableOutput("muts")),
        hr(),
        h4("Mutation Hotspot Plot"),
        withSpinner(g3LollipopOutput("llplot")),
        h4("Gene-Drug Interaction"),
        withSpinner(DT::dataTableOutput("drug"))
    )
  )
)
