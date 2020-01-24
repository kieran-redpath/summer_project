# Load Shiny
library(shiny)
# Load packages (outside the Shiny App)
library(limma)
library(edgeR)
library(data.table)
library(magrittr)
library(ggplot2)
library(CePa)
library(tidyverse)
library(ggbeeswarm)
library(dplyr)
library(ReactomePA)
library(org.Hs.eg.db)
library(reactome.db)
library(goseq)
library(gplots)
source('heatmap-mik.R')
library(knitr)
library(kableExtra)
library(devtools)
# Load Data (outside the Shiny App)
# Load CCLE
ccleExpData <- read.gct('data/CCLE_RNAseq_genes_counts_20180929.gct')
# Load GDSC
GDSC2 <- fread('data/GDSC2_fitted_dose_response_15Oct19.csv', sep= ';')










# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Drug Response in Cancer Tissues (CCLE Cell Lines)"),
  
  # Sidebar containing input options: Description, drug of interest, and selectable cell lines of interest
  sidebarLayout(
    sidebarPanel(
      helpText(em("Enter the name of a drug and a list of tissue types of interest."),
               textInput("drug_name",
                         label = h3("Drug of Interest (GDSC Name)"),
                         value = "Enter drug name..."
               ),
               
               textInput("cell_types",
                         label = h3("Enter cell line tissue type, separated by commas. A full list of supported names can be found",
                                  a("here.", href = "INSERT LINK HERE")),
                         value = "Enter CCLE names..."
               )
      )
    ),
    
    # Main panel containing the generated heatmap
    mainPanel(
      # plotOutput(outputID = "heatmap"),
      p("words.")
    )
  )
)










# Server that actually processes the data, with render calls etc.
server <- function(input, output) {
  output$heatmap <- renderPlot({
    # Should use: input$drug_name, input$cell_types_CCLE, input$cell_types_GDSC
    # Most of the code you have should go here.
    # It can be tidied up later so only the important parts are reactive, not everything.
  })  
  
}










# Run the application 
shinyApp(ui = ui, server = server)