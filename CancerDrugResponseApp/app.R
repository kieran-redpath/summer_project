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
logExpDat <- log(ccleExpData + 0.5)
dge <- DGEList(counts=ccleExpData)
dge <- calcNormFactors(dge)
dge_voom <- voom(dge, plot=FALSE)
expDat <- dge_voom$E
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
                         label = h3("Enter cell line tissue type, separated by commas. A full list of supported tissues can be found",
                                    a("here.", href = "INSERT LINK HERE")),
                         value = "Enter tissue types..."
               )
      )
    ),
    
    # Main panel containing the generated heatmap
    mainPanel(
      tabsetPanel(
        tabPanel("Enriched Pathways", tableOutput("pathways")),
        tabPanel("Overlap Between Pathways", plotOutput("overlap")),
        tabPanel("Heatmaps of Top 5 Pathways", plotOutput("heatmap"))
      )
    )
  )
)










# Server that actually processes the data, with render calls etc.
server <- function(input, output) {
  
  # Should use: input$drug_name (it does now), input$cell_types
  # Most of the code you have should go here.
  # It can be tidied up later so only the important parts are reactive, not everything.
  # This expression must return an object (the plot)
  
  # TAKE INPUT$DRUG_NAME and create the drug data file
  drug <- GDSC2 %>% filter(grepl(input$drug_name, DRUG_NAME, ignore.case=TRUE))
  expDatsplitnames <- strsplit(colnames(expDat), "_") %>% 
    lapply(., function(x) x[1]) %>% 
    unlist()
  expDatsplit <- expDat
  colnames(expDat) <- expDatsplitnames
  drug$CELL_LINE_NAME <- gsub("-","",drug$CELL_LINE_NAME, fixed=TRUE)
  commonSamples <- intersect(drug$CELL_LINE_NAME,colnames(expDat))
  expDat_match <- match(commonSamples, colnames(expDat))
  expDat_sort <- expDat[ , na.omit(expDat_match)]
  drug_match <- match(commonSamples, drug$CELL_LINE_NAME)
  drug_sort <- drug[na.omit(drug_match),]
  colnames(drug_sort)[5] <- "CCLE_Name"
  
  # MAKE TISSUETOOL, AND ADD A COLUMN TO IT SO IT CAN CONVERT TO CELL TYPE IDENTIFIER
  lng <- length(colnames(expDatsplit))
  tissuetool <-data.frame(Cell_Line=character(lng),Tissue_Type=character(lng))
  tissuetool$Cell_Line <- strsplit(colnames(expDatsplit), "_") %>% 
    lapply(., function(x) x[1]) %>% 
    unlist()
  tissuetool$Tissue_Type <- strsplit(colnames(expDatsplit), "_") %>% 
    lapply(., function(x) x[-1])%>%
    paste(., sep="_", collapse = NULL) %>% 
    unlist()
  tissuetool$Tissue_Type <- gsub('c', '', tissuetool$Tissue_Type)
  tissuetool$Tissue_Type <- gsub('[[:punct:]]', '', tissuetool$Tissue_Type)
  tissuetool$Tissue_Type <- gsub(' ', '_', tissuetool$Tissue_Type)
  drug_help <- drug_sort %>% subset(., select=c("CCLE_Name", "TCGA_DESC"))
  drug_help <- drug_help %>% dplyr::rename(., Cell_Line = CCLE_Name)
  tissuetool <- inner_join(tissuetool, drug_help, by="Cell_Line")
  
  #LOOPS ARE MISSING! ON THE OTHER DOCUMENT THOUGH
  
  
  
  
  
  
  
  
  output$pathways <- renderTable({
    # goseqPathways code goes here
  })
  
  output$overlap <- renderPlot({
    # overlap of pathways heatmap goes here
  })
  
  output$heatmap <- renderPlot({
    # heatmap generating code goes here. Only want the top 5-10
  })  
  
}










# Run the application 
shinyApp(ui = ui, server = server)