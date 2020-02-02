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
                         label = h3("Enter cell line tissue type, separated by commas. A full list of supported tissues (case-insensitive) can be found",
                                    a("here.", href = "https://raw.githubusercontent.com/kieran-redpath/summer_project/master/SupportedTissues_CancerDrugResponseApp.txt?token=ANXQUX7B3Y7NBE73E5OQRRK6G4VOK")),
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
  # NEED TO MAKE SURE IT RECOGNISES THE DRUG NAME EVEN IF IT HAS SPACES!!!!
  drug <- reactive({
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
  })
  # NO IDEA IF THIS WORKS 
  # LOOPS NEED FIXING
  # A filter for "tissuetoolsort", called "tissue_filter"
  tissue_filter <- reactive({
    tissue_filter <- strsplit(input$cell_types, ", ")
    tissue_filter[[1]] <- gsub(' ', '_', tissue_filter[[1]])
    tissue_filter <- as.vector(tissue_filter[[1]])
    tissue_filter <- toupper(tissue_filter)
    tissue_filter <- as.data.frame(tissue_filter)
    colnames(tissue_filter)<- "Tissue_Type"
    
    # A filter for "drug_sort", called "tcga_filter"
    tcga_filter <- inner_join(tissue_filter, tissuetool, by="Tissue_Type")
    tcga_filter <- subset(tcga_filter, TCGA_DESC!="UNCLASSIFIED")
    tcga_filter <- unique(unlist(tcga_filter$TCGA_DESC))
    tcga_filter <- as.data.frame(tcga_filter)
    colnames(tcga_filter) <- "TCGA_DESC"
    # Loops 1
    for(d in 1:length(tcga_filter)){
      drug_sort <- drug_sort %>% filter(., TCGA_DESC==tcga_filter$TCGA_DESC) 
    }
    drug_sort <- drug_sort[ drug_sort$AUC < quantile(drug_sort$AUC , 0.25 ) | drug_sort$AUC > quantile(drug_sort$AUC, 0.75), ]
    namelist <- drug_sort$CCLE_Name
    expDat_sort <- expDat_sort %>% subset(., select=which(colnames(expDat_sort) %in% namelist))
    # Loops 2
    for(e in 1:length(tissue_filter)){
      tissuetoolsort <- tissuetool %>% filter(., Tissue_Type==tissue_filter$Tissue_Type)
    }
    int <- intersect(colnames(expDat_sort), tissuetoolsort$Cell_Line)
    tissuetool_match <- match(int, tissuetoolsort$Cell_Line)
    tissuetoolsort <- tissuetoolsort[na.omit(tissuetool_match) , ]
    # NO IDEA IF THIS WORKS ^
  })
  
  # -log IC50 Dif Exp Analysis
  group <- reactive({
    group <- ifelse(drug_sort$LN_IC50 > median(drug_sort$LN_IC50), "High", "Low")
    colnames(design) = c("Mean"
                         ,"HighVsLow")
    fit = lmFit(expDat_sort, design)
    fit = eBayes(fit)
    tt = topTable(fit, coef="HighVsLow", adjust="BH",n=nrow(expDat_sort))
    options(digits=4)
    sigFC = (tt$adj.P.Val < 0.01)  & (abs(tt$logFC) > 1)
    split <- strsplit(rownames(tt),".", fixed=T) %>% lapply(., function(x) x[1]) %>% unlist()
    geneNames <- AnnotationDbi::select(org.Hs.eg.db, keys = split, column = c("SYMBOL","GENENAME"), key="ENSEMBL")
    tt$symbol <- geneNames$SYMBOL[match(split, geneNames$ENSEMBL)]
    invisible(setDT(tt, keep.rownames = TRUE)[])
    topExp <- expDat_sort[match(tt$rn[1], rownames(expDat_sort)),]
    df <- data.frame(topGene=topExp, ic50=group)
  })
  # AUC Dif Exp Analysis
  group2 <- reactive({
    group2 <- ifelse(drug_sort$AUC > median(drug_sort$AUC), "High", "Low")
    colnames(design2) = c("Mean"
                          ,"HighVsLow"
    )
    fit2 = lmFit(expDat_sort, design2)
    fit2 = eBayes(fit2)
    tt2 = topTable(fit2, coef="HighVsLow", adjust="BH",n=nrow(expDat_sort))
    options(digits=4)
    sigFC2 = (tt2$adj.P.Val < 0.01)  & (abs(tt2$logFC) > 1)
    split2 <- strsplit(rownames(tt2),".", fixed=T) %>% lapply(., function(x) x[1]) %>% unlist()
    geneNames2 <- AnnotationDbi::select(org.Hs.eg.db, keys = split2, column = c("SYMBOL","GENENAME"), key="ENSEMBL")
    tt2$symbol <- geneNames2$SYMBOL[match(split2, geneNames2$ENSEMBL)]
    invisible(setDT(tt2, keep.rownames = TRUE)[])
    topExp2 <- expDat_sort[match(tt2$rn[1], rownames(expDat_sort)),]
    df2 <- data.frame(topGene=topExp2, AUC=group)
  })
  # Combining -log IC50 and AUC
  tt3 <- reactive({
    tt3 <- full_join(tt, tt2, by= "rn")
    tt3 <- dplyr::select(tt3, -c("symbol.x", "P.Value.y", "P.Value.x", "t.y", "t.x", "B.x", "B.y", "AveExpr.x"))
    tt3 <- dplyr::rename(tt3, "Gene_Symbol" = "symbol.y", "AUC_logFC" = "logFC.y", "Avg_Exp" = "AveExpr.y", "AUC_Adj_PVal" = "adj.P.Val.y", "IC50_logFC" = "logFC.x", "IC50_Adj_PVal" = "adj.P.Val.x", "Ensembl_ID" = "rn")
    tt3 <- tt3[c(7,1,5,2,3,4,6)]
    tt3 <- mutate(tt3, sign(tt3$IC50_logFC), sign(tt3$AUC_logFC))
    tt3 <- mutate(tt3, sign(tt3$IC50_logFC)*sign(tt3$AUC_logFC))
    tt3 <- dplyr::rename(tt3, "logFC_Sign" = "sign(tt3$IC50_logFC) * sign(tt3$AUC_logFC)", "IC50_Sign" = "sign(tt3$IC50_logFC)", "AUC_Sign" = "sign(tt3$AUC_logFC)")
    SigSamples <- tt3[1:500,]
    SigSamples <- 
      SigSamples %>% mutate(., rank_ic50=rank(IC50_Adj_PVal)) %>%
      mutate(., rank_auc=rank(AUC_Adj_PVal)) %>% 
      mutate(., avg_rank=0.5*(rank_auc + rank_ic50)) %>% 
      arrange(., avg_rank)
    tt3 <- 
      tt3 %>% mutate(., rank_ic50=rank(IC50_Adj_PVal)) %>%
      mutate(., rank_auc=rank(AUC_Adj_PVal)) %>% 
      mutate(., avg_rank=0.5*(rank_auc + rank_ic50)) %>% 
      arrange(., avg_rank)
  })
  
  # GoSeq Analysis
  Sigsamples <- reactive({
    SigSamples <- filter(tt3, tt3$logFC_Sign == 1, abs(tt3$IC50_logFC) > log2(2), abs(tt3$AUC_logFC) > log2(2))
    SigSamples <- 
      SigSamples %>% mutate(., rank_ic50=rank(IC50_Adj_PVal)) %>%
      mutate(., rank_auc=rank(AUC_Adj_PVal)) %>% 
      mutate(., avg_rank=0.5*(rank_auc + rank_ic50)) %>% 
      arrange(., avg_rank)
    SigSamples <- SigSamples[1:500,]
    TopXGenes <- as.character(na.omit(SigSamples$Gene_Symbol))
    TopXGenesEntrez <- AnnotationDbi::select(hs, 
                                             keys = TopXGenes,
                                             columns = c("ENTREZID", "SYMBOL"),
                                             keytype = "SYMBOL")
    tt3Genes <- as.character(na.omit(tt3$Gene_Symbol))
    tt3GenesEntrez <- AnnotationDbi::select(hs, 
                                            keys = tt3Genes,
                                            columns = c("ENTREZID", "SYMBOL"),
                                            keytype = "SYMBOL")
    rName <- as.list(reactomePATHNAME2ID)
    rName <- rName[grep("Homo sapiens", names(rName))]
    rGenes <- as.list(reactomePATHID2EXTID)
    rGenesPath <- rGenes[match(rName, names(rGenes))]
    rGenesPath <- lapply(rGenesPath, unique)
    rGeneByPath <- as.list(reactomeEXTID2PATHID)
    allGenes <- intersect( tt3GenesEntrez$ENTREZID, unique(unlist(rGenesPath)) )
    sigGenes <- intersect( TopXGenesEntrez$ENTREZID, unique(unlist(rGenesPath)) )
    plotGenes <- rep(0, length(allGenes))
    names(plotGenes) <- allGenes
    plotGenes[match(sigGenes, names(plotGenes))] <- 1
    mt <- match(allGenes, names(rGeneByPath))
    rGeneByPath <- lapply(rGeneByPath[mt], function(x) intersect(x, names(rGenesPath)))
    invisible(pwf <- nullp(plotGenes, 'hg19', id = "knownGene", plot.fit=TRUE))
    goseqReactome <- goseq(pwf, gene2cat = rGeneByPath)
    hyperReactome <- goseq(pwf, gene2cat = rGeneByPath, method="Hypergeometric")
    goseqReactome$adjP <- p.adjust(goseqReactome$over_represented_pvalue, method="fdr")
    hyperReactome$adjP <- p.adjust(hyperReactome$over_represented_pvalue, method="fdr")
    goseqPathways <- filter(goseqReactome, goseqReactome$adjP <1)
    rPathName <- as.list(reactomePATHID2NAME)
    goseqPathways$Pathway <- gsub("Homo sapiens: ", "", rPathName[match(goseqPathways$category, names(rPathName))])
    SiggoseqPathways <- goseqPathways
  })
  
  # Extract Genes From GoSeq Analysis
  GeneLabelTool <- reactive({
    GeneLabelTool <- dplyr::pull(tt3, Gene_Symbol)
    GeneLabelTool <- AnnotationDbi::select(hs,
                                           keys = GeneLabelTool,
                                           columns = c("ENSEMBL", "ENTREZID", "SYMBOL"),
                                           keytype = "SYMBOL")
    SiggenesinPaths <-list()
    for(i in 1:nrow(goseqPathways)){
      SiggenesinPaths[[i]] <- rGenesPath[match(goseqPathways$category[i], names(rGenesPath))] %>% 
        .[[1]] %>% 
        intersect(., rownames(pwf)[pwf$DEgenes==1])
    }
    genesinPaths <- list()
    for(i in 1:nrow(goseqPathways)){
      genesinPaths[[i]] <- rGenesPath[match(goseqPathways$category[i], names(rGenesPath))] %>% 
        .[[1]] %>% 
        intersect(., rownames(pwf)[pwf$DEgenes==1 | pwf$DEgenes==0])
    }
    SigsymbolsinPaths <- lapply(SiggenesinPaths, function(x) GeneLabelTool$SYMBOL[na.omit(match(x, GeneLabelTool$ENTREZID))] )
    symbolsinPaths <- lapply(genesinPaths, function(x) GeneLabelTool$SYMBOL[na.omit(match(x, GeneLabelTool$ENTREZID))] )
    genesStick <- lapply(symbolsinPaths, function(x) paste0(x, collapse="::", sep="")) %>% unlist()
    goseqPathways$DEgenesInCat <- genesStick
    goseqPathways <- lapply(goseqPathways, gsub, pattern='/', replacement=' ') %>% as.data.frame()
    goseqPathways$Pathway <- lapply(goseqPathways$Pathway, gsub, pattern=':', replacement='-') %>% as.data.frame()
    SiggenesStick <- lapply(SigsymbolsinPaths, function(x) paste0(x, collapse="::", sep="")) %>% unlist()
  })
  SiggoseqPathways <- reactive({
    SiggoseqPathways$DEgenesInCat <- SiggenesStick
    SiggoseqPathways <- lapply(SiggoseqPathways, gsub, pattern='/', replacement=' ') %>% as.data.frame()
    SiggoseqPathways$Pathway <- lapply(SiggoseqPathways$Pathway, gsub, pattern=':', replacement='-') %>% as.data.frame()
  })
  
  # GoSeqPathways renderTable
  output$pathways <- renderTable(
    SiggoseqPathways(), striped=TRUE, bordered=TRUE, digits=4
  )
  
  
  
  # Heatmap Setup
  tissuetoollsort <- reactive({
    tissuetoolsort$IC50_Group <- ifelse(drug_sort$LN_IC50 > median(drug_sort$LN_IC50), "High", "Low")
    tissuetoolsort$AUC_Group <- ifelse(drug_sort$AUC > median(drug_sort$AUC), "High", "Low")
    tis <- tissuetoolsort$Tissue_Type
    cc <- rbind((as.factor(tis) %>% as.numeric() %>% rainbow(length(table(.)))[.]),
                c("yellow", "blue")[as.numeric(as.factor(tissuetoolsort$IC50_Group))],
                c("yellow", "blue")[as.numeric(as.factor(tissuetoolsort$AUC_Group))])
    rownames(cc) <- c("Tissue Type", "-log IC50 Group", "AUC Group")
    colord = order(cc["AUC Group",])
    ensginPaths <- lapply(genesinPaths, function(x) GeneLabelTool$ENSEMBL[na.omit(match(x, GeneLabelTool$ENTREZID))] )
    SigensginPaths <- lapply(SiggenesinPaths, function(x) GeneLabelTool$ENSEMBL[na.omit(match(x, GeneLabelTool$ENTREZID))] )
    rownames(expDat_sort) <- strsplit(rownames(expDat_sort),".", fixed=T) %>% lapply(., function(x) x[1]) %>% unlist()
    # Overlap Heatmap Setup
    sigPathCor <- matrix(0, length(SigsymbolsinPaths), length(SigsymbolsinPaths))
    rownames(sigPathCor) <- colnames(sigPathCor) <- paste0(goseqPathways$Pathway," (", goseqPathways$numDEInCat, ")")
    for(a in 1:length(SigsymbolsinPaths)){
      for(b in 1:length(SigsymbolsinPaths)){
        if(a >  b) sigPathCor[a,b] <- length(intersect(SigsymbolsinPaths[[a]], SigsymbolsinPaths[[b]])) / min(length(SigsymbolsinPaths[[a]]), length(SigsymbolsinPaths[[b]]))
        if(a <= b) sigPathCor[a,b] <- length(intersect(SigsymbolsinPaths[[a]], SigsymbolsinPaths[[b]])) / max(length(SigsymbolsinPaths[[a]]), length(SigsymbolsinPaths[[b]]))
      }
    }
    cols <- colorRampPalette(c("white", "red"))(n = 50)
    oo <- as.dendrogram(hclust(dist(sigPathCor)))
  })
  
  # Overlap renderPlot
  output$overlap <- renderPlot(
    heatmap.2(sigPathCor, scale='none', trace='none', col=cols, key=TRUE, 
              keysize=0.85, mar=c(26,26),cexRow = 0.8, cexCol=0.8,
              main="Pathway overlap (genes): drug associated pathways", 
              Rowv = oo, Colv=oo)
  )
  
  
  
  # Heatmap renderPlot
  output$heatmap <- renderPlot({
    # heatmap generating code goes here. Only want the top 5-10
    # NOT ALL OF THIS NEEDS TO BE IN HERE, PUT EVERYTHING BUT THE LOOP OUTSIDE. ALSO NEED TO FIGURE OUT HOW TO ALLOW MULTIPLE OUTPUTS, BUT ONLY 5, WHAT TO DO IF THERE'S LESS THAN 5 PATHWAYS,
    # AND HOW TO SKIP PATHWAYS THAT ONLY HAVE ONE GENE IN THEM, CAUSE THESE CAN'T BE TURNED INTO HEATMAPS (THAT WOULD JUST BE SILLY)
    # PROBABLY ALSO SOME OTHER IMPORTANT STUFF THAT NEEDS SORTING OUT
    # for(k in 1:length(ensginPaths)){
    for(k in 1){
      # Fetches expression data
      pathway_expDat_sort <- subset(expDat_sort, rownames(expDat_sort) %in% ensginPaths[[k]])
      rownames(pathway_expDat_sort) <- match(rownames(pathway_expDat_sort), GeneLabelTool$ENSEMBL) %>% 
        GeneLabelTool$SYMBOL[.]
      zz <- apply(pathway_expDat_sort, 1, scale, scale=TRUE) %>% t()
      colnames(zz) <- colnames(pathway_expDat_sort)
      # Order zz in the same order as ensginPaths[[k]]
      orderer <- as.vector(symbolsinPaths[[k]])
      new_order <- sapply(orderer, function(x,zz){which(rownames(zz) == x)}, zz=zz)
      zz <- zz[new_order,]
      # Orders rows and sorts out colours
      rowcc <- rep("red", length(ensginPaths[[k]]))
      rowcc[which(ensginPaths[[k]] %in% SigensginPaths[[k]])] = "green"
      roword <- as.factor(rowcc) %>% as.numeric(.) %>% order(.)
      # Changes row label font size based on how many rows there are, to avoid cluttering
      fonthelp <- length(ensginPaths[[k]])
      fontdefault <- 0.2 + 1/log10(fonthelp)
      l <- ifelse(fonthelp>200, 0.4, fontdefault)
      # Standardises extreme values for better imaging, and ensures heat maps aren't created for pathways with only one enriched gene
      if(nrow(zz)>1){
        zz[zz > 3] <- 3
        zz[zz < -3] <- -3
        # Creates the heatmaps
        pn <- gsub(" ","_",goseqPathways$Pathway[k])
        png(paste0("Heatmaps/", pn, "_heatmap.png"), width = 20, height = 20, units = 'in', res = 500)
        heatmap.mik(zz[roword,colord], trace='none', scale='none', col=bluered(50), ColSideColors=cc[,colord], Colv=FALSE, RowSideColors=rowcc[roword], Rowv=FALSE, mar=c(4,12), keysize=1, cexRow=l, main=(goseqPathways$Pathway[k]))
        legend(0.4,0.97,
               c("Breast", "High -log IC50", "High AUC", "Significant", "Large Intestine", "Low -log IC50", "Low AUC", "Not Significant", "Stomach", "", "", ""),
               fill=c("red","yellow","yellow","green","green","blue","blue","red","blue","white","white","white"), border="white", ncol=3)
        dev.off()
      }
    }
  })  
  
  
}



# Run the application 
shinyApp(ui = ui, server = server)