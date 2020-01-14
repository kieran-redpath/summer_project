#' ---
#' title: Heatmap example
#' author: Mik Black
#' date: 9 January 2020
#' output: html_document
#' ---

## Load gplots packages for heatmap.2 code (and associated helper functions)
library(gplots)

## Load custom code for added additional rows via ColSideColors
source('heatmap-mik.R')

## Load GSE62254 data set
load('gse62254_gastric_cancer.RData')

## list of genes to plot
emtGenes <- c("SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "CDH1", "CLDN4",
              "CLDN7", "TJP3", "MUC1")

## Extract expression data for those genes
gDat <- gse62254_expDat[match(emtGenes, rownames(gse62254_expDat)), ]
gDat[1:5,1:5]

## Standardize expression values for plotting
zz <- t(apply(gDat, 1, scale, scale=T))
zz[zz< -3] <- -3
zz[zz> 3] <- 3

## Molecular subtype data
table(gse62254_clinDat$molSub)

## Color by molecular subtype
cc <- c("black", "blue", "red", "green")[as.numeric(as.factor(gse62254_clinDat$molSub))]

## Check colours vs subtype
table(gse62254_clinDat$molSub, cc)

## standard heatmap plot
heatmap.2(zz, scale='none', trace='none', col='bluered', ColSideColors = cc)

## sort by molecular subtype
ord <- order(cc)
heatmap.2(zz[,ord], scale='none', trace='none', col='bluered', ColSideColors = cc[ord], Colv=FALSE)

## add additional row to cc (Lauren Classification)
dd <- rbind(c("black", "blue", "red", "green")[as.numeric(as.factor(gse62254_clinDat$molSub))],
            c("brown", "purple2", "orange")[as.numeric(as.factor(gse62254_clinDat$lc))])

## can now add row names as labels
rownames(dd) <- c("MolSub", "LC")

## Check colours
table(gse62254_clinDat$molSub, dd["MolSub",])
table(gse62254_clinDat$lc, dd["LC",])

## standard heatmap.mik plot
heatmap.mik(zz, scale='none', trace='none', col='bluered', ColSideColors = dd)

## with ordering by molecular subtype
ord <- order(dd["MolSub",])
heatmap.mik(zz[,ord], scale='none', trace='none', col='bluered', ColSideColors = dd[,ord], Colv=FALSE)

## with ordering by Lauren Classification
ord <- order(dd["LC",])
heatmap.mik(zz[,ord], scale='none', trace='none', col='bluered', ColSideColors = dd[,ord], Colv=FALSE)

## with ordering by Lauren Classification and then Moleclar Subtype (a bit kludgy...)
ord <- order( 100*as.numeric(as.factor(gse62254_clinDat$lc)) + 
                as.numeric(as.factor(gse62254_clinDat$molSub)) )

heatmap.mik(zz[,ord], scale='none', trace='none', col='bluered', ColSideColors = dd[,ord], Colv=FALSE)
                            