Breast\_Colorectal\_Stomach\_Heatmaps\_06\_01\_20
================
Kieran Redpath
6 January 2020

NOTE: Don’t run this at the same time as
“dasatinib\_dif\_exp\_analysis\_3\_tissues”, cause I’m lazy and didn’t
change the names of objects

limma differential expression analysis (lecture 3) look at high and low
ic50 groups, and then find out which genes are differentially expressed
between the two groups data frame with all the samples that have
dasatinib data, and all the genes. These will have high/low ic50
associated with them (50/50 split) Will show you which genes are
involved in dasatinib sensitivity fit model with lmfit and ebayes, gives
you a table with p values, showing significance of differential
expression volcano plot tells you how big and significant the difference
is, stuff that comes out the top will be the most important Pathway
analysis to see what it’s actually changing

``` r
library(limma)
library(edgeR)
library(data.table)
library(magrittr)
library(ggplot2)
library(CePa)
library(tidyverse)
```

    ## -- Attaching packages ------------------------------------------------------------------------------------------------- tidyverse 1.3.0 --

    ## v tibble  2.1.3     v dplyr   0.8.3
    ## v tidyr   1.0.0     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.4.0
    ## v purrr   0.3.3

    ## -- Conflicts ---------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::between()   masks data.table::between()
    ## x tidyr::extract()   masks magrittr::extract()
    ## x dplyr::filter()    masks stats::filter()
    ## x dplyr::first()     masks data.table::first()
    ## x dplyr::lag()       masks stats::lag()
    ## x dplyr::last()      masks data.table::last()
    ## x purrr::set_names() masks magrittr::set_names()
    ## x tidyr::spread()    masks CePa::spread()
    ## x purrr::transpose() masks data.table::transpose()

``` r
library(ggbeeswarm)
library(org.Hs.eg.db)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     first, second

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## The following object is masked from 'package:data.table':
    ## 
    ##     shift

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 

``` r
library(dplyr)
```

``` r
# Output suppressed cause setDT prints a whole lot of useless stuff
# expDat Processing
ccleExpData <- read.gct('CCLE_GDSC/Data/CCLE_RNAseq_genes_counts_20180929.gct')
logExpDat <- log(ccleExpData + 0.5)
dge <- DGEList(counts=ccleExpData)
dge <- calcNormFactors(dge)
dge_voom <- voom(dge, plot=TRUE)
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Create%20the%20first%20big%20file%20to%20use-1.png)<!-- -->

``` r
expDat <- dge_voom$E
# This keeps the tissue names for if you need them later
expDatTissues <- dge_voom$E

# GDSC2 Processing (from old document)
GDSC2 <- fread('CCLE_GDSC/Data/GDSC2_fitted_dose_response_15Oct19.csv', sep= ';')
dim(GDSC2)
dasatinib <- GDSC2 %>% filter(., DRUG_NAME=="Dasatinib")
expDatsplitnames <- strsplit(colnames(expDat), "_") %>% 
  lapply(., function(x) x[1]) %>% 
  unlist()
expDatsplit <- expDat
colnames(expDat) <- expDatsplitnames

dasatinib$CELL_LINE_NAME <- gsub("-","",dasatinib$CELL_LINE_NAME, fixed=TRUE)
commonSamples <- intersect(dasatinib$CELL_LINE_NAME,colnames(expDat))
expDat_match <- match(commonSamples, colnames(expDat))
expDat_sort <- expDat[ , na.omit(expDat_match)]
dasatinib_match <- match(commonSamples, dasatinib$CELL_LINE_NAME)
dasatinib_sort <- dasatinib[na.omit(dasatinib_match),]
colnames(dasatinib_sort)[5] <- "CCLE_Name"
dasatinib_sort$AUC_Level <- ifelse(dasatinib_sort$AUC < quantile(dasatinib_sort$AUC, 0.5, na.rm=TRUE), "Low","High")

# Filter "dasatinib_sort" so that you're only looking at 3 tissues for the rest of the document
dasatinib_sort <- dasatinib_sort %>% filter(., TCGA_DESC=="COREAD" | TCGA_DESC=="BRCA" | TCGA_DESC=="STAD")  

# expDat Processing
expDat_all <- t(expDat_sort)
expDat_all <- as.data.frame(expDat_all)
setDT(expDat_all, keep.rownames = "CCLE_Name")[]
# "dasatinib_sort" can now be combined with "expDat_sort"
dasatinib_all = full_join(x=dasatinib_sort, y= expDat_all, by=c("CCLE_Name" = "CCLE_Name"))

# Put together some sort of tool that will tell you what tissue a cell line is, to end this confusion. using expDatsplit, as it's an exact copy of expDat before removing the tissue info
tissuetool <-data.frame(Cell_Line=character(1019),Tissue_Type=character(1019))
tissuetool$Cell_Line <- strsplit(colnames(expDatsplit), "_") %>% 
  lapply(., function(x) x[1]) %>% 
  unlist()
tissuetool$Tissue_Type <- strsplit(colnames(expDatsplit), "_") %>% 
  lapply(., function(x) x[-1])%>%
  paste(., sep="_", collapse = NULL) %>% 
  unlist()
# Now want to remove weird characters and make it prettier
tissuetool$Tissue_Type <- gsub('c', '', tissuetool$Tissue_Type)
tissuetool$Tissue_Type <- gsub('[[:punct:]]', '', tissuetool$Tissue_Type)
tissuetool$Tissue_Type <- gsub(' ', '_', tissuetool$Tissue_Type)

# Create a cleaner version, that's only the tissues of interest and also only cell lines present in commonSamples
tissuetoolsort <- filter(tissuetool, Tissue_Type=="BREAST" | Tissue_Type=="STOMACH" | Tissue_Type=="LARGE_INTESTINE")

int <- intersect(dasatinib_sort$CCLE_Name, tissuetoolsort$Cell_Line)
tissuetool_match <- match(int, tissuetoolsort$Cell_Line)
tissuetoolsort <- tissuetoolsort[na.omit(tissuetool_match) , ]

# Now use this to filter expDat_sort
expDat_sort <- expDat_sort %>% subset(., select=which(tissuetoolsort$Cell_Line %in% colnames(expDat_sort)))
# Cool, now dasatinib_sort and expDat_sort both only refer to 89 cell lines
```

``` r
# Check all this stuff
dim(expDat_sort)
```

    ## [1] 56202    89

``` r
dim(dasatinib_sort)
```

    ## [1] 89 20

``` r
head(expDat_sort)[,1:5]
```

    ##                        SW48    SW1710     SW620     CAL51   TOV112D
    ## ENSG00000223972.4 -2.487965 -2.424805 -4.965088 -5.718232 -4.656524
    ## ENSG00000227232.4  3.255119  2.689605  3.393930  3.233053  3.080646
    ## ENSG00000243485.2 -3.224930 -3.462279 -4.965088 -5.718232 -4.946031
    ## ENSG00000237613.2 -2.608259 -3.066351 -4.434573 -5.232805 -5.308601
    ## ENSG00000268020.2 -3.961896 -4.250775 -5.813085 -5.718232 -5.308601
    ## ENSG00000240361.1 -2.177624 -4.540282 -4.434573 -6.455197 -5.308601

``` r
head(dasatinib_sort)[,1:5]
```

    ##   DATASET NLME_RESULT_ID NLME_CURVE_ID COSMIC_ID   CCLE_Name
    ## 1   GDSC2            290      14783021    909751        SW48
    ## 2   GDSC2            290      14783023    905962       SW620
    ## 3   GDSC2            290      14783024    910927       CAL51
    ## 4   GDSC2            290      14783037    908445        SNU5
    ## 5   GDSC2            290      14783041    909718       SKCO1
    ## 6   GDSC2            290      14783049    908120 MDAMB175VII

``` r
sum(colnames(expDat_sort)==dasatinib_sort$CCLE_Name)
```

    ## [1] 1

``` r
names(dasatinib_sort)
```

    ##  [1] "DATASET"         "NLME_RESULT_ID"  "NLME_CURVE_ID"   "COSMIC_ID"      
    ##  [5] "CCLE_Name"       "SANGER_MODEL_ID" "TCGA_DESC"       "DRUG_ID"        
    ##  [9] "DRUG_NAME"       "PUTATIVE_TARGET" "PATHWAY_NAME"    "COMPANY_ID"     
    ## [13] "WEBRELEASE"      "MIN_CONC"        "MAX_CONC"        "LN_IC50"        
    ## [17] "AUC"             "RMSE"            "Z_SCORE"         "AUC_Level"

``` r
# lmFit
group <- ifelse(dasatinib_sort$LN_IC50 > median(dasatinib_sort$LN_IC50), "High", "Low")
table(group)
```

    ## group
    ## High  Low 
    ##   44   45

``` r
boxplot(dasatinib_sort$LN_IC50 ~ group)
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20IC50-1.png)<!-- -->

``` r
design = model.matrix(~group);
design %>% head()
```

    ##   (Intercept) groupLow
    ## 1           1        1
    ## 2           1        0
    ## 3           1        1
    ## 4           1        0
    ## 5           1        1
    ## 6           1        0

``` r
colnames(design) = c("Mean"
,"HighVsLow"
)

# Fit the data to the model
fit = lmFit(expDat_sort, design)
fit = eBayes(fit)
tt = topTable(fit, coef="HighVsLow", adjust="BH",n=nrow(expDat_sort))
options(digits=4)
tt[1:20,]
```

    ##                     logFC AveExpr      t   P.Value adj.P.Val      B
    ## ENSG00000178734.4  1.5040 -6.2460  4.432 2.654e-05    0.4564 -1.930
    ## ENSG00000241739.1 -1.0466 -6.1297 -4.278 4.742e-05    0.4564 -2.098
    ## ENSG00000260689.2  0.9576 -3.6336  4.151 7.575e-05    0.4564 -2.235
    ## ENSG00000249242.3  2.0343 -0.0662  4.051 1.086e-04    0.4564 -2.340
    ## ENSG00000241088.2 -1.0050 -6.9632 -3.983 1.389e-04    0.4564 -2.412
    ## ENSG00000197157.6 -0.4814  8.1291 -3.960 1.506e-04    0.4564 -2.436
    ## ENSG00000213309.3 -1.0145 -4.0815 -3.943 1.601e-04    0.4564 -2.454
    ## ENSG00000134046.7 -0.5946  5.7289 -3.934 1.651e-04    0.4564 -2.463
    ## ENSG00000142513.4  1.0871 -4.4066  3.893 1.912e-04    0.4564 -2.506
    ## ENSG00000230721.1 -0.8954 -3.2931 -3.882 1.983e-04    0.4564 -2.517
    ## ENSG00000239917.3 -0.8015 -3.6744 -3.856 2.179e-04    0.4564 -2.545
    ## ENSG00000244050.2  1.4653 -5.7474  3.854 2.191e-04    0.4564 -2.546
    ## ENSG00000186973.6  1.9809 -3.4341  3.842 2.288e-04    0.4564 -2.559
    ## ENSG00000256304.1  1.3702 -0.6076  3.809 2.565e-04    0.4564 -2.593
    ## ENSG00000260582.1  1.1203 -6.3847  3.770 2.929e-04    0.4564 -2.632
    ## ENSG00000213873.4 -0.7457 -7.4602 -3.760 3.034e-04    0.4564 -2.642
    ## ENSG00000261313.1 -0.4783 -7.8201 -3.737 3.282e-04    0.4564 -2.665
    ## ENSG00000074211.9  2.6775  0.9198  3.722 3.460e-04    0.4564 -2.681
    ## ENSG00000213676.6 -0.4520  6.5023 -3.714 3.552e-04    0.4564 -2.689
    ## ENSG00000237357.1  1.7050 -3.5188  3.670 4.134e-04    0.4564 -2.734

``` r
# Take a quick look at the first gene for the hell of it
plot(density(expDat_sort[,1]))
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20IC50-2.png)<!-- -->

``` r
# Count number of significant genes
sum(tt$adj.P.Val<0.01)
```

    ## [1] 0

``` r
## Plot log fold-change _versus_ -log(P-value). This graph is the important one, things up the top and spread out are useful associations
## (i.e., higher number = lower p-value):
volcanoplot(fit, coef="HighVsLow")

# Fancier Volcano plot, with lines for adjusted p-value, fold change, and highlighted important samples.
sigFC = (tt$adj.P.Val < 0.01)  & (abs(tt$logFC) > 1)

volcanoplot(fit, coef="HighVsLow")
points(tt$logFC[which(sigFC)], 
       -log10(tt$P.Value[which(sigFC)]), 
       cex=0.6, col='red', pch=16)
abline(h = min(-log10(tt$P.Value[which(sigFC)])), lty=2, col='blue')
```

    ## Warning in min(-log10(tt$P.Value[which(sigFC)])): no non-missing arguments to
    ## min; returning Inf

``` r
abline(v = c(-1,1), lty=2, col='blue')
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20IC50-3.png)<!-- -->

``` r
split <- strsplit(rownames(tt),".", fixed=T) %>% lapply(., function(x) x[1]) %>% unlist()
geneNames <- AnnotationDbi::select(org.Hs.eg.db, keys = split, column = c("SYMBOL","GENENAME"), key="ENSEMBL")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
head(split)
```

    ## [1] "ENSG00000178734" "ENSG00000241739" "ENSG00000260689" "ENSG00000249242"
    ## [5] "ENSG00000241088" "ENSG00000197157"

``` r
# What does the top table look like?
dim(tt)
```

    ## [1] 56202     6

``` r
tt$symbol <- geneNames$SYMBOL[match(split, geneNames$ENSEMBL)]
setDT(tt, keep.rownames = TRUE)[]
```

    ##                        rn      logFC AveExpr          t   P.Value adj.P.Val
    ##     1:  ENSG00000178734.4  1.504e+00 -6.2460  4.432e+00 2.654e-05    0.4564
    ##     2:  ENSG00000241739.1 -1.047e+00 -6.1297 -4.278e+00 4.742e-05    0.4564
    ##     3:  ENSG00000260689.2  9.576e-01 -3.6336  4.151e+00 7.575e-05    0.4564
    ##     4:  ENSG00000249242.3  2.034e+00 -0.0662  4.051e+00 1.086e-04    0.4564
    ##     5:  ENSG00000241088.2 -1.005e+00 -6.9632 -3.983e+00 1.389e-04    0.4564
    ##    ---                                                                     
    ## 56198:  ENSG00000154237.8 -3.005e-05  3.2806 -5.512e-05 1.000e+00    1.0000
    ## 56199: ENSG00000125414.14  1.772e-05 -6.6127  4.522e-05 1.000e+00    1.0000
    ## 56200: ENSG00000158014.10 -2.552e-05 -3.7510 -4.494e-05 1.000e+00    1.0000
    ## 56201: ENSG00000138386.12  8.728e-06  5.4153  4.290e-05 1.000e+00    1.0000
    ## 56202:  ENSG00000170954.7 -2.743e-05 -0.3447 -3.635e-05 1.000e+00    1.0000
    ##             B   symbol
    ##     1: -1.930   LMO7DN
    ##     2: -2.098     <NA>
    ##     3: -2.235     <NA>
    ##     4: -2.340 TMEM150C
    ##     5: -2.412     <NA>
    ##    ---                
    ## 56198: -4.805    LRRK1
    ## 56199: -4.805     MYH2
    ## 56200: -4.805  SLC30A2
    ## 56201: -4.805     NAB1
    ## 56202: -4.805   ZNF415

``` r
# Check the data: Output shows that some ENSEMBL ID's have multiple regular gene names
length(split)
```

    ## [1] 56202

``` r
dim(geneNames)
```

    ## [1] 56436     3

``` r
# Print the top 300 associated genes. Can then look them up on Enrichr, GeneSetDB, etc
cat(na.omit(tt$symbol[1:300]),sep="\n")
```

    ## LMO7DN
    ## TMEM150C
    ## SND1
    ## MBD2
    ## ACP4
    ## FAM183A
    ## PPP2R2C
    ## ATF6B
    ## GYS2
    ## KRTAP17-1
    ## IMPDH2
    ## PIGN
    ## NADK2
    ## KIAA1328
    ## ADRA2C
    ## RPL17
    ## TP73
    ## PRH1
    ## NAA80
    ## OR14A16
    ## WASF2
    ## LINC01563
    ## ACSBG1
    ## LOC107986453
    ## MMP3
    ## TXNDC5
    ## FLG-AS1
    ## DOCK4-AS1
    ## CCDC173
    ## SNORA2A
    ## MIR377
    ## MPP2
    ## NMRK1
    ## SNORD48
    ## LOC102724265
    ## DYM
    ## ADGRE3
    ## MIR3200
    ## NIPAL4
    ## ADGRB2
    ## FAM66C
    ## CTB-178M22.2
    ## EDARADD
    ## CACNA1D
    ## LINC02652
    ## KCNIP3
    ## MMP17
    ## SHB
    ## MOGS
    ## STMN3
    ## WSB2
    ## MIR181B2
    ## RIOK3
    ## JAKMIP3
    ## FZD10-AS1
    ## PRLR
    ## ZNF217
    ## MIR874
    ## PCCB
    ## LOC100996447
    ## DPEP3
    ## APBA2
    ## DISC1-IT1
    ## ANKRD34C-AS1
    ## LOC101928516
    ## MRO
    ## DNAH9
    ## MSX2
    ## OR4S1
    ## CHN2
    ## LOC100506551
    ## LINC01701
    ## INHA
    ## ENO2
    ## CHRNA1
    ## HADHA
    ## SRRM3
    ## C2CD4C
    ## C18orf21
    ## VAX1
    ## NRIP2
    ## PSME1
    ## ARHGEF11
    ## PPP2R2B-IT1
    ## PRKAA2
    ## PSMG2
    ## OR14J1
    ## HCN2
    ## ELP2
    ## OBSL1
    ## TUBB2B
    ## ZCCHC12
    ## CITED1
    ## SLC35G5
    ## NOA1
    ## GABRR1
    ## SORCS3
    ## ZFYVE16
    ## HIC1
    ## ZNF703
    ## GNRHR2
    ## MGC27382
    ## RPN1
    ## FAM174B
    ## IRF1
    ## MEF2C
    ## BAZ2B
    ## DPYSL5
    ## DRD3
    ## AKNAD1
    ## SCAMP5
    ## LOC101927040

``` r
#Check log fold change
sign(tt$logFC)
```

    ##     [1]  1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1 -1  1
    ##    [25] -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1 -1
    ##    [49]  1  1 -1  1  1  1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1
    ##    [73]  1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1
    ##    [97] -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1  1 -1
    ##   [121] -1  1 -1  1  1  1 -1 -1  1 -1  1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1
    ##   [145]  1  1  1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1
    ##   [169] -1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1  1  1  1 -1 -1 -1  1  1  1 -1
    ##   [193]  1 -1 -1  1  1  1  1  1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1  1 -1
    ##   [217]  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1
    ##   [241] -1  1  1 -1  1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1
    ##   [265] -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1
    ##   [289] -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1  1 -1  1  1 -1  1  1 -1  1  1  1 -1
    ##   [313] -1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1
    ##   [337]  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1
    ##   [361] -1  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1
    ##   [385] -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1
    ##   [409]  1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1
    ##   [433] -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1  1  1
    ##   [457]  1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ##   [481]  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1
    ##   [505] -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ##   [529]  1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1
    ##   [553] -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1  1 -1 -1  1 -1  1
    ##   [577]  1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1  1  1  1  1 -1  1
    ##   [601]  1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1
    ##   [625]  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1
    ##   [649]  1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1
    ##   [673]  1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1
    ##   [697] -1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1  1  1  1  1 -1 -1  1  1 -1 -1 -1
    ##   [721]  1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1
    ##   [745] -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1
    ##   [769]  1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1
    ##   [793] -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ##   [817] -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1
    ##   [841]  1  1 -1 -1 -1  1 -1  1 -1  1 -1  1  1  1 -1  1  1 -1 -1  1  1 -1 -1 -1
    ##   [865] -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1
    ##   [889]  1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1
    ##   [913] -1  1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1
    ##   [937]  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1
    ##   [961] -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1
    ##   [985] -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1  1 -1
    ##  [1009] -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1
    ##  [1033] -1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [1057] -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1
    ##  [1081] -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1
    ##  [1105]  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1
    ##  [1129]  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1
    ##  [1153] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1  1
    ##  [1177] -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1 -1 -1  1
    ##  [1201] -1 -1 -1 -1  1 -1  1 -1  1  1 -1  1 -1  1 -1 -1  1 -1  1 -1  1  1 -1  1
    ##  [1225] -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [1249] -1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1
    ##  [1273] -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1
    ##  [1297] -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1
    ##  [1321] -1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1
    ##  [1345]  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [1369] -1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1
    ##  [1393] -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1
    ##  [1417] -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1
    ##  [1441] -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1
    ##  [1465]  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1
    ##  [1489] -1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1
    ##  [1513] -1 -1 -1  1  1 -1  1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1
    ##  [1537]  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1  1  1
    ##  [1561] -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1
    ##  [1585] -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ##  [1609]  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1
    ##  [1633]  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [1657] -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1
    ##  [1681] -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1
    ##  [1705]  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1  1 -1  1
    ##  [1729] -1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1
    ##  [1753]  1 -1  1  1 -1  1  1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1
    ##  [1777] -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1
    ##  [1801] -1  1  1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1
    ##  [1825] -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1  1
    ##  [1849] -1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1
    ##  [1873] -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1
    ##  [1897] -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1
    ##  [1921]  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1
    ##  [1945] -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1
    ##  [1969]  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1
    ##  [1993] -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1
    ##  [2017] -1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1
    ##  [2041] -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ##  [2065]  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1
    ##  [2089] -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1
    ##  [2113] -1 -1  1 -1  1 -1  1  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [2137] -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1  1  1 -1
    ##  [2161] -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1
    ##  [2185] -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1
    ##  [2209] -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1  1  1 -1 -1  1 -1  1 -1
    ##  [2233] -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1  1
    ##  [2257] -1 -1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1  1 -1
    ##  [2281] -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1
    ##  [2305] -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1
    ##  [2329] -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1
    ##  [2353] -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1
    ##  [2377]  1  1  1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [2401]  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1
    ##  [2425] -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1
    ##  [2449] -1 -1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1
    ##  [2473] -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1
    ##  [2497] -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1
    ##  [2521] -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1
    ##  [2545] -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1
    ##  [2569]  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1
    ##  [2593]  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1
    ##  [2617] -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1
    ##  [2641] -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1
    ##  [2665] -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1
    ##  [2689]  1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1  1  1 -1
    ##  [2713] -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1
    ##  [2737] -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1
    ##  [2761]  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [2785]  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1
    ##  [2809]  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ##  [2833]  1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1  1  1 -1  1 -1 -1  1  1 -1  1 -1
    ##  [2857] -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1
    ##  [2881] -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1  1 -1
    ##  [2905] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [2929] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1
    ##  [2953]  1 -1  1  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1
    ##  [2977] -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1
    ##  [3001] -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1
    ##  [3025]  1 -1  1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1 -1 -1  1  1  1
    ##  [3049] -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1  1
    ##  [3073] -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1
    ##  [3097] -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1
    ##  [3121] -1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1  1 -1
    ##  [3145] -1 -1  1 -1  1  1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1
    ##  [3169] -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3193] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1
    ##  [3217]  1 -1  1 -1 -1 -1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1
    ##  [3241] -1 -1 -1  1  1  1 -1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1
    ##  [3265]  1  1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1
    ##  [3289] -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1
    ##  [3313] -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1
    ##  [3337] -1 -1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1
    ##  [3361] -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3385] -1 -1  1 -1  1  1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1
    ##  [3409] -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ##  [3433] -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3457] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1
    ##  [3481] -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1  1  1 -1 -1  1 -1  1 -1
    ##  [3505] -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1  1 -1  1 -1
    ##  [3529]  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ##  [3553]  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1  1 -1
    ##  [3577] -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1
    ##  [3601] -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ##  [3625] -1 -1 -1  1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1
    ##  [3649] -1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1
    ##  [3673] -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1
    ##  [3697]  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3721] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3745] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3769] -1 -1 -1 -1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1
    ##  [3793]  1 -1 -1  1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3817] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3841] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3865] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3889] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3913] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3937] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [3961] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1  1
    ##  [3985] -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4009] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1  1  1 -1 -1  1  1
    ##  [4033]  1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1
    ##  [4057] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4081] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1
    ##  [4105] -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1
    ##  [4129]  1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1  1  1  1  1
    ##  [4153]  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4177] -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4201] -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1
    ##  [4225]  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1
    ##  [4249] -1  1 -1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1
    ##  [4273] -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4297] -1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1  1 -1  1 -1 -1  1 -1 -1  1 -1  1  1
    ##  [4321] -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1
    ##  [4345] -1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1
    ##  [4369] -1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1  1  1 -1 -1 -1
    ##  [4393]  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1
    ##  [4417] -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1
    ##  [4441]  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4465] -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4489] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1
    ##  [4513] -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1
    ##  [4537] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1
    ##  [4561] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4585] -1  1  1  1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1
    ##  [4609] -1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1
    ##  [4633] -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1  1
    ##  [4657] -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1  1  1  1 -1 -1  1 -1 -1 -1
    ##  [4681] -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ##  [4705]  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4729] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1
    ##  [4753]  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1
    ##  [4777] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4801] -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1  1 -1  1
    ##  [4825]  1  1  1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1  1  1  1  1  1  1  1 -1  1 -1
    ##  [4849] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ##  [4873] -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [4897] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1
    ##  [4921]  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1  1 -1
    ##  [4945]  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1
    ##  [4969] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1
    ##  [4993] -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1
    ##  [5017]  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1  1  1  1 -1 -1  1 -1
    ##  [5041] -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1
    ##  [5065] -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1
    ##  [5089]  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [5113] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1 -1
    ##  [5137]  1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1
    ##  [5161] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1
    ##  [5185]  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [5209] -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [5233] -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1  1  1  1 -1  1  1 -1 -1
    ##  [5257] -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1
    ##  [5281] -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1
    ##  [5305] -1  1  1  1  1  1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ##  [5329] -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1
    ##  [5353] -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ##  [5377] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1
    ##  [5401]  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [5425] -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1
    ##  [5449]  1  1 -1 -1  1  1  1  1  1 -1  1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1
    ##  [5473] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1
    ##  [5497] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [5521] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [5545] -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1  1  1 -1 -1 -1  1  1  1
    ##  [5569] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1
    ##  [5593] -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1  1  1 -1  1
    ##  [5617] -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [5641] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1
    ##  [5665] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1
    ##  [5689] -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [5713] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1
    ##  [5737] -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ##  [5761]  1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1  1
    ##  [5785] -1 -1 -1  1  1  1  1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1
    ##  [5809]  1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [5833] -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1
    ##  [5857] -1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1
    ##  [5881] -1  1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1
    ##  [5905] -1  1  1  1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1  1 -1 -1 -1  1 -1  1  1
    ##  [5929]  1 -1 -1  1  1  1  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1
    ##  [5953] -1  1 -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1
    ##  [5977] -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1 -1 -1  1  1 -1
    ##  [6001] -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1
    ##  [6025] -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1
    ##  [6049]  1  1  1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1
    ##  [6073] -1 -1  1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1  1  1  1  1 -1 -1
    ##  [6097] -1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ##  [6121]  1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1
    ##  [6145] -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1  1 -1 -1  1 -1 -1  1  1 -1  1 -1
    ##  [6169]  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1  1  1
    ##  [6193]  1 -1  1  1 -1  1  1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ##  [6217]  1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1
    ##  [6241] -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [6265] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1
    ##  [6289] -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1
    ##  [6313]  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1
    ##  [6337] -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1
    ##  [6361]  1  1  1  1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1  1 -1  1  1
    ##  [6385] -1 -1  1  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1
    ##  [6409] -1 -1 -1 -1  1  1 -1  1  1 -1  1  1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1
    ##  [6433]  1  1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1
    ##  [6457]  1 -1 -1 -1 -1  1  1  1  1 -1 -1  1  1  1  1 -1  1 -1 -1  1 -1 -1  1  1
    ##  [6481] -1 -1 -1  1  1 -1 -1 -1  1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1
    ##  [6505]  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1
    ##  [6529]  1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1
    ##  [6553]  1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1
    ##  [6577]  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1  1
    ##  [6601] -1 -1 -1  1  1  1 -1 -1  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1
    ##  [6625]  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1
    ##  [6649]  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1
    ##  [6673]  1 -1 -1  1  1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1
    ##  [6697] -1 -1 -1  1  1  1 -1  1 -1 -1  1  1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1
    ##  [6721]  1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1
    ##  [6745] -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1
    ##  [6769] -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ##  [6793] -1 -1 -1 -1  1  1  1  1  1 -1 -1  1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1
    ##  [6817]  1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1
    ##  [6841]  1  1  1  1 -1 -1  1  1  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1
    ##  [6865]  1  1 -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1
    ##  [6889] -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1
    ##  [6913]  1  1  1 -1  1 -1  1 -1  1  1  1  1  1 -1 -1  1  1  1  1 -1  1 -1  1 -1
    ##  [6937] -1  1 -1 -1 -1  1  1 -1  1  1 -1 -1  1  1  1 -1  1  1  1  1 -1 -1 -1 -1
    ##  [6961] -1 -1  1  1  1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1
    ##  [6985]  1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [7009] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ##  [7033] -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1
    ##  [7057] -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1  1  1 -1 -1 -1
    ##  [7081] -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1
    ##  [7105] -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1
    ##  [7129] -1  1  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1
    ##  [7153]  1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1  1  1  1  1 -1
    ##  [7177] -1 -1  1 -1  1 -1 -1  1  1 -1  1  1  1 -1 -1  1 -1  1  1 -1 -1  1  1 -1
    ##  [7201] -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1
    ##  [7225] -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [7249]  1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1  1  1 -1 -1 -1  1  1  1 -1 -1 -1  1
    ##  [7273]  1  1 -1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1
    ##  [7297]  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1
    ##  [7321] -1  1 -1  1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1
    ##  [7345] -1 -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1  1 -1  1  1  1 -1 -1 -1 -1  1 -1
    ##  [7369] -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ##  [7393] -1 -1  1  1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1
    ##  [7417] -1  1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1
    ##  [7441]  1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ##  [7465] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [7489] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1
    ##  [7513]  1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1 -1
    ##  [7537] -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1  1  1  1 -1 -1  1  1 -1
    ##  [7561] -1 -1  1  1 -1  1 -1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1
    ##  [7585]  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1
    ##  [7609] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1
    ##  [7633]  1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1
    ##  [7657] -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1
    ##  [7681] -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1
    ##  [7705]  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ##  [7729] -1 -1  1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1
    ##  [7753]  1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ##  [7777]  1  1 -1  1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1
    ##  [7801]  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1
    ##  [7825]  1 -1  1  1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1  1  1 -1  1 -1 -1  1 -1  1
    ##  [7849] -1 -1  1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1
    ##  [7873] -1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1
    ##  [7897] -1 -1  1  1 -1  1 -1  1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1
    ##  [7921]  1 -1 -1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1  1 -1  1  1 -1 -1  1 -1  1
    ##  [7945] -1 -1 -1  1 -1  1  1  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1
    ##  [7969]  1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1
    ##  [7993] -1 -1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1
    ##  [8017] -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1
    ##  [8041] -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1
    ##  [8065] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8089] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8113] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8137] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8161] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8185] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8209] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8233] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8257] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8281] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8305] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8329] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8353] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8377] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8401] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8425] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8449] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8473] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8497] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8521] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8545] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8569] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8593] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8617] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8641] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8665] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8689] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8713] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8737] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8761] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8785] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8809] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8833] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8857] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8881] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8905] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8929] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8953] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [8977] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9001] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9025] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9049] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9073] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9097] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9121] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9145] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9169] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9193] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9217] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9241] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9265] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9289] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9313] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9337] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9361] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9385] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9409] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9433] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9457] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9481] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9505] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9529] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9553] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9577] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9601] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9625] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9649] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9673] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9697] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9721] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9745] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9769] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9793] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9817] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9841] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9865] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9889] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9913] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9937] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9961] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ##  [9985] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10009] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10033] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10057] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10081] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10105] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10129] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10153] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10177] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10201] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10225] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10249] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10273] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10297] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10321] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10345] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10369] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10393] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10417] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10441] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10465] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10489] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10513] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10537] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10561] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10585] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10609] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10633] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10657] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10681] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10705] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10729] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10753] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10777] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10801] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10825] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10849] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10873] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10897] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10921] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10945] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10969] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [10993] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11017] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11041] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11065] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11089] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11113] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11137] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11161] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11185] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11209] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11233] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11257] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11281] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11305] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11329] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11353] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11377] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11401] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11425] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11449] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11473] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11497] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11521] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11545] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11569] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11593] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11617] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11641] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11665] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11689] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11713] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11737] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11761] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11785] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11809] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11833] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11857] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11881] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11905] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1
    ## [11929] -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [11953] -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1 -1
    ## [11977] -1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1  1 -1 -1  1 -1  1  1 -1
    ## [12001] -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1
    ## [12025]  1  1 -1 -1  1 -1  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1  1  1  1  1  1
    ## [12049] -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1
    ## [12073] -1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1  1  1 -1  1  1
    ## [12097] -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1
    ## [12121]  1 -1  1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1
    ## [12145]  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1  1  1 -1  1 -1 -1  1  1
    ## [12169] -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1
    ## [12193]  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1  1
    ## [12217] -1 -1 -1  1 -1  1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1
    ## [12241] -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1  1
    ## [12265] -1  1 -1 -1 -1  1  1 -1  1  1  1  1  1  1  1 -1 -1 -1 -1 -1  1 -1  1  1
    ## [12289]  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1
    ## [12313] -1  1  1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1
    ## [12337] -1  1  1 -1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1
    ## [12361] -1 -1  1  1 -1  1 -1  1 -1 -1  1 -1  1  1 -1  1 -1  1  1 -1  1 -1  1 -1
    ## [12385]  1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1
    ## [12409] -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1 -1
    ## [12433] -1 -1 -1 -1  1  1  1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1  1
    ## [12457]  1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [12481] -1  1  1  1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1
    ## [12505] -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1  1 -1 -1
    ## [12529] -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1  1  1  1 -1 -1  1
    ## [12553] -1  1  1 -1 -1  1 -1 -1  1  1  1  1 -1 -1  1  1  1 -1 -1  1 -1  1 -1  1
    ## [12577] -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1
    ## [12601]  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1  1  1 -1  1  1 -1 -1
    ## [12625] -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1  1  1 -1  1  1 -1 -1
    ## [12649] -1  1  1  1  1  1  1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [12673]  1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1
    ## [12697] -1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1
    ## [12721]  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1
    ## [12745]  1  1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1  1  1  1 -1
    ## [12769]  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1  1 -1 -1  1  1  1 -1 -1 -1
    ## [12793]  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1
    ## [12817] -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1
    ## [12841] -1 -1  1  1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1  1  1 -1
    ## [12865]  1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1
    ## [12889]  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1 -1  1  1  1 -1 -1  1
    ## [12913] -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1
    ## [12937] -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1
    ## [12961]  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1
    ## [12985] -1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1
    ## [13009]  1  1  1 -1  1 -1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [13033] -1 -1  1  1  1 -1  1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1
    ## [13057] -1 -1  1 -1  1  1 -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1  1  1
    ## [13081]  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1
    ## [13105]  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1
    ## [13129]  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1  1
    ## [13153] -1  1 -1 -1  1 -1  1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1
    ## [13177] -1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1
    ## [13201] -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1  1 -1  1  1 -1  1 -1  1  1 -1  1 -1
    ## [13225] -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1  1  1  1 -1 -1
    ## [13249] -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1
    ## [13273]  1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1 -1 -1  1 -1  1 -1
    ## [13297] -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1
    ## [13321] -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1
    ## [13345]  1 -1  1  1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1
    ## [13369]  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1
    ## [13393]  1  1 -1  1  1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1
    ## [13417] -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1
    ## [13441]  1  1 -1  1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1
    ## [13465] -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1  1
    ## [13489] -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [13513] -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1
    ## [13537] -1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1  1 -1 -1  1
    ## [13561]  1 -1  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1  1 -1
    ## [13585]  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1
    ## [13609]  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1
    ## [13633] -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [13657]  1  1  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1
    ## [13681]  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1  1  1 -1
    ## [13705] -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1  1
    ## [13729]  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1
    ## [13753] -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1
    ## [13777] -1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1
    ## [13801] -1  1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1  1
    ## [13825]  1  1  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1
    ## [13849] -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1  1  1
    ## [13873] -1 -1 -1  1 -1  1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1
    ## [13897] -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1
    ## [13921]  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1
    ## [13945]  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [13969] -1 -1 -1 -1 -1 -1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1
    ## [13993] -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1
    ## [14017] -1 -1  1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1
    ## [14041] -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [14065]  1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [14089] -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1
    ## [14113]  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1  1  1 -1
    ## [14137]  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1
    ## [14161] -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ## [14185]  1  1  1  1 -1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1 -1 -1  1 -1 -1 -1  1
    ## [14209]  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [14233] -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1
    ## [14257]  1 -1 -1  1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [14281]  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1
    ## [14305]  1  1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1
    ## [14329] -1 -1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [14353] -1  1 -1  1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1
    ## [14377] -1 -1 -1  1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1  1
    ## [14401] -1 -1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1
    ## [14425] -1 -1 -1  1 -1 -1  1  1  1 -1  1  1 -1  1  1 -1  1 -1 -1 -1  1 -1  1 -1
    ## [14449]  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1  1 -1  1
    ## [14473]  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1
    ## [14497]  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [14521] -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1
    ## [14545] -1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1
    ## [14569] -1 -1 -1  1 -1 -1 -1  1  1  1 -1  1  1  1 -1  1  1 -1  1 -1 -1 -1  1 -1
    ## [14593]  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1
    ## [14617] -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1
    ## [14641]  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1
    ## [14665]  1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1 -1 -1  1 -1  1  1 -1
    ## [14689] -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1
    ## [14713] -1 -1 -1  1 -1 -1 -1  1  1  1  1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1
    ## [14737]  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1
    ## [14761]  1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1
    ## [14785] -1  1 -1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1
    ## [14809] -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1  1
    ## [14833] -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [14857]  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1
    ## [14881]  1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1  1
    ## [14905]  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1  1  1 -1  1
    ## [14929]  1  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ## [14953] -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1  1 -1  1 -1 -1  1  1 -1
    ## [14977] -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1
    ## [15001]  1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1
    ## [15025] -1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1  1  1 -1
    ## [15049] -1  1  1  1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1  1 -1 -1  1  1 -1  1 -1
    ## [15073]  1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1  1 -1 -1  1 -1
    ## [15097] -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [15121] -1  1 -1 -1  1  1 -1 -1  1  1  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1
    ## [15145]  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1  1
    ## [15169] -1  1  1 -1  1 -1  1 -1  1 -1 -1  1  1  1  1  1 -1  1  1 -1 -1  1 -1 -1
    ## [15193]  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1  1
    ## [15217]  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1 -1  1 -1
    ## [15241] -1  1 -1 -1  1 -1  1  1  1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1
    ## [15265]  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1
    ## [15289] -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1
    ## [15313]  1 -1 -1  1 -1  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1
    ## [15337] -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1  1
    ## [15361] -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1
    ## [15385] -1 -1 -1  1 -1  1 -1  1 -1  1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [15409]  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1
    ## [15433] -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1
    ## [15457]  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1
    ## [15481] -1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1
    ## [15505] -1 -1 -1 -1 -1  1 -1 -1  1  1  1  1  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1
    ## [15529]  1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1  1
    ## [15553] -1  1  1 -1  1 -1  1 -1  1  1 -1 -1  1  1  1 -1  1  1 -1  1 -1 -1  1 -1
    ## [15577]  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1
    ## [15601] -1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1  1 -1 -1 -1  1 -1 -1  1
    ## [15625] -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1  1
    ## [15649] -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [15673] -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1
    ## [15697] -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1
    ## [15721]  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1
    ## [15745]  1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1
    ## [15769] -1 -1  1 -1 -1  1  1  1  1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1 -1 -1  1
    ## [15793]  1 -1 -1 -1 -1  1  1  1  1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1
    ## [15817] -1  1  1  1 -1  1 -1  1  1  1 -1  1 -1 -1  1  1  1 -1  1 -1 -1  1 -1 -1
    ## [15841] -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1
    ## [15865]  1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1  1 -1  1  1  1  1 -1 -1  1  1
    ## [15889] -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1
    ## [15913]  1 -1 -1  1 -1  1 -1 -1  1  1  1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1  1
    ## [15937]  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1
    ## [15961] -1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1
    ## [15985] -1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1
    ## [16009] -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1
    ## [16033] -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [16057]  1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1
    ## [16081] -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1  1 -1  1  1 -1
    ## [16105]  1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1
    ## [16129]  1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1
    ## [16153] -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1 -1
    ## [16177]  1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [16201] -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1
    ## [16225]  1 -1 -1  1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1
    ## [16249] -1 -1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1  1 -1 -1  1  1
    ## [16273] -1  1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1
    ## [16297]  1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1
    ## [16321]  1  1  1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1
    ## [16345] -1 -1  1 -1 -1 -1 -1  1  1 -1  1  1  1 -1  1  1 -1 -1  1  1  1  1 -1 -1
    ## [16369]  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1  1  1 -1
    ## [16393]  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1 -1  1  1 -1  1  1 -1  1 -1  1
    ## [16417] -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1 -1  1 -1  1  1
    ## [16441]  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1
    ## [16465]  1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1  1 -1
    ## [16489]  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1  1  1  1  1  1 -1  1  1
    ## [16513]  1  1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1  1  1  1 -1
    ## [16537] -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1  1  1  1  1 -1  1  1 -1  1
    ## [16561] -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [16585] -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1
    ## [16609]  1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1  1  1 -1 -1 -1 -1  1 -1  1
    ## [16633]  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1  1
    ## [16657] -1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1  1  1  1 -1 -1 -1
    ## [16681]  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1
    ## [16705] -1  1  1  1  1  1  1 -1  1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1
    ## [16729] -1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1
    ## [16753] -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1
    ## [16777]  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1
    ## [16801] -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1
    ## [16825]  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1 -1 -1  1 -1  1
    ## [16849] -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1 -1
    ## [16873] -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1
    ## [16897]  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1  1  1  1  1
    ## [16921]  1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1
    ## [16945] -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1
    ## [16969] -1  1 -1  1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1
    ## [16993]  1  1  1 -1  1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1  1
    ## [17017] -1  1  1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1
    ## [17041]  1  1  1 -1  1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1
    ## [17065]  1  1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1
    ## [17089]  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1  1  1 -1 -1  1 -1
    ## [17113]  1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [17137] -1 -1  1 -1  1  1 -1  1 -1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1
    ## [17161]  1  1  1  1  1 -1  1  1 -1  1  1  1  1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1
    ## [17185] -1  1 -1 -1  1  1 -1 -1  1  1  1 -1  1 -1  1  1  1  1 -1  1 -1 -1  1  1
    ## [17209] -1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1  1
    ## [17233]  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1  1 -1  1 -1  1
    ## [17257]  1  1  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1
    ## [17281]  1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1
    ## [17305] -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1  1 -1 -1
    ## [17329] -1 -1  1  1 -1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1
    ## [17353] -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1
    ## [17377] -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [17401]  1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [17425]  1 -1 -1  1  1  1  1 -1  1  1 -1  1  1 -1  1  1  1  1 -1  1 -1 -1 -1  1
    ## [17449]  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1
    ## [17473]  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1
    ## [17497] -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1  1  1
    ## [17521]  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ## [17545] -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [17569] -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1
    ## [17593] -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1  1 -1 -1  1 -1 -1
    ## [17617] -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [17641]  1  1  1 -1  1 -1 -1  1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1
    ## [17665] -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1
    ## [17689] -1  1 -1  1 -1 -1  1 -1 -1  1  1  1  1 -1  1  1  1  1 -1 -1 -1  1 -1 -1
    ## [17713] -1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1  1 -1  1  1 -1
    ## [17737]  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1 -1 -1
    ## [17761] -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1  1 -1  1
    ## [17785]  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1
    ## [17809] -1 -1  1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1  1  1 -1
    ## [17833] -1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1
    ## [17857] -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1
    ## [17881] -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1
    ## [17905] -1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1
    ## [17929]  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1
    ## [17953] -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1
    ## [17977] -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ## [18001] -1  1 -1 -1  1 -1  1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1 -1
    ## [18025] -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1 -1 -1
    ## [18049] -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1  1 -1 -1 -1  1
    ## [18073] -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1
    ## [18097] -1  1  1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1
    ## [18121]  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1  1
    ## [18145]  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1  1 -1 -1
    ## [18169]  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1  1  1 -1 -1  1 -1  1  1  1
    ## [18193] -1  1 -1 -1  1  1  1 -1  1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1
    ## [18217] -1 -1  1  1  1  1  1  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1  1
    ## [18241] -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1
    ## [18265] -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [18289] -1 -1  1  1  1 -1  1  1  1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1
    ## [18313]  1  1  1 -1  1  1 -1 -1  1  1  1  1 -1  1  1 -1 -1  1 -1  1 -1  1  1 -1
    ## [18337]  1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1
    ## [18361] -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1
    ## [18385]  1  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1
    ## [18409]  1  1  1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1
    ## [18433]  1  1 -1  1 -1 -1  1 -1  1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1
    ## [18457] -1  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1
    ## [18481]  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1  1 -1 -1  1 -1 -1  1
    ## [18505]  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1  1 -1  1 -1 -1
    ## [18529]  1  1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1  1
    ## [18553] -1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1
    ## [18577]  1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1 -1 -1  1  1
    ## [18601] -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1
    ## [18625]  1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1
    ## [18649] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [18673] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1
    ## [18697]  1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [18721] -1  1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1
    ## [18745] -1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1  1  1 -1 -1  1 -1 -1  1  1 -1  1 -1
    ## [18769]  1 -1  1  1 -1 -1  1 -1 -1  1 -1  1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1
    ## [18793]  1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1  1
    ## [18817]  1  1  1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1  1  1  1 -1
    ## [18841]  1  1  1  1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1  1  1  1  1 -1  1
    ## [18865]  1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1  1 -1  1
    ## [18889] -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1  1  1  1  1  1  1  1 -1
    ## [18913] -1  1 -1  1  1 -1 -1 -1  1 -1  1  1 -1  1 -1  1  1  1  1  1 -1  1 -1  1
    ## [18937] -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1  1  1 -1 -1 -1 -1  1 -1
    ## [18961] -1  1  1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1  1  1
    ## [18985] -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1
    ## [19009]  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1  1  1
    ## [19033] -1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1
    ## [19057]  1  1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1  1  1  1  1 -1  1  1  1 -1 -1
    ## [19081] -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1
    ## [19105]  1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1  1
    ## [19129] -1  1 -1  1  1 -1  1  1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1  1 -1
    ## [19153]  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1
    ## [19177] -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1
    ## [19201]  1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1
    ## [19225] -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1  1  1  1  1 -1 -1
    ## [19249] -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1
    ## [19273] -1 -1  1  1  1 -1  1 -1 -1  1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1 -1  1
    ## [19297]  1  1  1  1  1 -1 -1 -1 -1  1  1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1
    ## [19321] -1 -1  1  1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1
    ## [19345] -1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1
    ## [19369] -1  1  1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1
    ## [19393] -1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1
    ## [19417] -1 -1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1
    ## [19441]  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1  1 -1 -1  1
    ## [19465]  1  1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1
    ## [19489] -1  1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1
    ## [19513] -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1
    ## [19537]  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1 -1 -1  1  1  1  1
    ## [19561]  1  1 -1  1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1
    ## [19585]  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1
    ## [19609] -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1
    ## [19633] -1  1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1  1
    ## [19657] -1  1  1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1
    ## [19681]  1  1  1 -1  1  1  1 -1  1 -1 -1  1 -1 -1  1  1  1  1  1 -1 -1  1 -1  1
    ## [19705]  1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [19729] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [19753] -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [19777] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1
    ## [19801]  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1  1 -1  1 -1
    ## [19825]  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [19849]  1  1 -1 -1 -1 -1 -1  1  1  1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1
    ## [19873]  1  1 -1  1 -1 -1  1  1  1 -1  1 -1 -1  1  1  1 -1 -1  1  1  1  1  1 -1
    ## [19897] -1  1 -1  1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1
    ## [19921] -1  1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1
    ## [19945] -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1
    ## [19969]  1 -1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1
    ## [19993] -1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1 -1  1  1 -1 -1  1 -1 -1  1
    ## [20017]  1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [20041] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [20065] -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1
    ## [20089] -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [20113] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1
    ## [20137] -1  1  1  1 -1  1  1 -1  1  1  1 -1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1
    ## [20161] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [20185] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1  1  1 -1  1  1  1
    ## [20209]  1  1  1 -1 -1  1  1  1 -1  1  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1  1
    ## [20233] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [20257] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [20281]  1 -1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1
    ## [20305]  1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1  1  1 -1
    ## [20329] -1 -1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1  1  1 -1 -1  1 -1 -1  1 -1
    ## [20353] -1 -1  1 -1  1 -1 -1  1 -1 -1  1  1  1  1 -1  1  1  1 -1  1  1  1  1 -1
    ## [20377] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1
    ## [20401]  1 -1  1 -1 -1 -1  1  1 -1  1  1  1 -1 -1  1  1  1 -1  1  1  1  1  1 -1
    ## [20425] -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1  1  1 -1 -1  1 -1
    ## [20449] -1  1 -1 -1  1 -1  1  1 -1  1 -1  1  1  1  1 -1  1 -1 -1 -1  1 -1  1  1
    ## [20473] -1 -1  1  1 -1  1  1  1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1
    ## [20497] -1  1  1  1 -1  1 -1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1  1
    ## [20521] -1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1  1  1  1 -1  1 -1 -1  1  1  1 -1  1
    ## [20545] -1  1 -1 -1 -1  1  1  1  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1
    ## [20569] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [20593]  1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1
    ## [20617]  1  1  1  1 -1  1  1 -1 -1  1 -1  1  1  1 -1 -1  1  1  1 -1 -1  1 -1  1
    ## [20641] -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1
    ## [20665] -1  1  1  1  1 -1 -1  1 -1 -1 -1  1  1  1  1 -1  1  1 -1  1  1  1  1 -1
    ## [20689] -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1  1 -1  1  1  1 -1 -1 -1  1  1
    ## [20713]  1  1  1 -1  1 -1  1  1 -1  1 -1  1  1  1 -1  1  1 -1 -1  1 -1 -1  1  1
    ## [20737] -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [20761] -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1  1  1  1 -1 -1  1  1 -1 -1
    ## [20785]  1  1 -1  1 -1 -1  1 -1  1 -1  1 -1  1 -1  1  1 -1  1  1  1 -1  1 -1  1
    ## [20809] -1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1 -1  1  1 -1 -1
    ## [20833]  1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1 -1  1  1
    ## [20857] -1 -1 -1  1 -1  1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1  1  1  1 -1 -1
    ## [20881] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [20905] -1 -1 -1  1  1  1  1 -1  1 -1 -1  1  1 -1  1  1 -1  1  1  1 -1 -1 -1 -1
    ## [20929] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1  1 -1
    ## [20953] -1 -1  1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1
    ## [20977] -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1 -1  1  1
    ## [21001] -1 -1 -1 -1 -1 -1  1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21025] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21049] -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1  1 -1
    ## [21073]  1 -1 -1 -1 -1  1  1  1  1  1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [21097] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1  1  1  1 -1  1  1 -1
    ## [21121] -1  1 -1  1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1 -1  1 -1  1  1  1  1
    ## [21145] -1 -1 -1  1  1 -1  1  1 -1 -1  1 -1  1 -1  1  1  1 -1 -1  1 -1 -1  1  1
    ## [21169]  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21193] -1 -1 -1 -1 -1 -1  1  1 -1  1  1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ## [21217] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1  1
    ## [21241] -1 -1 -1  1  1  1 -1  1 -1 -1 -1  1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1
    ## [21265]  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21289] -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21313] -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1  1  1 -1
    ## [21337] -1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21361] -1 -1 -1 -1  1  1  1  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21385] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1
    ## [21409] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1  1
    ## [21433] -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [21457] -1 -1 -1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1  1 -1  1  1 -1  1 -1 -1 -1
    ## [21481] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21505] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1
    ## [21529] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1  1 -1 -1
    ## [21553] -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1
    ## [21577] -1 -1  1  1 -1  1  1 -1  1  1  1 -1  1 -1  1  1  1 -1  1  1 -1  1 -1  1
    ## [21601] -1 -1  1  1  1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1
    ## [21625]  1  1 -1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1
    ## [21649] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [21673] -1 -1  1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1  1  1
    ## [21697] -1 -1 -1  1 -1  1  1 -1  1  1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1  1 -1
    ## [21721] -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [21745]  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21769] -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1  1  1
    ## [21793] -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1
    ## [21817]  1  1  1  1 -1  1  1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1
    ## [21841] -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [21865] -1 -1 -1  1  1 -1  1  1  1  1 -1  1  1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1
    ## [21889] -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1
    ## [21913] -1 -1 -1  1 -1  1  1 -1  1 -1  1  1 -1  1 -1  1  1 -1 -1  1  1  1 -1 -1
    ## [21937] -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1  1
    ## [21961]  1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1  1  1 -1  1 -1  1  1  1  1 -1  1  1
    ## [21985]  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1
    ## [22009] -1 -1 -1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1
    ## [22033] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1
    ## [22057] -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [22081] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1  1  1  1 -1 -1
    ## [22105] -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1  1  1  1 -1
    ## [22129]  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1
    ## [22153] -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1
    ## [22177]  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1
    ## [22201] -1 -1  1 -1 -1  1  1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1
    ## [22225] -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1  1  1  1  1  1
    ## [22249]  1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1
    ## [22273] -1 -1 -1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1  1
    ## [22297] -1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1
    ## [22321] -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1  1
    ## [22345] -1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1  1 -1  1  1  1 -1 -1
    ## [22369] -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1  1
    ## [22393] -1 -1  1 -1  1  1 -1  1  1 -1  1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1 -1
    ## [22417]  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [22441] -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1
    ## [22465] -1  1  1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [22489]  1  1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [22513] -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1  1  1  1 -1  1 -1  1 -1  1 -1
    ## [22537] -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1  1  1  1  1 -1 -1
    ## [22561] -1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1
    ## [22585] -1 -1  1 -1  1  1 -1  1  1  1 -1 -1  1 -1  1  1 -1  1  1  1 -1  1 -1  1
    ## [22609] -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [22633] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1
    ## [22657] -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1
    ## [22681]  1  1  1  1 -1  1  1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1  1  1 -1 -1 -1
    ## [22705]  1  1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1  1 -1  1 -1  1 -1  1 -1
    ## [22729] -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [22753]  1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [22777] -1  1 -1 -1 -1 -1  1  1  1  1 -1  1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1
    ## [22801] -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [22825]  1 -1  1  1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1
    ## [22849] -1  1 -1 -1  1 -1 -1  1 -1  1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1
    ## [22873] -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1  1 -1  1 -1  1  1  1 -1 -1  1
    ## [22897] -1  1 -1  1  1 -1 -1  1  1  1 -1  1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1
    ## [22921]  1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [22945]  1 -1  1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1
    ## [22969] -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1  1 -1 -1  1
    ## [22993] -1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [23017] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1
    ## [23041]  1 -1  1  1 -1  1 -1  1  1  1  1 -1  1 -1 -1  1  1  1 -1  1 -1  1  1  1
    ## [23065]  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [23089]  1  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1 -1
    ## [23113] -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1
    ## [23137] -1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1  1
    ## [23161] -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1
    ## [23185]  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1  1
    ## [23209]  1  1 -1 -1  1 -1  1  1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1
    ## [23233] -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1
    ## [23257]  1 -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1
    ## [23281] -1  1  1  1  1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1
    ## [23305] -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [23329] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1
    ## [23353] -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1
    ## [23377]  1 -1  1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1
    ## [23401] -1  1  1  1 -1  1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1
    ## [23425] -1  1 -1  1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1  1  1 -1 -1  1 -1 -1
    ## [23449]  1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1
    ## [23473] -1 -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1
    ## [23497]  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1  1  1  1 -1 -1 -1  1 -1 -1
    ## [23521] -1 -1  1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1
    ## [23545]  1 -1  1  1  1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1 -1
    ## [23569] -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1
    ## [23593]  1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1  1
    ## [23617] -1 -1 -1  1  1  1 -1  1  1  1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1  1
    ## [23641] -1 -1  1  1  1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1
    ## [23665] -1  1  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1 -1
    ## [23689]  1 -1  1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1
    ## [23713] -1 -1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1  1 -1 -1  1  1  1 -1  1
    ## [23737] -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1  1  1
    ## [23761] -1 -1 -1  1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1  1  1  1 -1  1
    ## [23785] -1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1  1 -1  1
    ## [23809]  1  1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1  1
    ## [23833]  1 -1  1  1 -1  1 -1  1 -1  1 -1  1  1 -1  1 -1  1 -1  1  1  1 -1  1  1
    ## [23857]  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1
    ## [23881] -1 -1 -1 -1 -1 -1  1  1 -1  1  1  1  1  1 -1 -1  1 -1 -1  1  1  1  1 -1
    ## [23905]  1  1 -1 -1 -1  1 -1  1  1  1 -1  1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1
    ## [23929]  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1
    ## [23953]  1  1 -1  1  1 -1  1  1 -1 -1  1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1
    ## [23977]  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1
    ## [24001]  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [24025] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [24049] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1
    ## [24073]  1  1  1  1 -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1  1  1 -1
    ## [24097]  1 -1 -1  1  1  1  1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [24121] -1  1 -1  1  1 -1  1 -1 -1  1  1  1  1 -1  1  1 -1  1  1 -1 -1 -1  1 -1
    ## [24145] -1  1  1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [24169] -1  1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1  1
    ## [24193] -1 -1 -1  1  1  1  1  1  1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1
    ## [24217]  1  1 -1 -1  1  1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1
    ## [24241] -1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1
    ## [24265]  1 -1 -1 -1 -1  1 -1  1  1  1  1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1
    ## [24289]  1  1  1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1  1  1  1  1  1 -1 -1 -1  1
    ## [24313]  1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1
    ## [24337] -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1  1
    ## [24361]  1  1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1  1  1 -1  1  1  1  1 -1
    ## [24385] -1  1 -1  1  1 -1 -1  1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1
    ## [24409] -1 -1 -1 -1  1  1  1  1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [24433] -1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1 -1 -1 -1  1  1
    ## [24457] -1  1  1  1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1  1  1 -1
    ## [24481]  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1
    ## [24505]  1 -1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1
    ## [24529] -1  1 -1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1  1 -1  1 -1
    ## [24553] -1  1 -1  1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1
    ## [24577]  1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1
    ## [24601]  1  1  1 -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1  1
    ## [24625] -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1 -1  1  1  1
    ## [24649] -1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1
    ## [24673]  1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1 -1 -1 -1  1
    ## [24697] -1  1  1 -1  1  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [24721] -1 -1 -1  1  1  1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1
    ## [24745] -1  1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1  1 -1 -1  1  1  1 -1  1  1
    ## [24769]  1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1
    ## [24793] -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1
    ## [24817] -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1  1  1  1 -1 -1 -1  1
    ## [24841] -1 -1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ## [24865] -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1  1
    ## [24889]  1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1  1  1 -1  1 -1  1 -1  1  1  1
    ## [24913]  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [24937] -1  1 -1 -1  1  1  1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1
    ## [24961] -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1  1 -1  1 -1  1
    ## [24985]  1  1  1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1
    ## [25009]  1 -1 -1  1  1 -1 -1  1  1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1  1
    ## [25033] -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1 -1
    ## [25057] -1  1 -1 -1  1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1
    ## [25081] -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1
    ## [25105] -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1  1  1 -1  1 -1  1  1  1 -1 -1
    ## [25129] -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1  1 -1 -1  1 -1
    ## [25153] -1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1
    ## [25177] -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1
    ## [25201] -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1
    ## [25225] -1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1
    ## [25249] -1 -1 -1  1  1 -1  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1
    ## [25273] -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1  1
    ## [25297] -1  1  1  1 -1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1
    ## [25321] -1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1
    ## [25345] -1  1  1  1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1
    ## [25369] -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1
    ## [25393] -1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1  1
    ## [25417] -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1  1 -1  1 -1 -1  1 -1
    ## [25441]  1 -1 -1  1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [25465] -1 -1  1  1  1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1  1  1 -1  1  1
    ## [25489]  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1 -1  1
    ## [25513]  1 -1  1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1
    ## [25537] -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1
    ## [25561]  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1  1 -1 -1 -1 -1
    ## [25585] -1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1  1 -1  1  1  1 -1  1 -1 -1  1  1
    ## [25609] -1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1
    ## [25633] -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1  1  1  1 -1  1 -1  1  1 -1 -1  1  1
    ## [25657]  1  1 -1  1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1
    ## [25681] -1  1 -1  1  1  1 -1 -1  1  1 -1 -1  1  1  1 -1  1  1  1  1  1 -1  1  1
    ## [25705] -1  1  1  1  1  1 -1  1 -1  1  1  1 -1  1  1 -1  1 -1  1 -1 -1  1 -1  1
    ## [25729] -1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ## [25753] -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1  1  1  1 -1 -1  1 -1
    ## [25777]  1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1 -1
    ## [25801]  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1  1
    ## [25825] -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1  1  1 -1  1
    ## [25849] -1  1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1  1  1
    ## [25873] -1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1
    ## [25897] -1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1
    ## [25921] -1  1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1
    ## [25945]  1 -1  1  1  1 -1 -1 -1 -1  1  1 -1  1 -1  1  1  1 -1  1 -1  1  1  1 -1
    ## [25969]  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1
    ## [25993] -1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1  1  1  1  1 -1 -1  1 -1  1 -1  1  1
    ## [26017] -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1  1  1  1 -1  1 -1 -1  1 -1 -1
    ## [26041] -1  1 -1  1  1  1 -1  1 -1  1  1  1  1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1
    ## [26065] -1  1  1  1 -1 -1  1  1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [26089]  1  1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1 -1  1  1  1 -1
    ## [26113] -1 -1  1  1  1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1
    ## [26137]  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [26161] -1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1 -1 -1
    ## [26185] -1  1 -1  1  1 -1 -1 -1 -1  1  1 -1  1  1 -1  1 -1  1 -1 -1  1 -1  1 -1
    ## [26209] -1  1  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1
    ## [26233] -1  1  1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1
    ## [26257] -1  1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1  1  1 -1  1  1
    ## [26281] -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1  1 -1
    ## [26305]  1 -1 -1  1 -1 -1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1
    ## [26329] -1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1
    ## [26353]  1 -1 -1 -1  1 -1  1  1 -1  1 -1  1  1  1 -1 -1  1 -1  1  1 -1  1 -1 -1
    ## [26377]  1 -1  1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1  1 -1  1  1  1  1  1 -1 -1
    ## [26401]  1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1 -1 -1
    ## [26425]  1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1
    ## [26449]  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ## [26473] -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1
    ## [26497]  1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1
    ## [26521] -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [26545]  1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1  1  1
    ## [26569] -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1
    ## [26593] -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1
    ## [26617] -1  1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1
    ## [26641]  1 -1  1  1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1  1  1  1  1  1 -1
    ## [26665] -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1  1  1  1 -1  1  1  1 -1  1  1  1
    ## [26689] -1  1  1  1  1 -1 -1  1  1  1  1 -1  1  1  1  1 -1 -1  1  1 -1 -1  1 -1
    ## [26713] -1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1  1  1 -1 -1  1 -1 -1 -1
    ## [26737]  1 -1  1 -1  1  1  1 -1 -1  1  1  1  1  1  1 -1 -1  1  1  1 -1  1 -1 -1
    ## [26761]  1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1
    ## [26785] -1  1 -1  1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1
    ## [26809] -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1  1
    ## [26833] -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1  1  1 -1
    ## [26857] -1 -1 -1 -1  1  1 -1  1  1 -1  1 -1 -1  1  1  1  1 -1 -1 -1 -1  1  1  1
    ## [26881]  1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1
    ## [26905] -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1
    ## [26929] -1 -1  1  1  1  1 -1 -1 -1  1  1  1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1
    ## [26953] -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1  1 -1 -1 -1
    ## [26977]  1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1
    ## [27001] -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1  1 -1 -1
    ## [27025]  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1  1  1  1  1  1 -1
    ## [27049] -1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1  1  1  1 -1
    ## [27073]  1 -1  1  1  1 -1  1  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1
    ## [27097]  1  1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1
    ## [27121] -1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1  1 -1
    ## [27145] -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1  1  1
    ## [27169]  1  1 -1  1  1  1  1 -1  1 -1  1  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1
    ## [27193] -1  1  1 -1  1  1  1  1  1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1
    ## [27217] -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1  1  1 -1 -1  1 -1  1 -1  1 -1 -1  1
    ## [27241] -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1
    ## [27265] -1  1 -1  1 -1 -1  1  1  1  1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1
    ## [27289]  1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1
    ## [27313] -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1 -1  1 -1 -1  1
    ## [27337] -1  1  1  1 -1 -1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1
    ## [27361]  1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1
    ## [27385] -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1  1  1  1  1
    ## [27409] -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1
    ## [27433]  1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1  1 -1 -1
    ## [27457] -1 -1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [27481] -1  1  1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1  1 -1  1  1  1  1
    ## [27505] -1  1  1 -1 -1  1  1 -1  1 -1  1  1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1
    ## [27529] -1 -1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1 -1  1 -1  1  1 -1 -1 -1  1 -1
    ## [27553]  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1
    ## [27577] -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1  1 -1  1  1  1  1
    ## [27601]  1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1
    ## [27625]  1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1  1 -1 -1  1  1
    ## [27649] -1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1  1  1 -1  1  1  1 -1  1 -1  1
    ## [27673] -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1  1  1  1
    ## [27697] -1 -1  1  1  1 -1  1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1
    ## [27721]  1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1
    ## [27745]  1 -1 -1  1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1
    ## [27769]  1 -1  1 -1 -1  1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ## [27793]  1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [27817]  1  1 -1  1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1  1
    ## [27841] -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1
    ## [27865]  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1
    ## [27889]  1  1  1  1  1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1 -1
    ## [27913]  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1
    ## [27937]  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1  1 -1  1  1
    ## [27961]  1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1  1  1  1  1 -1 -1 -1  1 -1  1
    ## [27985] -1  1 -1  1  1  1 -1  1  1  1  1  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1 -1
    ## [28009]  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1 -1
    ## [28033] -1  1  1 -1  1  1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1
    ## [28057]  1  1 -1  1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1  1  1
    ## [28081] -1  1 -1 -1  1 -1 -1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1 -1  1 -1  1 -1
    ## [28105]  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1  1  1  1 -1  1  1 -1  1 -1 -1 -1
    ## [28129]  1  1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1
    ## [28153] -1  1  1  1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1
    ## [28177]  1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1  1  1  1
    ## [28201]  1  1 -1  1  1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1
    ## [28225] -1 -1  1  1 -1  1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1 -1 -1
    ## [28249]  1 -1 -1  1  1  1 -1 -1 -1  1  1  1  1  1  1 -1  1  1 -1  1 -1  1  1  1
    ## [28273] -1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1 -1  1  1  1 -1 -1  1 -1  1 -1  1
    ## [28297] -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1
    ## [28321] -1  1  1  1  1 -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1  1 -1  1 -1 -1  1  1
    ## [28345] -1 -1  1 -1 -1  1  1  1 -1  1 -1  1  1  1  1  1  1 -1 -1 -1 -1 -1  1 -1
    ## [28369] -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1 -1 -1
    ## [28393] -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1  1  1  1  1  1 -1 -1  1  1  1  1 -1
    ## [28417] -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1  1 -1 -1
    ## [28441]  1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1  1 -1  1 -1  1 -1 -1
    ## [28465] -1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1
    ## [28489]  1 -1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1
    ## [28513]  1  1  1 -1  1 -1 -1  1  1  1  1  1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1
    ## [28537]  1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1  1  1 -1  1  1  1  1 -1
    ## [28561] -1  1  1  1 -1  1  1 -1 -1  1 -1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1
    ## [28585] -1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1
    ## [28609] -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1  1 -1  1  1 -1 -1 -1  1 -1
    ## [28633]  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1
    ## [28657] -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1 -1  1  1
    ## [28681]  1  1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1  1  1  1  1  1 -1 -1 -1  1  1 -1
    ## [28705] -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1
    ## [28729] -1 -1 -1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1  1  1  1  1  1  1  1  1  1
    ## [28753]  1  1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1
    ## [28777]  1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1
    ## [28801]  1 -1  1  1  1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1
    ## [28825] -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1  1  1  1  1  1
    ## [28849] -1  1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1
    ## [28873]  1  1  1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1 -1 -1 -1
    ## [28897]  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1
    ## [28921] -1  1  1  1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1 -1
    ## [28945] -1 -1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1
    ## [28969]  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1
    ## [28993] -1  1  1  1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1  1
    ## [29017]  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1
    ## [29041]  1 -1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1
    ## [29065]  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1 -1  1 -1  1 -1  1 -1  1
    ## [29089] -1 -1  1  1 -1  1  1 -1  1 -1 -1 -1  1  1 -1  1  1 -1  1 -1 -1  1 -1 -1
    ## [29113] -1  1 -1  1 -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1
    ## [29137]  1  1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1  1
    ## [29161] -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1
    ## [29185] -1  1 -1  1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1  1
    ## [29209] -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1
    ## [29233] -1  1 -1  1  1  1 -1  1 -1  1  1  1  1 -1 -1 -1  1 -1 -1  1  1  1  1 -1
    ## [29257] -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1
    ## [29281]  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1
    ## [29305] -1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1
    ## [29329] -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1
    ## [29353] -1 -1 -1  1  1  1  1  1  1  1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [29377]  1 -1  1  1  1  1 -1  1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1
    ## [29401] -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1  1  1  1
    ## [29425]  1  1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1  1
    ## [29449] -1  1 -1  1  1 -1  1  1  1 -1 -1  1 -1  1  1  1 -1 -1 -1  1 -1  1 -1 -1
    ## [29473] -1  1 -1  1 -1 -1 -1 -1  1 -1  1  1  1  1  1  1  1  1  1 -1  1 -1 -1 -1
    ## [29497]  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1  1  1 -1 -1
    ## [29521] -1 -1  1  1 -1 -1  1  1 -1  1  1  1  1  1  1  1 -1 -1  1  1 -1  1 -1 -1
    ## [29545]  1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1 -1  1 -1 -1 -1  1  1  1
    ## [29569]  1  1 -1 -1 -1  1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1
    ## [29593] -1  1  1  1  1 -1  1  1 -1  1 -1  1  1  1  1 -1 -1  1  1 -1  1 -1 -1  1
    ## [29617]  1  1  1  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1
    ## [29641]  1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1  1  1  1  1 -1 -1  1 -1  1 -1  1
    ## [29665]  1  1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1
    ## [29689] -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [29713] -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1  1 -1  1 -1
    ## [29737] -1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1
    ## [29761]  1  1 -1 -1  1  1  1 -1  1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1
    ## [29785] -1 -1 -1  1 -1  1 -1  1 -1  1  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1
    ## [29809]  1  1  1 -1  1  1  1  1 -1  1  1 -1  1  1 -1  1  1  1 -1  1 -1  1  1  1
    ## [29833]  1 -1 -1 -1 -1  1  1  1 -1  1  1 -1  1 -1 -1 -1  1  1  1  1 -1  1 -1 -1
    ## [29857]  1  1 -1  1  1  1 -1  1  1  1  1 -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1  1
    ## [29881]  1 -1  1  1 -1  1 -1  1  1 -1  1  1  1  1 -1  1 -1  1  1 -1  1 -1 -1 -1
    ## [29905]  1  1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [29929]  1 -1  1  1 -1  1 -1  1  1 -1  1  1  1 -1 -1  1  1  1  1 -1 -1  1  1 -1
    ## [29953] -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1  1  1 -1  1 -1
    ## [29977]  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1  1 -1
    ## [30001]  1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1  1  1  1  1 -1
    ## [30025] -1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1  1
    ## [30049] -1  1  1  1 -1  1  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1  1  1  1
    ## [30073] -1  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1
    ## [30097] -1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1
    ## [30121] -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1
    ## [30145]  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1
    ## [30169]  1  1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1  1 -1
    ## [30193] -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1  1  1 -1 -1 -1  1 -1 -1
    ## [30217]  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1
    ## [30241] -1  1 -1  1  1 -1  1 -1 -1  1  1 -1  1  1 -1  1 -1  1  1  1  1 -1 -1 -1
    ## [30265] -1  1  1  1  1  1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1
    ## [30289]  1 -1 -1 -1  1  1  1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1
    ## [30313] -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1
    ## [30337] -1  1  1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1
    ## [30361] -1 -1 -1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1  1
    ## [30385] -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1
    ## [30409] -1  1 -1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1  1  1  1  1  1
    ## [30433]  1 -1 -1  1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1
    ## [30457]  1 -1  1  1  1  1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1  1  1  1
    ## [30481] -1 -1  1  1 -1  1 -1  1  1 -1  1 -1 -1  1  1 -1  1  1 -1  1 -1  1  1 -1
    ## [30505]  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [30529] -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [30553]  1  1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1  1  1  1 -1  1 -1 -1  1  1
    ## [30577] -1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1  1 -1  1  1 -1  1 -1  1  1 -1 -1  1
    ## [30601] -1 -1 -1  1  1 -1 -1  1 -1 -1  1  1  1 -1  1  1  1  1 -1  1  1 -1  1 -1
    ## [30625] -1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1  1  1  1  1 -1  1  1
    ## [30649] -1  1 -1  1 -1  1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1
    ## [30673] -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1  1
    ## [30697] -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1  1  1 -1  1 -1  1
    ## [30721]  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1
    ## [30745] -1 -1  1  1  1  1 -1  1  1  1 -1  1  1 -1  1  1 -1  1 -1 -1  1  1  1  1
    ## [30769]  1  1 -1  1  1  1  1 -1 -1  1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1
    ## [30793] -1 -1 -1 -1  1  1  1 -1  1  1 -1  1 -1 -1  1  1  1  1 -1  1 -1 -1 -1 -1
    ## [30817] -1  1 -1  1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1  1
    ## [30841] -1  1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1  1  1  1 -1 -1  1  1 -1  1 -1
    ## [30865]  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1
    ## [30889] -1 -1  1  1  1 -1 -1  1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1
    ## [30913] -1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1  1  1  1  1 -1 -1  1 -1 -1 -1  1
    ## [30937] -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1  1  1 -1 -1  1  1 -1
    ## [30961]  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1 -1
    ## [30985]  1  1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1  1  1  1 -1 -1  1  1 -1 -1  1
    ## [31009] -1 -1  1  1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1  1 -1  1 -1 -1  1 -1
    ## [31033] -1  1  1  1  1  1  1  1 -1  1  1 -1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1
    ## [31057]  1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1
    ## [31081] -1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1  1  1  1  1 -1 -1 -1  1 -1 -1 -1
    ## [31105] -1  1  1 -1  1  1 -1  1 -1  1  1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1  1
    ## [31129] -1  1  1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1  1  1  1  1  1  1  1
    ## [31153] -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1
    ## [31177]  1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1  1  1 -1  1  1 -1  1 -1  1  1 -1
    ## [31201]  1 -1 -1 -1  1  1 -1  1 -1  1  1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1  1  1
    ## [31225] -1  1 -1  1  1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1
    ## [31249] -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1  1  1  1 -1
    ## [31273] -1  1 -1  1 -1  1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1
    ## [31297] -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1  1 -1  1  1 -1 -1 -1  1 -1
    ## [31321]  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1
    ## [31345] -1 -1  1 -1  1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1
    ## [31369] -1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1
    ## [31393] -1  1  1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1
    ## [31417]  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1  1  1  1 -1
    ## [31441] -1 -1  1  1  1  1  1 -1  1 -1  1  1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1
    ## [31465] -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1
    ## [31489]  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1
    ## [31513] -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1
    ## [31537] -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1
    ## [31561] -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1  1  1
    ## [31585]  1  1  1  1 -1 -1  1  1 -1  1  1  1  1 -1  1 -1  1  1  1  1  1  1 -1 -1
    ## [31609]  1 -1  1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1
    ## [31633]  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1
    ## [31657]  1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1  1  1  1 -1 -1  1  1  1 -1 -1 -1 -1
    ## [31681] -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1  1 -1
    ## [31705]  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1 -1 -1 -1  1  1 -1
    ## [31729]  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1  1  1 -1 -1
    ## [31753] -1 -1 -1  1 -1 -1  1  1 -1  1  1 -1 -1  1  1 -1 -1 -1  1  1  1  1 -1 -1
    ## [31777] -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1  1  1  1  1  1  1  1
    ## [31801]  1  1  1  1  1  1  1 -1  1  1 -1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1 -1
    ## [31825]  1 -1 -1 -1  1 -1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [31849] -1  1  1  1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1
    ## [31873] -1  1 -1  1  1 -1 -1  1  1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1
    ## [31897] -1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1
    ## [31921]  1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [31945]  1 -1  1 -1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1  1  1  1 -1  1
    ## [31969]  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1
    ## [31993] -1  1  1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1
    ## [32017] -1 -1  1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1
    ## [32041] -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1
    ## [32065]  1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1
    ## [32089] -1  1  1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1
    ## [32113]  1  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1
    ## [32137]  1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1
    ## [32161] -1 -1  1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1
    ## [32185]  1  1 -1  1 -1  1  1  1  1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1 -1
    ## [32209]  1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1  1  1 -1  1  1 -1  1  1 -1
    ## [32233] -1 -1 -1  1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1
    ## [32257]  1  1  1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1  1  1  1
    ## [32281] -1  1 -1 -1  1 -1 -1  1  1  1  1 -1  1 -1  1  1 -1  1  1  1 -1  1  1 -1
    ## [32305]  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1
    ## [32329] -1  1  1 -1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1 -1  1  1 -1  1 -1 -1  1
    ## [32353] -1 -1  1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1
    ## [32377]  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1
    ## [32401] -1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1
    ## [32425]  1  1 -1 -1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1  1 -1 -1  1 -1 -1
    ## [32449]  1  1  1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1
    ## [32473]  1  1 -1  1  1 -1 -1  1  1  1 -1  1 -1  1  1  1  1  1  1 -1  1 -1 -1 -1
    ## [32497] -1  1  1 -1 -1 -1  1  1  1  1 -1 -1  1  1  1 -1  1 -1 -1 -1  1  1  1  1
    ## [32521]  1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1  1  1
    ## [32545] -1  1  1 -1 -1  1  1 -1  1  1  1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1  1
    ## [32569] -1  1  1  1 -1  1  1  1 -1  1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1 -1  1
    ## [32593]  1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1  1  1 -1 -1
    ## [32617]  1  1 -1  1 -1 -1 -1  1  1  1  1 -1  1  1  1  1  1 -1  1  1  1  1 -1  1
    ## [32641]  1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1
    ## [32665] -1 -1  1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1
    ## [32689] -1 -1  1  1  1  1  1  1 -1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1  1  1
    ## [32713] -1 -1  1 -1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1
    ## [32737] -1  1  1  1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1
    ## [32761]  1 -1  1  1  1 -1  1  1  1 -1 -1  1 -1 -1  1 -1  1 -1  1  1  1  1  1 -1
    ## [32785] -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1 -1  1 -1
    ## [32809] -1  1 -1  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1
    ## [32833] -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1  1 -1
    ## [32857] -1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1
    ## [32881]  1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1  1
    ## [32905] -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1
    ## [32929]  1  1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1  1 -1  1 -1
    ## [32953]  1  1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1 -1  1  1  1  1
    ## [32977]  1  1 -1 -1  1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1
    ## [33001] -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1  1  1  1  1
    ## [33025] -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1
    ## [33049]  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1
    ## [33073]  1  1  1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1
    ## [33097]  1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1
    ## [33121] -1 -1 -1 -1  1  1  1  1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1
    ## [33145] -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1  1 -1  1  1  1 -1  1  1  1 -1 -1  1
    ## [33169]  1  1  1  1  1 -1 -1  1  1 -1 -1 -1  1 -1  1  1  1  1  1  1  1  1  1  1
    ## [33193]  1 -1  1  1  1 -1  1  1  1 -1  1  1  1 -1  1  1 -1  1 -1  1  1 -1  1 -1
    ## [33217]  1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1 -1  1  1  1
    ## [33241] -1 -1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1  1
    ## [33265] -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1
    ## [33289] -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1
    ## [33313]  1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1
    ## [33337] -1  1  1  1 -1  1 -1  1 -1  1 -1 -1  1  1 -1  1  1 -1  1 -1 -1 -1  1  1
    ## [33361] -1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1  1  1  1 -1
    ## [33385] -1 -1 -1 -1  1  1  1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1  1 -1  1  1 -1
    ## [33409] -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [33433]  1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1
    ## [33457]  1  1 -1  1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1  1
    ## [33481] -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1
    ## [33505] -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1 -1  1  1 -1
    ## [33529]  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1
    ## [33553]  1 -1 -1  1  1 -1  1  1 -1 -1  1 -1 -1  1  1  1 -1 -1  1 -1  1  1 -1 -1
    ## [33577]  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1 -1 -1  1
    ## [33601] -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1
    ## [33625] -1  1 -1 -1  1 -1  1  1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1  1 -1  1 -1  1
    ## [33649] -1 -1  1 -1 -1  1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1  1  1 -1  1  1  1
    ## [33673]  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1
    ## [33697]  1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1
    ## [33721]  1  1 -1  1  1  1 -1 -1  1 -1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1  1  1
    ## [33745] -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1  1  1
    ## [33769]  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1
    ## [33793] -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1 -1  1 -1  1  1 -1  1 -1 -1  1  1  1
    ## [33817]  1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1 -1  1  1  1  1 -1
    ## [33841]  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1
    ## [33865]  1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1  1
    ## [33889]  1 -1 -1  1 -1  1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1
    ## [33913]  1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1
    ## [33937]  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1
    ## [33961] -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1
    ## [33985] -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1
    ## [34009] -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1
    ## [34033] -1  1  1  1 -1  1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1  1  1  1  1 -1 -1  1
    ## [34057] -1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1  1  1  1  1 -1  1  1  1 -1  1
    ## [34081] -1  1 -1 -1 -1 -1  1  1  1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1
    ## [34105]  1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1  1
    ## [34129]  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1
    ## [34153] -1  1  1  1 -1  1 -1  1 -1 -1 -1  1  1 -1  1  1  1 -1 -1 -1  1 -1  1  1
    ## [34177] -1  1  1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1
    ## [34201]  1  1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1
    ## [34225]  1 -1  1  1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1
    ## [34249]  1 -1  1  1 -1  1 -1  1 -1 -1 -1  1  1  1  1  1 -1 -1  1  1 -1  1  1 -1
    ## [34273]  1  1 -1 -1  1  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [34297] -1 -1 -1 -1  1  1 -1  1  1 -1  1  1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ## [34321] -1 -1  1 -1  1  1  1 -1 -1 -1  1 -1  1  1  1 -1 -1  1  1 -1 -1  1 -1 -1
    ## [34345] -1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1
    ## [34369]  1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1
    ## [34393]  1 -1 -1 -1  1 -1  1  1  1  1  1 -1  1  1 -1 -1 -1  1  1 -1  1  1 -1 -1
    ## [34417]  1  1 -1 -1 -1  1  1  1 -1  1  1  1  1  1 -1  1  1  1  1 -1  1 -1  1  1
    ## [34441] -1  1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1  1  1  1  1 -1 -1 -1  1 -1  1  1
    ## [34465]  1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1 -1  1  1
    ## [34489]  1 -1  1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1 -1 -1
    ## [34513] -1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1  1  1 -1  1 -1  1  1  1
    ## [34537] -1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1 -1  1 -1
    ## [34561] -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1
    ## [34585] -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1
    ## [34609]  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1
    ## [34633]  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ## [34657] -1 -1 -1 -1  1  1  1  1 -1  1 -1  1  1 -1 -1 -1  1  1 -1  1 -1  1  1  1
    ## [34681]  1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1
    ## [34705]  1  1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1  1  1  1
    ## [34729] -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1  1  1 -1 -1  1 -1 -1 -1 -1
    ## [34753]  1  1 -1  1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1
    ## [34777] -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1
    ## [34801] -1  1  1  1 -1  1 -1  1 -1  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1
    ## [34825]  1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1  1
    ## [34849] -1  1  1  1  1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1
    ## [34873]  1  1  1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1 -1  1
    ## [34897] -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1
    ## [34921] -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1
    ## [34945]  1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1  1
    ## [34969] -1 -1 -1  1  1 -1  1  1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1
    ## [34993] -1  1 -1  1  1  1  1 -1  1  1 -1  1  1  1  1 -1  1  1  1  1 -1  1 -1 -1
    ## [35017]  1  1  1 -1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1
    ## [35041] -1 -1 -1  1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1  1 -1 -1  1 -1  1  1 -1
    ## [35065] -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1
    ## [35089] -1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [35113] -1  1  1  1 -1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1
    ## [35137] -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1  1  1  1 -1 -1 -1 -1 -1  1 -1
    ## [35161] -1  1  1  1 -1  1  1 -1  1 -1 -1  1 -1  1  1 -1  1 -1  1 -1 -1  1  1  1
    ## [35185]  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [35209]  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1
    ## [35233] -1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1
    ## [35257] -1 -1 -1 -1  1 -1  1  1 -1  1  1  1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1
    ## [35281] -1 -1 -1  1 -1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1
    ## [35305]  1 -1  1  1 -1  1 -1  1  1  1 -1  1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1
    ## [35329]  1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1
    ## [35353] -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1  1  1 -1
    ## [35377]  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1 -1  1 -1  1  1
    ## [35401] -1 -1  1 -1  1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1
    ## [35425] -1  1 -1 -1  1  1  1  1 -1 -1  1  1  1 -1  1  1 -1  1  1 -1 -1  1  1 -1
    ## [35449]  1  1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1  1  1 -1
    ## [35473]  1  1  1  1  1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1 -1 -1  1  1
    ## [35497] -1  1 -1 -1  1 -1  1 -1 -1 -1  1  1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1
    ## [35521] -1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1
    ## [35545]  1  1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1  1  1 -1 -1  1 -1 -1 -1
    ## [35569] -1  1  1 -1  1  1 -1  1  1 -1  1  1 -1  1  1 -1 -1  1 -1 -1 -1  1  1  1
    ## [35593]  1 -1 -1  1 -1  1 -1 -1  1  1  1  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1
    ## [35617] -1  1 -1  1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1  1  1  1  1  1  1 -1 -1
    ## [35641] -1  1 -1 -1  1  1 -1  1 -1 -1  1  1  1  1  1 -1 -1  1  1 -1  1  1  1 -1
    ## [35665] -1 -1  1  1  1  1 -1  1  1  1 -1  1 -1  1 -1 -1 -1  1  1 -1  1 -1  1 -1
    ## [35689] -1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1  1
    ## [35713] -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1  1  1
    ## [35737] -1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1  1
    ## [35761] -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1
    ## [35785] -1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1
    ## [35809] -1  1  1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1
    ## [35833]  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1  1  1 -1  1  1 -1  1 -1  1
    ## [35857] -1 -1 -1  1  1 -1 -1  1  1  1 -1  1 -1 -1  1 -1  1 -1  1  1 -1 -1  1  1
    ## [35881] -1 -1  1 -1  1 -1 -1  1  1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1
    ## [35905]  1  1 -1  1  1  1 -1 -1  1  1  1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1  1
    ## [35929] -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1
    ## [35953] -1  1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1  1  1  1  1 -1 -1  1 -1 -1
    ## [35977]  1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1
    ## [36001] -1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1  1  1 -1 -1  1 -1 -1  1  1  1  1
    ## [36025] -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1  1  1  1
    ## [36049]  1  1 -1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1
    ## [36073] -1  1  1 -1 -1  1 -1 -1  1  1  1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1
    ## [36097] -1 -1 -1  1  1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1
    ## [36121] -1  1  1 -1  1  1  1  1 -1  1  1 -1 -1  1  1  1 -1 -1  1 -1  1 -1 -1  1
    ## [36145] -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1  1
    ## [36169] -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1  1  1 -1 -1
    ## [36193]  1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1
    ## [36217]  1  1  1 -1  1  1  1 -1 -1 -1 -1  1 -1 -1  1  1  1  1 -1  1  1  1  1 -1
    ## [36241]  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1
    ## [36265] -1  1  1 -1  1 -1  1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1  1  1  1  1  1
    ## [36289] -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1
    ## [36313]  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1  1 -1  1 -1  1  1  1  1 -1
    ## [36337]  1 -1  1  1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1
    ## [36361]  1  1 -1 -1  1 -1  1  1  1 -1  1 -1  1 -1  1 -1 -1  1  1 -1  1  1  1 -1
    ## [36385] -1  1  1 -1  1  1  1  1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [36409] -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1
    ## [36433] -1 -1 -1 -1  1  1  1  1 -1 -1  1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1
    ## [36457]  1 -1  1  1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1  1  1  1  1  1  1  1  1
    ## [36481] -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1
    ## [36505]  1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1  1  1
    ## [36529] -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1
    ## [36553]  1  1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1 -1 -1  1 -1
    ## [36577] -1 -1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1  1
    ## [36601]  1  1  1  1  1  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1  1
    ## [36625] -1  1 -1  1 -1  1  1  1  1  1  1  1 -1  1 -1  1 -1  1  1 -1  1 -1 -1  1
    ## [36649]  1 -1  1 -1 -1 -1  1  1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1
    ## [36673]  1 -1  1  1  1 -1 -1  1 -1  1 -1  1 -1  1  1  1  1 -1  1 -1 -1  1 -1 -1
    ## [36697] -1 -1  1  1 -1  1  1 -1 -1  1  1  1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1
    ## [36721] -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1  1  1  1 -1 -1
    ## [36745] -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1 -1  1  1
    ## [36769]  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1  1 -1
    ## [36793] -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1 -1 -1  1  1 -1  1 -1
    ## [36817]  1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1  1  1  1  1  1  1
    ## [36841] -1 -1 -1  1  1  1 -1  1 -1  1 -1  1 -1  1  1  1 -1 -1  1  1 -1  1 -1 -1
    ## [36865] -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1  1  1 -1 -1  1 -1 -1  1
    ## [36889]  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1  1
    ## [36913]  1  1  1  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1
    ## [36937]  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1  1 -1  1  1  1  1 -1  1 -1 -1 -1 -1
    ## [36961] -1  1  1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1
    ## [36985] -1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1 -1  1 -1 -1 -1  1  1
    ## [37009] -1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1  1 -1 -1  1  1  1  1 -1 -1 -1 -1  1
    ## [37033]  1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1
    ## [37057]  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1  1
    ## [37081] -1  1  1  1 -1  1  1  1 -1  1 -1 -1  1 -1  1 -1  1  1 -1 -1  1 -1 -1  1
    ## [37105]  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1
    ## [37129] -1 -1 -1 -1  1 -1  1 -1  1  1 -1  1  1  1  1  1 -1 -1 -1 -1 -1  1  1 -1
    ## [37153] -1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [37177]  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1 -1  1
    ## [37201]  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1  1  1  1 -1  1 -1  1  1  1
    ## [37225] -1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1  1  1  1  1  1  1  1 -1  1 -1 -1 -1
    ## [37249] -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1  1 -1  1  1  1  1 -1 -1  1  1
    ## [37273]  1  1  1  1 -1  1 -1 -1  1  1  1  1  1 -1  1  1  1  1  1  1  1 -1 -1  1
    ## [37297] -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1
    ## [37321] -1  1 -1 -1  1 -1 -1 -1  1 -1  1  1  1  1  1  1  1  1  1  1  1 -1  1 -1
    ## [37345] -1  1  1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1 -1  1  1  1
    ## [37369] -1 -1 -1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1
    ## [37393] -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1 -1  1
    ## [37417]  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1
    ## [37441]  1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1  1
    ## [37465]  1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1  1  1 -1 -1
    ## [37489] -1  1  1 -1 -1  1  1  1 -1 -1  1  1  1  1 -1 -1  1 -1 -1  1  1  1 -1 -1
    ## [37513] -1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1  1
    ## [37537]  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1 -1  1
    ## [37561]  1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1  1  1  1 -1
    ## [37585]  1  1 -1 -1 -1  1  1  1  1  1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1
    ## [37609] -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1
    ## [37633] -1 -1 -1 -1 -1  1  1  1 -1  1 -1  1  1  1  1  1 -1 -1  1 -1 -1  1  1 -1
    ## [37657]  1  1  1 -1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1
    ## [37681] -1  1  1  1 -1  1  1  1  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1
    ## [37705]  1  1  1  1  1  1 -1  1 -1  1  1  1  1 -1  1  1  1  1 -1  1  1 -1 -1  1
    ## [37729]  1 -1  1  1  1 -1 -1  1  1 -1  1 -1  1  1  1  1  1 -1 -1 -1  1  1  1 -1
    ## [37753] -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1
    ## [37777]  1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1  1 -1  1
    ## [37801] -1  1 -1  1 -1 -1  1  1  1 -1 -1  1 -1  1 -1 -1  1  1  1  1 -1 -1 -1  1
    ## [37825] -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1  1  1
    ## [37849] -1 -1  1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [37873]  1  1 -1 -1 -1  1  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1  1
    ## [37897] -1 -1  1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1
    ## [37921]  1 -1 -1  1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1  1  1  1  1
    ## [37945]  1  1 -1  1 -1 -1 -1  1  1  1 -1 -1  1 -1  1  1  1  1  1 -1 -1 -1  1  1
    ## [37969] -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1  1 -1
    ## [37993]  1  1 -1  1  1 -1  1  1 -1  1 -1  1  1 -1 -1  1  1  1  1  1  1  1 -1 -1
    ## [38017]  1 -1  1  1 -1  1  1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1  1  1  1
    ## [38041]  1 -1  1  1 -1  1  1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1  1  1
    ## [38065] -1  1 -1 -1 -1 -1 -1  1  1  1 -1  1  1  1  1 -1  1  1 -1 -1  1  1 -1 -1
    ## [38089]  1  1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1
    ## [38113] -1  1  1 -1  1 -1  1  1  1 -1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1
    ## [38137]  1 -1  1  1 -1 -1 -1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1
    ## [38161]  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1  1 -1 -1  1
    ## [38185]  1 -1  1  1  1 -1  1  1  1 -1  1 -1  1 -1  1  1  1  1  1  1  1  1 -1  1
    ## [38209] -1 -1  1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1  1  1  1  1  1 -1 -1  1  1
    ## [38233] -1  1 -1  1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1  1 -1
    ## [38257]  1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1  1  1  1 -1 -1 -1 -1
    ## [38281] -1 -1 -1  1  1 -1 -1 -1  1  1  1  1  1  1 -1 -1 -1  1 -1 -1  1  1  1 -1
    ## [38305]  1  1  1  1 -1  1 -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1  1  1  1  1 -1  1
    ## [38329] -1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1
    ## [38353]  1 -1  1 -1  1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1
    ## [38377]  1  1  1  1 -1 -1  1  1  1  1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1
    ## [38401] -1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1  1 -1 -1  1 -1
    ## [38425] -1 -1  1  1  1  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1  1  1
    ## [38449] -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1  1  1 -1 -1  1  1  1  1  1 -1 -1
    ## [38473]  1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1
    ## [38497]  1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1
    ## [38521] -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1  1  1  1  1
    ## [38545] -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1  1  1  1 -1  1  1  1  1  1 -1
    ## [38569] -1  1  1 -1 -1  1  1  1  1  1  1  1 -1  1  1  1 -1 -1  1 -1 -1 -1  1 -1
    ## [38593]  1 -1  1  1 -1  1  1  1 -1  1  1  1  1  1  1  1  1  1  1 -1  1  1  1 -1
    ## [38617]  1  1  1  1  1 -1  1  1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1
    ## [38641]  1  1 -1  1 -1 -1 -1  1  1 -1 -1  1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1
    ## [38665]  1  1 -1 -1 -1  1  1  1 -1 -1  1 -1  1  1  1 -1 -1 -1  1 -1  1  1  1 -1
    ## [38689] -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1
    ## [38713] -1  1  1 -1 -1 -1  1  1 -1  1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1
    ## [38737]  1 -1  1 -1  1  1 -1  1 -1  1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1
    ## [38761]  1  1 -1  1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1 -1  1 -1  1
    ## [38785] -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1
    ## [38809] -1 -1 -1 -1 -1 -1  1  1  1  1  1  1  1  1 -1  1  1  1 -1  1  1  1 -1  1
    ## [38833] -1  1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1
    ## [38857] -1 -1  1  1 -1 -1  1  1  1  1 -1 -1 -1  1  1  1 -1  1 -1  1  1 -1  1  1
    ## [38881] -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1  1
    ## [38905]  1  1 -1  1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1
    ## [38929]  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1  1  1  1 -1 -1  1 -1 -1  1
    ## [38953] -1  1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1
    ## [38977] -1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1  1 -1 -1
    ## [39001]  1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1  1 -1  1 -1  1  1  1  1 -1 -1  1
    ## [39025] -1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [39049] -1 -1  1  1  1  1  1 -1  1 -1  1 -1  1 -1  1  1 -1 -1 -1  1  1  1 -1 -1
    ## [39073]  1 -1 -1  1  1  1 -1  1  1 -1  1  1  1  1 -1 -1  1  1 -1  1 -1 -1  1  1
    ## [39097]  1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1  1 -1
    ## [39121]  1 -1  1  1 -1  1 -1 -1  1  1 -1  1  1  1 -1 -1  1 -1  1  1  1  1 -1 -1
    ## [39145]  1  1 -1  1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1  1  1
    ## [39169]  1 -1 -1 -1  1 -1  1  1  1  1 -1 -1  1 -1  1 -1  1 -1  1  1 -1  1 -1  1
    ## [39193] -1 -1 -1  1 -1  1  1 -1  1  1 -1  1  1  1 -1 -1 -1  1  1 -1 -1  1  1  1
    ## [39217] -1  1  1  1  1 -1  1 -1  1  1  1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1 -1
    ## [39241] -1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1
    ## [39265]  1 -1  1  1  1 -1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [39289] -1  1  1  1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1  1
    ## [39313]  1 -1  1 -1 -1  1  1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1
    ## [39337] -1 -1 -1 -1 -1  1  1 -1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1  1
    ## [39361]  1 -1  1 -1  1 -1  1  1  1 -1  1 -1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1
    ## [39385] -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1
    ## [39409]  1  1 -1 -1  1  1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1 -1  1 -1  1
    ## [39433] -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1
    ## [39457] -1  1 -1  1  1 -1 -1  1  1  1  1  1  1 -1  1  1 -1 -1 -1  1  1  1  1  1
    ## [39481]  1  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1 -1  1  1  1 -1  1 -1 -1  1
    ## [39505] -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1  1 -1 -1  1
    ## [39529] -1 -1  1 -1  1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1 -1 -1  1  1  1 -1 -1
    ## [39553] -1  1  1  1 -1 -1  1  1  1  1  1  1  1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1
    ## [39577]  1 -1 -1 -1  1  1 -1 -1  1  1  1 -1  1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1
    ## [39601] -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1 -1  1  1  1  1 -1  1  1  1 -1 -1  1
    ## [39625]  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1  1
    ## [39649] -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1  1
    ## [39673]  1 -1 -1  1  1 -1 -1 -1  1  1 -1  1 -1  1  1 -1  1 -1 -1 -1  1  1 -1 -1
    ## [39697]  1  1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1
    ## [39721]  1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1
    ## [39745] -1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1 -1 -1 -1
    ## [39769]  1 -1  1  1 -1  1  1  1 -1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1
    ## [39793]  1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1  1  1  1 -1 -1  1  1  1 -1 -1 -1 -1
    ## [39817]  1  1  1  1  1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1
    ## [39841] -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1
    ## [39865] -1  1  1 -1 -1  1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1  1
    ## [39889] -1 -1  1 -1  1  1  1  1  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1
    ## [39913]  1  1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1  1  1  1  1  1 -1  1
    ## [39937] -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1
    ## [39961]  1  1 -1 -1  1  1 -1  1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1
    ## [39985]  1  1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1 -1 -1
    ## [40009]  1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1  1
    ## [40033] -1  1  1  1  1 -1  1  1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1
    ## [40057] -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1
    ## [40081] -1 -1 -1  1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1
    ## [40105]  1  1  1  1  1 -1  1  1  1 -1 -1 -1  1 -1  1 -1  1  1 -1  1  1 -1 -1  1
    ## [40129] -1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1
    ## [40153]  1 -1  1 -1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1
    ## [40177] -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1  1 -1
    ## [40201] -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1
    ## [40225]  1 -1 -1  1 -1  1 -1  1  1  1  1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1
    ## [40249] -1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1  1  1  1 -1 -1  1 -1 -1 -1 -1
    ## [40273] -1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1 -1  1  1  1  1 -1  1  1  1 -1  1
    ## [40297]  1  1  1 -1  1  1  1 -1 -1  1 -1  1  1  1  1  1 -1 -1  1 -1 -1 -1 -1  1
    ## [40321]  1  1 -1  1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1
    ## [40345] -1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1
    ## [40369]  1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1  1  1
    ## [40393]  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1  1 -1 -1  1 -1
    ## [40417] -1 -1 -1  1 -1  1  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1  1
    ## [40441] -1  1  1 -1 -1  1  1  1  1  1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1  1 -1
    ## [40465] -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1
    ## [40489] -1  1  1  1  1  1 -1 -1  1 -1  1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1
    ## [40513] -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1  1 -1  1 -1
    ## [40537] -1 -1  1  1  1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1
    ## [40561] -1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1
    ## [40585]  1 -1  1  1 -1  1  1 -1 -1 -1 -1  1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1
    ## [40609]  1 -1 -1 -1  1 -1  1 -1  1  1  1  1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1
    ## [40633]  1 -1  1  1  1 -1  1  1 -1 -1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1  1  1
    ## [40657] -1  1 -1  1  1  1 -1 -1  1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1  1  1  1  1
    ## [40681] -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1  1
    ## [40705] -1  1  1 -1 -1  1 -1  1 -1  1  1  1  1 -1 -1 -1 -1  1  1  1  1  1  1  1
    ## [40729]  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1  1  1 -1 -1  1 -1 -1 -1
    ## [40753] -1 -1 -1  1  1  1  1  1  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1  1  1  1 -1
    ## [40777]  1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1  1  1  1 -1  1 -1 -1 -1  1  1  1
    ## [40801] -1 -1 -1  1 -1  1 -1  1  1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1 -1  1
    ## [40825]  1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1  1  1  1 -1  1 -1 -1
    ## [40849]  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1
    ## [40873] -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1
    ## [40897] -1  1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1  1
    ## [40921] -1  1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1 -1  1  1  1 -1 -1  1 -1  1 -1
    ## [40945] -1  1  1 -1 -1  1 -1 -1  1  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1
    ## [40969]  1  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1
    ## [40993]  1  1  1 -1  1 -1 -1  1 -1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1
    ## [41017] -1  1 -1  1 -1  1 -1 -1  1  1  1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1
    ## [41041] -1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1 -1  1  1  1  1 -1  1
    ## [41065] -1  1  1 -1 -1  1  1  1  1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1  1
    ## [41089] -1 -1  1  1 -1  1  1 -1  1 -1 -1 -1 -1  1  1  1  1 -1  1  1 -1  1  1 -1
    ## [41113] -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1  1 -1  1  1  1 -1 -1  1 -1
    ## [41137]  1  1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1
    ## [41161]  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1 -1
    ## [41185] -1 -1  1 -1 -1  1  1 -1  1  1  1  1  1 -1  1 -1  1  1  1 -1 -1  1  1 -1
    ## [41209]  1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1  1 -1
    ## [41233]  1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1  1 -1 -1 -1  1  1
    ## [41257] -1 -1 -1 -1  1  1 -1 -1 -1  1  1  1  1  1  1 -1 -1  1  1 -1  1  1  1 -1
    ## [41281]  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1
    ## [41305]  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1  1  1  1 -1  1  1
    ## [41329] -1  1 -1  1  1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1
    ## [41353] -1  1 -1 -1 -1  1  1  1  1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1
    ## [41377] -1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1  1  1  1  1  1  1  1  1
    ## [41401] -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1
    ## [41425]  1  1 -1 -1  1 -1 -1  1  1  1  1  1  1  1 -1  1  1 -1 -1  1  1  1 -1 -1
    ## [41449]  1 -1 -1 -1  1 -1  1 -1  1  1  1  1  1  1  1  1 -1  1 -1 -1 -1  1 -1 -1
    ## [41473]  1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1  1  1  1 -1  1 -1  1 -1 -1  1
    ## [41497]  1 -1 -1  1 -1  1  1  1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1
    ## [41521] -1 -1  1  1  1 -1 -1 -1  1  1  1  1 -1  1  1  1  1 -1 -1 -1  1  1  1 -1
    ## [41545]  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1  1 -1  1  1
    ## [41569]  1 -1  1  1  1  1 -1 -1 -1  1  1 -1  1  1  1  1  1 -1  1  1  1 -1 -1 -1
    ## [41593]  1  1  1  1 -1 -1  1  1 -1  1  1  1  1 -1  1 -1 -1  1  1  1  1  1  1 -1
    ## [41617]  1 -1  1  1  1 -1 -1  1  1  1 -1 -1 -1 -1  1  1 -1  1  1  1 -1 -1  1  1
    ## [41641]  1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1  1
    ## [41665]  1  1  1 -1  1 -1  1  1 -1 -1  1  1  1 -1  1  1 -1  1 -1  1 -1  1 -1 -1
    ## [41689] -1  1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1  1  1  1 -1 -1  1  1 -1
    ## [41713] -1  1  1  1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1  1  1  1 -1  1  1
    ## [41737]  1 -1  1 -1 -1  1  1  1  1 -1  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1 -1 -1
    ## [41761]  1  1  1  1  1  1  1  1 -1 -1  1 -1  1 -1  1  1 -1  1  1 -1  1  1 -1 -1
    ## [41785]  1  1  1  1  1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1  1  1  1  1 -1  1  1
    ## [41809] -1  1 -1  1 -1  1  1 -1 -1  1 -1  1  1  1  1 -1  1 -1  1  1 -1  1 -1 -1
    ## [41833]  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1  1  1 -1 -1  1  1
    ## [41857] -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1  1  1  1  1
    ## [41881] -1 -1 -1  1  1 -1  1  1  1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1
    ## [41905]  1  1  1  1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1  1
    ## [41929] -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1  1  1 -1 -1 -1  1 -1  1
    ## [41953]  1 -1  1  1  1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1
    ## [41977]  1  1 -1  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1
    ## [42001]  1  1  1  1 -1  1  1  1 -1  1  1  1  1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1
    ## [42025]  1 -1 -1  1 -1  1 -1 -1  1  1  1  1 -1  1 -1 -1 -1  1  1  1 -1  1  1 -1
    ## [42049]  1  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1  1 -1  1
    ## [42073] -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1
    ## [42097] -1 -1 -1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1
    ## [42121] -1  1  1  1 -1  1  1  1 -1  1  1  1  1 -1  1  1 -1  1  1 -1 -1  1 -1 -1
    ## [42145] -1 -1  1  1 -1  1  1  1 -1 -1  1  1  1  1  1 -1 -1  1 -1  1  1  1  1 -1
    ## [42169]  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1 -1 -1
    ## [42193]  1  1  1  1  1  1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1 -1  1  1  1 -1  1
    ## [42217]  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1
    ## [42241] -1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1
    ## [42265]  1  1 -1  1 -1 -1  1  1  1  1  1  1  1 -1  1  1  1  1 -1 -1  1  1  1 -1
    ## [42289]  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1  1  1  1 -1
    ## [42313]  1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1
    ## [42337] -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1 -1  1 -1  1  1  1 -1 -1  1  1  1 -1
    ## [42361]  1  1  1  1  1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1  1  1  1 -1  1  1  1
    ## [42385] -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1  1  1 -1 -1 -1  1
    ## [42409]  1  1 -1  1  1  1  1  1  1  1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1  1  1
    ## [42433] -1  1  1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1  1 -1  1  1 -1  1 -1 -1  1
    ## [42457]  1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1
    ## [42481] -1  1  1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1
    ## [42505]  1  1  1  1 -1  1 -1  1 -1  1  1 -1 -1  1 -1 -1  1  1  1  1 -1 -1  1 -1
    ## [42529]  1 -1  1  1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1
    ## [42553]  1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1  1  1  1 -1
    ## [42577] -1  1  1 -1  1 -1  1  1  1  1  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1 -1
    ## [42601]  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1
    ## [42625]  1  1 -1 -1  1 -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1
    ## [42649] -1  1  1  1 -1  1  1  1  1  1  1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1  1 -1
    ## [42673]  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1  1  1  1
    ## [42697]  1 -1  1  1 -1  1  1 -1  1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1
    ## [42721]  1  1  1  1 -1  1 -1  1 -1 -1  1  1 -1  1  1  1  1 -1  1 -1  1  1 -1 -1
    ## [42745] -1  1 -1 -1 -1  1  1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1  1  1  1 -1 -1
    ## [42769]  1  1  1  1  1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1
    ## [42793] -1 -1 -1  1 -1  1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1
    ## [42817] -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1
    ## [42841]  1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1  1 -1 -1
    ## [42865]  1 -1  1 -1 -1  1  1 -1  1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1
    ## [42889] -1  1  1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1  1
    ## [42913]  1  1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1
    ## [42937]  1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1  1  1
    ## [42961] -1  1  1  1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1  1  1 -1  1  1 -1 -1 -1
    ## [42985] -1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1  1 -1  1  1  1  1
    ## [43009] -1 -1 -1  1  1  1  1  1  1  1  1  1 -1  1  1  1  1  1 -1 -1 -1  1 -1 -1
    ## [43033] -1  1  1  1 -1 -1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1  1  1
    ## [43057]  1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1 -1
    ## [43081]  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1  1  1 -1  1  1 -1  1  1  1
    ## [43105] -1  1 -1  1  1  1  1 -1  1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1  1  1  1
    ## [43129] -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1
    ## [43153]  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1 -1 -1  1
    ## [43177]  1  1 -1  1 -1  1 -1 -1  1  1  1  1 -1  1  1 -1  1  1 -1  1 -1 -1 -1 -1
    ## [43201]  1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1 -1  1  1 -1
    ## [43225] -1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [43249] -1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1  1  1
    ## [43273]  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1
    ## [43297]  1  1  1 -1  1 -1 -1 -1 -1  1  1  1 -1  1  1  1 -1 -1  1 -1 -1 -1  1  1
    ## [43321]  1  1  1 -1 -1 -1  1 -1  1  1  1  1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1
    ## [43345] -1 -1  1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1
    ## [43369]  1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1
    ## [43393] -1  1 -1  1 -1  1 -1 -1  1  1 -1  1 -1 -1  1  1  1 -1  1 -1  1 -1 -1 -1
    ## [43417]  1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1  1 -1
    ## [43441] -1 -1  1  1  1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1  1
    ## [43465]  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1  1  1 -1  1 -1  1 -1  1 -1 -1
    ## [43489]  1 -1 -1  1  1 -1 -1  1 -1  1  1  1 -1 -1  1 -1  1  1  1  1 -1  1 -1 -1
    ## [43513]  1 -1 -1  1 -1  1  1  1  1  1  1  1  1  1 -1  1  1  1 -1  1 -1  1 -1 -1
    ## [43537] -1 -1 -1  1 -1  1  1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1
    ## [43561]  1 -1 -1  1  1  1 -1  1  1  1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1
    ## [43585] -1 -1 -1  1 -1  1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1
    ## [43609] -1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1
    ## [43633]  1  1 -1 -1  1  1  1  1 -1 -1 -1  1 -1  1  1  1 -1  1 -1  1 -1 -1 -1  1
    ## [43657]  1  1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1  1
    ## [43681]  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1
    ## [43705] -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1  1 -1  1
    ## [43729] -1 -1  1  1  1 -1 -1  1  1  1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1
    ## [43753] -1  1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1  1 -1  1  1  1 -1  1 -1 -1 -1
    ## [43777]  1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1  1 -1 -1  1
    ## [43801]  1 -1  1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1
    ## [43825] -1 -1 -1 -1  1  1  1 -1  1  1  1  1 -1  1 -1 -1  1  1  1 -1  1 -1  1  1
    ## [43849]  1 -1  1 -1  1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1  1 -1 -1  1 -1
    ## [43873] -1  1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1
    ## [43897] -1  1  1 -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1  1  1 -1 -1
    ## [43921]  1  1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1  1  1 -1  1 -1  1  1  1  1 -1
    ## [43945] -1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1
    ## [43969] -1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [43993]  1 -1 -1  1  1 -1  1  1 -1  1 -1 -1  1 -1  1  1 -1  1  1 -1  1  1 -1  1
    ## [44017] -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1
    ## [44041] -1 -1  1  1  1 -1 -1  1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1
    ## [44065] -1 -1  1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1
    ## [44089] -1  1  1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1
    ## [44113] -1  1 -1  1  1  1 -1 -1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ## [44137] -1  1 -1  1  1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1
    ## [44161] -1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1 -1 -1  1 -1 -1 -1  1 -1
    ## [44185] -1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1  1
    ## [44209] -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1  1  1  1  1 -1  1 -1  1  1 -1 -1  1
    ## [44233] -1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1  1
    ## [44257]  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1  1 -1  1
    ## [44281]  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1  1
    ## [44305] -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1
    ## [44329]  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1
    ## [44353]  1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1  1
    ## [44377] -1 -1  1  1  1  1  1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1
    ## [44401]  1 -1  1  1 -1 -1  1  1  1  1 -1 -1 -1  1  1 -1  1  1  1 -1 -1  1 -1 -1
    ## [44425]  1  1 -1  1  1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1  1  1  1
    ## [44449] -1  1 -1  1 -1  1  1  1  1  1  1 -1  1  1 -1  1 -1 -1  1  1  1 -1 -1  1
    ## [44473]  1 -1  1 -1 -1 -1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1
    ## [44497]  1  1  1 -1  1  1 -1 -1  1  1  1 -1 -1  1  1 -1  1  1  1  1  1 -1  1  1
    ## [44521]  1 -1 -1  1  1  1 -1  1  1  1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1  1
    ## [44545]  1  1 -1 -1 -1  1  1  1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1  1
    ## [44569]  1 -1 -1  1 -1  1  1  1  1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [44593] -1  1  1 -1 -1  1 -1 -1  1 -1  1  1  1 -1  1  1 -1 -1  1  1 -1  1 -1  1
    ## [44617]  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1  1  1  1 -1
    ## [44641] -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1  1  1  1
    ## [44665]  1  1  1 -1  1  1  1 -1 -1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1  1  1  1
    ## [44689]  1 -1  1  1 -1  1  1  1 -1  1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1
    ## [44713] -1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1
    ## [44737]  1  1  1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1  1 -1  1 -1  1 -1  1 -1 -1
    ## [44761]  1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1 -1  1
    ## [44785] -1 -1 -1  1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1
    ## [44809] -1 -1  1  1 -1  1  1  1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1  1
    ## [44833]  1 -1 -1  1 -1  1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1
    ## [44857]  1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1  1  1  1  1  1  1 -1
    ## [44881] -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1
    ## [44905]  1  1 -1 -1  1 -1  1  1  1  1  1 -1 -1  1  1 -1  1  1  1 -1 -1 -1  1  1
    ## [44929] -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1  1 -1
    ## [44953]  1  1  1 -1  1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1 -1
    ## [44977] -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1 -1  1  1  1
    ## [45001]  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1
    ## [45025]  1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    ## [45049]  1  1 -1  1 -1  1 -1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1
    ## [45073] -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [45097]  1  1 -1  1  1 -1  1  1 -1  1 -1  1  1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1
    ## [45121]  1  1 -1  1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1
    ## [45145]  1 -1 -1  1 -1  1 -1  1  1 -1  1 -1  1 -1  1 -1 -1  1  1  1  1 -1 -1  1
    ## [45169]  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1 -1
    ## [45193] -1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1
    ## [45217]  1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1  1  1  1  1  1  1
    ## [45241] -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1
    ## [45265]  1  1  1  1 -1 -1  1 -1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1
    ## [45289] -1  1  1  1  1  1  1  1 -1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1  1 -1  1
    ## [45313] -1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1
    ## [45337] -1 -1 -1 -1 -1  1  1 -1  1 -1  1  1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1
    ## [45361] -1  1  1 -1 -1  1  1  1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1
    ## [45385] -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1  1  1  1 -1  1 -1  1  1  1 -1  1
    ## [45409] -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1  1  1 -1 -1  1  1  1 -1
    ## [45433] -1  1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1
    ## [45457] -1 -1  1  1  1  1 -1  1  1  1 -1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1  1
    ## [45481]  1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1  1  1 -1
    ## [45505] -1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1
    ## [45529] -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1  1 -1 -1 -1  1  1
    ## [45553] -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1  1 -1 -1  1 -1  1 -1
    ## [45577] -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1  1 -1  1 -1  1  1
    ## [45601]  1  1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1
    ## [45625] -1 -1  1  1  1 -1  1  1  1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1
    ## [45649]  1  1  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1
    ## [45673] -1 -1 -1  1  1 -1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1  1  1 -1 -1 -1  1
    ## [45697] -1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [45721]  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1  1  1 -1  1 -1  1  1 -1 -1 -1
    ## [45745]  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1
    ## [45769] -1 -1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1 -1  1  1 -1  1 -1  1 -1  1
    ## [45793] -1  1 -1  1  1 -1  1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1
    ## [45817] -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1
    ## [45841]  1 -1  1  1  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1  1  1 -1  1  1 -1 -1
    ## [45865] -1  1  1  1 -1 -1  1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1
    ## [45889] -1 -1  1  1 -1 -1  1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1
    ## [45913]  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1  1  1 -1  1  1 -1 -1
    ## [45937]  1  1  1  1 -1  1  1 -1 -1  1  1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1
    ## [45961]  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1
    ## [45985] -1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1
    ## [46009]  1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1
    ## [46033] -1  1  1  1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [46057]  1 -1  1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1
    ## [46081]  1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1  1  1 -1
    ## [46105]  1 -1  1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1  1 -1
    ## [46129] -1 -1  1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1  1 -1  1
    ## [46153] -1  1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1  1  1 -1 -1  1  1  1 -1 -1  1
    ## [46177] -1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1
    ## [46201] -1  1  1  1 -1  1 -1  1  1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1  1
    ## [46225] -1 -1 -1 -1 -1  1 -1 -1  1  1  1  1  1 -1  1 -1 -1  1  1  1  1  1 -1 -1
    ## [46249]  1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1  1
    ## [46273]  1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [46297] -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1
    ## [46321]  1 -1 -1  1 -1  1  1  1 -1 -1 -1  1  1  1 -1 -1  1  1  1 -1 -1  1  1  1
    ## [46345] -1  1  1 -1  1  1 -1  1  1 -1  1  1  1 -1  1  1  1 -1 -1  1 -1  1  1  1
    ## [46369] -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1
    ## [46393]  1 -1  1  1  1  1 -1 -1  1 -1 -1  1  1  1  1 -1 -1  1  1 -1  1  1 -1  1
    ## [46417] -1  1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1  1
    ## [46441]  1  1  1 -1  1  1  1 -1  1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1
    ## [46465]  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1  1  1
    ## [46489] -1  1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1 -1  1 -1  1 -1 -1  1  1  1  1
    ## [46513] -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1  1 -1  1
    ## [46537]  1  1  1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1
    ## [46561] -1  1  1  1  1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1  1  1  1  1 -1 -1
    ## [46585]  1  1 -1  1  1 -1  1 -1  1  1  1 -1 -1  1 -1  1 -1  1  1 -1  1 -1  1 -1
    ## [46609]  1 -1  1  1 -1  1 -1 -1  1  1 -1  1  1  1 -1  1 -1  1  1 -1 -1  1  1  1
    ## [46633] -1 -1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1  1  1 -1 -1 -1  1 -1  1 -1
    ## [46657] -1  1  1  1 -1  1 -1  1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1
    ## [46681]  1 -1  1  1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1
    ## [46705]  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1
    ## [46729]  1  1  1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1  1  1 -1  1  1  1  1  1  1  1
    ## [46753] -1  1  1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1  1 -1  1  1  1  1  1  1
    ## [46777] -1 -1 -1 -1  1  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1
    ## [46801]  1 -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1  1 -1
    ## [46825] -1 -1 -1  1 -1  1 -1  1  1 -1  1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1
    ## [46849]  1  1  1  1  1 -1 -1  1 -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1
    ## [46873] -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1
    ## [46897] -1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1
    ## [46921]  1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1  1 -1  1 -1 -1 -1  1 -1  1
    ## [46945] -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1 -1  1  1
    ## [46969]  1 -1  1 -1  1 -1 -1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1
    ## [46993]  1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1 -1  1
    ## [47017] -1 -1 -1  1 -1  1  1  1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1  1  1  1  1
    ## [47041] -1  1  1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1  1  1  1  1  1 -1  1  1  1
    ## [47065]  1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1  1  1  1  1  1 -1  1  1 -1 -1
    ## [47089] -1 -1  1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1  1 -1  1  1 -1  1
    ## [47113] -1 -1  1  1 -1  1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1 -1 -1  1  1  1 -1
    ## [47137] -1 -1  1 -1  1 -1  1  1 -1 -1 -1  1 -1  1 -1  1 -1  1  1  1  1 -1 -1 -1
    ## [47161] -1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1
    ## [47185] -1  1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1  1  1  1
    ## [47209] -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1 -1
    ## [47233]  1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1 -1
    ## [47257] -1  1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1  1  1 -1  1  1 -1
    ## [47281]  1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1  1  1 -1  1
    ## [47305]  1  1 -1 -1 -1 -1  1  1  1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1 -1  1 -1
    ## [47329] -1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1  1  1  1 -1 -1  1  1 -1  1
    ## [47353] -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1  1  1 -1  1 -1  1 -1  1  1
    ## [47377] -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1 -1  1  1  1  1 -1  1  1 -1
    ## [47401] -1  1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1
    ## [47425] -1  1  1 -1 -1  1  1 -1  1 -1  1 -1  1  1  1  1  1 -1 -1  1  1 -1 -1  1
    ## [47449]  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1
    ## [47473]  1  1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1  1  1  1  1  1 -1  1 -1  1 -1
    ## [47497]  1  1  1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1
    ## [47521]  1  1  1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1  1
    ## [47545]  1 -1  1 -1  1  1  1  1 -1 -1 -1  1 -1  1  1  1  1  1 -1  1  1  1 -1 -1
    ## [47569]  1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1  1  1  1 -1  1  1
    ## [47593]  1 -1  1  1  1  1  1 -1  1  1  1  1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1
    ## [47617]  1  1 -1  1 -1 -1  1  1  1  1 -1  1 -1 -1  1 -1 -1  1 -1  1  1  1  1  1
    ## [47641]  1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1
    ## [47665]  1 -1  1  1  1  1  1 -1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1  1  1  1 -1
    ## [47689]  1 -1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1
    ## [47713] -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1  1 -1  1  1  1 -1 -1  1  1
    ## [47737]  1  1  1 -1  1  1  1  1 -1  1 -1 -1 -1  1  1 -1  1 -1  1 -1  1  1 -1 -1
    ## [47761]  1 -1  1 -1  1 -1  1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1  1 -1  1
    ## [47785]  1 -1 -1 -1  1  1  1 -1  1 -1  1  1 -1 -1 -1  1  1  1 -1 -1  1  1  1  1
    ## [47809]  1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1
    ## [47833] -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1  1 -1  1  1  1  1  1  1  1  1
    ## [47857]  1  1  1 -1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1
    ## [47881]  1  1  1  1  1  1  1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1 -1
    ## [47905]  1 -1 -1  1  1  1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1 -1  1  1  1
    ## [47929]  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ## [47953]  1 -1 -1  1 -1  1  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1
    ## [47977] -1  1  1 -1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1 -1
    ## [48001] -1 -1 -1  1 -1  1 -1  1  1  1 -1  1  1  1  1 -1 -1 -1  1  1 -1  1  1  1
    ## [48025] -1 -1 -1 -1 -1  1  1 -1  1  1  1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1
    ## [48049] -1  1 -1 -1  1  1 -1 -1  1  1  1  1  1  1  1 -1  1 -1  1  1 -1 -1  1  1
    ## [48073] -1 -1  1 -1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1 -1  1
    ## [48097]  1  1 -1  1 -1  1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1
    ## [48121] -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1  1 -1 -1  1
    ## [48145]  1  1  1 -1  1 -1 -1 -1 -1  1  1 -1  1  1  1  1  1  1  1 -1  1 -1  1  1
    ## [48169]  1 -1 -1  1 -1  1  1  1  1  1  1 -1  1 -1  1 -1  1 -1  1  1 -1  1  1  1
    ## [48193] -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1 -1  1  1  1
    ## [48217] -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1
    ## [48241]  1 -1 -1 -1  1  1  1  1 -1  1  1 -1  1  1  1  1  1  1 -1  1  1 -1 -1 -1
    ## [48265]  1  1  1 -1 -1  1  1  1  1  1 -1  1 -1  1  1  1 -1  1 -1  1  1  1 -1  1
    ## [48289] -1  1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1  1  1  1 -1  1 -1
    ## [48313] -1  1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1 -1  1 -1 -1  1 -1  1  1 -1 -1
    ## [48337] -1 -1 -1  1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1  1  1  1  1  1 -1 -1 -1
    ## [48361] -1  1  1 -1  1 -1  1 -1 -1  1  1 -1  1 -1  1 -1 -1  1 -1  1  1 -1 -1  1
    ## [48385] -1  1  1  1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1 -1  1
    ## [48409] -1  1  1  1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1
    ## [48433] -1  1  1  1  1  1  1  1 -1  1  1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1
    ## [48457] -1  1  1  1  1  1 -1  1 -1  1  1  1  1  1 -1  1 -1 -1 -1 -1  1  1  1 -1
    ## [48481] -1  1 -1 -1  1  1  1  1  1  1 -1  1 -1 -1  1  1  1  1 -1 -1 -1 -1  1 -1
    ## [48505] -1  1 -1  1 -1  1 -1 -1 -1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1  1  1  1
    ## [48529]  1 -1  1 -1  1 -1  1  1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1
    ## [48553] -1  1 -1  1  1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1
    ## [48577] -1 -1  1  1  1  1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1  1 -1 -1  1  1
    ## [48601] -1  1  1 -1  1  1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1  1  1 -1 -1  1  1
    ## [48625] -1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1  1  1 -1 -1  1  1 -1  1 -1 -1 -1
    ## [48649]  1  1  1 -1  1  1  1  1  1  1  1 -1  1  1 -1 -1 -1  1  1 -1 -1  1  1  1
    ## [48673] -1  1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1
    ## [48697] -1  1 -1 -1  1  1  1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1
    ## [48721] -1 -1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1
    ## [48745]  1  1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1  1  1 -1  1  1  1 -1  1 -1  1 -1
    ## [48769] -1 -1  1 -1 -1  1  1  1  1 -1  1  1  1 -1  1 -1  1  1 -1  1  1  1  1 -1
    ## [48793]  1  1 -1 -1  1  1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1
    ## [48817]  1  1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1 -1 -1
    ## [48841] -1 -1  1 -1  1 -1 -1  1 -1 -1  1  1  1  1  1 -1  1 -1  1 -1  1  1 -1 -1
    ## [48865]  1  1 -1  1 -1  1  1  1 -1  1 -1  1  1 -1  1  1  1 -1  1  1  1 -1  1 -1
    ## [48889]  1 -1  1  1  1  1  1  1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1
    ## [48913]  1 -1 -1 -1  1 -1 -1  1  1  1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1  1
    ## [48937] -1 -1 -1  1 -1  1 -1  1  1  1 -1  1 -1  1  1 -1  1 -1  1 -1 -1  1 -1 -1
    ## [48961] -1  1  1  1  1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1
    ## [48985]  1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1
    ## [49009] -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1
    ## [49033]  1  1  1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1  1  1 -1  1 -1  1  1 -1
    ## [49057] -1 -1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1  1
    ## [49081] -1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1  1 -1
    ## [49105] -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1 -1
    ## [49129]  1 -1  1  1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1
    ## [49153] -1  1 -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1  1  1 -1 -1  1  1  1 -1  1
    ## [49177] -1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1
    ## [49201]  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1  1 -1 -1  1  1  1  1  1
    ## [49225] -1  1  1 -1  1  1  1  1 -1  1  1 -1 -1 -1  1  1 -1  1  1 -1  1  1 -1  1
    ## [49249]  1 -1  1  1 -1  1 -1  1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1  1  1  1 -1
    ## [49273]  1  1 -1  1 -1  1  1  1  1 -1  1 -1  1  1  1  1  1 -1  1  1  1 -1  1 -1
    ## [49297]  1  1  1  1  1  1  1  1 -1  1 -1 -1  1  1  1  1  1  1  1  1  1 -1 -1 -1
    ## [49321]  1  1  1  1  1 -1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1
    ## [49345] -1 -1  1  1  1 -1  1  1 -1  1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1  1  1
    ## [49369]  1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1 -1  1 -1  1  1  1  1 -1  1 -1  1
    ## [49393] -1  1 -1 -1  1 -1 -1  1  1 -1 -1  1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1
    ## [49417]  1 -1 -1  1  1 -1 -1 -1  1  1 -1  1 -1  1  1  1  1 -1 -1  1  1 -1  1 -1
    ## [49441]  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1
    ## [49465] -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1  1 -1 -1 -1  1
    ## [49489]  1  1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1  1  1  1  1  1  1  1  1 -1 -1
    ## [49513]  1  1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1  1
    ## [49537] -1 -1  1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1
    ## [49561]  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1  1 -1  1 -1  1 -1  1
    ## [49585] -1 -1  1  1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1  1
    ## [49609] -1  1  1  1 -1 -1  1  1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1
    ## [49633]  1  1  1  1  1 -1  1 -1  1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1  1
    ## [49657]  1 -1  1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1
    ## [49681] -1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1
    ## [49705] -1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1  1  1  1
    ## [49729] -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1
    ## [49753]  1  1 -1  1  1 -1  1 -1 -1  1  1  1 -1 -1  1  1  1  1 -1  1 -1  1 -1  1
    ## [49777]  1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1 -1  1 -1 -1  1 -1  1
    ## [49801]  1  1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1 -1
    ## [49825] -1  1  1 -1 -1  1  1 -1 -1  1  1  1  1  1 -1 -1  1  1  1 -1  1  1  1 -1
    ## [49849]  1  1 -1 -1 -1  1  1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1  1  1 -1
    ## [49873] -1  1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1
    ## [49897] -1  1 -1  1 -1  1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [49921]  1  1 -1 -1 -1  1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1
    ## [49945]  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1  1
    ## [49969]  1 -1  1  1 -1  1  1 -1  1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1
    ## [49993]  1  1  1  1  1 -1  1  1 -1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1
    ## [50017] -1  1  1  1 -1  1  1 -1  1  1  1  1 -1 -1  1  1  1  1  1 -1  1 -1  1 -1
    ## [50041] -1  1  1 -1 -1  1 -1  1 -1  1  1  1  1  1 -1  1 -1 -1  1  1 -1  1  1 -1
    ## [50065] -1  1  1  1  1 -1  1 -1 -1  1 -1 -1  1  1 -1  1 -1  1  1 -1  1 -1  1  1
    ## [50089] -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1
    ## [50113]  1 -1  1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1  1 -1  1  1 -1  1 -1  1  1  1
    ## [50137] -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1  1 -1 -1  1
    ## [50161]  1 -1  1 -1  1  1  1  1 -1  1 -1 -1  1 -1 -1 -1  1  1  1 -1  1  1  1  1
    ## [50185]  1  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1  1 -1  1
    ## [50209] -1 -1  1  1 -1 -1 -1  1 -1  1 -1  1  1 -1  1 -1 -1  1 -1  1  1  1  1  1
    ## [50233] -1  1  1  1 -1  1  1  1  1  1 -1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1
    ## [50257] -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ## [50281]  1 -1  1 -1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1
    ## [50305]  1 -1  1  1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1 -1  1  1  1  1 -1 -1
    ## [50329]  1 -1 -1  1 -1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1  1  1  1
    ## [50353] -1 -1  1  1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1 -1  1 -1
    ## [50377] -1 -1  1  1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1  1  1  1  1  1 -1 -1 -1
    ## [50401]  1  1 -1  1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1
    ## [50425]  1 -1  1  1  1 -1 -1 -1 -1  1  1  1  1 -1  1  1  1  1  1 -1  1 -1  1 -1
    ## [50449]  1  1 -1  1  1  1  1  1  1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1
    ## [50473] -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1  1 -1  1  1  1  1  1
    ## [50497] -1  1  1  1  1  1 -1 -1 -1  1  1  1  1  1 -1  1  1  1  1  1  1 -1 -1 -1
    ## [50521]  1  1  1  1  1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1  1  1
    ## [50545] -1  1 -1  1  1 -1  1  1  1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1  1  1 -1
    ## [50569]  1 -1  1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1  1 -1  1  1  1  1 -1
    ## [50593] -1  1 -1 -1  1  1  1 -1 -1  1  1  1 -1 -1 -1  1 -1  1 -1 -1  1  1  1  1
    ## [50617] -1  1  1 -1  1  1  1  1  1 -1 -1 -1  1  1 -1 -1  1 -1  1  1 -1  1 -1  1
    ## [50641]  1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1
    ## [50665] -1  1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1  1 -1  1  1 -1 -1  1  1  1 -1  1
    ## [50689]  1 -1  1 -1  1  1  1  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1 -1
    ## [50713]  1 -1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1 -1  1  1  1  1 -1  1  1 -1  1
    ## [50737] -1 -1 -1  1  1  1  1  1 -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1  1
    ## [50761] -1  1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1  1 -1  1 -1  1  1 -1  1  1
    ## [50785]  1  1  1 -1  1  1  1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1 -1  1
    ## [50809] -1  1 -1  1 -1  1  1  1 -1 -1 -1  1  1  1  1 -1 -1 -1  1 -1 -1  1  1  1
    ## [50833] -1  1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1  1  1  1 -1  1  1  1  1 -1 -1
    ## [50857] -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1  1  1  1  1  1
    ## [50881] -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1  1  1  1  1  1 -1 -1
    ## [50905] -1  1  1 -1  1 -1 -1  1 -1  1  1 -1  1  1  1  1  1  1  1  1  1 -1 -1 -1
    ## [50929] -1  1 -1 -1  1  1 -1  1 -1 -1  1  1  1 -1 -1  1  1  1  1  1 -1  1 -1 -1
    ## [50953]  1  1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1 -1  1 -1  1  1 -1
    ## [50977]  1  1  1  1 -1 -1  1 -1  1 -1  1  1  1 -1 -1  1 -1 -1 -1  1  1  1 -1  1
    ## [51001] -1 -1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1
    ## [51025] -1  1  1  1 -1  1 -1  1  1  1  1 -1  1  1  1  1  1 -1  1 -1  1  1 -1 -1
    ## [51049] -1  1 -1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1  1  1  1  1  1  1
    ## [51073]  1 -1  1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1
    ## [51097]  1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1  1  1 -1  1  1  1 -1  1 -1  1
    ## [51121] -1 -1 -1  1 -1  1  1  1  1  1 -1  1 -1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1
    ## [51145]  1 -1  1  1 -1 -1  1  1  1 -1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1
    ## [51169] -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1
    ## [51193] -1  1  1  1  1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1  1
    ## [51217] -1 -1 -1  1  1 -1  1  1 -1  1 -1  1  1 -1  1 -1  1  1 -1 -1 -1  1  1 -1
    ## [51241]  1 -1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1  1 -1  1  1
    ## [51265] -1  1  1 -1  1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1  1  1
    ## [51289]  1  1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1
    ## [51313]  1 -1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1  1  1  1  1
    ## [51337] -1 -1  1  1  1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1 -1
    ## [51361]  1  1 -1  1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1  1  1
    ## [51385]  1  1  1  1 -1 -1 -1  1  1  1  1  1 -1  1 -1  1  1 -1 -1  1  1  1 -1 -1
    ## [51409] -1  1 -1 -1  1  1  1  1  1 -1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1
    ## [51433]  1  1  1  1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1  1  1  1 -1 -1  1 -1  1
    ## [51457]  1  1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1  1  1  1  1 -1  1 -1 -1 -1 -1  1
    ## [51481]  1  1  1 -1  1 -1  1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1  1  1 -1
    ## [51505]  1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1 -1  1  1  1 -1  1 -1 -1 -1  1  1
    ## [51529]  1  1  1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1 -1  1 -1 -1 -1 -1
    ## [51553]  1  1 -1 -1 -1  1 -1  1  1  1 -1  1 -1  1 -1  1 -1  1  1 -1 -1  1  1 -1
    ## [51577]  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1 -1
    ## [51601] -1  1  1 -1 -1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1  1
    ## [51625]  1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1  1  1  1  1  1  1 -1 -1 -1 -1 -1  1
    ## [51649] -1 -1 -1 -1 -1  1 -1  1  1 -1  1 -1  1 -1  1 -1  1  1 -1 -1  1 -1  1 -1
    ## [51673]  1 -1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1 -1  1  1  1  1
    ## [51697] -1  1  1  1 -1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1  1  1 -1 -1
    ## [51721] -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1
    ## [51745]  1 -1  1  1 -1  1  1 -1 -1  1  1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1  1  1
    ## [51769]  1 -1  1 -1 -1  1 -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1
    ## [51793] -1 -1 -1  1  1  1  1 -1  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1 -1  1  1  1
    ## [51817] -1 -1  1  1  1  1  1  1 -1  1  1 -1  1  1 -1  1  1  1 -1  1 -1 -1 -1  1
    ## [51841]  1  1  1  1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1
    ## [51865] -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1
    ## [51889] -1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1  1  1 -1
    ## [51913]  1 -1 -1 -1  1 -1 -1  1  1  1  1 -1  1 -1  1 -1 -1  1 -1 -1 -1 -1 -1  1
    ## [51937]  1  1  1  1  1 -1  1  1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1
    ## [51961] -1 -1 -1 -1  1 -1 -1  1 -1  1 -1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1  1
    ## [51985]  1 -1 -1 -1 -1  1 -1  1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1
    ## [52009]  1 -1 -1  1  1 -1 -1  1  1  1 -1 -1  1  1  1  1 -1  1  1  1  1  1 -1 -1
    ## [52033] -1  1 -1  1  1  1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1  1  1  1 -1  1
    ## [52057] -1  1  1 -1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1
    ## [52081] -1  1 -1 -1  1  1 -1 -1  1  1  1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1
    ## [52105] -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1  1  1 -1 -1  1  1  1 -1  1 -1
    ## [52129] -1 -1  1  1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1  1
    ## [52153]  1 -1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1 -1 -1  1 -1 -1  1  1  1  1 -1
    ## [52177] -1 -1  1 -1 -1 -1 -1  1  1  1  1  1  1  1 -1 -1  1 -1 -1  1  1 -1 -1  1
    ## [52201] -1 -1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1  1  1 -1
    ## [52225]  1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1 -1  1  1  1  1  1  1 -1 -1  1
    ## [52249] -1 -1  1  1  1  1 -1 -1  1  1  1 -1 -1  1  1  1 -1  1  1 -1  1  1 -1 -1
    ## [52273]  1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1 -1  1  1 -1  1  1  1 -1
    ## [52297]  1  1 -1  1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1  1  1  1 -1  1 -1 -1  1
    ## [52321]  1  1  1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1  1 -1  1  1  1  1  1 -1 -1
    ## [52345]  1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1  1 -1 -1  1 -1 -1 -1  1 -1
    ## [52369]  1 -1  1  1  1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1  1  1  1  1  1
    ## [52393] -1 -1 -1 -1 -1  1  1  1  1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1 -1  1
    ## [52417]  1 -1 -1 -1  1  1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1
    ## [52441] -1  1  1 -1  1  1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1 -1  1
    ## [52465] -1  1  1  1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1  1  1
    ## [52489] -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1
    ## [52513]  1 -1 -1  1 -1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1  1  1 -1  1
    ## [52537]  1  1  1  1  1 -1  1  1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1 -1 -1 -1
    ## [52561]  1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1 -1  1
    ## [52585]  1 -1  1  1 -1 -1  1 -1 -1  1 -1  1 -1  1  1 -1  1  1  1  1 -1 -1 -1  1
    ## [52609] -1 -1  1  1  1  1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1 -1  1  1  1 -1 -1
    ## [52633]  1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1
    ## [52657] -1 -1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1 -1 -1 -1 -1
    ## [52681] -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1  1  1  1 -1  1  1  1 -1
    ## [52705] -1  1 -1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1  1  1 -1 -1 -1 -1  1 -1
    ## [52729] -1 -1  1  1 -1  1  1  1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1  1  1  1 -1
    ## [52753] -1  1 -1  1  1  1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1  1  1  1 -1  1  1  1
    ## [52777] -1  1  1  1  1 -1 -1  1 -1  1 -1  1  1  1  1  1  1  1 -1  1  1 -1 -1  1
    ## [52801]  1  1  1 -1  1 -1  1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1  1  1 -1  1
    ## [52825]  1 -1  1  1  1 -1  1  1  1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1 -1 -1  1
    ## [52849] -1  1 -1 -1 -1 -1  1  1  1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1
    ## [52873]  1  1 -1  1  1 -1  1  1 -1  1  1  1  1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1
    ## [52897] -1  1  1 -1  1 -1  1  1  1  1  1 -1  1  1 -1 -1  1  1  1 -1  1 -1 -1 -1
    ## [52921]  1  1  1  1 -1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1  1 -1  1  1  1 -1  1 -1
    ## [52945] -1  1 -1  1  1  1  1 -1 -1  1  1  1  1 -1  1  1  1 -1  1  1  1  1  1  1
    ## [52969] -1 -1 -1 -1 -1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1  1 -1  1  1  1 -1  1
    ## [52993]  1  1  1  1  1  1  1  1  1  1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1  1 -1  1
    ## [53017]  1  1  1 -1  1 -1  1  1  1 -1 -1  1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1
    ## [53041] -1  1 -1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1 -1 -1
    ## [53065]  1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1  1  1
    ## [53089] -1 -1 -1 -1  1  1 -1 -1  1  1  1 -1  1  1  1 -1 -1 -1  1  1 -1  1  1  1
    ## [53113]  1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1
    ## [53137] -1  1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1
    ## [53161]  1  1  1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1  1 -1 -1
    ## [53185] -1 -1  1  1  1  1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1  1  1  1  1
    ## [53209] -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1
    ## [53233]  1 -1 -1  1 -1  1 -1 -1  1  1  1 -1  1  1  1  1  1  1 -1  1  1  1  1 -1
    ## [53257]  1 -1 -1 -1  1  1  1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1 -1  1 -1  1
    ## [53281]  1  1  1  1  1 -1 -1 -1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1  1 -1  1
    ## [53305] -1 -1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1  1 -1 -1  1 -1 -1  1  1  1 -1
    ## [53329] -1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1  1  1  1 -1 -1  1  1  1  1
    ## [53353] -1 -1  1 -1 -1  1 -1  1  1 -1  1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1
    ## [53377]  1 -1  1  1 -1 -1  1 -1  1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1 -1
    ## [53401] -1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1 -1  1  1  1 -1 -1
    ## [53425] -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1  1  1  1  1 -1  1  1 -1  1 -1 -1 -1
    ## [53449]  1 -1  1  1  1 -1  1 -1  1  1 -1  1 -1  1  1  1  1  1  1 -1  1 -1  1 -1
    ## [53473]  1 -1 -1 -1  1  1  1  1 -1  1 -1  1 -1 -1  1  1 -1 -1 -1  1  1  1  1 -1
    ## [53497] -1  1  1  1 -1  1  1  1 -1  1  1 -1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1
    ## [53521] -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1  1  1 -1 -1  1
    ## [53545]  1 -1  1 -1  1 -1  1  1  1  1  1  1  1 -1  1  1 -1  1 -1 -1  1  1 -1 -1
    ## [53569] -1  1  1  1  1  1  1 -1  1  1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1  1 -1
    ## [53593]  1 -1  1  1 -1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1  1  1 -1 -1
    ## [53617] -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1  1  1 -1  1 -1  1
    ## [53641] -1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1
    ## [53665] -1  1  1 -1  1 -1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1  1  1 -1
    ## [53689] -1  1 -1 -1 -1  1  1  1  1 -1  1 -1  1 -1 -1  1 -1 -1  1  1 -1  1  1  1
    ## [53713]  1 -1  1  1 -1  1 -1  1 -1 -1  1 -1 -1  1  1  1  1  1 -1  1 -1  1 -1 -1
    ## [53737] -1  1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1  1  1 -1 -1 -1  1  1 -1  1
    ## [53761] -1 -1 -1  1 -1 -1 -1 -1  1 -1  1  1  1  1  1  1  1 -1 -1 -1 -1  1  1  1
    ## [53785]  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1  1  1  1  1  1
    ## [53809] -1 -1 -1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1  1
    ## [53833]  1 -1  1 -1  1  1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1  1
    ## [53857] -1  1  1  1 -1  1  1 -1 -1  1 -1  1  1  1  1 -1  1 -1 -1 -1 -1 -1 -1 -1
    ## [53881] -1  1  1  1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1  1  1 -1  1  1 -1 -1
    ## [53905]  1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1  1 -1 -1  1  1  1 -1
    ## [53929]  1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1  1 -1 -1 -1
    ## [53953] -1  1 -1  1  1  1 -1  1  1  1 -1  1 -1 -1  1 -1  1  1  1 -1  1  1 -1 -1
    ## [53977] -1 -1 -1  1  1  1  1 -1  1  1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1  1 -1  1
    ## [54001] -1 -1  1 -1  1  1 -1  1 -1  1  1 -1  1  1 -1  1  1  1 -1  1  1  1  1  1
    ## [54025]  1 -1  1 -1 -1  1  1 -1  1 -1  1  1 -1  1  1  1  1  1  1  1 -1 -1 -1  1
    ## [54049]  1  1  1  1 -1  1 -1  1 -1 -1  1  1  1  1  1  1 -1  1 -1 -1 -1  1  1 -1
    ## [54073] -1 -1 -1  1  1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1
    ## [54097]  1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1  1  1  1  1 -1  1
    ## [54121]  1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1
    ## [54145]  1  1  1  1  1 -1  1  1 -1  1 -1  1  1 -1 -1  1 -1 -1  1  1  1 -1  1  1
    ## [54169]  1 -1  1 -1  1 -1  1  1  1  1  1 -1  1 -1  1 -1 -1 -1  1  1  1 -1 -1  1
    ## [54193] -1  1 -1  1 -1 -1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1  1  1 -1  1 -1 -1
    ## [54217]  1 -1 -1  1 -1  1  1 -1 -1  1 -1  1  1  1  1  1  1  1 -1 -1 -1  1 -1 -1
    ## [54241] -1 -1  1 -1  1  1 -1 -1  1  1  1 -1  1  1  1  1  1  1  1 -1 -1  1 -1 -1
    ## [54265]  1 -1 -1  1 -1  1 -1  1 -1  1  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1
    ## [54289] -1  1  1 -1  1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1 -1  1
    ## [54313]  1  1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1  1 -1 -1 -1  1  1 -1  1 -1
    ## [54337] -1  1  1  1  1  1 -1  1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1  1 -1 -1
    ## [54361] -1  1  1  1  1 -1  1 -1  1 -1  1  1  1 -1  1  1 -1 -1  1 -1 -1 -1  1 -1
    ## [54385]  1  1 -1  1 -1  1 -1 -1  1  1 -1  1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1
    ## [54409] -1  1 -1 -1 -1  1  1  1 -1 -1 -1  1 -1  1 -1  1 -1  1 -1 -1 -1  1 -1  1
    ## [54433] -1  1  1 -1  1 -1 -1 -1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1  1 -1
    ## [54457] -1  1  1  1  1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1
    ## [54481]  1 -1 -1  1  1  1 -1  1 -1 -1 -1  1 -1  1  1 -1  1  1  1  1 -1 -1 -1 -1
    ## [54505] -1 -1 -1 -1 -1  1  1 -1  1  1 -1  1  1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1
    ## [54529] -1  1 -1  1 -1  1 -1  1  1  1  1 -1  1  1  1 -1  1  1 -1 -1  1 -1 -1  1
    ## [54553]  1 -1 -1  1  1 -1 -1  1  1  1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1  1  1  1
    ## [54577]  1  1  1 -1  1  1 -1  1  1  1 -1  1 -1  1  1  1 -1 -1 -1 -1  1  1 -1 -1
    ## [54601]  1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1 -1 -1 -1  1  1  1 -1 -1  1
    ## [54625] -1  1 -1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1  1  1  1  1  1  1 -1  1 -1
    ## [54649] -1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1  1 -1  1  1 -1  1 -1
    ## [54673] -1  1  1  1 -1  1 -1  1  1  1  1  1 -1 -1 -1 -1  1  1  1  1  1  1 -1  1
    ## [54697]  1 -1  1  1 -1 -1 -1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1  1  1 -1  1 -1  1
    ## [54721] -1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1 -1  1  1  1 -1  1 -1  1 -1  1 -1
    ## [54745]  1  1  1 -1  1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1  1 -1  1  1  1  1 -1  1
    ## [54769]  1 -1 -1  1  1  1  1  1 -1  1  1 -1  1  1 -1  1 -1 -1  1 -1 -1  1  1  1
    ## [54793]  1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1  1 -1 -1  1 -1
    ## [54817]  1  1  1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1  1 -1  1  1 -1  1  1  1
    ## [54841] -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [54865] -1 -1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1
    ## [54889]  1  1  1 -1  1 -1 -1 -1  1  1  1 -1  1  1 -1  1 -1 -1  1  1 -1  1 -1  1
    ## [54913] -1  1 -1  1 -1  1  1 -1 -1 -1  1  1  1 -1  1 -1  1 -1 -1  1  1 -1  1  1
    ## [54937]  1  1  1  1  1 -1  1  1 -1  1  1 -1  1  1  1 -1 -1  1 -1 -1  1  1 -1 -1
    ## [54961] -1 -1 -1 -1  1 -1 -1 -1 -1  1 -1 -1  1  1  1  1 -1 -1  1  1 -1  1 -1  1
    ## [54985]  1  1  1 -1  1 -1  1  1  1  1 -1  1  1 -1  1 -1  1  1  1 -1  1  1 -1 -1
    ## [55009] -1  1  1  1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1 -1  1  1  1
    ## [55033] -1 -1  1  1  1  1 -1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1 -1  1
    ## [55057]  1  1  1  1 -1 -1 -1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1 -1 -1
    ## [55081]  1  1 -1 -1 -1 -1  1 -1  1  1  1  1  1  1  1 -1 -1 -1  1  1  1 -1  1  1
    ## [55105] -1  1 -1 -1 -1 -1  1  1  1  1 -1  1  1  1  1  1 -1  1 -1 -1 -1  1  1  1
    ## [55129] -1 -1  1  1 -1  1  1  1  1  1 -1  1  1  1  1 -1 -1 -1  1  1  1  1 -1  1
    ## [55153]  1 -1  1 -1 -1  1  1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1 -1
    ## [55177] -1 -1  1  1  1 -1 -1  1 -1  1  1  1  1 -1 -1 -1  1 -1  1  1 -1  1 -1 -1
    ## [55201] -1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1  1 -1  1  1  1 -1  1  1  1  1 -1
    ## [55225]  1  1  1  1  1  1  1 -1  1  1  1 -1  1  1  1  1 -1 -1 -1  1 -1 -1 -1 -1
    ## [55249]  1  1 -1  1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1  1  1 -1  1 -1  1  1  1
    ## [55273] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1  1  1  1 -1 -1 -1  1 -1  1
    ## [55297] -1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 -1  1  1  1
    ## [55321] -1  1  1  1  1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1  1
    ## [55345]  1  1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1 -1  1 -1  1  1
    ## [55369] -1  1 -1 -1  1 -1 -1 -1  1  1  1  1  1 -1  1 -1  1 -1  1  1  1 -1 -1  1
    ## [55393] -1 -1 -1  1  1 -1  1 -1 -1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1
    ## [55417]  1  1 -1 -1  1  1  1  1  1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1  1 -1
    ## [55441] -1  1  1 -1  1 -1  1 -1  1 -1  1  1  1 -1 -1  1  1 -1  1  1  1 -1  1 -1
    ## [55465]  1  1  1  1 -1  1 -1  1 -1 -1  1  1  1 -1  1 -1 -1  1 -1 -1  1 -1 -1 -1
    ## [55489]  1 -1 -1  1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 -1 -1 -1 -1  1
    ## [55513]  1 -1  1  1  1 -1  1  1 -1  1  1  1 -1 -1  1  1 -1 -1 -1 -1  1  1 -1  1
    ## [55537] -1 -1  1 -1  1 -1  1  1  1 -1 -1  1  1 -1  1 -1  1  1 -1 -1  1 -1 -1 -1
    ## [55561] -1  1 -1 -1 -1  1 -1  1  1  1 -1  1 -1  1 -1  1 -1  1  1  1 -1 -1  1 -1
    ## [55585] -1  1 -1 -1  1 -1 -1  1 -1  1  1 -1 -1  1  1  1 -1 -1  1 -1 -1 -1 -1 -1
    ## [55609] -1  1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1  1  1  1 -1 -1 -1  1  1
    ## [55633]  1  1  1 -1  1  1  1 -1 -1  1  1  1  1 -1  1 -1 -1 -1 -1  1  1 -1  1 -1
    ## [55657] -1  1  1 -1 -1  1  1  1  1 -1  1 -1  1  1 -1 -1 -1 -1  1 -1  1  1  1  1
    ## [55681] -1 -1  1  1  1  1 -1 -1 -1  1  1  1  1  1 -1 -1 -1 -1  1 -1 -1 -1 -1  1
    ## [55705]  1  1 -1  1  1 -1  1  1 -1  1 -1  1  1  1 -1 -1  1 -1  1  1  1 -1 -1 -1
    ## [55729] -1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1  1
    ## [55753] -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1  1  1 -1  1 -1  1 -1 -1 -1  1
    ## [55777]  1  1  1  1 -1  1  1  1  1  1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1
    ## [55801]  1  1 -1  1 -1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1
    ## [55825]  1 -1  1  1  1 -1  1 -1  1 -1  1  1 -1  1  1 -1 -1  1  1 -1  1 -1 -1 -1
    ## [55849]  1  1  1 -1  1  1  1  1 -1  1 -1 -1  1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1
    ## [55873]  1  1  1  1  1 -1  1  1 -1 -1 -1  1 -1  1 -1  1  1  1  1 -1 -1 -1  1 -1
    ## [55897] -1 -1 -1  1  1  1  1  1 -1  1 -1  1  1  1  1 -1  1  1 -1 -1 -1 -1 -1  1
    ## [55921]  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1  1  1  1 -1
    ## [55945]  1  1  1 -1 -1  1  1 -1 -1 -1  1  1  1 -1  1 -1 -1 -1 -1 -1  1 -1  1 -1
    ## [55969]  1  1 -1 -1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1  1  1  1  1  1 -1 -1
    ## [55993] -1  1  1 -1 -1  1 -1 -1 -1 -1  1  1  1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1
    ## [56017] -1  1 -1 -1  1  1  1  1  1  1  1  1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1  1
    ## [56041]  1  1  1 -1  1 -1  1 -1  1 -1 -1 -1  1  1  1 -1 -1  1 -1  1  1 -1 -1  1
    ## [56065] -1 -1  1  1  1 -1 -1 -1  1  1 -1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1  1  1
    ## [56089]  1 -1 -1  1  1 -1 -1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1 -1 -1
    ## [56113]  1  1  1 -1 -1  1  1 -1 -1  1  1  1  1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1
    ## [56137]  1 -1 -1  1  1  1 -1  1  1 -1  1 -1  1 -1  1  1 -1 -1  1 -1  1 -1  1 -1
    ## [56161] -1  1  1  1 -1  1 -1  1  1  1 -1 -1  1 -1  1 -1 -1 -1 -1  1 -1 -1 -1 -1
    ## [56185]  1 -1  1 -1 -1 -1  1  1  1 -1  1  1  1 -1  1 -1  1 -1

``` r
# Create a histogram of IC50 values to visualise them
hist(dasatinib_sort$LN_IC50,50)
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20IC50-4.png)<!-- -->

``` r
# Take a quick look at the median, just out of interest
median(dasatinib_sort$LN_IC50)
```

    ## [1] 2.544

``` r
head(tt)
```

    ##                   rn   logFC AveExpr      t   P.Value adj.P.Val      B   symbol
    ## 1: ENSG00000178734.4  1.5040 -6.2460  4.432 2.654e-05    0.4564 -1.930   LMO7DN
    ## 2: ENSG00000241739.1 -1.0466 -6.1297 -4.278 4.742e-05    0.4564 -2.098     <NA>
    ## 3: ENSG00000260689.2  0.9576 -3.6336  4.151 7.575e-05    0.4564 -2.235     <NA>
    ## 4: ENSG00000249242.3  2.0343 -0.0662  4.051 1.086e-04    0.4564 -2.340 TMEM150C
    ## 5: ENSG00000241088.2 -1.0050 -6.9632 -3.983 1.389e-04    0.4564 -2.412     <NA>
    ## 6: ENSG00000197157.6 -0.4814  8.1291 -3.960 1.506e-04    0.4564 -2.436     SND1

``` r
#Find the top expression values. Shows that high expression correlates to high ic50, and V.V (e.g. this gene is protective against dasatinib)
topExp <- expDat_sort[match(tt$rn[1], rownames(expDat_sort)),]
df <- data.frame(topGene=topExp, ic50=group)
ggplot(df, aes(x=ic50, y=topGene)) + geom_boxplot()
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20IC50-5.png)<!-- -->

``` r
# Check all this stuff
dim(expDat_sort)
```

    ## [1] 56202    89

``` r
dim(dasatinib_sort)
```

    ## [1] 89 20

``` r
head(expDat_sort)[,1:5]
```

    ##                     SW48 SW1710  SW620  CAL51 TOV112D
    ## ENSG00000223972.4 -2.488 -2.425 -4.965 -5.718  -4.657
    ## ENSG00000227232.4  3.255  2.690  3.394  3.233   3.081
    ## ENSG00000243485.2 -3.225 -3.462 -4.965 -5.718  -4.946
    ## ENSG00000237613.2 -2.608 -3.066 -4.435 -5.233  -5.309
    ## ENSG00000268020.2 -3.962 -4.251 -5.813 -5.718  -5.309
    ## ENSG00000240361.1 -2.178 -4.540 -4.435 -6.455  -5.309

``` r
head(dasatinib_sort)[,1:5]
```

    ##   DATASET NLME_RESULT_ID NLME_CURVE_ID COSMIC_ID   CCLE_Name
    ## 1   GDSC2            290      14783021    909751        SW48
    ## 2   GDSC2            290      14783023    905962       SW620
    ## 3   GDSC2            290      14783024    910927       CAL51
    ## 4   GDSC2            290      14783037    908445        SNU5
    ## 5   GDSC2            290      14783041    909718       SKCO1
    ## 6   GDSC2            290      14783049    908120 MDAMB175VII

``` r
sum(colnames(expDat_sort)==dasatinib_sort$CCLE_Name)
```

    ## [1] 1

``` r
names(dasatinib_sort)
```

    ##  [1] "DATASET"         "NLME_RESULT_ID"  "NLME_CURVE_ID"   "COSMIC_ID"      
    ##  [5] "CCLE_Name"       "SANGER_MODEL_ID" "TCGA_DESC"       "DRUG_ID"        
    ##  [9] "DRUG_NAME"       "PUTATIVE_TARGET" "PATHWAY_NAME"    "COMPANY_ID"     
    ## [13] "WEBRELEASE"      "MIN_CONC"        "MAX_CONC"        "LN_IC50"        
    ## [17] "AUC"             "RMSE"            "Z_SCORE"         "AUC_Level"

``` r
# lmFit
group2 <- ifelse(dasatinib_sort$AUC > median(dasatinib_sort$AUC), "High", "Low")
table(group2)
```

    ## group2
    ## High  Low 
    ##   44   45

``` r
boxplot(dasatinib_sort$AUC ~ group2)
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20AUC-1.png)<!-- -->

``` r
design2 = model.matrix(~group2);
design2 %>% head()
```

    ##   (Intercept) group2Low
    ## 1           1         0
    ## 2           1         0
    ## 3           1         1
    ## 4           1         0
    ## 5           1         1
    ## 6           1         0

``` r
colnames(design2) = c("Mean"
,"HighVsLow"
)
# Fit the data to the model
fit2 = lmFit(expDat_sort, design2)
fit2 = eBayes(fit2)
tt2 = topTable(fit2, coef="HighVsLow", adjust="BH",n=nrow(expDat_sort))
options(digits=4)
tt2[1:20,]
```

    ##                      logFC AveExpr      t   P.Value adj.P.Val      B
    ## ENSG00000154479.8   1.4857 -1.2359  4.862 4.970e-06    0.1231 2.8392
    ## ENSG00000186860.3  -1.0677 -6.8270 -4.847 5.288e-06    0.1231 2.7936
    ## ENSG00000213309.3  -1.1914 -4.0815 -4.792 6.572e-06    0.1231 2.6337
    ## ENSG00000241739.1  -1.1063 -6.1297 -4.577 1.521e-05    0.1890 2.0157
    ## ENSG00000197157.6  -0.5404  8.1291 -4.551 1.682e-05    0.1890 1.9417
    ## ENSG00000152620.8  -0.5668  5.1901 -4.498 2.066e-05    0.1918 1.7902
    ## ENSG00000186973.6   2.2447 -3.4341  4.460 2.389e-05    0.1918 1.6829
    ## ENSG00000249242.3   2.1887 -0.0662  4.424 2.737e-05    0.1923 1.5829
    ## ENSG00000239264.4  -0.8616  8.1796 -4.295 4.442e-05    0.2699 1.2259
    ## ENSG00000175564.8   0.8985  0.3289  4.274 4.803e-05    0.2699 1.1683
    ## ENSG00000178035.7  -0.6662  7.6231 -4.211 6.079e-05    0.3106 0.9947
    ## ENSG00000224346.1   1.1469 -6.1748  4.167 7.141e-05    0.3345 0.8761
    ## ENSG00000153046.13 -0.4763  5.1241 -4.123 8.394e-05    0.3613 0.7570
    ## ENSG00000232450.1   1.1997 -1.2796  4.098 9.179e-05    0.3613 0.6911
    ## ENSG00000270714.1  -1.0875 -5.4596 -4.084 9.642e-05    0.3613 0.6549
    ## ENSG00000231335.1  -1.1336 -6.9190 -4.013 1.246e-04    0.3902 0.4662
    ## ENSG00000111752.6   1.3810  3.9467  4.008 1.270e-04    0.3902 0.4521
    ## ENSG00000250471.2  -0.6051 -2.0419 -3.983 1.388e-04    0.3902 0.3867
    ## ENSG00000261355.1   0.8207  1.6939  3.980 1.404e-04    0.3902 0.3780
    ## ENSG00000214145.2   2.1684 -2.0940  3.962 1.494e-04    0.3902 0.3324

``` r
# Take a quick look at the first gene for the hell of it
plot(density(expDat_sort[,1]))
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20AUC-2.png)<!-- -->

``` r
# Count number of significant genes
sum(tt2$adj.P.Val<0.01)
```

    ## [1] 0

``` r
## Plot log fold-change _versus_ -log(P-value). This graph is the important one, things up the top and spread out are useful associations
## (i.e., higher number = lower p-value):
volcanoplot(fit2, coef="HighVsLow")

# Fancier Volcano plot, with lines for adjusted p-value, fold change, and highlighted important samples.
sigFC2 = (tt2$adj.P.Val < 0.01)  & (abs(tt2$logFC) > 1)
volcanoplot(fit2, coef="HighVsLow")
points(tt2$logFC[which(sigFC2)], 
       -log10(tt2$P.Value[which(sigFC2)]), 
       cex=0.6, col='red', pch=16)
abline(h = min(-log10(tt2$P.Value[which(sigFC2)])), lty=2, col='blue')
```

    ## Warning in min(-log10(tt2$P.Value[which(sigFC2)])): no non-missing arguments to
    ## min; returning Inf

``` r
abline(v = c(-1,1), lty=2, col='blue')
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20AUC-3.png)<!-- -->

``` r
split2 <- strsplit(rownames(tt2),".", fixed=T) %>% lapply(., function(x) x[1]) %>% unlist()
geneNames2 <- AnnotationDbi::select(org.Hs.eg.db, keys = split2, column = c("SYMBOL","GENENAME"), key="ENSEMBL")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
head(split2)
```

    ## [1] "ENSG00000154479" "ENSG00000186860" "ENSG00000213309" "ENSG00000241739"
    ## [5] "ENSG00000197157" "ENSG00000152620"

``` r
# What does the top table look like?
dim(tt2)
```

    ## [1] 56202     6

``` r
tt2$symbol <- geneNames2$SYMBOL[match(split2, geneNames2$ENSEMBL)]
setDT(tt2, keep.rownames = TRUE)[]
```

    ##                        rn      logFC AveExpr          t   P.Value adj.P.Val
    ##     1:  ENSG00000154479.8  1.486e+00  -1.236  4.862e+00 4.970e-06    0.1231
    ##     2:  ENSG00000186860.3 -1.068e+00  -6.827 -4.847e+00 5.288e-06    0.1231
    ##     3:  ENSG00000213309.3 -1.191e+00  -4.081 -4.792e+00 6.572e-06    0.1231
    ##     4:  ENSG00000241739.1 -1.106e+00  -6.130 -4.577e+00 1.521e-05    0.1890
    ##     5:  ENSG00000197157.6 -5.404e-01   8.129 -4.551e+00 1.682e-05    0.1890
    ##    ---                                                                     
    ## 56198:  ENSG00000247844.1 -3.293e-04  -2.058 -3.587e-04 9.997e-01    0.9998
    ## 56199:  ENSG00000253623.1 -5.410e-05  -7.767 -2.973e-04 9.998e-01    0.9998
    ## 56200: ENSG00000061273.13 -2.949e-05   6.094 -1.551e-04 9.999e-01    0.9999
    ## 56201: ENSG00000115107.15 -3.192e-05   4.868 -7.708e-05 9.999e-01    0.9999
    ## 56202:  ENSG00000250441.1  1.997e-05  -6.875  7.659e-05 9.999e-01    0.9999
    ##             B    symbol
    ##     1:  2.839   CCDC173
    ##     2:  2.794 KRTAP17-1
    ##     3:  2.634      <NA>
    ##     4:  2.016      <NA>
    ##     5:  1.942      SND1
    ##    ---                 
    ## 56198: -5.391      <NA>
    ## 56199: -5.391      <NA>
    ## 56200: -5.391     HDAC7
    ## 56201: -5.391    STEAP3
    ## 56202: -5.391      <NA>

``` r
# Check the data: Output shows that some ENSEMBL ID's have multiple regular gene names
length(split2)
```

    ## [1] 56202

``` r
dim(geneNames2)
```

    ## [1] 56436     3

``` r
# Print the top 300 associated genes. Can then look them up on Enrichr, GeneSetDB, etc
cat(na.omit(tt2$symbol[1:300]),sep="\n")
```

    ## CCDC173
    ## KRTAP17-1
    ## SND1
    ## NADK2
    ## FAM183A
    ## TMEM150C
    ## TXNDC5
    ## UCP3
    ## IMPDH2
    ## CDYL
    ## PHC1
    ## LINC00887
    ## RIMBP2
    ## FAM66C
    ## MAPK14
    ## FAM104B
    ## MLLT11
    ## PPP2R2C
    ## RNF213-AS1
    ## GNG3
    ## STMN3
    ## DCTN1-AS1
    ## SHANK2-AS1
    ## NOA1
    ## SCAMP5
    ## TTC9B
    ## CACNA1D
    ## UBOX5-AS1
    ## MBD2
    ## PCDH12
    ## SDHA
    ## LINC01563
    ## SLC26A3
    ## KCNAB3
    ## LINC01772
    ## GJC2
    ## KIAA2013
    ## SYP
    ## ZNF704
    ## CPLX2
    ## LINC02290
    ## IGFBP2
    ## HLA-E
    ## SHISA2
    ## LINC00641
    ## ZNF536
    ## RAB41
    ## ADGRB2
    ## ADRA2C
    ## LINCR-0001
    ## TMEM169
    ## ATL1
    ## ULK2
    ## NDUFAF3
    ## LIN7A
    ## PRR36
    ## NHLH1
    ## MDH1B
    ## NPDC1
    ## KCNH3
    ## CACHD1
    ## UBE2E1
    ## LINC01620
    ## CHRNB2
    ## ACP4
    ## FSIP1
    ## LINC01431
    ## LINC01107
    ## ZNF713
    ## CHGB
    ## SEZ6
    ## QARS1
    ## TCAP
    ## SNURF
    ## TNFRSF1B
    ## LRRC71
    ## SEH1L
    ## ZNF513
    ## PROX1-AS1
    ## SNORA47
    ## IFI30
    ## ZBTB20
    ## CCDC65
    ## LOC112268389
    ## SMAD1-AS2
    ## SPIN2A
    ## B4GALNT4
    ## LAG3
    ## LINC02570
    ## SNORD101
    ## NAAA
    ## BCL2L12
    ## ADGRE3
    ## CBLL1
    ## PLCB1
    ## ESRRB
    ## FAM169A
    ## MACROH2A2
    ## SP7
    ## JAM3
    ## CAVIN4
    ## ENO2
    ## EIF3CL
    ## ZNF286A
    ## SCN3B
    ## SUCLG2
    ## USP17L7
    ## PI4K2B
    ## LINC01854
    ## IRGC
    ## KDM1B
    ## CTB-178M22.2
    ## AMH
    ## RALBP1
    ## CHN2
    ## LYRM4
    ## MIR33B
    ## SSR1
    ## C12orf57
    ## TTLL7
    ## LINC01701
    ## MICOS13
    ## BAZ2B

``` r
# Create a histogram of AUC values to visualise them
hist(dasatinib_sort$AUC,50)
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20AUC-4.png)<!-- -->

``` r
# Take a quick look at the median, just out of interest
median(dasatinib_sort$AUC)
```

    ## [1] 0.9541

``` r
head(tt2)
```

    ##                   rn   logFC AveExpr      t   P.Value adj.P.Val     B    symbol
    ## 1: ENSG00000154479.8  1.4857  -1.236  4.862 4.970e-06    0.1231 2.839   CCDC173
    ## 2: ENSG00000186860.3 -1.0677  -6.827 -4.847 5.288e-06    0.1231 2.794 KRTAP17-1
    ## 3: ENSG00000213309.3 -1.1914  -4.081 -4.792 6.572e-06    0.1231 2.634      <NA>
    ## 4: ENSG00000241739.1 -1.1063  -6.130 -4.577 1.521e-05    0.1890 2.016      <NA>
    ## 5: ENSG00000197157.6 -0.5404   8.129 -4.551 1.682e-05    0.1890 1.942      SND1
    ## 6: ENSG00000152620.8 -0.5668   5.190 -4.498 2.066e-05    0.1918 1.790     NADK2

``` r
#Find the top expression values. Shows that high expression correlates to low AUC, and V.V (e.g. this gene is protective against dasatinib)
topExp2 <- expDat_sort[match(tt2$rn[1], rownames(expDat_sort)),]
df2 <- data.frame(topGene=topExp2, AUC=group)
ggplot(df2, aes(x=AUC, y=topGene)) + geom_boxplot()
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Differential%20Expression%20Analysis%20With%20AUC-5.png)<!-- -->

``` r
library(gplots)
```

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     space

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     space

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
## top 500 from each:
venn(list(ic50=tt$rn[1:500], auc=tt2$rn[1:500]))
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/look%20at%20overlap%20between%20IC50%20and%20AUC%20results-1.png)<!-- -->

``` r
overlap <- c()
sz <- c(50, 100, 150, 200, 300, 500, 750, 1000, 2000, 5000, 10000, 20000, 40000, 50000)
for(i in 1:length(sz)) overlap[i] <- length(intersect(tt$rn[1:sz[i]], tt2$rn[1:sz[i]]))
expect <- sz^2/nrow(tt) 
plot(sz, expect, type='b', ylim=c(0,max(overlap)))
points(sz, overlap, type='b', col='red')
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/look%20at%20overlap%20between%20IC50%20and%20AUC%20results-2.png)<!-- -->

``` r
## This shows that we're still getting more overlap between the two lists (tt/IC50, tt2/AUC) than we would expect to see by chance, so the relationship is still important, even if nothing appears to be significant. Instead, we'll just use the top 500 genes from the final list (tt3).
```

``` r
# Create a table with ensembl ID, gene symbol, adj. p value, and sign for AUC and IC50
tt3 <- full_join(tt, tt2, by= "rn")
tt3 <- dplyr::select(tt3, -c("symbol.x", "P.Value.y", "P.Value.x", "t.y", "t.x", "B.x", "B.y", "AveExpr.x"))
tt3 <- dplyr::rename(tt3, "Gene_Symbol" = "symbol.y", "AUC_logFC" = "logFC.y", "Avg_Exp" = "AveExpr.y", "AUC_Adj_PVal" = "adj.P.Val.y", "IC50_logFC" = "logFC.x", "IC50_Adj_PVal" = "adj.P.Val.x", "Ensembl_ID" = "rn")
tt3 <- tt3[c(7,1,5,2,3,4,6)]
# Find out if the log fold change sign is consistent with significance. Remember, -logFC relates to resistance
tt3 <- mutate(tt3, sign(tt3$IC50_logFC), sign(tt3$AUC_logFC))
tt3 <- mutate(tt3, sign(tt3$IC50_logFC)*sign(tt3$AUC_logFC))
# Multiply signs together, keep genes where the result is negative (cause you want two different signs)
tt3 <- dplyr::rename(tt3, "logFC_Sign" = "sign(tt3$IC50_logFC) * sign(tt3$AUC_logFC)", "IC50_Sign" = "sign(tt3$IC50_logFC)", "AUC_Sign" = "sign(tt3$AUC_logFC)")
# Look at the samples that are potentially important- don't get good results :(
###################### IMPORTANT!!! ######################
### IF ALL THE DATA IS ACTUALLY -LOG, NOT LN, THEN THE RELATIONSHIP BETWEEN AUC AND IC50(-log) SHOULD BE +VE. THIS MEANS THAT THE SIGN WE'RE LOOKING FOR WILL BE +VE, NOT -VE! FROM HERE, WE'RE ASSUMING THEIR LABELLING IS WRONG! Adding in +ve AUC sign now for SigSamples, genes that are related to susceptibility to dasatinib will have a larger AUC (and theoretically larger "-log" IC50)
# Coming back to this, 05/12/19
# Filter based on significance, increased susceptibility to dasatinib, larger log fold change and positive correlation between AUC and "-log" IC50
SigSamples <- tt3[1:500,]
# Rank the samples based on average rank between the two p-values
SigSamples <- 
  SigSamples %>% mutate(., rank_ic50=rank(IC50_Adj_PVal)) %>%
  mutate(., rank_auc=rank(AUC_Adj_PVal)) %>% 
  mutate(., avg_rank=0.5*(rank_auc + rank_ic50)) %>% 
  arrange(., avg_rank)
# Find the top 400 ranked genes, and look them up
cat(na.omit(SigSamples$Gene_Symbol[1:400]),sep="\n")
```

    ## KRTAP17-1
    ## CCDC173
    ## SND1
    ## FAM183A
    ## NADK2
    ## TMEM150C
    ## TXNDC5
    ## IMPDH2
    ## PPP2R2C
    ## FAM66C
    ## CACNA1D
    ## STMN3
    ## NOA1
    ## SCAMP5
    ## SHANK2-AS1
    ## MBD2
    ## LINC01563
    ## KIAA2013
    ## ACP4
    ## ADRA2C
    ## ADGRB2
    ## HLA-E
    ## QARS1
    ## KCNH3
    ## SNORA47
    ## LAG3
    ## SMAD1-AS2
    ## ADGRE3
    ## ENO2
    ## CAVIN4
    ## CTB-178M22.2
    ## CHN2
    ## LINC01701
    ## BAZ2B
    ## SNORA2A
    ## PSMG2
    ## RFTN2
    ## RPL17
    ## HADHA
    ## JAKMIP3
    ## FZD10-AS1
    ## PCCB
    ## ZNF620
    ## RPS10
    ## RPL29P2
    ## DPYSL5
    ## NSG1
    ## FOXO6
    ## MIR9-3HG
    ## MEF2C
    ## ATF6B
    ## PRH1
    ## NAA80
    ## APBA2
    ## ARHGEF11
    ## LINC01506
    ## RPL27A
    ## QPCTL
    ## TMOD2
    ## TMEM81
    ## SORCS3
    ## TP73
    ## MGC27382
    ## OR14J1
    ## RNF157-AS1
    ## LOC101929586
    ## ZCCHC12
    ## LOC541473
    ## PRLR
    ## DNMT3A
    ## CLEC17A
    ## WASF2
    ## SHB
    ## OR13C3
    ## PSME1
    ## FZD10
    ## LOC101928222
    ## SRRM3
    ## FAM160B2
    ## DYM
    ## DNAH9
    ## GNRHR2
    ## RPN1
    ## FAM174B
    ## DDOST
    ## SNORA28
    ## MMP3
    ## PRKAA2
    ## MEF2C-AS1
    ## SLC26A8
    ## C18orf21
    ## NMRK1
    ## MADD-AS1
    ## ALG3
    ## EZR-AS1
    ## GLRA1
    ## IRF1
    ## MCEE
    ## ARHGEF10
    ## LOC100996447
    ## LOC107986453
    ## PEX11B
    ## MPP2
    ## NIPAL4
    ## WSB2
    ## HCN2
    ## ZFYVE16
    ## AMDHD1
    ## ERCC3
    ## GRN
    ## GNG4
    ## OR52B4
    ## MIR4505
    ## MOGS
    ## TUBB2B
    ## RGS12
    ## RPL10A
    ## ATP5F1A
    ## TDH
    ## CD47
    ## GBX2
    ## RIOK3
    ## PSMB11
    ## LINC01487
    ## TNFRSF17
    ## CTD-2194D22.4
    ## DPEP3
    ## MIR377
    ## KCNIP3
    ## C2CD4C
    ## NRIP2
    ## OBSL1
    ## MRO
    ## HAUS4
    ## LMO7DN
    ## GYS2
    ## PIGN
    ## KIAA1328
    ## OR14A16
    ## ACSBG1
    ## FLG-AS1
    ## DOCK4-AS1
    ## SNORD48
    ## LOC102724265
    ## MIR3200
    ## EDARADD
    ## LINC02652
    ## MMP17
    ## MIR181B2
    ## ZNF217
    ## MIR874
    ## DISC1-IT1
    ## ANKRD34C-AS1
    ## LOC101928516
    ## MSX2
    ## OR4S1
    ## LOC100506551
    ## INHA
    ## CHRNA1
    ## VAX1
    ## PPP2R2B-IT1
    ## ELP2
    ## CITED1
    ## SLC35G5
    ## GABRR1

``` r
# Sort based on average significance rank across AUC and IC50 (filter based on significant adj. p values, and consistent signs between the two analysis)
tt3 <- 
  tt3 %>% mutate(., rank_ic50=rank(IC50_Adj_PVal)) %>%
  mutate(., rank_auc=rank(AUC_Adj_PVal)) %>% 
  mutate(., avg_rank=0.5*(rank_auc + rank_ic50)) %>% 
  arrange(., avg_rank)

# When you actually get a good list of genes, analyse with enrichr and genesetDB (Use more strict requirements for these two (e.g. Adj_P_Val < 0.01, abs(logFC > 1) than you do for the following RectomePA and GoSeq stuff)
```

``` r
# Gene expression vs IC50
# Boxplot
ggplot(data=tt3, mapping=aes(x=IC50_Sign,y=Avg_Exp, group=IC50_Sign, fill=IC50_Sign)) +
  geom_boxplot()
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Plot%20gene%20expression%20vs%20AUC%20and%20gene%20expression%20vs%20IC50,%20on%20a%20scatter%20plot%20and%20boxplot-1.png)<!-- -->

``` r
ggplot(data=tt3, mapping=aes(x=IC50_logFC,y=Avg_Exp)) +geom_point()
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Plot%20gene%20expression%20vs%20AUC%20and%20gene%20expression%20vs%20IC50,%20on%20a%20scatter%20plot%20and%20boxplot-2.png)<!-- -->

``` r
# Gene expression vs AUC
# Boxplot
ggplot(data=tt3, mapping=aes(x=AUC_Sign,y=Avg_Exp, group=AUC_Sign, fill=AUC_Sign)) +
  geom_boxplot()
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Plot%20gene%20expression%20vs%20AUC%20and%20gene%20expression%20vs%20IC50,%20on%20a%20scatter%20plot%20and%20boxplot-3.png)<!-- -->

``` r
# Scatter plot
ggplot(data=tt3, mapping=aes(x=AUC_logFC,y=Avg_Exp)) +geom_point()
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Plot%20gene%20expression%20vs%20AUC%20and%20gene%20expression%20vs%20IC50,%20on%20a%20scatter%20plot%20and%20boxplot-4.png)<!-- -->

Looks like either it’s all -logged, or the relationship between ic50 and
auc is not what we think it is.

``` r
# Load packages
library(ReactomePA)
```

    ## 

    ## Registered S3 method overwritten by 'enrichplot':
    ##   method               from
    ##   fortify.enrichResult DOSE

    ## ReactomePA v1.30.0  For help: https://guangchuangyu.github.io/ReactomePA
    ## 
    ## If you use ReactomePA in published research, please cite:
    ## Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular BioSystems 2016, 12(2):477-479

``` r
library(org.Hs.eg.db)

# Create list of top associated genes
Top400Genes <- as.character(na.omit(SigSamples$Gene_Symbol[1:400]))
# Add Entrez ID's
hs <- org.Hs.eg.db
Top400GenesEntrez <- AnnotationDbi::select(hs, 
                            keys = Top400Genes,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
# Perform pathway analysis via ReactomePA
# Gives you a nice table with some useful pathways
dasat_rPAoverrep <- enrichPathway(gene=Top400GenesEntrez$ENTREZID, organism = "human", pvalueCutoff = 10, readable=T)
View(as.data.frame(dasat_rPAoverrep))
# Visualise these as plots (note that they look weird in the inline output, have to be opened properly)
barplot(dasat_rPAoverrep, showCategory=11, font.size=8)
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/reactome%20PA%20Analysis%20from%20Miks%20Slides-1.png)<!-- -->

``` r
dotplot(dasat_rPAoverrep, showCategory=11, font.size=8)
```

    ## wrong orderBy parameter; set to default `orderBy = "x"`

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/reactome%20PA%20Analysis%20from%20Miks%20Slides-2.png)<!-- -->

# Find all the genes, note which ones are significant. tt3 is all genes, assign those that are significant as 1’s, everything else as 0

``` r
## Load packages
library(reactome.db)
library(goseq)
```

    ## Loading required package: BiasedUrn

    ## Loading required package: geneLenDataBase

    ## 
    ## Attaching package: 'geneLenDataBase'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     unfactor

``` r
library(gplots)
hs <- org.Hs.eg.db

## Play around with the parameters for what's significant
SigSamples <- filter(tt3, tt3$logFC_Sign == 1, abs(tt3$IC50_logFC) > log2(2), abs(tt3$AUC_logFC) > log2(2))
# Rank the samples based on average rank between the two p-values
SigSamples <- 
  SigSamples %>% mutate(., rank_ic50=rank(IC50_Adj_PVal)) %>%
  mutate(., rank_auc=rank(AUC_Adj_PVal)) %>% 
  mutate(., avg_rank=0.5*(rank_auc + rank_ic50)) %>% 
  arrange(., avg_rank)
SigSamples <- SigSamples[1:500,]
# Add Entrez ID's to the SigSamples samples. These will become the 1's.
TopXGenes <- as.character(na.omit(SigSamples$Gene_Symbol))
TopXGenesEntrez <- AnnotationDbi::select(hs, 
                            keys = TopXGenes,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")
```

    ## 'select()' returned many:1 mapping between keys and columns

``` r
# SAME # Do something similar to all the genes, from tt3. This will give you the 0's (those that aren't 1's).
tt3Genes <- as.character(na.omit(tt3$Gene_Symbol))
tt3GenesEntrez <- AnnotationDbi::select(hs, 
                            keys = tt3Genes,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")
```

    ## 'select()' returned many:many mapping between keys and columns

``` r
# SAME # Get pathway names, and extract human-only pathways
rName <- as.list(reactomePATHNAME2ID)
rName <- rName[grep("Homo sapiens", names(rName))]
rGenes <- as.list(reactomePATHID2EXTID)
rGenesPath <- rGenes[match(rName, names(rGenes))]
rGenesPath <- lapply(rGenesPath, unique)
rGeneByPath <- as.list(reactomeEXTID2PATHID)
allGenes <- intersect( tt3GenesEntrez$ENTREZID, unique(unlist(rGenesPath)) )
length(allGenes)
```

    ## [1] 10420

``` r
## sigGenesX is the intersection of all the significant genes from the dataset with the human reactome genes
sigGenes <- intersect( TopXGenesEntrez$ENTREZID, unique(unlist(rGenesPath)) )
length(sigGenes)
```

    ## [1] 184

``` r
## sigGenes is the binary list of whether or not the gene was significant.
## same length and order as allGenes
plotGenes <- rep(0, length(allGenes))
names(plotGenes) <- allGenes

# Now match the significant and total gene lists
plotGenes[match(sigGenes, names(plotGenes))] <- 1
table(plotGenes)
```

    ## plotGenes
    ##     0     1 
    ## 10236   184

``` r
# SAME ## Figure out which genes are relevant to your data set (some may not be present)
mt <- match(allGenes, names(rGeneByPath))
# SAME ## Count how many of the reactome genes aren't in your data set
sum(is.na(mt))
```

    ## [1] 0

``` r
# SAME ## Get genes and pathways relevant to our analysis
rGeneByPath <- lapply(rGeneByPath[mt], function(x) intersect(x, names(rGenesPath)))
head(rGeneByPath)
```

    ## $`83902`
    ## [1] "R-HSA-1266738" "R-HSA-6805567"
    ## 
    ## $`27044`
    ## [1] "R-HSA-1643685" "R-HSA-5663202" "R-HSA-6802952" "R-HSA-6802957"
    ## 
    ## $`133686`
    ## [1] "R-HSA-1430728" "R-HSA-196807"  "R-HSA-196849"  "R-HSA-196854" 
    ## 
    ## $`81567`
    ## [1] "R-HSA-168249"  "R-HSA-168256"  "R-HSA-199991"  "R-HSA-199992" 
    ## [5] "R-HSA-432720"  "R-HSA-432722"  "R-HSA-5653656" "R-HSA-6798695"
    ## 
    ## $`7352`
    ## [1] "R-HSA-1428517" "R-HSA-1430728" "R-HSA-163200"  "R-HSA-166187" 
    ## [5] "R-HSA-167826"  "R-HSA-167827" 
    ## 
    ## $`3615`
    ## [1] "R-HSA-1430728" "R-HSA-15869"   "R-HSA-168249"  "R-HSA-168256" 
    ## [5] "R-HSA-6798695" "R-HSA-73817"   "R-HSA-8956320"

``` r
## GOseq analysis
pwf <- nullp(plotGenes, 'hg19', id = "knownGene", plot.fit=TRUE)
```

    ## Loading hg19 length data...

    ## Warning in pcls(G): initial point very close to some inequality constraints

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/GoSeq%20Analysis%20for%20top%20X%20genes-1.png)<!-- -->

``` r
## Length-corrected
goseqReactome <- goseq(pwf, gene2cat = rGeneByPath)
```

    ## Using manually entered categories.

    ## Calculating the p-values...

``` r
## Non-length-corrected
hyperReactome <- goseq(pwf, gene2cat = rGeneByPath, method="Hypergeometric")
```

    ## Using manually entered categories.
    ## Calculating the p-values...

``` r
## Adjust p-values
goseqReactome$adjP <- p.adjust(goseqReactome$over_represented_pvalue, method="fdr")
hyperReactome$adjP <- p.adjust(hyperReactome$over_represented_pvalue, method="fdr")


venn(list(goseq=goseqReactome$category[goseqReactome$adjP <1],
          hyper=hyperReactome$category[hyperReactome$adjP <1]))
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/GoSeq%20Analysis%20for%20top%20X%20genes-2.png)<!-- -->

``` r
# Look at the genes in the venn diagram
goseqPathways <- filter(goseqReactome, goseqReactome$adjP <1)
View(goseqPathways)

# Add Pathway info
rPathName <- as.list(reactomePATHID2NAME)
goseqPathways$Pathway <- gsub("Homo sapiens: ", "", rPathName[match(goseqPathways$category, names(rPathName))])
```

# Extract genes from goseqPathways

``` r
# Make heatmaps for the top pathways and the genes involved from your data (e.g. 45 in ECM organisation)
# To do this:
# match entrez ID's from "goseqPathwaysX" stuff to ensembl ID's from expDat, to get the expression data
# extract the genes, their names, and match them to expDat.
# Then make a heatmap (STAT 435 slides)

# Create a dataframe with gene symbol, entrez ID, and ensembl ID, for easy matching of all genes
GeneLabelTool <- dplyr::pull(tt3, Gene_Symbol)
# Add Entrez ID's
GeneLabelTool <- AnnotationDbi::select(hs,
                                       keys = GeneLabelTool,
                                       columns = c("ENSEMBL", "ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL")
```

    ## 'select()' returned many:many mapping between keys and columns

``` r
# Create a list of all the significant genes for each goseq pathway
genesinPaths <-list()
for(i in 1:nrow(goseqPathways)){
  genesinPaths[[i]] <- rGenesPath[match(goseqPathways$category[i], names(rGenesPath))] %>% 
    .[[1]] %>% 
    intersect(., rownames(pwf)[pwf$DEgenes==1])
}

# convert first list element from entrez ids to symbols like this (as an example):
# GeneLabelTool$SYMBOL[na.omit(match(genesinPaths[[1]], GeneLabelTool$ENTREZID))]

# This converts *all* the list elements from entrez ID to symbols, and you can then concenate it seperately
symbolsinPaths <- lapply(genesinPaths, function(x) GeneLabelTool$SYMBOL[na.omit(match(x, GeneLabelTool$ENTREZID))] )
# Concenate these for ease of reading
genesStick <- lapply(symbolsinPaths, function(x) paste0(x, collapse="::", sep="")) %>% unlist()
# Add these genes as a column on goseqPathways
goseqPathways$DEgenesInCat <- genesStick
# Alter the "goseqPathways" data.frame so that no pathway names contain "/" (for file names later)
goseqPathways <- lapply(goseqPathways, gsub, pattern='/', replacement=' ') %>% as.data.frame()
```

``` r
## SETUP STUFF
# THE RAINBOW IS CURRENTLY CODED TO ALL TISSUES. WE WANT IT ONLY CODED TO 3.

# Look at the tissue names, colour coded the same as the rainbow in the heatmap. Create a makeshift legend

# SPARE I THINK
# tis <- lapply(strsplit(colnames(expDat_sort), "_"), function(x) paste0(x[-1], collapse="_", sep="")) %>%
#  unlist()
tis <- tissuetoolsort$Tissue_Type
cc <- as.factor(tis) %>% as.numeric() %>% rainbow(length(table(.)))[.]
tisTab <- table(tis)
plot(1,1,type='p',col="white")
legend('topleft', paste0(names(tisTab), " (", tisTab, ")"), fill=rainbow(length(table(tis))), cex=0.6)
```

![](Breast_Colorectal_Stomach_Heatmaps_06_01_20_files/figure-gfm/Find%20expression%20data%20and%20make%20a%20heatmap-1.png)<!-- -->

``` r
# Create an object that shows you which ensembl genes are in which paths
ensginPaths <- lapply(genesinPaths, function(x) GeneLabelTool$ENSEMBL[na.omit(match(x, GeneLabelTool$ENTREZID))] )

# Standardise rownames between expDat and GeneLabelTool
rownames(expDat_sort) <- strsplit(rownames(expDat_sort),".", fixed=T) %>% lapply(., function(x) x[1]) %>% unlist()

#SPARE I THINK
# tisx <- lapply(strsplit(colnames(expDat_sort), "_"), function(x) paste0(x[-1], collapse="_", sep="")) %>%
#   unlist()
# ccx <- as.factor(tisx) %>% as.numeric() %>% rainbow(length(table(.)))[.]
# table(tisx,ccx)
```

``` r
## The Version That Definitely Works
for(k in 1:length(ensginPaths)){
  pathway_expDat_sort <- subset(expDat_sort, rownames(expDat_sort) %in% ensginPaths[[k]])
  rownames(pathway_expDat_sort) <- match(rownames(pathway_expDat_sort), GeneLabelTool$ENSEMBL) %>% 
    GeneLabelTool$SYMBOL[.]
  zz <- apply(pathway_expDat_sort, 1, scale, scale=TRUE) %>% t()
  colnames(zz) <- colnames(pathway_expDat_sort)
  zz[zz > 3] <- 3
  zz[zz < -3] <- -3
  pn <- gsub(" ","_",goseqPathways$Pathway[k])
  pdf(paste0("BCG-", pn, "-heatmap.pdf"), height=12, width=20)
  heatmap.2(zz, trace='none', scale='none', col=bluered(50), ColSideColors=cc)
  dev.off()
}
```

# Trying to order the heatmap columns by IC50 log fold change

for(k in 1:length(ensginPaths)){ pathway\_expDat\_sort \<-
subset(expDat\_sort, rownames(expDat\_sort) %in% ensginPaths\[\[k\]\])
rownames(pathway\_expDat\_sort) \<-
match(rownames(pathway\_expDat\_sort),
GeneLabelTool\(ENSEMBL) %>%  GeneLabelTool\)SYMBOL\[.\] zz \<-
apply(pathway\_expDat\_sort, 1, scale, scale=TRUE) %\>% t() colnames(zz)
\<- colnames(pathway\_expDat\_sort) zz\[zz \> 3\] \<- 3 zz\[zz \< -3\]
\<- -3 pn \<- gsub("
“,”\_",goseqPathways\(Pathway[k])  pdf(paste0("BCG-", pn, "-heatmap.pdf"), height=12, width=20)  ord = order(SigSamples\)IC50\_logFC)
heatmap.2(zz\[ ,ord\], trace=‘none’, scale=‘none’, col=bluered(50),
ColSideColors=cc, Colv = FALSE) dev.off() }
