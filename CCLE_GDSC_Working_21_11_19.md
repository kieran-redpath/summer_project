CCLE\_GDSC\_Analysis\_Working\_21\_11\_19
================
Kieran Redpath
21 November 2019

## Shortcuts

Ctrl Shift M = %\>% Alt - = \<- Ctrl Alt I = Insert new R chunk Ctrl Alt
M = Git commit

## What exactly are these files?

DepMap-2018q3-celllines.csv: tells you more information about the cell
lines CCLE\_mutations.csv: information on what mutations are present in
CCLE cell lines CCLE\_RNAseq\_genes\_counts\_20180929.gct: expression
data for a bunch of genes in CCLE cell lines

gdsc\_codes\_1.csv: tells you about drugs (GDSC ID, actual name, etc.)
GDSC\_AUC.csv: AUC data for GDSC cell lines, treated with a bunch of
different drugs GDSC\_IC50.csv: IC50 data for GDSC cell lines, treated
with a bunch of different drugs

# 18\_11\_19

``` r
setwd("~/Documents/R Data/summer_project")
# Packages
library(limma)
library(edgeR)
library(data.table)
```

    ## Warning: package 'data.table' was built under R version 3.5.3

``` r
library(magrittr)
```

    ## Warning: package 'magrittr' was built under R version 3.5.3

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 3.5.3

``` r
library(CePa)
```

    ## Warning: package 'CePa' was built under R version 3.5.3

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 3.5.3

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     between, first, last

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 3.5.3

    ## -- Attaching packages -------------------------------------------------------------------------------- tidyverse 1.2.1 --

    ## v tibble  2.1.3     v purrr   0.3.3
    ## v tidyr   1.0.0     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.4.0

    ## Warning: package 'tibble' was built under R version 3.5.3

    ## Warning: package 'tidyr' was built under R version 3.5.3

    ## Warning: package 'readr' was built under R version 3.5.3

    ## Warning: package 'purrr' was built under R version 3.5.3

    ## Warning: package 'stringr' was built under R version 3.5.3

    ## Warning: package 'forcats' was built under R version 3.5.3

    ## -- Conflicts ----------------------------------------------------------------------------------- tidyverse_conflicts() --
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
# Read Data
depMapCCLEmut <- fread('CCLE_GDSC/Data/CCLE_mutations.csv', sep=',')
dim(depMapCCLEmut)
```

    ## [1] 1239235      34

``` r
depMapCCLElines <- fread('CCLE_GDSC/Data/DepMap-2018q3-celllines.csv', sep=',')
dim(depMapCCLElines)
```

    ## [1] 1673    9

``` r
# Check format of cell line info:
head(depMapCCLElines)
```

    ##      Broad_ID                                 CCLE_Name Aliases COSMIC_ID
    ## 1: ACH-000557 AML193_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE AML-193        NA
    ## 2: ACH-001000             1321N1_CENTRAL_NERVOUS_SYSTEM                NA
    ## 3: ACH-000198   EOL1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE    EOL1    906856
    ## 4: ACH-000956                            22RV1_PROSTATE  22-RV1    924100
    ## 5: ACH-000948                           2313287_STOMACH 2313287    910924
    ## 6: ACH-000011                        253J_URINARY_TRACT    253J        NA
    ##    Sanger ID Primary Disease         Subtype Disease Gender Source
    ## 1:        NA        Leukemia                         Female   ATCC
    ## 2:        NA    Brain Cancer                  glioma              
    ## 3:       426        Leukemia haematopoietic_neoplasm   Male   DSMZ
    ## 4:      1027 Prostate Cancer               carcinoma   Male   ATCC
    ## 5:       558  Gastric Cancer          adenocarcinoma   Male   DSMZ
    ## 6:        NA  Bladder Cancer               carcinoma          KCLB

``` r
# Get cell line barcodes:
bc <- depMapCCLEmut$Tumor_Sample_Barcode
head(bc)
```

    ## [1] "ACH-000001" "ACH-000001" "ACH-000001" "ACH-000001" "ACH-000001"
    ## [6] "ACH-000001"

``` r
# Match barcodes to Broad IDs
mt <- match(bc, depMapCCLElines$Broad_ID)
sum(is.na(mt))
```

    ## [1] 25914

``` r
# Cell lines that don't match to Broad IDs
is.na(mt) %>%  which() %>% bc[.] %>% table()
```

    ## .
    ## ACH-001335 ACH-001366 ACH-001403 ACH-001421 ACH-001422 ACH-001451 ACH-001484 
    ##        310        496        232        192        186        264       1073 
    ## ACH-001485 ACH-001495 ACH-001509 ACH-001513 ACH-001515 ACH-001520 ACH-001538 
    ##       1113        306        250        836        684        192        437 
    ## ACH-001540 ACH-001541 ACH-001548 ACH-001573 ACH-001607 ACH-001608 ACH-001609 
    ##        310        440        256        273        421        595        317 
    ## ACH-001618 ACH-001622 ACH-001628 ACH-001632 ACH-001634 ACH-001647 ACH-001651 
    ##        229        426        246        249        208        332        259 
    ## ACH-001664 ACH-001668 ACH-001677 ACH-001735 ACH-001736 ACH-001791 ACH-001793 
    ##       2416        271        245       1059        252        221        214 
    ## ACH-001794 ACH-001796 ACH-001814 ACH-001818 ACH-001836 ACH-001838 ACH-001839 
    ##        218        336        241        165        285        208        230 
    ## ACH-001841 ACH-001842 ACH-001843 ACH-001848 ACH-001849 ACH-001852 ACH-001853 
    ##        372        243        286        181        248        401       1443 
    ## ACH-001856 ACH-001857 ACH-001959 ACH-001960 ACH-001997 ACH-001999 ACH-002001 
    ##        309        451        295        366        224        170        682 
    ## ACH-002002 ACH-002003 ACH-002004 ACH-002005 
    ##        691        696        531        832

``` r
# Extract disease and gene info
dis <- depMapCCLElines$`Primary Disease`[mt]
gn <- depMapCCLEmut$Hugo_Symbol
```

``` r
# Number of CDH1 mutations
sum(gn == "CDH1")
```

    ## [1] 107

``` r
# Lines with CDH1 mutations (ordered by mutation counts)
which(gn == "CDH1") %>% dis[.] %>% table() %>% sort(., decreasing=TRUE) 
```

    ## .
    ##                    Colon Cancer                  Gastric Cancer 
    ##                              13                              13 
    ##                      Lung NSCLC                   Breast Cancer 
    ##                              11                               9 
    ##              Endometrial Cancer                        Leukemia 
    ##                               5                               5 
    ##                 Cervical Cancer                       Lung SCLC 
    ##                               3                               3 
    ##                        Melanoma                  Ovarian Cancer 
    ##                               3                               3 
    ##                 Prostate Cancer                  Thyroid Cancer 
    ##                               3                               3 
    ##                  Bladder Cancer                    Brain Cancer 
    ##                               2                               2 
    ##               Esophageal Cancer                   Head and Neck 
    ##                               2                               2 
    ##            Head and Neck Cancer                   Kidney Cancer 
    ##                               2                               2 
    ##                        Lymphoma                     Skin Cancer 
    ##                               2                               2 
    ##        T-Lymphoblastic Leukemia B-Cell Non-Hodgkin\\'s Lymphoma 
    ##                               2                               1 
    ##        B-Lymphoblastic Leukemia          Bile Duct/Liver Cancer 
    ##                               1                               1 
    ##                Breast Carcinoma                Kidney Carcinoma 
    ##                               1                               1 
    ##                     lung cancer                     Lung Cancer 
    ##                               1                               1 
    ##                    Mesothelioma                Multiple Myeloma 
    ##                               1                               1 
    ##               Pancreatic Cancer                        Rhabdoid 
    ##                               1                               1 
    ##       Small Cell Lung Carcinoma    Squamous Cell Lung Carcinoma 
    ##                               1                               1 
    ##                         unknown 
    ##                               1

``` r
# Ignore NSC Lung Cancer for now
diseases = c("Colon Cancer", "Gastric Cancer", "Breast Cancer")
kp <- bc[which(gn == "CDH1" & dis %in% diseases)]
# Colon, breast and gastric cell lines with CDH1 mutations
cdh1Lines <- depMapCCLElines %>% dplyr::filter(., Broad_ID %in% kp)
cdh1Lines %>% 
  dplyr::arrange(., `Primary Disease`) 
```

    ##      Broad_ID               CCLE_Name      Aliases COSMIC_ID Sanger ID
    ## 1  ACH-000536             BT20_BREAST        BT-20    906801      2264
    ## 2  ACH-000902           CAL148_BREAST      CAL-148    924106      1610
    ## 3  ACH-000783            CAMA1_BREAST       CAMA-1    946382       363
    ## 4  ACH-001065            EVSAT_BREAST                 906862      1863
    ## 5  ACH-000643            HDQP1_BREAST HDQP1_BREAST   1290922       743
    ## 6  ACH-000910         MDAMB453_BREAST   MDA-MB-453    908122       101
    ## 7  ACH-000568          UACC812_BREAST     UACC-812    910910      1754
    ## 8  ACH-000554          UACC893_BREAST     UACC-893    909778      1090
    ## 9  ACH-000828           ZR7530_BREAST     ZR-75-30    909907       429
    ## 10 ACH-000998     CW2_LARGE_INTESTINE          CW2    910554      1833
    ## 11 ACH-000982    GP2D_LARGE_INTESTINE         GP2d        NA        NA
    ## 12 ACH-001345    GP5D_LARGE_INTESTINE         GP5d    907291       741
    ## 13 ACH-001081 HCC2998_LARGE_INTESTINE      HCC2998    905971        43
    ## 14 ACH-000971  HCT116_LARGE_INTESTINE       HCT116    905936      1920
    ## 15 ACH-000986   HT115_LARGE_INTESTINE                 907289      2099
    ## 16 ACH-000999 SNU1040_LARGE_INTESTINE      SNU1040   1659823      2191
    ## 17 ACH-000955  SNU407_LARGE_INTESTINE                1660034      1907
    ## 18 ACH-000532   SNU61_LARGE_INTESTINE        SNU61   1660035       151
    ## 19 ACH-000967  SNUC2A_LARGE_INTESTINE      SNU-C2A        NA        NA
    ## 20 ACH-000948         2313287_STOMACH      2313287    910924       558
    ## 21 ACH-000880             AGS_STOMACH          AGS    906790      1011
    ## 22 ACH-000047            GCIY_STOMACH         GCIY    906869      1414
    ## 23 ACH-000919            IM95_STOMACH         IM95   1240155       405
    ## 24 ACH-000793         KATOIII_STOMACH     KATO III    907276       406
    ## 25 ACH-000507            KE39_STOMACH         KE39        NA        NA
    ## 26 ACH-000356           MKN45_STOMACH       MKN-45    925340      1542
    ## 27 ACH-000110      NCCSTCK140_STOMACH   NCCSTCK140        NA        NA
    ## 28 ACH-000674           NUGC4_STOMACH        NUGC4   1298357       941
    ## 29 ACH-000247           OCUM1_STOMACH        OCUM1   1298358      1712
    ## 30 ACH-000144        RERFGC1B_STOMACH     RERFGC1B   1240209       285
    ## 31 ACH-000303            SNU5_STOMACH         SNU5    908445        28
    ##    Primary Disease Subtype Disease Gender Source
    ## 1    Breast Cancer           Basal Female  NIBRI
    ## 2    Breast Cancer         Luminal Female   DSMZ
    ## 3    Breast Cancer         Luminal Female   ATCC
    ## 4    Breast Cancer                 Female       
    ## 5    Breast Cancer           Basal Female   DSMZ
    ## 6    Breast Cancer         Luminal Female  NIBRI
    ## 7    Breast Cancer         Luminal Female   ATCC
    ## 8    Breast Cancer         Luminal Female   ATCC
    ## 9    Breast Cancer         Luminal Female    GNF
    ## 10    Colon Cancer  adenocarcinoma Female  RIKEN
    ## 11    Colon Cancer  adenocarcinoma Female  ECACC
    ## 12    Colon Cancer  adenocarcinoma Female       
    ## 13    Colon Cancer                              
    ## 14    Colon Cancer  adenocarcinoma   Male  NIBRI
    ## 15    Colon Cancer  adenocarcinoma         ECACC
    ## 16    Colon Cancer  adenocarcinoma   Male   KCLB
    ## 17    Colon Cancer  adenocarcinoma   Male   KCLB
    ## 18    Colon Cancer  adenocarcinoma   Male   KCLB
    ## 19    Colon Cancer  adenocarcinoma Female  NIBRI
    ## 20  Gastric Cancer  adenocarcinoma   Male   DSMZ
    ## 21  Gastric Cancer  adenocarcinoma Female  NIBRI
    ## 22  Gastric Cancer  adenocarcinoma Female    GNF
    ## 23  Gastric Cancer  adenocarcinoma   Male  HSRRB
    ## 24  Gastric Cancer  adenocarcinoma   Male   ATCC
    ## 25  Gastric Cancer  adenocarcinoma   Male  RIKEN
    ## 26  Gastric Cancer  adenocarcinoma Female  HSRRB
    ## 27  Gastric Cancer                 Female  RIKEN
    ## 28  Gastric Cancer  adenocarcinoma Female    GNF
    ## 29  Gastric Cancer  adenocarcinoma Female  NIBRI
    ## 30  Gastric Cancer                   Male  NIBRI
    ## 31  Gastric Cancer  adenocarcinoma Female  NIBRI

``` r
col_order_cdh1 <- c("CCLE_Name","Broad_ID","Aliases","COSMIC_ID","Sanger ID","Primary Disease","Subtype Disease","Gender","Source")
cdh1Lines <- cdh1Lines[, col_order_cdh1]
```

``` r
# Load Data
gdsc_auc <- read.csv('CCLE_GDSC/Data/GDSC_AUC.csv', sep=',')
gdsc_ic50 <- read.csv('CCLE_GDSC/Data/GDSC_IC50.csv', sep=',')
dim(gdsc_auc)
```

    ## [1] 266 970

``` r
dim(gdsc_ic50)
```

    ## [1] 266 970

``` r
gdsc_auc[1:5,1:5]
```

    ##           X ACH.002137 ACH.000474 ACH.002089 ACH.000956
    ## 1    GDSC:1         NA         NA         NA         NA
    ## 2 GDSC:1001   0.817796   0.943611   0.971663   0.899492
    ## 3 GDSC:1004   0.681053   0.409378   0.581949   0.600888
    ## 4 GDSC:1005   0.956814   0.966637   0.792002   0.913204
    ## 5 GDSC:1006   0.973314   0.509397   0.537315   0.802140

``` r
gdsc_ic50[1:5,1:5]
```

    ##           X ACH.002137 ACH.000474 ACH.002089 ACH.000956
    ## 1    GDSC:1         NA         NA         NA         NA
    ## 2 GDSC:1001   7.258918   9.131374  10.182594   8.332992
    ## 3 GDSC:1004  -3.802467  -5.702659  -4.499864  -4.366928
    ## 4 GDSC:1005   4.146364   4.551663   1.772586   3.263056
    ## 5 GDSC:1006   3.171367  -2.014854  -1.818771   0.252367

``` r
# GDSC drug codes
drugCodes <- read.csv('CCLE_GDSC/Data/gdsc_codes_1.csv', sep= ',')
head(drugCodes)
```

    ##   drug_id  drug_name                   synonyms                  pathway_name
    ## 1    1050   ZM447439       ZM-447439, ZM 447439                       Mitosis
    ## 2    1058 Pictilisib GDC-0941, GDC0941, RG-7621           PI3K/MTOR signaling
    ## 3     199  Pazopanib                   Votrient                 RTK signaling
    ## 4     200 Dacinostat         NVP-LAQ824, LAQ824 Chromatin histone acetylation
    ## 5     219    AT-7519                     AT7519                    Cell cycle
    ## 6     255   CP724714                  CP-724714                 RTK signaling
    ##                        targets  pubchem
    ## 1                 AURKA, AURKB  9914412
    ## 2               PI3K (class 1) 17755052
    ## 3  CSF1R, KIT,  PDGFRA, PDGFRB 10113978
    ## 4                        HDAC1  6445533
    ## 5 CDK1, CDK2, CDK4, CDK6, CDK9 11338033
    ## 6                        ERBB2  9874913

``` r
# Find Dasatinib (drug_id is "51")
drugCodes %>% dplyr::filter(., drug_name == "Dasatinib")
```

    ##   drug_id drug_name                           synonyms  pathway_name
    ## 1      51 Dasatinib BMS-354825-03, BMS-354825, Sprycel RTK signaling
    ##                         targets pubchem
    ## 1 ABL, SRC, Ephrins, PDGFR, KIT 3062316

``` r
# Extract Dasatinib data
dasatinib_ic50 <- gdsc_ic50 %>% dplyr::filter(., X=="GDSC:51")
dasatinib_auc <- gdsc_auc %>% dplyr::filter(., X=="GDSC:51")
# Count non-NA values
sum(!is.na(dasatinib_ic50))
```

    ## [1] 386

``` r
sum(!is.na(dasatinib_auc))
```

    ## [1] 386

``` r
# GDSC cell lines with Dasatinib AUC data
names(dasatinib_auc)[!is.na(dasatinib_auc)][-1]
```

    ##   [1] "ACH.000070" "ACH.000137" "ACH.000008" "ACH.000740" "ACH.001442"
    ##   [6] "ACH.000697" "ACH.000555" "ACH.000157" "ACH.000429" "ACH.002209"
    ##  [11] "ACH.000105" "ACH.000269" "ACH.000838" "ACH.002210" "ACH.002212"
    ##  [16] "ACH.002213" "ACH.002090" "ACH.002214" "ACH.002215" "ACH.002216"
    ##  [21] "ACH.002217" "ACH.001016" "ACH.000927" "ACH.000432" "ACH.000009"
    ##  [26] "ACH.000440" "ACH.000902" "ACH.000264" "ACH.000464" "ACH.001738"
    ##  [31] "ACH.002218" "ACH.001023" "ACH.000136" "ACH.000249" "ACH.000641"
    ##  [36] "ACH.002219" "ACH.000803" "ACH.000864" "ACH.000401" "ACH.000644"
    ##  [41] "ACH.000257" "ACH.000508" "ACH.000743" "ACH.002097" "ACH.002098"
    ##  [46] "ACH.001049" "ACH.002221" "ACH.002222" "ACH.000998" "ACH.002224"
    ##  [51] "ACH.002225" "ACH.000055" "ACH.002226" "ACH.002229" "ACH.000786"
    ##  [56] "ACH.000334" "ACH.000233" "ACH.002232" "ACH.002100" "ACH.000530"
    ##  [61] "ACH.000703" "ACH.000981" "ACH.000056" "ACH.002235" "ACH.000258"
    ##  [66] "ACH.001064" "ACH.002236" "ACH.000305" "ACH.000038" "ACH.000821"
    ##  [71] "ACH.000706" "ACH.000295" "ACH.000198" "ACH.002103" "ACH.000906"
    ##  [76] "ACH.002104" "ACH.002105" "ACH.002106" "ACH.002107" "ACH.002108"
    ##  [81] "ACH.002109" "ACH.002237" "ACH.001065" "ACH.002110" "ACH.002111"
    ##  [86] "ACH.002112" "ACH.002113" "ACH.002114" "ACH.002115" "ACH.002117"
    ##  [91] "ACH.002118" "ACH.001498" "ACH.000162" "ACH.002120" "ACH.000738"
    ##  [96] "ACH.000047" "ACH.000081" "ACH.000756" "ACH.001344" "ACH.001716"
    ## [101] "ACH.000073" "ACH.002238" "ACH.002240" "ACH.002139" "ACH.002241"
    ## [106] "ACH.002242" "ACH.000111" "ACH.000196" "ACH.000691" "ACH.000755"
    ## [111] "ACH.000515" "ACH.002243" "ACH.000267" "ACH.000190" "ACH.000004"
    ## [116] "ACH.001519" "ACH.000061" "ACH.000002" "ACH.000861" "ACH.002245"
    ## [121] "ACH.002141" "ACH.000914" "ACH.000322" "ACH.000538" "ACH.000672"
    ## [126] "ACH.002247" "ACH.002248" "ACH.002143" "ACH.000569" "ACH.002144"
    ## [131] "ACH.002145" "ACH.001529" "ACH.002252" "ACH.000653" "ACH.000151"
    ## [136] "ACH.002253" "ACH.000106" "ACH.000346" "ACH.002148" "ACH.000551"
    ## [141] "ACH.000231" "ACH.002255" "ACH.000053" "ACH.000315" "ACH.002256"
    ## [146] "ACH.000193" "ACH.000263" "ACH.000983" "ACH.000101" "ACH.000259"
    ## [151] "ACH.000386" "ACH.002149" "ACH.002257" "ACH.000293" "ACH.000969"
    ## [156] "ACH.000815" "ACH.002258" "ACH.000380" "ACH.000622" "ACH.002259"
    ## [161] "ACH.001106" "ACH.000345" "ACH.000227" "ACH.002261" "ACH.000631"
    ## [166] "ACH.000074" "ACH.000524" "ACH.002262" "ACH.000702" "ACH.000183"
    ## [171] "ACH.000754" "ACH.000806" "ACH.000301" "ACH.001355" "ACH.002150"
    ## [176] "ACH.002151" "ACH.002152" "ACH.002153" "ACH.002154" "ACH.002265"
    ## [181] "ACH.002266" "ACH.002155" "ACH.002157" "ACH.002267" "ACH.000977"
    ## [186] "ACH.000104" "ACH.000750" "ACH.000204" "ACH.000252" "ACH.000501"
    ## [191] "ACH.000985" "ACH.000007" "ACH.000438" "ACH.000787" "ACH.002270"
    ## [196] "ACH.000935" "ACH.000439" "ACH.000072" "ACH.002271" "ACH.000078"
    ## [201] "ACH.002272" "ACH.002273" "ACH.002274" "ACH.002162" "ACH.002275"
    ## [206] "ACH.000362" "ACH.000745" "ACH.001127" "ACH.000006" "ACH.000319"
    ## [211] "ACH.002163" "ACH.001131" "ACH.000335" "ACH.000045" "ACH.002276"
    ## [216] "ACH.002164" "ACH.002165" "ACH.002166" "ACH.000938" "ACH.000944"
    ## [221] "ACH.000804" "ACH.002278" "ACH.002280" "ACH.002281" "ACH.002282"
    ## [226] "ACH.000294" "ACH.002283" "ACH.002284" "ACH.002285" "ACH.002277"
    ## [231] "ACH.000514" "ACH.000980" "ACH.002286" "ACH.002169" "ACH.000666"
    ## [236] "ACH.000327" "ACH.001591" "ACH.000830" "ACH.000015" "ACH.000766"
    ## [241] "ACH.000448" "ACH.000431" "ACH.001362" "ACH.000559" "ACH.000733"
    ## [246] "ACH.000894" "ACH.001136" "ACH.000586" "ACH.000729" "ACH.000587"
    ## [251] "ACH.000394" "ACH.000290" "ACH.000639" "ACH.001363" "ACH.001138"
    ## [256] "ACH.000525" "ACH.000399" "ACH.000367" "ACH.000900" "ACH.002172"
    ## [261] "ACH.001364" "ACH.000800" "ACH.000871" "ACH.000816" "ACH.000767"
    ## [266] "ACH.001599" "ACH.000358" "ACH.000491" "ACH.002174" "ACH.000403"
    ## [271] "ACH.002176" "ACH.000355" "ACH.001365" "ACH.002288" "ACH.001603"
    ## [276] "ACH.002289" "ACH.002290" "ACH.000200" "ACH.001605" "ACH.001606"
    ## [281] "ACH.000168" "ACH.002291" "ACH.000287" "ACH.000113" "ACH.000336"
    ## [286] "ACH.000065" "ACH.000124" "ACH.000751" "ACH.002179" "ACH.002180"
    ## [291] "ACH.000776" "ACH.000024" "ACH.000159" "ACH.000617" "ACH.000770"
    ## [296] "ACH.002293" "ACH.000937" "ACH.000218" "ACH.000320" "ACH.002299"
    ## [301] "ACH.000654" "ACH.002300" "ACH.000189" "ACH.000960" "ACH.002301"
    ## [306] "ACH.002194" "ACH.000943" "ACH.000371" "ACH.000965" "ACH.001639"
    ## [311] "ACH.001182" "ACH.000817" "ACH.000636" "ACH.002302" "ACH.000874"
    ## [316] "ACH.002195" "ACH.002196" "ACH.000254" "ACH.001642" "ACH.002198"
    ## [321] "ACH.000609" "ACH.000655" "ACH.000273" "ACH.000441" "ACH.000790"
    ## [326] "ACH.000112" "ACH.000099" "ACH.000748" "ACH.000087" "ACH.000145"
    ## [331] "ACH.000373" "ACH.000465" "ACH.001190" "ACH.000363" "ACH.000366"
    ## [336] "ACH.000341" "ACH.001193" "ACH.000939" "ACH.000504" "ACH.000932"
    ## [341] "ACH.000581" "ACH.000303" "ACH.000722" "ACH.001199" "ACH.001665"
    ## [346] "ACH.000567" "ACH.000664" "ACH.002307" "ACH.000365" "ACH.000611"
    ## [351] "ACH.000656" "ACH.000059" "ACH.001203" "ACH.000226" "ACH.000953"
    ## [356] "ACH.002309" "ACH.002310" "ACH.001400" "ACH.002200" "ACH.001274"
    ## [361] "ACH.000197" "ACH.000424" "ACH.000647" "ACH.000318" "ACH.000488"
    ## [366] "ACH.002311" "ACH.000353" "ACH.000772" "ACH.000408" "ACH.000605"
    ## [371] "ACH.000452" "ACH.000694" "ACH.002312" "ACH.001674" "ACH.000146"
    ## [376] "ACH.001208" "ACH.002314" "ACH.001680" "ACH.000075" "ACH.000579"
    ## [381] "ACH.001702" "ACH.000115" "ACH.000534" "ACH.001709" "ACH.002317"

``` r
head(dasatinib_ic50)[,1:5]
```

    ##         X ACH.002137 ACH.000474 ACH.002089 ACH.000956
    ## 1 GDSC:51         NA         NA         NA         NA

``` r
head(dasatinib_auc)[,1:5]
```

    ##         X ACH.002137 ACH.000474 ACH.002089 ACH.000956
    ## 1 GDSC:51         NA         NA         NA         NA

``` r
# Replace names with CCLE names
ic50_names <- match(names(dasatinib_ic50), 
                  gsub("-",".",depMapCCLElines$Broad_ID, fixed=TRUE))
names(dasatinib_ic50) <- depMapCCLElines$CCLE_Name[ic50_names]
head(dasatinib_ic50)[,1:5]
```

    ##        NA [MERGED_TO_ACH-000109]H3255_LUNG [MERGED_TO_ACH-000474]NCIH292_LUNG
    ## 1 GDSC:51                               NA                                 NA
    ##   201T_LUNG 22RV1_PROSTATE
    ## 1        NA             NA

``` r
auc_names <- match(names(dasatinib_auc), 
                  gsub("-",".",depMapCCLElines$Broad_ID, fixed=TRUE))
names(dasatinib_auc) <- depMapCCLElines$CCLE_Name[auc_names]
head(dasatinib_auc)[,1:5]
```

    ##        NA [MERGED_TO_ACH-000109]H3255_LUNG [MERGED_TO_ACH-000474]NCIH292_LUNG
    ## 1 GDSC:51                               NA                                 NA
    ##   201T_LUNG 22RV1_PROSTATE
    ## 1        NA             NA

``` r
# Load data, convert to table
ccleExpData <- read.gct('CCLE_GDSC/Data/CCLE_RNAseq_genes_counts_20180929.gct')
dim(ccleExpData)
```

    ## [1] 56202  1019

``` r
# Data are read counts and ensembl gene IDs
ccleExpData[1:5,1:5]
```

    ##                   X22RV1_PROSTATE X2313287_STOMACH X253JBV_URINARY_TRACT
    ## ENSG00000223972.4              12                8                     8
    ## ENSG00000227232.4            1340              821                   678
    ## ENSG00000243485.2               4                1                     3
    ## ENSG00000237613.2               6                3                     2
    ## ENSG00000268020.2               0                2                     3
    ##                   X253J_URINARY_TRACT X42MGBA_CENTRAL_NERVOUS_SYSTEM
    ## ENSG00000223972.4                   6                              2
    ## ENSG00000227232.4                 677                            884
    ## ENSG00000243485.2                   3                              2
    ## ENSG00000237613.2                   4                              1
    ## ENSG00000268020.2                   1                              1

``` r
# Log the data (add 0.5 to avoid log(0))
logExpDat <- log(ccleExpData + 0.5)
# Density plots to assess consistency-> data are consistent with one another -> will normalise this next
plot(density(logExpDat[,1]))
for(i in 2:ncol(logExpDat)) lines(density(logExpDat[,i]))
```

![](CCLE_GDSC_Working_21_11_19_files/figure-gfm/CCLE%20Expression%20Data%20Setup/Test-1.png)<!-- -->

``` r
# Normalise via limma::voom (see page 70 of Limma usersguide)
# Note - I'm not filtering the genes
dge <- DGEList(counts=ccleExpData)
dge <- calcNormFactors(dge)

v <- voom(dge, plot=TRUE)
```

![](CCLE_GDSC_Working_21_11_19_files/figure-gfm/CCLE%20Expression%20Data%20Normalised/Unfiltered-1.png)<!-- -->

``` r
# Extract the normalised expression data
# Note that the voom procedure produces logged data
expDat <- v$E

# Plot the post-normalised data -> the important part
# Note there is now less variation in the "bump" above 5 - instead the 
# variation has been "moved" to the very low abundance genes (below -5).
plot(density(expDat[,1]))
for(i in 2:ncol(expDat)) lines(density(expDat[,i]))
```

![](CCLE_GDSC_Working_21_11_19_files/figure-gfm/CCLE%20Expression%20Data%20Normalised/Unfiltered-2.png)<!-- -->

``` r
expDat[1:5,1:5]
```

    ##                   X22RV1_PROSTATE X2313287_STOMACH X253JBV_URINARY_TRACT
    ## ENSG00000223972.4       -3.610620        -3.683642             -3.887225
    ## ENSG00000227232.4        3.134080         2.911012              2.431517
    ## ENSG00000243485.2       -5.084551        -6.186143             -5.167333
    ## ENSG00000237613.2       -4.554036        -4.963750             -5.652760
    ## ENSG00000268020.2       -8.254476        -5.449177             -5.167333
    ##                   X253J_URINARY_TRACT X42MGBA_CENTRAL_NERVOUS_SYSTEM
    ## ENSG00000223972.4           -4.146037                      -5.401171
    ## ENSG00000227232.4            2.557601                       3.065620
    ## ENSG00000243485.2           -5.039122                      -5.401171
    ## ENSG00000237613.2           -4.676551                      -6.138136
    ## ENSG00000268020.2           -6.261514                      -6.138136

## Now you can start having a play with the expression data in combination with the Dasatinib AUC and IC50 information. Can also add in CDH1 mutation data as needed.

# 19\_11\_19

``` r
#Create a vector containing the shared sample names (the intersect of identical data)
commonSamples <- intersect(names(dasatinib_ic50),colnames(expDat))
head(commonSamples)
```

    ## [1] "A101D_SKIN"                  "A172_CENTRAL_NERVOUS_SYSTEM"
    ## [3] "A204_SOFT_TISSUE"            "A2058_SKIN"                 
    ## [5] "A253_SALIVARY_GLAND"         "A2780_OVARY"

``` r
# Figure out which IC50 observations are within commonSamples (match)
ic50_match <- match(commonSamples, names(dasatinib_ic50))
# Interpret this (e.g. the 18th entry in dasatinib_ic50 matches the first entry in commonSamples)
head(ic50_match)
```

    ## [1] 18 19 20 21 22 23

``` r
# Figure out which expDat observations are within commonSamples (match)
expDat_match <- match(commonSamples, colnames(expDat))
head(expDat_match)
```

    ## [1] 16 18 19 20 21 22

``` r
# Create a new IC50 variable containing samples from commonSamples (na.omit) and sorted in the same order (also due to na.omit)
dasatinib_ic50_sort <- dasatinib_ic50[na.omit(ic50_match)]
# Create a new expDat variable containing samples from commonSamples (na.omit) and sorted in the same order (also due to na.omit)
expDat_sort <- expDat[ , na.omit(expDat_match)]

# strsplit takes each "names" data, splitting it at each "_", creating a list.
# lapply takes each list entry from this, dropping the first value (in this case, cell-line), then pastes the rest back together with "_".
# unlist converts the data from a list back to a vector.
# Define "tissue_sort"
tissue_sort <- strsplit(names(dasatinib_ic50_sort), "_") %>% 
  lapply(., function(x) paste0(x[-1], collapse='_')) %>%  
  unlist()
# Check numbers for each tissue type
table(tissue_sort)
```

    ## tissue_sort
    ##                  AUTONOMIC_GANGLIA                      BILIARY_TRACT 
    ##                                 12                                  1 
    ##                               BONE                             BREAST 
    ##                                 13                                 45 
    ##             CENTRAL_NERVOUS_SYSTEM                             CERVIX 
    ##                                 34                                  2 
    ##                        ENDOMETRIUM HAEMATOPOIETIC_AND_LYMPHOID_TISSUE 
    ##                                 10                                106 
    ##                             KIDNEY                    LARGE_INTESTINE 
    ##                                 12                                 42 
    ##                              LIVER                               LUNG 
    ##                                 15                                133 
    ##                         OESOPHAGUS                              OVARY 
    ##                                 24                                 29 
    ##                           PANCREAS                             PLEURA 
    ##                                 28                                  6 
    ##                           PROSTATE                     SALIVARY_GLAND 
    ##                                  4                                  1 
    ##                               SKIN                    SMALL_INTESTINE 
    ##                                 33                                  1 
    ##                        SOFT_TISSUE                            STOMACH 
    ##                                 14                                 21 
    ##                            THYROID          UPPER_AERODIGESTIVE_TRACT 
    ##                                  6                                 15 
    ##                      URINARY_TRACT 
    ##                                 14

``` r
# Create a vector containing the shared sample names (the intersect of identical data)
commonSamples <- intersect(names(dasatinib_auc),colnames(expDat))
head(commonSamples)
```

    ## [1] "A101D_SKIN"                  "A172_CENTRAL_NERVOUS_SYSTEM"
    ## [3] "A204_SOFT_TISSUE"            "A2058_SKIN"                 
    ## [5] "A253_SALIVARY_GLAND"         "A2780_OVARY"

``` r
# Figure out which IC50 observations are within commonSamples (match)
auc_match <- match(commonSamples, names(dasatinib_auc))
# Interpret this 
head(auc_match)
```

    ## [1] 18 19 20 21 22 23

``` r
# Figure out which expDat observations are within commonSamples (match)
expDat_auc_match <- match(commonSamples, colnames(expDat))
head(expDat_match)
```

    ## [1] 16 18 19 20 21 22

``` r
# Create a new auc variable containing samples from commonSamples (na.omit) and sorted in the same order (also due to na.omit)
dasatinib_auc_sort <- dasatinib_auc[na.omit(auc_match)]
# Create a new expDat variable containing samples from commonSamples (na.omit) and sorted in the same order (also due to na.omit)
expDat_auc_sort <- expDat[ , na.omit(expDat_auc_match)]

# Note that the order is now the same
expDat_auc_sort[1:5,1:5]
```

    ##                   A101D_SKIN A172_CENTRAL_NERVOUS_SYSTEM A204_SOFT_TISSUE
    ## ENSG00000223972.4  -4.891718                   -4.635942        -5.287784
    ## ENSG00000227232.4   2.749732                    2.480690         2.186328
    ## ENSG00000243485.2  -4.361203                   -5.483939        -4.802357
    ## ENSG00000237613.2  -6.476680                   -4.635942        -6.024749
    ## ENSG00000268020.2  -6.476680                   -7.805867        -6.024749
    ##                   A2058_SKIN A253_SALIVARY_GLAND
    ## ENSG00000223972.4  -4.248102           -7.983744
    ## ENSG00000227232.4   2.699097            2.687797
    ## ENSG00000243485.2  -5.833065           -5.661816
    ## ENSG00000237613.2  -5.347638           -5.661816
    ## ENSG00000268020.2  -8.154993           -6.398781

``` r
dasatinib_auc_sort[1:5]
```

    ##   A101D_SKIN A172_CENTRAL_NERVOUS_SYSTEM A204_SOFT_TISSUE A2058_SKIN
    ## 1   0.731988                          NA               NA         NA
    ##   A253_SALIVARY_GLAND
    ## 1            0.909876

``` r
dasatinib_ic50_sort[1:5]
```

    ##   A101D_SKIN A172_CENTRAL_NERVOUS_SYSTEM A204_SOFT_TISSUE A2058_SKIN
    ## 1  -0.829648                          NA               NA         NA
    ##   A253_SALIVARY_GLAND
    ## 1            0.619424

# 20\_11\_19

## CCLE\_dasat\_BCCCGC is *all* BC, CC, and GC lines with GDSC data

``` r
# Create the filtered file, removing cell line info. Leaves you with IC50 for BC, CC, GC
tissue_sort_BCCCGC <- strsplit(names(dasatinib_ic50_sort), "_") %>% 
  lapply(., function(x) paste0(x[-1], collapse='_')) %>%  
  unlist()
dasat_ic50_BCCCGC <- dasatinib_ic50_sort[tissue_sort == "STOMACH" | tissue_sort == "BREAST" | tissue_sort == "LARGE_INTESTINE"]
# Create the filtered file, removing cell line info. Leaves you with AUC for BC, CC, GC
tissue_sort_BCCCGC2 <- strsplit(names(dasatinib_auc_sort), "_") %>% 
  lapply(., function(x) paste0(x[-1], collapse='_')) %>%  
  unlist()
dasat_auc_BCCCGC <- dasatinib_auc_sort[tissue_sort_BCCCGC2 == "STOMACH" | tissue_sort_BCCCGC2 == "BREAST" | tissue_sort_BCCCGC2 == "LARGE_INTESTINE"]
# Combine these two files
dasat_data_BCCCGC <- rbind(dasat_auc_BCCCGC, dasat_ic50_BCCCGC)
rownames(dasat_data_BCCCGC) <- c("AUC","IC50")
dasat_data_BCCCGC <- t(dasat_data_BCCCGC)

# Create a file with cell line info to combine with "dasat_data_BCCCGC"
CCLElinesBCCCGC <- depMapCCLElines %>%
  filter(`Primary Disease`=="Gastric Cancer"|`Primary Disease`=="Breast Cancer" |`Primary Disease`=="Colon Cancer")
# Rename CCLElinesBCCCGC rows
l <- as.list(CCLElinesBCCCGC$CCLE_Name)
rownames(CCLElinesBCCCGC) <- l
CCLElinesBCCCGC <- select (CCLElinesBCCCGC, -c(CCLE_Name))

#Combine dasatinib data and CCLElinesBCCCGC
CCLE_dasat_BCCCGC <- merge(CCLElinesBCCCGC,dasat_data_BCCCGC,by=0, all=T)
names(CCLE_dasat_BCCCGC)[names(CCLE_dasat_BCCCGC) == 'Row.names'] <- 'CCLE_Name'
View(CCLE_dasat_BCCCGC)
```

# 21\_11\_19

## Create a file with CDH1 non-mutants

``` r
# Number of lines without CDH1 mutations
sum(gn != "CDH1")
```

    ## [1] 1239128

``` r
# Lines without CDH1 mutations (ordered by mutation counts)
which(gn != "CDH1") %>% dis[.] %>% table() %>% sort(., descending = TRUE)
```

    ## .
    ##                      Lymphoma/Leukemia                          Primary Cells 
    ##                                     36                                     95 
    ##                               leukemia                    T-cell_ALL/Leukemia 
    ##                                    159                                    172 
    ##                        Prostate cancer              Esophageal Adenocarcinoma 
    ##                                    175                                    186 
    ##                     Epitheliod Sarcoma                         Retinoblastoma 
    ##                                    201                                    202 
    ##                      Colorectal Cancer                      Hodgkins Lymphoma 
    ##                                    207                                    208 
    ##             Neuroblastoma/Brain Cancer                        T-cell lymphoma 
    ##                                    208                                    208 
    ##                             fibroblast                         Rhabdoid tumor 
    ##                                    216                                    223 
    ##       Multiple Myeloma/ Gastric Cancer                   Melanoma/Lung Cancer 
    ##                                    226                                    228 
    ##          Medulloblastoma/ Brain Cancer                       Lymphoma_Hodgkin 
    ##                                    261                                    280 
    ##                       Rhabdoid/Sarcoma                        T-cell Lymphoma 
    ##                                    285                                    303 
    ##                Thyroid Gland Carcinoma                       Embryonal Cancer 
    ##                                    303                                    321 
    ##                           Immortalized                         Rhabdoid Tumor 
    ##                                    338                                    352 
    ##                   Melanoma/Skin Cancer             Chondrosarcoma/Bone Cancer 
    ##                                    361                                    381 
    ##         B-cell_lymphoma_other/Leukemia                    Glioma/Brain Cancer 
    ##                                    432                                    448 
    ##        T-Cell Non-Hodgkin\\'s Lymphoma                          neuroblastoma 
    ##                                    448                                    455 
    ##                    B-cell_ALL/Leukemia         B-cell_lymphoma_other/Lymphoma 
    ##                                    461                                    465 
    ##          Lung Mesothelioma/Lung Cancer           Chronic Myelogenous Leukemia 
    ##                                    557                                    600 
    ##           Squamous Cell Lung Carcinoma                               melanoma 
    ##                                    605                                    634 
    ##                     Gallbladder Cancer                                Sarcoma 
    ##                                    665                                    689 
    ##                     Cervical Carcinoma                  Melanoma/ Skin Cancer 
    ##                                    711                                    712 
    ##                           Osteosarcoma                      Bladder Carcinoma 
    ##                                    722                                    725 
    ##                         Ewings Sarcoma          Non-Small Cell Lung Carcinoma 
    ##                                    777                                    794 
    ##                           CML/Leukemia                Biliary Tract Carcinoma 
    ##                                    846                                    875 
    ##                        Medulloblastoma                    Upper Head and Neck 
    ##                                    879                                    939 
    ##                   Pancreatic Carcinoma                Soft Tissue/Bone Cancer 
    ##                                    960                                   1034 
    ##                    Lung Adenocarcinoma                          Non-Cancerous 
    ##                                   1085                                   1100 
    ##                       Low Grade Glioma                         Chondrosarcoma 
    ##                                   1106                                   1119 
    ##                       GBM/Brain Cancer                             Fibroblast 
    ##                                   1126                                   1333 
    ## Multiple Myeloma/Leukemia, Plasma Cell                        Choriocarcinoma 
    ##                                   1412                                   1450 
    ##                                Control                       Bile Duct Cancer 
    ##                                   1472                                   1626 
    ##                    Soft Tissue/Sarcoma                      Ovarian Carcinoma 
    ##                                   1737                                   1757 
    ##                 Other Blood Carcinomas     Esophageal Squamous Cell Carcinoma 
    ##                                   1888                                   2056 
    ##                      Gastric Carcinoma               Mesothelioma/Lung Cancer 
    ##                                   2121                                   2280 
    ##                     Soft Tissue Cancer                 Bile Duct/Liver Cancer 
    ##                                   2307                                   2361 
    ##                       Breast Carcinoma                           Mesothelioma 
    ##                                   2386                                   2452 
    ##                           Glioblastoma                      Leukemia/Lymphoma 
    ##                                   2516                                   3152 
    ##                   Colorectal Carcinoma            Soft Tissue/ Thyroid Cancer 
    ##                                   3208                                   3305 
    ##                            lung cancer                   Burkitt\\'s Lymphoma 
    ##                                   3326                                   3424 
    ##                               Rhabdoid                       Rhabdomyosarcoma 
    ##                                   3920                                   4229 
    ##                          Head and Neck                       Kidney Carcinoma 
    ##                                   4332                                   4453 
    ##                    Soft Tissue Sarcoma                Head and Neck Carcinoma 
    ##                                   4461                                   5295 
    ##                             Meningioma                      Ewing\\'s Sarcoma 
    ##                                   5630                                   6643 
    ##               B-Lymphoblastic Leukemia        B-Cell Non-Hodgkin\\'s Lymphoma 
    ##                                   7015                                   7374 
    ##              Small Cell Lung Carcinoma                 Acute Myeloid Leukemia 
    ##                                   7815                                   7823 
    ##               T-Lymphoblastic Leukemia                         Thyroid Cancer 
    ##                                   8390                                   8649 
    ##                                unknown                            Bone Cancer 
    ##                                   8991                                  11570 
    ##                 Other Solid Carcinomas                            Lung Cancer 
    ##                                  11608                                  11851 
    ##                          Kidney Cancer                        Cervical Cancer 
    ##                                  12335                                  12375 
    ##                               Melanoma                           Liver Cancer 
    ##                                  12686                                  12919 
    ##                          Neuroblastoma                       Multiple Myeloma 
    ##                                  14963                                  15218 
    ##                        Prostate Cancer                      Pancreatic Cancer 
    ##                                  16049                                  17721 
    ##                         Bladder Cancer                      Esophageal Cancer 
    ##                                  18888                                  21171 
    ##                           Brain Cancer                         Gastric Cancer 
    ##                                  29276                                  30956 
    ##                          Breast Cancer                              Lung SCLC 
    ##                                  33737                                  36860 
    ##                         Ovarian Cancer                            Skin Cancer 
    ##                                  38222                                  64162 
    ##                               Leukemia                   Head and Neck Cancer 
    ##                                  73802                                  80494 
    ##                               Lymphoma                             Lung NSCLC 
    ##                                  83140                                 102520 
    ##                     Endometrial Cancer                           Colon Cancer 
    ##                                 108683                                 171316

``` r
# Ignore NSC Lung Cancer for now
diseases = c("Colon Cancer", "Gastric Cancer", "Breast Cancer")
kk <- bc[which(gn != "CDH1" & dis %in% diseases)]
# Colon, breast and gastric cell lines without CDH1 mutations
Noncdh1Lines <- depMapCCLElines %>% dplyr::filter(., Broad_ID %in% kk)
Noncdh1Lines %>% 
  dplyr::arrange(., `Primary Disease`)
```

    ##       Broad_ID               CCLE_Name               Aliases COSMIC_ID
    ## 1   ACH-000248            AU565_BREAST                 AU565    910704
    ## 2   ACH-000536             BT20_BREAST                 BT-20    906801
    ## 3   ACH-000927            BT474_BREAST                BT-474    946359
    ## 4   ACH-000818            BT483_BREAST                BT-483    949093
    ## 5   ACH-000288            BT549_BREAST                BT-549    905951
    ## 6   ACH-000212           CAL120_BREAST               CAL-120    906826
    ## 7   ACH-000902           CAL148_BREAST               CAL-148    924106
    ## 8   ACH-000856            CAL51_BREAST                CAL-51    910927
    ## 9   ACH-000857           CAL851_BREAST              CAL-85-1    910852
    ## 10  ACH-000783            CAMA1_BREAST                CAMA-1    946382
    ## 11  ACH-001820          COLO824_BREAST              COLO-824        NA
    ## 12  ACH-000258           DU4475_BREAST                Du4475    906844
    ## 13  ACH-000330            EFM19_BREAST                EFM-19    906851
    ## 14  ACH-000117          EFM192A_BREAST        EFM192A_BREAST   1290798
    ## 15  ACH-001065            EVSAT_BREAST                          906862
    ## 16  ACH-000374          HCC1143_BREAST               HCC1143    749710
    ## 17  ACH-000111          HCC1187_BREAST               HCC1187    749711
    ## 18  ACH-000699          HCC1395_BREAST               HCC1395    749712
    ## 19  ACH-000277          HCC1419_BREAST               HCC1419    907045
    ## 20  ACH-000352          HCC1428_BREAST               HCC1428   1290905
    ## 21  ACH-000349          HCC1500_BREAST               HCC1500   1303900
    ## 22  ACH-000930          HCC1569_BREAST               HCC1569    907046
    ## 23  ACH-000196          HCC1599_BREAST               HCC1599    749713
    ## 24  ACH-000624          HCC1806_BREAST               HCC1806    907047
    ## 25  ACH-000223          HCC1937_BREAST               HCC1937    749714
    ## 26  ACH-000859          HCC1954_BREAST               HCC1954    749709
    ## 27  ACH-000725           HCC202_BREAST                HCC202   1290906
    ## 28  ACH-000691          HCC2157_BREAST               HCC2157    749715
    ## 29  ACH-000755          HCC2218_BREAST               HCC2218    749716
    ## 30  ACH-000276            HCC38_BREAST                 HCC38    749717
    ## 31  ACH-000668            HCC70_BREAST                 HCC70    907048
    ## 32  ACH-000643            HDQP1_BREAST          HDQP1_BREAST   1290922
    ## 33  ACH-000721            HMC18_BREAST                 HMC18        NA
    ## 34  ACH-000134           HS274T_BREAST                HS274T        NA
    ## 35  ACH-000184           HS343T_BREAST                              NA
    ## 36  ACH-000148           HS578T_BREAST               Hs 578T    905957
    ## 37  ACH-000540           HS606T_BREAST                              NA
    ## 38  ACH-000413           HS739T_BREAST                              NA
    ## 39  ACH-000230           HS742T_BREAST                              NA
    ## 40  ACH-000711            JIMT1_BREAST                 JIMT1   1298157
    ## 41  ACH-000028             KPL1_BREAST                 KPL-1        NA
    ## 42  ACH-000019             MCF7_BREAST                  MCF7    905946
    ## 43  ACH-000044       MDAMB134VI_BREAST         MDA-MB-134-VI        NA
    ## 44  ACH-000621         MDAMB157_BREAST            MDA-MB-157    925338
    ## 45  ACH-000759      MDAMB175VII_BREAST        MDA-MB-175-VII    908120
    ## 46  ACH-000768         MDAMB231_BREAST            MDA-MB-231    905960
    ## 47  ACH-001358         MDAMB330_BREAST            MDA-MB-330   1330941
    ## 48  ACH-000934         MDAMB361_BREAST            MDA-MB-361    908121
    ## 49  ACH-000876         MDAMB415_BREAST            MDA-MB-415    924240
    ## 50  ACH-000573         MDAMB436_BREAST            MDA-MB-436   1240172
    ## 51  ACH-000910         MDAMB453_BREAST            MDA-MB-453    908122
    ## 52  ACH-000849         MDAMB468_BREAST            MDA-MB-468    908123
    ## 53  ACH-001819           MFM223_BREAST               MFM-223        NA
    ## 54  ACH-000017            SKBR3_BREAST               SK-BR-3        NA
    ## 55  ACH-001390         SUM149PT_BREAST             SUM-149PT        NA
    ## 56  ACH-001391         SUM159PT_BREAST             SUM-159PT        NA
    ## 57  ACH-001392         SUM185PE_BREAST             SUM-185PE        NA
    ## 58  ACH-001394         SUM229PE_BREAST             SUM-229PE        NA
    ## 59  ACH-001396          SUM52PE_BREAST        SUM-52PE;SUM52        NA
    ## 60  ACH-000147             T47D_BREAST                 T-47D    905945
    ## 61  ACH-000568          UACC812_BREAST              UACC-812    910910
    ## 62  ACH-000554          UACC893_BREAST              UACC-893    909778
    ## 63  ACH-000097            ZR751_BREAST               ZR-75-1        NA
    ## 64  ACH-000828           ZR7530_BREAST              ZR-75-30    909907
    ## 65  ACH-001454     C10_LARGE_INTESTINE                   C10        NA
    ## 66  ACH-001456  C125PM_LARGE_INTESTINE                C125PM        NA
    ## 67  ACH-000009  C2BBE1_LARGE_INTESTINE                C2BBe1    910700
    ## 68  ACH-001458     C75_LARGE_INTESTINE                   C75        NA
    ## 69  ACH-001459     C80_LARGE_INTESTINE                   C80        NA
    ## 70  ACH-001460     C84_LARGE_INTESTINE                   C84        NA
    ## 71  ACH-001461     C99_LARGE_INTESTINE                   C99        NA
    ## 72  ACH-000963   CCK81_LARGE_INTESTINE                         1240123
    ## 73  ACH-000249    CL11_LARGE_INTESTINE                  CL11   1290769
    ## 74  ACH-000342    CL14_LARGE_INTESTINE                  CL14        NA
    ## 75  ACH-000895    CL34_LARGE_INTESTINE                         1290771
    ## 76  ACH-000798    CL40_LARGE_INTESTINE                  CL40   1240124
    ## 77  ACH-000253 COLO201_LARGE_INTESTINE               COLO201        NA
    ## 78  ACH-001039 COLO205_LARGE_INTESTINE              COLO 205    905961
    ## 79  ACH-000202 COLO320_LARGE_INTESTINE              COLO-320        NA
    ## 80  ACH-000350 COLO678_LARGE_INTESTINE              COLO-678    910689
    ## 81  ACH-000998     CW2_LARGE_INTESTINE                   CW2    910554
    ## 82  ACH-000982    GP2D_LARGE_INTESTINE                  GP2d        NA
    ## 83  ACH-001345    GP5D_LARGE_INTESTINE                  GP5d    907291
    ## 84  ACH-001081 HCC2998_LARGE_INTESTINE               HCC2998    905971
    ## 85  ACH-000467   HCC56_LARGE_INTESTINE                 HCC56   1290907
    ## 86  ACH-000971  HCT116_LARGE_INTESTINE                HCT116    905936
    ## 87  ACH-000997   HCT15_LARGE_INTESTINE                HCT-15    905937
    ## 88  ACH-001091   HRT18_LARGE_INTESTINE                              NA
    ## 89  ACH-000199  HS255T_LARGE_INTESTINE                              NA
    ## 90  ACH-000214  HS675T_LARGE_INTESTINE                              NA
    ## 91  ACH-000850  HS698T_LARGE_INTESTINE                              NA
    ## 92  ACH-000986   HT115_LARGE_INTESTINE                          907289
    ## 93  ACH-000552    HT29_LARGE_INTESTINE                 HT-29    905939
    ## 94  ACH-000926    HT55_LARGE_INTESTINE                  HT55    907287
    ## 95  ACH-000538  HUTU80_SMALL_INTESTINE                HuTu80    907073
    ## 96  ACH-000969    KM12_LARGE_INTESTINE                  KM12    905989
    ## 97  ACH-001546 LIM1215_LARGE_INTESTINE               LIM1215        NA
    ## 98  ACH-000950    LOVO_LARGE_INTESTINE                  LoVo    907790
    ## 99  ACH-000252  LS1034_LARGE_INTESTINE                          917486
    ## 100 ACH-000501   LS123_LARGE_INTESTINE                          907792
    ## 101 ACH-000957   LS180_LARGE_INTESTINE                 LS180    998189
    ## 102 ACH-000985  LS411N_LARGE_INTESTINE                LS411N    907794
    ## 103 ACH-000007   LS513_LARGE_INTESTINE                 LS513    907795
    ## 104 ACH-000935   MDST8_LARGE_INTESTINE                         1240173
    ## 105 ACH-000360 NCIH508_LARGE_INTESTINE              NCI-H508    908442
    ## 106 ACH-000491 NCIH716_LARGE_INTESTINE              NCI-H716    908458
    ## 107 ACH-000403 NCIH747_LARGE_INTESTINE              NCI-H747    908457
    ## 108 ACH-000296  OUMS23_LARGE_INTESTINE                              NA
    ## 109 ACH-000565    RCM1_LARGE_INTESTINE                  RCM1    909263
    ## 110 ACH-000943     RKO_LARGE_INTESTINE                   RKO    909698
    ## 111 ACH-000400   SKCO1_LARGE_INTESTINE               SK-CO-1    909718
    ## 112 ACH-000286 SNU1033_LARGE_INTESTINE               SNU1033        NA
    ## 113 ACH-000999 SNU1040_LARGE_INTESTINE               SNU1040   1659823
    ## 114 ACH-000412 SNU1197_LARGE_INTESTINE              SNU-1197        NA
    ## 115 ACH-000989  SNU175_LARGE_INTESTINE                SNU175   1659928
    ## 116 ACH-000708  SNU283_LARGE_INTESTINE                SNU283   1659929
    ## 117 ACH-000955  SNU407_LARGE_INTESTINE                         1660034
    ## 118 ACH-000683  SNU503_LARGE_INTESTINE                              NA
    ## 119 ACH-000532   SNU61_LARGE_INTESTINE                 SNU61   1660035
    ## 120 ACH-000991   SNU81_LARGE_INTESTINE                         1660036
    ## 121 ACH-000722   SNUC1_LARGE_INTESTINE                SNU-C1    910905
    ## 122 ACH-000967  SNUC2A_LARGE_INTESTINE               SNU-C2A        NA
    ## 123 ACH-001199  SNUC2B_LARGE_INTESTINE                SNUC2B    909740
    ## 124 ACH-000959   SNUC4_LARGE_INTESTINE                              NA
    ## 125 ACH-000970   SNUC5_LARGE_INTESTINE                 SNUC5   1674021
    ## 126 ACH-000489  SW1116_LARGE_INTESTINE                SW1116    909746
    ## 127 ACH-000236  SW1417_LARGE_INTESTINE                SW1417    909747
    ## 128 ACH-000470  SW1463_LARGE_INTESTINE                          909748
    ## 129 ACH-000820   SW403_LARGE_INTESTINE                 SW403        NA
    ## 130 ACH-000958    SW48_LARGE_INTESTINE                 SW 48    909751
    ## 131 ACH-000842   SW480_LARGE_INTESTINE                SW 480        NA
    ## 132 ACH-000651   SW620_LARGE_INTESTINE                SW 620    905962
    ## 133 ACH-001399             SW626_OVARY                SW 626    909753
    ## 134 ACH-000421   SW837_LARGE_INTESTINE                SW 837    909755
    ## 135 ACH-000680   SW948_LARGE_INTESTINE SW948_LARGE_INTESTINE    909757
    ## 136 ACH-000381     T84_LARGE_INTESTINE                          909761
    ## 137 ACH-002394     GEO_LARGE_INTESTINE                              NA
    ## 138 ACH-000948         2313287_STOMACH               2313287    910924
    ## 139 ACH-000880             AGS_STOMACH                   AGS    906790
    ## 140 ACH-000560           ECC10_STOMACH                 ECC10    906848
    ## 141 ACH-000225           ECC12_STOMACH                 ECC12    906849
    ## 142 ACH-000633            FU97_STOMACH                  FU97   1290806
    ## 143 ACH-000047            GCIY_STOMACH                  GCIY    906869
    ## 144 ACH-000746             GSS_STOMACH                   GSS        NA
    ## 145 ACH-000485             GSU_STOMACH                   GSU        NA
    ## 146 ACH-000847           HGC27_STOMACH                 HGC27    907055
    ## 147 ACH-000616          HS746T_STOMACH                HS746T   1240151
    ## 148 ACH-000239           HUG1N_STOMACH                Hu-G1N        NA
    ## 149 ACH-000919            IM95_STOMACH                  IM95   1240155
    ## 150 ACH-000793         KATOIII_STOMACH              KATO III    907276
    ## 151 ACH-000507            KE39_STOMACH                  KE39        NA
    ## 152 ACH-000255            LMSU_STOMACH                  LMSU        NA
    ## 153 ACH-000351            MKN1_STOMACH                  MKN1    908138
    ## 154 ACH-000356           MKN45_STOMACH                MKN-45    925340
    ## 155 ACH-000678            MKN7_STOMACH                  MKN7    924250
    ## 156 ACH-000758           MKN74_STOMACH                 MKN74        NA
    ## 157 ACH-000110      NCCSTCK140_STOMACH            NCCSTCK140        NA
    ## 158 ACH-000427          NCIN87_STOMACH               NCI-N87    908461
    ## 159 ACH-000761           NUGC2_STOMACH                 NUGC2        NA
    ## 160 ACH-000911           NUGC3_STOMACH                NUGC-3    908455
    ## 161 ACH-000674           NUGC4_STOMACH                 NUGC4   1298357
    ## 162 ACH-000247           OCUM1_STOMACH                 OCUM1   1298358
    ## 163 ACH-000144        RERFGC1B_STOMACH              RERFGC1B   1240209
    ## 164 ACH-000764          SH10TC_STOMACH                SH10TC        NA
    ## 165 ACH-001653           SKGT2_STOMACH               SK-GT-2   1503364
    ## 166 ACH-000932            SNU1_STOMACH          SNU1_STOMACH    908444
    ## 167 ACH-000581           SNU16_STOMACH                 SNU16    908446
    ## 168 ACH-000466          SNU216_STOMACH                SNU216        NA
    ## 169 ACH-000303            SNU5_STOMACH                  SNU5    908445
    ## 170 ACH-000908          SNU520_STOMACH                SNU520        NA
    ## 171 ACH-000736          SNU601_STOMACH                SNU601        NA
    ## 172 ACH-000325          SNU620_STOMACH                SNU620        NA
    ## 173 ACH-000344          SNU668_STOMACH                SNU668        NA
    ## 174 ACH-000898          SNU719_STOMACH                SNU719        NA
    ## 175 ACH-000949       TGBC11TKB_STOMACH             TGBC11TKB    909770
    ##     Sanger ID Primary Disease          Subtype Disease Gender Source
    ## 1         726   Breast Cancer                  Luminal Female   ATCC
    ## 2        2264   Breast Cancer                    Basal Female  NIBRI
    ## 3         412   Breast Cancer                  Luminal Female  NIBRI
    ## 4        1231   Breast Cancer                  Luminal Female   ATCC
    ## 5        1835   Breast Cancer                    Basal Female   ATCC
    ## 6        1748   Breast Cancer                  Luminal Female   DSMZ
    ## 7        1610   Breast Cancer                  Luminal Female   DSMZ
    ## 8           4   Breast Cancer                    Basal Female   DSMZ
    ## 9         113   Breast Cancer                    Basal Female   DSMZ
    ## 10        363   Breast Cancer                  Luminal Female   ATCC
    ## 11         NA   Breast Cancer                                       
    ## 12        302   Breast Cancer                    Basal Female   ATCC
    ## 13       1786   Breast Cancer                  Luminal Female   DSMZ
    ## 14       1187   Breast Cancer                          Female  NIBRI
    ## 15       1863   Breast Cancer                          Female       
    ## 16       1316   Breast Cancer                    Basal Female   ATCC
    ## 17       1644   Breast Cancer                    Basal Female   ATCC
    ## 18         94   Breast Cancer                    Basal Female   ATCC
    ## 19       1551   Breast Cancer                  Luminal Female    GNF
    ## 20       1567   Breast Cancer                  Luminal Female   ATCC
    ## 21        653   Breast Cancer            Luminal/Basal Female   ATCC
    ## 22       2127   Breast Cancer                    Basal Female       
    ## 23       1264   Breast Cancer                    Basal Female   ATCC
    ## 24       1317   Breast Cancer                    Basal Female  NIBRI
    ## 25        936   Breast Cancer                    Basal Female   ATCC
    ## 26       2143   Breast Cancer                    Basal Female       
    ## 27       1552   Breast Cancer                  Luminal Female   ATCC
    ## 28        351   Breast Cancer                    Basal Female  NIBRI
    ## 29        352   Breast Cancer                  Luminal Female   ATCC
    ## 30       1318   Breast Cancer                    Basal Female   ATCC
    ## 31        937   Breast Cancer                    Basal Female  NIBRI
    ## 32        743   Breast Cancer                    Basal Female   DSMZ
    ## 33         NA   Breast Cancer                          Female  HSSRB
    ## 34         NA   Breast Cancer                          Female   ATCC
    ## 35         NA   Breast Cancer                          Female    GNF
    ## 36        655   Breast Cancer                    Basal Female   ATCC
    ## 37         NA   Breast Cancer                          Female   ATCC
    ## 38         NA   Breast Cancer                          Female       
    ## 39         NA   Breast Cancer                          Female   ATCC
    ## 40        479   Breast Cancer                          Female   DSMZ
    ## 41         NA   Breast Cancer                  Luminal Female   DSMZ
    ## 42        588   Breast Cancer                  Luminal Female   ATCC
    ## 43         NA   Breast Cancer                  Luminal Female       
    ## 44       1796   Breast Cancer                    Basal Female   ATCC
    ## 45         50   Breast Cancer                  Luminal Female   ATCC
    ## 46       1013   Breast Cancer                    Basal Female   ATCC
    ## 47       1989   Breast Cancer                          Female       
    ## 48        880   Breast Cancer                          Female   ATCC
    ## 49        344   Breast Cancer                  Luminal Female   ATCC
    ## 50        927   Breast Cancer                    Basal Female   ATCC
    ## 51        101   Breast Cancer                  Luminal Female  NIBRI
    ## 52        415   Breast Cancer                    Basal Female   ATCC
    ## 53         NA   Breast Cancer Ductal Mammary Carcinoma              
    ## 54         NA   Breast Cancer                  Luminal Female  NIBRI
    ## 55         NA   Breast Cancer                    Basal Female       
    ## 56         NA   Breast Cancer                    Basal Female       
    ## 57         NA   Breast Cancer                  Luminal Female       
    ## 58         NA   Breast Cancer                    Basal Female       
    ## 59         NA   Breast Cancer                  Luminal Female       
    ## 60       1286   Breast Cancer                  Luminal Female   ATCC
    ## 61       1754   Breast Cancer                  Luminal Female   ATCC
    ## 62       1090   Breast Cancer                  Luminal Female   ATCC
    ## 63         NA   Breast Cancer                  Luminal Female   ATCC
    ## 64        429   Breast Cancer                  Luminal Female    GNF
    ## 65         NA    Colon Cancer           adenocarcinoma   Male       
    ## 66         NA    Colon Cancer           adenocarcinoma              
    ## 67       2104    Colon Cancer           adenocarcinoma   Male   ATCC
    ## 68         NA    Colon Cancer           adenocarcinoma   Male       
    ## 69         NA    Colon Cancer           adenocarcinoma   Male       
    ## 70         NA    Colon Cancer           adenocarcinoma   Male       
    ## 71         NA    Colon Cancer           adenocarcinoma   Male       
    ## 72       1915    Colon Cancer           adenocarcinoma              
    ## 73       1153    Colon Cancer           adenocarcinoma          DSMZ
    ## 74         NA    Colon Cancer           adenocarcinoma   Male   DSMZ
    ## 75       1005    Colon Cancer           adenocarcinoma Female   DSMZ
    ## 76        555    Colon Cancer           adenocarcinoma Female   DSMZ
    ## 77         NA    Colon Cancer           adenocarcinoma   Male       
    ## 78       1687    Colon Cancer           adenocarcinoma   Male       
    ## 79         NA    Colon Cancer           adenocarcinoma Female   DSMZ
    ## 80       2148    Colon Cancer           adenocarcinoma   Male   DSMZ
    ## 81       1833    Colon Cancer           adenocarcinoma Female  RIKEN
    ## 82         NA    Colon Cancer           adenocarcinoma Female  ECACC
    ## 83        741    Colon Cancer           adenocarcinoma Female       
    ## 84         43    Colon Cancer                                       
    ## 85       1036    Colon Cancer           adenocarcinoma         NIBRI
    ## 86       1920    Colon Cancer           adenocarcinoma   Male  NIBRI
    ## 87       1553    Colon Cancer           adenocarcinoma   Male  NIBRI
    ## 88         NA    Colon Cancer           adenocarcinoma   Male       
    ## 89         NA    Colon Cancer           adenocarcinoma Female  NIBRI
    ## 90         NA    Colon Cancer           adenocarcinoma   Male  NIBRI
    ## 91         NA    Colon Cancer           adenocarcinoma   Male   ATCC
    ## 92       2099    Colon Cancer           adenocarcinoma         ECACC
    ## 93       1212    Colon Cancer           adenocarcinoma Female  NIBRI
    ## 94       1688    Colon Cancer           adenocarcinoma         ECACC
    ## 95        305    Colon Cancer  duodenal_adenocarcinoma   Male  NIBRI
    ## 96       1990    Colon Cancer           adenocarcinoma         NIBRI
    ## 97         NA    Colon Cancer                                       
    ## 98        306    Colon Cancer           adenocarcinoma   Male  NIBRI
    ## 99       1756    Colon Cancer           adenocarcinoma   Male   ATCC
    ## 100      2266    Colon Cancer           adenocarcinoma Female       
    ## 101      2081    Colon Cancer           adenocarcinoma Female       
    ## 102       724    Colon Cancer           adenocarcinoma   Male   ATCC
    ## 103       569    Colon Cancer           adenocarcinoma   Male       
    ## 104       928    Colon Cancer           adenocarcinoma         ECACC
    ## 105        35    Colon Cancer           adenocarcinoma   Male  NIBRI
    ## 106       649    Colon Cancer           adenocarcinoma   Male  NIBRI
    ## 107       573    Colon Cancer           adenocarcinoma   Male  NIBRI
    ## 108        NA    Colon Cancer           adenocarcinoma         HSSRB
    ## 109      1271    Colon Cancer           adenocarcinoma Female  HSSRB
    ## 110       607    Colon Cancer           adenocarcinoma          ATCC
    ## 111        40    Colon Cancer           adenocarcinoma   Male   ATCC
    ## 112        NA    Colon Cancer           adenocarcinoma Female   KCLB
    ## 113      2191    Colon Cancer           adenocarcinoma   Male   KCLB
    ## 114        NA    Colon Cancer           adenocarcinoma   Male   KCLB
    ## 115       154    Colon Cancer           adenocarcinoma Female   KCLB
    ## 116      2189    Colon Cancer           adenocarcinoma Female   KCLB
    ## 117      1907    Colon Cancer           adenocarcinoma   Male   KCLB
    ## 118        NA    Colon Cancer           adenocarcinoma   Male   KCLB
    ## 119       151    Colon Cancer           adenocarcinoma   Male   KCLB
    ## 120       697    Colon Cancer           adenocarcinoma   Male   KCLB
    ## 121      1066    Colon Cancer           adenocarcinoma   Male  NIBRI
    ## 122        NA    Colon Cancer           adenocarcinoma Female  NIBRI
    ## 123        37    Colon Cancer           adenocarcinoma Female       
    ## 124        NA    Colon Cancer           adenocarcinoma   Male   KCLB
    ## 125       152    Colon Cancer           adenocarcinoma Female   KCLB
    ## 126       480    Colon Cancer           adenocarcinoma   Male   ATCC
    ## 127       974    Colon Cancer           adenocarcinoma Female  NIBRI
    ## 128       271    Colon Cancer           adenocarcinoma Female   ATCC
    ## 129        NA    Colon Cancer           adenocarcinoma Female  NIBRI
    ## 130         1    Colon Cancer           adenocarcinoma Female  NIBRI
    ## 131        NA    Colon Cancer           adenocarcinoma   Male   ATCC
    ## 132         3    Colon Cancer           adenocarcinoma   Male   ATCC
    ## 133      2038    Colon Cancer                          Female       
    ## 134      1141    Colon Cancer           adenocarcinoma   Male   ATCC
    ## 135      1689    Colon Cancer           adenocarcinoma Female  NIBRI
    ## 136       997    Colon Cancer           adenocarcinoma   Male  NIBRI
    ## 137        NA    Colon Cancer                                 Sanger
    ## 138       558  Gastric Cancer           adenocarcinoma   Male   DSMZ
    ## 139      1011  Gastric Cancer           adenocarcinoma Female  NIBRI
    ## 140       799  Gastric Cancer          carcinoid_tumor         RIKEN
    ## 141      2094  Gastric Cancer          carcinoid_tumor   Male  Riken
    ## 142      2062  Gastric Cancer           adenocarcinoma Female  HSRRB
    ## 143      1414  Gastric Cancer           adenocarcinoma Female    GNF
    ## 144        NA  Gastric Cancer           adenocarcinoma         RIKEN
    ## 145        NA  Gastric Cancer           adenocarcinoma         RIKEN
    ## 146      1102  Gastric Cancer                                    GNF
    ## 147      1850  Gastric Cancer           adenocarcinoma   Male  NIBRI
    ## 148        NA  Gastric Cancer           adenocarcinoma   Male    GNF
    ## 149       405  Gastric Cancer           adenocarcinoma   Male  HSRRB
    ## 150       406  Gastric Cancer           adenocarcinoma   Male   ATCC
    ## 151        NA  Gastric Cancer           adenocarcinoma   Male  RIKEN
    ## 152        NA  Gastric Cancer                                  RIKEN
    ## 153      1647  Gastric Cancer            adenosquamous   Male  HSSRB
    ## 154      1542  Gastric Cancer           adenocarcinoma Female  HSRRB
    ## 155      2113  Gastric Cancer           adenocarcinoma   Male    GNF
    ## 156        NA  Gastric Cancer           adenocarcinoma   Male  HSSRB
    ## 157        NA  Gastric Cancer                          Female  RIKEN
    ## 158        88  Gastric Cancer           adenocarcinoma   Male       
    ## 159        NA  Gastric Cancer           adenocarcinoma Female  NIBRI
    ## 160      1251  Gastric Cancer           adenocarcinoma   Male   JHSF
    ## 161       941  Gastric Cancer           adenocarcinoma Female    GNF
    ## 162      1712  Gastric Cancer           adenocarcinoma Female  NIBRI
    ## 163       285  Gastric Cancer                            Male  NIBRI
    ## 164        NA  Gastric Cancer           adenocarcinoma         RIKEN
    ## 165      1868  Gastric Cancer           adenocarcinoma   Male       
    ## 166       564  Gastric Cancer                            Male  NIBRI
    ## 167       291  Gastric Cancer                          Female  NIBRI
    ## 168        NA  Gastric Cancer           adenocarcinoma Female   KCLB
    ## 169        28  Gastric Cancer           adenocarcinoma Female  NIBRI
    ## 170        NA  Gastric Cancer           adenocarcinoma Female   KCLB
    ## 171        NA  Gastric Cancer           adenocarcinoma   Male   KCLB
    ## 172        NA  Gastric Cancer           adenocarcinoma Female   KCLB
    ## 173        NA  Gastric Cancer           adenocarcinoma   Male   KCLB
    ## 174        NA  Gastric Cancer           adenocarcinoma   Male   KCLB
    ## 175      1917  Gastric Cancer           adenocarcinoma Female  RIKEN

``` r
col_order_Noncdh1 <- c("CCLE_Name","Broad_ID","Aliases","COSMIC_ID","Sanger ID","Primary Disease","Subtype Disease","Gender","Source")
Noncdh1Lines <- Noncdh1Lines[, col_order_Noncdh1]

# Open the new data
View(Noncdh1Lines)
View(cdh1Lines)
```

# 25\_11\_19

## Add CDH1 mutation data to the existing table, then create a boxplot with CDH1 mutants and non-CDH1 mutants next to each other, telling you about IC50 values.

``` r
CCLE_dasat_BCCCGC$CDH1_Mutations=FALSE
# Add AUC, IC50, and CDH1_Mutations columns to CCLE_dasat_BCCCGC
cdh1Lines$AUC=NA
cdh1Lines$IC50=NA
cdh1Lines$CDH1_Mutations=TRUE
# Use match to combine the two
Boxplot_data <- match(cdh1Lines$CCLE_Name, CCLE_dasat_BCCCGC$CCLE_Name)
CCLE_dasat_BCCCGC$CDH1_Mutations[Boxplot_data] <-  TRUE

ggplot(data = CCLE_dasat_BCCCGC, mapping =aes(x=CDH1_Mutations, y=AUC)) +
  geom_boxplot()
```

    ## Warning: Removed 178 rows containing non-finite values (stat_boxplot).

![](CCLE_GDSC_Working_21_11_19_files/figure-gfm/AUC%20and%20IC50%20Plots-1.png)<!-- -->

``` r
ggplot(data = CCLE_dasat_BCCCGC, mapping =aes(x=CDH1_Mutations, y=IC50)) +
  geom_boxplot()
```

    ## Warning: Removed 178 rows containing non-finite values (stat_boxplot).

![](CCLE_GDSC_Working_21_11_19_files/figure-gfm/AUC%20and%20IC50%20Plots-2.png)<!-- -->
