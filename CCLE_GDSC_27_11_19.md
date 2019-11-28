CCLE\_GDSC\_27\_11\_19
================
Kieran Redpath
27 November 2019

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

    ## -- Attaching packages --------------------------------------------------------------------------------------------------------- tidyverse 1.2.1 --

    ## v tibble  2.1.3     v purrr   0.3.3
    ## v tidyr   1.0.0     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.4.0

    ## Warning: package 'tibble' was built under R version 3.5.3

    ## Warning: package 'tidyr' was built under R version 3.5.3

    ## Warning: package 'readr' was built under R version 3.5.3

    ## Warning: package 'purrr' was built under R version 3.5.3

    ## Warning: package 'stringr' was built under R version 3.5.3

    ## Warning: package 'forcats' was built under R version 3.5.3

    ## -- Conflicts ------------------------------------------------------------------------------------------------------------ tidyverse_conflicts() --
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
```

    ## Warning: package 'ggbeeswarm' was built under R version 3.5.3

``` r
CCLE_dasat_BCCCGC <- read.csv("CCLE_GDSC/Data/CCLE_dasat_BCCCGC.csv")
CDH1_Exp_Levels <- read.csv("CCLE_GDSC/Data/CDH1_Exp_Levels.csv")
```

\#27\_11\_19

``` r
# Remove the unneccesary "CDH1_Expression" column from "CCLE_dasat_BCCCGC"
CCLE_dasat_BCCCGC <- subset(CCLE_dasat_BCCCGC, select= -c(CDH1_Expression))
# Merge "CCLE_dasat_BCCCGC" and "CDH1_Exp_Levels", to add expression data to "CCLE_dasat_BCCCGC"
CCLE_dasat_CDH1 = full_join(x=CCLE_dasat_BCCCGC, y=CDH1_Exp_Levels,  by=c("CCLE_Name" = "CCLE_Name"))
```

    ## Warning: Column `CCLE_Name` joining factors with different levels, coercing to
    ## character vector

``` r
CCLE_dasat_CDH1 <- CCLE_dasat_CDH1 %>%
  filter(`Primary.Disease`=="Gastric Cancer"|`Primary.Disease`=="Breast Cancer" |`Primary.Disease`=="Colon Cancer")
# Note that "CCLE_dasat_BCCCGC" has 202 entries, but we only have CDH1 expression data from ~108 of these (however, all drug treated samples are within this).
#Tidy up column names a bit
CCLE_dasat_CDH1 <- CCLE_dasat_CDH1 %>% 
  rename(
    Normalised_CDH1_Expression = ENSG00000039068.14,
    CDH1_Expression_Levels = CDH1_Expression,
    Sanger_ID = Sanger.ID,
    Primary_Disease = Primary.Disease,
    Subtype_Disease = Subtype.Disease
    )
View(CCLE_dasat_CDH1)
# Save the file as a .csv
write.csv(CCLE_dasat_CDH1,"~/Documents/R Data/summer_project/CCLE_GDSC/Data/CCLE_dasat_CDH1.csv", row.names = FALSE)
```

``` r
# Run this if  you don't already have "CCLE_dasat_CDH1" available:
# CCLE_dasat_CDH1 <- read.csv("~/Documents/R Data/summer_project/CCLE_GDSC/Data/CCLE_dasat_CDH1.csv")
# AUC ggplot (with CDH1 expression levels)
ggplot(data=CCLE_dasat_CDH1, mapping=aes(x=CDH1_Expression_Levels, y=AUC, fill=CDH1_Expression_Levels)) +
  labs(x= "Level of CDH1 Expression (High= top 60%, Low= bottom 40%)", y="AUC", caption="(Data from GDSC and CCLE)") +
  ggtitle("Dasatinib Data", subtitle="AUC vs Level of CDH1 Expression") +
  geom_boxplot() +
  geom_beeswarm()
```

    ## Warning: Removed 178 rows containing non-finite values (stat_boxplot).

    ## Warning: Removed 178 rows containing missing values (position_beeswarm).

![](CCLE_GDSC_27_11_19_files/figure-gfm/Plot%20AUC/IC50%20vs.%20CDH1%20Expression%20Levels-1.png)<!-- -->

``` r
# IC50 ggplot (with CDH1 expression levels)
ggplot(data=CCLE_dasat_CDH1, mapping =aes(x=CDH1_Expression_Levels, y=IC50, fill=CDH1_Expression_Levels)) +
  labs(x= "Level of CDH1 Expression (High= top 60%, Low= bottom 40%)",y="IC50", caption="(Data from GDSC and CCLE)") +
  ggtitle("Dasatinib Data", subtitle="IC50 vs Level of CDH1 Expresion") +
  geom_boxplot() +
  geom_beeswarm()
```

    ## Warning: Removed 178 rows containing non-finite values (stat_boxplot).
    
    ## Warning: Removed 178 rows containing missing values (position_beeswarm).

![](CCLE_GDSC_27_11_19_files/figure-gfm/Plot%20AUC/IC50%20vs.%20CDH1%20Expression%20Levels-2.png)<!-- -->

``` r
# Hmmm, doesn't look good.
```

``` r
ggplot(data=CCLE_dasat_CDH1, mapping=aes(x=CDH1_Mutations, y=Normalised_CDH1_Expression, fill=CDH1_Mutations)) +
  labs(x= "Number of CDH1 Mutations", y="Normalised CDH1 Expression (logCPM)", caption="(Data from GDSC and CCLE)") +
  ggtitle("CDH1 Data", subtitle="Normalised CDH1 Expression vs Number of CDH1 Mutations") +
  geom_boxplot() +
  geom_beeswarm()
```

    ## Warning: Removed 93 rows containing non-finite values (stat_boxplot).

    ## Warning: Removed 93 rows containing missing values (position_beeswarm).

![](CCLE_GDSC_27_11_19_files/figure-gfm/Plot%20Expression%20vs.%20Mutation%20Status,%20to%20Try%20Find%20out%20Where%20the%20Problem%20Lies-1.png)<!-- -->

``` r
#Hmmmmm... This still doesn't look good, does it?
```

# 28\_11\_19

## Look at the other dataset

``` r
# Original Set
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
#GDSC Drug Codes
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
#Other Data we Need Later (so process it now)
depMapCCLEmut <- fread('CCLE_GDSC/Data/CCLE_mutations.csv', sep=',')
depMapCCLElines <- fread('CCLE_GDSC/Data/DepMap-2018q3-celllines.csv', sep=',')
bc <- depMapCCLEmut$Tumor_Sample_Barcode
mt <- match(bc, depMapCCLElines$Broad_ID)
dis <- depMapCCLElines$`Primary Disease`[mt]
gn <- depMapCCLEmut$Hugo_Symbol
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
diseases = c("Colon Cancer", "Gastric Cancer", "Breast Cancer")
kp <- bc[which(gn == "CDH1" & dis %in% diseases)]
kk <- bc[which(gn != "CDH1" & dis %in% diseases)]
col_order_cdh1 <- c("CCLE_Name","Broad_ID","Aliases","COSMIC_ID","Sanger ID","Primary Disease","Subtype Disease","Gender","Source")

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
cdh1Lines <- cdh1Lines[, col_order_cdh1]
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
Noncdh1Lines <- cdh1Lines[, col_order_cdh1]

ccleExpData <- read.gct('CCLE_GDSC/Data/CCLE_RNAseq_genes_counts_20180929.gct')
logExpDat <- log(ccleExpData + 0.5)
dge <- DGEList(counts=ccleExpData)
dge <- calcNormFactors(dge)
v <- voom(dge, plot=TRUE)
```

![](CCLE_GDSC_27_11_19_files/figure-gfm/Load%20the%20GDSC%20data-1.png)<!-- -->

``` r
expDat <- v$E
```

``` r
# Load GDSC2 Data
GDSC2 <- fread('CCLE_GDSC/Data/GDSC2_fitted_dose_response_15Oct19.csv', sep= ';')
dim(GDSC2)
```

    ## [1] 118428     19

``` r
# Extract Dasatinib Data
dasatinib <- GDSC2 %>% filter(., DRUG_NAME=="Dasatinib")
# Fix the names of the expDat columns, so they're just cell line identifiers instead of full names
expDatsplit <- strsplit(colnames(expDat), "_") %>% 
  lapply(., function(x) x[1]) %>% 
  unlist()
colnames(expDat) <- expDatsplit
View(expDat)
# Match "dasatinib$CELL_LINE_NAME" formatting to expDat names
dasatinib$CELL_LINE_NAME <- gsub("-","",dasatinib$CELL_LINE_NAME, fixed=TRUE)
# Now keep sorting
commonSamples <- intersect(dasatinib$CELL_LINE_NAME,colnames(expDat))
# A better way to do it would be: dasatinib_sort <- dasatinib %>% filter(CELL_LINE_NAME %in% commonSamples)
expDat_match <- match(commonSamples, colnames(expDat))
expDat_sort <- expDat[ , na.omit(expDat_match)]
dasatinib_match <- match(commonSamples, dasatinib$CELL_LINE_NAME)
dasatinib_sort <- dasatinib[na.omit(dasatinib_match),]
```

``` r
# Wrangle CCLE Info for CDH1
CDH1_Exp_Levels <- subset(expDat_sort, rownames(expDat_sort)=="ENSG00000039068.14")
#Reorient "CDH1_Exp_Levels" so it has the same orientation as "CCLE_dasat_BCCCGC"
CDH1_Exp_Levels <- t(CDH1_Exp_Levels)
# Move row names into the first column, and name it "CCLE_Name"
CDH1_Exp_Levels <- as.data.frame(CDH1_Exp_Levels)
setDT(CDH1_Exp_Levels, keep.rownames = "CCLE_Name")[]
```

    ##      CCLE_Name ENSG00000039068.14
    ##   1:      SW48          8.4608406
    ##   2:    SW1710          0.8172700
    ##   3:     SW620          5.3111658
    ##   4:     CAL51          6.0239935
    ##   5:   TOV112D         -0.1911435
    ##  ---                             
    ## 468:     LS513          9.1008931
    ## 469:      TE15          9.0483323
    ## 470:    CORL95          7.0271447
    ## 471:   JURLMK1          3.9995901
    ## 472:   OCIAML5          0.1849553

``` r
View(CDH1_Exp_Levels)
# Filter to define "High" and "Low"
CDH1_Exp_Levels$CDH1_Expression <- ifelse(CDH1_Exp_Levels$ENSG00000039068.14 < quantile(CDH1_Exp_Levels$ENSG00000039068.14, 0.4), "Low","High")
# Merge with dasatinib sort
dasatinib_sort <- dasatinib_sort %>% 
  rename(
    CCLE_Name = CELL_LINE_NAME
    )
dasatinib_CDH1 = full_join(x=dasatinib_sort, y=CDH1_Exp_Levels,  by=c("CCLE_Name" = "CCLE_Name"))
View(dasatinib_CDH1)
```

``` r
# Plot everything
ggplot(data=dasatinib_CDH1, mapping=aes(x=CDH1_Expression, y=AUC, fill=CDH1_Expression)) +
  labs(x= "Level of CDH1 Expression (High= top 60%, Low= bottom 40%)", y="AUC", caption="(Data from GDSC2 and CCLE)") +
  ggtitle("Dasatinib Data (all tissues)", subtitle="AUC vs Level of CDH1 Expression") +
  geom_boxplot() +
  geom_beeswarm()
```

![](CCLE_GDSC_27_11_19_files/figure-gfm/Create%20Some%20Plots-1.png)<!-- -->

``` r
ggplot(data=dasatinib_CDH1, mapping=aes(x=CDH1_Expression, y=LN_IC50, fill=CDH1_Expression)) +
  labs(x= "Level of CDH1 Expression (High= top 60%, Low= bottom 40%)", y="ln IC50", caption="(Data from GDSC2 and CCLE)") +
  ggtitle("Dasatinib Data (all tissues)", subtitle="ln IC50 vs Level of CDH1 Expression") +
  geom_boxplot() +
  geom_beeswarm()
```

![](CCLE_GDSC_27_11_19_files/figure-gfm/Create%20Some%20Plots-2.png)<!-- -->

``` r
# Look at just the three main tissues
dasatinib_CDH1_3 <- dasatinib_CDH1 %>%
  filter(`TCGA_DESC`=="STAD"|`TCGA_DESC`=="BRCA" |`TCGA_DESC`=="COREAD")
View(dasatinib_CDH1_3)
# Plot these
ggplot(data=dasatinib_CDH1_3, mapping=aes(x=CDH1_Expression, y=AUC, fill=CDH1_Expression)) +
  labs(x= "Level of CDH1 Expression (High= top 60%, Low= bottom 40%)", y="AUC", caption="(Data from GDSC2 and CCLE)") +
  ggtitle("Dasatinib Data (all tissues)", subtitle="AUC vs Level of CDH1 Expression") +
  geom_boxplot() +
  geom_beeswarm()
```

![](CCLE_GDSC_27_11_19_files/figure-gfm/Plot%203%20Tissues-1.png)<!-- -->

``` r
ggplot(data=dasatinib_CDH1_3, mapping=aes(x=CDH1_Expression, y=LN_IC50, fill=CDH1_Expression)) +
  labs(x= "Level of CDH1 Expression (High= top 60%, Low= bottom 40%)", y="ln IC50", caption="(Data from GDSC2 and CCLE)") +
  ggtitle("Dasatinib Data (all tissues)", subtitle="ln IC50 vs Level of CDH1 Expression") +
  geom_boxplot() +
  geom_beeswarm()
```

![](CCLE_GDSC_27_11_19_files/figure-gfm/Plot%203%20Tissues-2.png)<!-- -->

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
