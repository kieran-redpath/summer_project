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

    ## -- Attaching packages ------------------------------------------------------------------------------------------------ tidyverse 1.2.1 --

    ## v tibble  2.1.3     v purrr   0.3.3
    ## v tidyr   1.0.0     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.4.0

    ## Warning: package 'tibble' was built under R version 3.5.3

    ## Warning: package 'tidyr' was built under R version 3.5.3

    ## Warning: package 'readr' was built under R version 3.5.3

    ## Warning: package 'purrr' was built under R version 3.5.3

    ## Warning: package 'stringr' was built under R version 3.5.3

    ## Warning: package 'forcats' was built under R version 3.5.3

    ## -- Conflicts --------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
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
View(CCLE_dasat_BCCCGC)
CDH1_Exp_Levels <- read.csv("CCLE_GDSC/Data/CDH1_Exp_Levels.csv")
View(CDH1_Exp_Levels)
```

\#27\_11\_19

``` r
# Remove the unneccesary "CDH1_Expression" column from "CCLE_dasat_BCCCGC"
CCLE_dasat_BCCCGC <- subset(CCLE_dasat_BCCCGC, select= -c(CDH1_Expression))
View(CCLE_dasat_BCCCGC)
# Merge "CCLE_dasat_BCCCGC" and "CDH1_Exp_Levels", to add expression data to "CCLE_dasat_BCCCGC"
CCLE_dasat_CDH1 = full_join(x=CCLE_dasat_BCCCGC, y=CDH1_Exp_Levels,  by=c("CCLE_Name" = "CCLE_Name"))
```

    ## Warning: Column `CCLE_Name` joining factors with different levels, coercing to
    ## character vector

``` r
CCLE_dasat_CDH1 <- CCLE_dasat_CDH1 %>%
  filter(`Primary.Disease`=="Gastric Cancer"|`Primary.Disease`=="Breast Cancer" |`Primary.Disease`=="Colon Cancer")
View(CCLE_dasat_CDH1)
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
  labs(x= "Level of CDH1 Expression", y="AUC", caption="(Data from GDSC and CCLE)") +
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
  labs(x= "Level of CDH1 Expression",y="IC50", caption="(Data from GDSC and CCLE)") +
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
  labs(x= "Number of CDH1 Mutations", y="Normalised CDH1 Expression", caption="(Data from GDSC and CCLE)") +
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

Also do expression vs mutation status (I.e does expression level
correlate with mutation)
