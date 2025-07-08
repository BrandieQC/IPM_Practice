---
title: "IPM_from_PopSim"
author: "Brandie QC"
date: "2025-07-08"
output: 
  html_document: 
    keep_md: true
---



# IPM from Population Simulation

\*Simulations done in "Population_Sim_Trial.Rmd" file (currently using AlphaSimR)

**Goal:** Compare IPM created from just phenotypes from the simulation to an IPM created from phenotypes + genotype scores.

## Libraries


``` r
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
## ✔ purrr     1.0.2     
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(broom)
library(broom.mixed)
library(magrittr)
```

```
## 
## Attaching package: 'magrittr'
## 
## The following object is masked from 'package:purrr':
## 
##     set_names
## 
## The following object is masked from 'package:tidyr':
## 
##     extract
```

``` r
library(lmerTest)
```

```
## Loading required package: lme4
## Loading required package: Matrix
## 
## Attaching package: 'Matrix'
## 
## The following objects are masked from 'package:tidyr':
## 
##     expand, pack, unpack
## 
## 
## Attaching package: 'lmerTest'
## 
## The following object is masked from 'package:lme4':
## 
##     lmer
## 
## The following object is masked from 'package:stats':
## 
##     step
```

``` r
conflicted::conflicts_prefer(lmerTest::lmer)
```

```
## [conflicted] Will prefer lmerTest::lmer over any other package.
```

``` r
conflicted::conflicts_prefer(dplyr::filter)
```

```
## [conflicted] Will prefer dplyr::filter over any other package.
```

## Load the data

``` r
simpop <- read_csv("../output/PopSim_Phenos-Genos_v1.csv")
```

```
## Rows: 12000 Columns: 31
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## dbl (31): germination.logit_pheno, establishment.logit_pheno, y1surv.logit_p...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
names(simpop)
```

```
##  [1] "germination.logit_pheno"   "establishment.logit_pheno"
##  [3] "y1surv.logit_pheno"        "alpha_pheno"              
##  [5] "beta_pheno"                "k_pheno"                  
##  [7] "delta_pheno"               "flowering.logit_pheno"    
##  [9] "fruitPerPlant_pheno"       "germination.prob"         
## [11] "establishment.prob"        "y1surv.prob"              
## [13] "flowering.prob"            "germinated"               
## [15] "established"               "y1surv"                   
## [17] "flowered"                  "Indiv_ID"                 
## [19] "week"                      "week_next"                
## [21] "elapsed_weeks"             "size"                     
## [23] "germination.logit_geno"    "establishment.logit_geno" 
## [25] "y1surv.logit_geno"         "alpha_geno"               
## [27] "beta_geno"                 "k_geno"                   
## [29] "delta_geno"                "flowering.logit_geno"     
## [31] "fruitPerPlant_geno"
```

``` r
#logit_pheno = scaled probability/phenotype
#_geno = genotype score for that trait 
#size = height calculated from simulated Weibull parameters (alpha, beta, k, and delta)
####Note: we don't have a genotype score for this since it's not explicitly estimated by the pop sim...
#elapsed_weeks = interval between weeks - should be all be 1 in this data
```

### Add size next and weekly surv

``` r
unique(simpop$established) #NA = did not germinate, 1 = established (survived first 3 weeks in field); 0 = did not establish
```

```
## [1] NA  1  0
```

``` r
unique(simpop$y1surv) #NA = did not establish, 1 = survived to week 12; 0 = did not survive 
```

```
## [1] NA  0  1
```

``` r
simpop %>% select(Indiv_ID, germinated, established, y1surv, k_pheno, week, size) %>% filter(Indiv_ID<5)
```

```
## # A tibble: 48 × 7
##    Indiv_ID germinated established y1surv k_pheno  week  size
##       <dbl>      <dbl>       <dbl>  <dbl>   <dbl> <dbl> <dbl>
##  1        1          0          NA     NA      NA     1  5.01
##  2        1          0          NA     NA      NA     2 NA   
##  3        1          0          NA     NA      NA     3 NA   
##  4        1          0          NA     NA      NA     4 NA   
##  5        1          0          NA     NA      NA     5 NA   
##  6        1          0          NA     NA      NA     6 NA   
##  7        1          0          NA     NA      NA     7 NA   
##  8        1          0          NA     NA      NA     8 NA   
##  9        1          0          NA     NA      NA     9 NA   
## 10        1          0          NA     NA      NA    10 NA   
## # ℹ 38 more rows
```

``` r
simpop_timeprep <- simpop %>%
  select(-alpha_pheno, -alpha_geno, -beta_pheno, -beta_geno, -k_pheno, -k_geno, -delta_pheno, -delta_geno) %>% 
  group_by(Indiv_ID) %>% 
  mutate(size = if_else(germinated==0, NA, #if a plant did not germinate it cannot have a size at any time
                        size), #in popsim code already have week 1 size = beta (min size)
         size_next = lead(size, order_by = week), #next time point's size; *see note about not having a genotype score for this...
         surv=if_else(germinated==0, NA, #surv = NA if you did not germ 
                        if_else(week==1, 1, #if you did germinate, assume you survived to week 1 (equivalent to pre-transplant size time point)
                                if_else(established==0, 0, #surv = 0 after week 1 if you did not establish 
                                        if_else(week<5, 1, #if you did establish, your surv is equal to 1 for the first 3 weeks post-transplant 
                                                y1surv) #your survival for weeks 5-12 is then equal to your y1surv
                                        )))) %>% #note we only have establishment and y1surv survival from the pop sim --> little weekly variation 
  #also geno_score is for the pheno in logit scale... does that matter?
  ungroup() %>%
  mutate(weeks = as.numeric(week - 1), #Weeks since week 1 (equivalent to pre-transplant or initial size) 
         size = if_else(surv==0, NA, size) #no size if no surv 
         ) %>%   
  drop_na(surv) #for models

simpop_timeprep %>% filter(Indiv_ID<10) %>% select(Indiv_ID, week, weeks, germinated, established, y1surv, size, size_next, surv) #look at data for a few indivs to see if above worked
```

```
## # A tibble: 24 × 9
##    Indiv_ID  week weeks germinated established y1surv  size size_next  surv
##       <dbl> <dbl> <dbl>      <dbl>       <dbl>  <dbl> <dbl>     <dbl> <dbl>
##  1        2     1     0          1           1      0  7.74      16.6     1
##  2        2     2     1          1           1      0 16.6       17.8     1
##  3        2     3     2          1           1      0 17.8       18.8     1
##  4        2     4     3          1           1      0 18.8       19.6     1
##  5        2     5     4          1           1      0 NA         20.4     0
##  6        2     6     5          1           1      0 NA         21.0     0
##  7        2     7     6          1           1      0 NA         21.5     0
##  8        2     8     7          1           1      0 NA         22.1     0
##  9        2     9     8          1           1      0 NA         22.5     0
## 10        2    10     9          1           1      0 NA         23.0     0
## # ℹ 14 more rows
```

``` r
simpop %>% filter(Indiv_ID==4) #was removed from the time prep dataset since it did not germinate --> no size or surv 
```

```
## # A tibble: 12 × 31
##    germination.logit_pheno establishment.logit_…¹ y1surv.logit_pheno alpha_pheno
##                      <dbl>                  <dbl>              <dbl>       <dbl>
##  1                  -0.842                  0.839             -0.259        53.5
##  2                  -0.842                  0.839             -0.259        53.5
##  3                  -0.842                  0.839             -0.259        53.5
##  4                  -0.842                  0.839             -0.259        53.5
##  5                  -0.842                  0.839             -0.259        53.5
##  6                  -0.842                  0.839             -0.259        53.5
##  7                  -0.842                  0.839             -0.259        53.5
##  8                  -0.842                  0.839             -0.259        53.5
##  9                  -0.842                  0.839             -0.259        53.5
## 10                  -0.842                  0.839             -0.259        53.5
## 11                  -0.842                  0.839             -0.259        53.5
## 12                  -0.842                  0.839             -0.259        53.5
## # ℹ abbreviated name: ¹​establishment.logit_pheno
## # ℹ 27 more variables: beta_pheno <dbl>, k_pheno <dbl>, delta_pheno <dbl>,
## #   flowering.logit_pheno <dbl>, fruitPerPlant_pheno <dbl>,
## #   germination.prob <dbl>, establishment.prob <dbl>, y1surv.prob <dbl>,
## #   flowering.prob <dbl>, germinated <dbl>, established <dbl>, y1surv <dbl>,
## #   flowered <dbl>, Indiv_ID <dbl>, week <dbl>, week_next <dbl>,
## #   elapsed_weeks <dbl>, size <dbl>, germination.logit_geno <dbl>, …
```

``` r
simpop_timeprep_size <- simpop_timeprep %>% drop_na(size, size_next) #for size models 
```

### Scaling and transformations

``` r
simpop_timeprep %>%  #skewed
  ggplot(aes(x=size)) +
  geom_histogram()

simpop_timeprep_scaled <- simpop_timeprep %>% 
  mutate(logSize=log(size))  #transform height to be more normal - only for surv models 
#this doesn't work because if you take NAs for size out then you only have surv = 1 ...
```

## Survival Models

### Phenos Only

#### Predicted vs. Observed Survival

### Phenos + Genos

#### Predicted vs. Observed Survival

## Growth Models

### Observed Patterns

### Phenos Only

#### Predicted vs. Observed Growth

### Phenos + Genos

#### Predicted vs. Observed Growth

## P Matrix

### Make a dataframe to store the parameters

### Use broom:tidy to create a df with the coef from each model

### Define the functions

### Define the structure of the IPM

### Make the matrices (G, S, and P)

### Plot the matrix

### Check for eviction
