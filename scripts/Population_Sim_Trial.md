---
title: "Population_Sim_Trial"
author: "Brandie QC"
date: "2025-10-14"
output: 
  html_document: 
    keep_md: true
editor_options: 
  markdown: 
    wrap: 72
---



# Population Simulation Attempt

Rishav demoed AlphaSimR:
<https://cran.r-project.org/web/packages/AlphaSimR/index.html>

-   Rishav's demo is on his GitHub:
    <https://github.com/rishavray/IPM_practice/blob/main/IPM_genetic_simulation.Rmd>

## Libraries


``` r
library(AlphaSimR)
```

```
## Loading required package: R6
```

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
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ✖ dplyr::mutate() masks AlphaSimR::mutate()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
logit <- brms::logit_scaled
inv_logit <- brms::inv_logit_scaled
```

## Create Founder Population


``` r
founderPop = runMacs(
  nInd = 1000, #number of individuals to simulate 
  nChr = 10, #number of chromosomes to simulate - may not be important to match to Streps' 14
  segSites = NULL, #number of segregating sites to keep per chromosome.
  species = "GENERIC" #The GENERIC history is meant to be a reasonable all-purpose choice. It runs quickly and models a population with an effective populations size that has gone through several historic bottlenecks.
)
```

## Get germ prob estimate from Maya's data

``` r
germprob <- read_csv("../data/WL2_2025_germ_prob.csv")
```

```
## Rows: 66 Columns: 4
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: ","
## chr (2): pop.id, Cross_Type
## dbl (2): prob, SE
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
summary(germprob) #overall mean prob = 0.39
```

```
##     pop.id               prob                 SE             Cross_Type       
##  Length:66          Min.   :0.0000001   Min.   :0.0000847   Length:66         
##  Class :character   1st Qu.:0.1751701   1st Qu.:0.0328350   Class :character  
##  Mode  :character   Median :0.3541667   Median :0.0755195   Mode  :character  
##                     Mean   :0.3856414   Mean   :0.0866275                     
##                     3rd Qu.:0.5484085   3rd Qu.:0.1155918                     
##                     Max.   :0.9999999   Max.   :0.3535534
```

``` r
germprob %>% 
  group_by(Cross_Type) %>% 
  summarise(n=n(), meangerm=mean(prob)) #parents mean = 0.18 (lowest); F2s = highest 0.47
```

```
## # A tibble: 4 × 3
##   Cross_Type     n meangerm
##   <chr>      <int>    <dbl>
## 1 BC1           14    0.414
## 2 F1            17    0.359
## 3 F2            25    0.469
## 4 Parent        10    0.182
```

``` r
germprob %>% 
  filter(Cross_Type=="Parent") #WL2 prob = 0.27
```

```
## # A tibble: 10 × 4
##    pop.id   prob     SE Cross_Type
##    <chr>   <dbl>  <dbl> <chr>     
##  1 BH     0.190  0.0495 Parent    
##  2 CC     0.0635 0.0307 Parent    
##  3 DPR    0.222  0.0524 Parent    
##  4 LV1    0.0317 0.0221 Parent    
##  5 SQ3    0.254  0.0548 Parent    
##  6 TM2    0.170  0.0219 Parent    
##  7 WL1    0.333  0.0594 Parent    
##  8 WL2    0.272  0.0258 Parent    
##  9 WV     0.222  0.0524 Parent    
## 10 YO11   0.0635 0.0307 Parent
```

## Define Genetic Architecture for Traits

``` r
SP <- SimParam$new(founderPop)
#Starts the process of building a new simulation by creating a new SimParam object and assigning a founder population to the class.

traitMeans <- c(germination.logit = logit(0.18), #logit to scale it (from gene --> trait); used avg for all parents from Maya's data 
                establishment.logit = logit(0.67),
                y1surv.logit = logit(0.45),
                #alpha = 40, #Hmax = alpha -Sets the asymptote (max growth) for the Weibull model 
                beta = 7, #Hmin = beta -Changes the y intercept (the lower limit) for the Weibull model- recruit size?
                k100 = 1, # Divide by 100 to get K. Often described as the growth rate --\> related to the slope of the line for the Weibull model 
                #delta = 0.5, #Sets the inflection point - when you shift b/t concave to convex for the Weibull model 
                flowering.logit = logit(0.2),
                fruitPerPlant = 5)

SP$addTraitA(nQtlPerChr = 10, 
              mean = traitMeans, 
              var = c(germination.logit=1,
                                 establishment.logit=1,
                                 y1surv.logit=1,
                                 beta=2,
                                 k100 = 1,
                                 flowering.logit=1,
                                 fruitPerPlant=2)^2, #genetic Var.  Numbers are in SD
              gamma = TRUE, #a lot of small effect loci and a few large effect loci 
              shape = 1, #shape parameter for gamma distribution 
              name=names(traitMeans))
#addTraitAG - Randomly assigns eligible QTLs for one or more additive GxE traits. If simulating more than one trait, all traits will be pleiotropic with correlated effects.
#addTraitA - Randomly assigns eligible QTLs for one or more additive traits. If simulating more than one trait, all traits will be pleiotropic with correlated additive effects. (we will keep it simple for now)
##nQtlPerChr = number of QTLs per chromosome. Can be a single value or nChr values.
##mean = a vector of desired mean genetic values for one or more traits
##var = a vector of desired genetic variances for one or more traits
##varGxE = a vector of total genotype-by-environment variances for the traits
##name = optional name for trait(s)

SP$setVarE(h2=rep(0.7, length(traitMeans))) #set heritability
```

## Create Initial Population


``` r
pop = newPop(founderPop, simParam = SP)
pop1 <- setPheno(pop)
```


``` r
pheno <- pop1@pheno %>% as_tibble()
dim(pheno)
```

```
## [1] 1000    7
```

``` r
head(pheno)
```

```
## # A tibble: 6 × 7
##   germination.logit establishment.logit y1surv.logit  beta   k100
##               <dbl>               <dbl>        <dbl> <dbl>  <dbl>
## 1           -0.520            -0.000729        0.309  8.49 -1.30 
## 2           -0.0644            0.573           0.362  3.54  0.590
## 3            0.163             0.158          -0.323  9.01  1.32 
## 4           -0.408             2.72           -0.565  9.55  1.05 
## 5           -1.60             -0.807          -1.88   5.56  1.87 
## 6           -2.52              2.48            0.208  7.36  0.863
## # ℹ 2 more variables: flowering.logit <dbl>, fruitPerPlant <dbl>
```

``` r
summary(pheno)
```

```
##  germination.logit establishment.logit  y1surv.logit          beta       
##  Min.   :-6.0997   Min.   :-2.95863    Min.   :-3.6815   Min.   :-2.783  
##  1st Qu.:-2.3015   1st Qu.:-0.03985    1st Qu.:-0.9592   1st Qu.: 5.397  
##  Median :-1.4774   Median : 0.63469    Median :-0.1940   Median : 7.091  
##  Mean   :-1.5036   Mean   : 0.68785    Mean   :-0.1907   Mean   : 7.004  
##  3rd Qu.:-0.6356   3rd Qu.: 1.47446    3rd Qu.: 0.5776   3rd Qu.: 8.605  
##  Max.   : 1.5168   Max.   : 5.03885    Max.   : 3.4960   Max.   :14.474  
##       k100         flowering.logit   fruitPerPlant   
##  Min.   :-2.9164   Min.   :-4.8794   Min.   :-1.778  
##  1st Qu.: 0.2195   1st Qu.:-2.2190   1st Qu.: 3.450  
##  Median : 1.0343   Median :-1.4758   Median : 5.001  
##  Mean   : 1.0421   Mean   :-1.3757   Mean   : 4.996  
##  3rd Qu.: 1.8400   3rd Qu.:-0.5779   3rd Qu.: 6.594  
##  Max.   : 5.2876   Max.   : 3.3413   Max.   :12.776
```

### Convert pheno logits to probabilities and phenotypes

``` r
pheno <- pheno %>%
  mutate(across(ends_with(".logit"), .fns = inv_logit, 
         .names = "{.col}.prob")) %>%
  rename_with(.fn = \(n) str_replace(n, "\\.logit\\.prob", "\\.prob")) %>%
  mutate(germinated = rbinom(n(), size=1, prob=germination.prob),
         established = ifelse(germinated, rbinom(n(), size=1, prob=establishment.prob), NA),
         beta,
         k = k100/100,
         k = ifelse(k <=0, min(k[k>0], na.rm = TRUE), k), # no negative k
         k = ifelse(established, k, NA), #can only have a growth rate if you establish
         fruitPerPlant = round(fruitPerPlant), # no fractional fruit
         fruitPerPlant = ifelse(fruitPerPlant < 0, 0, fruitPerPlant) # no negative fruit
  ) %>% select(-k100)
```


``` r
pheno %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  facet_wrap(~name, scales="free")
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

```
## Warning: Removed 1595 rows containing non-finite outside the scale range
## (`stat_bin()`).
```

![](Population_Sim_Trial_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

## Try to get weekly growth data

### Model growth with the Weibull model


``` r
weibull <- function (t, alpha, beta, k, delta) {
  result <- alpha - (alpha - beta) * exp(-(k * t)^delta)
  return(result)
}

##Example:
growth <- tibble(week = seq(0,12, 1)) %>%
  mutate(size = weibull(t = week,
                          alpha = 40, #Hmax = alpha -  Sets the asymptote (max growth)
                          beta = 0.5, #Hmin = beta - Changes the y intercept (the lower limit)
                          k = 0.01, #k - Often described as the growth rate --\> related to the slope of the line
                          delta = 0.5)) #delta - Sets the inflection point - when you shift b/t concave to convex

growth %>%  ggplot(aes(x=week, y=size)) +
  geom_line()
```

![](Population_Sim_Trial_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

### Function for weekly survival that includes size

``` r
weekly_surv <- function(d, weeks, estab_weeks, size.coef=.25, size2.coef=0) {
  # function to calculate weekly survival as a function of
  # year1 survival and size
  d %>%
    mutate(weekly.surv.prob = ifelse(week <= estab_weeks, established,
                                     min(1, # make sure it is not > 1
                                     y1surv.prob^(1/(weeks-estab_weeks)) + # basic weekly survival
                                       (size.coef * size.scale + size2.coef * size.scale^2) #*  # the proportion of the growth rate is affected by size
                                       #y1surv.prob^(1/(weeks-estab_weeks))  # scaling the size adjustment to the surv.prob scale
                                     ) # min
                                     )
           ) %>% 
           pull(weekly.surv.prob)
}
```

### Pheno weekly

``` r
weeks <- 12
estab_weeks <- 3

pheno_weekly <- pheno %>% 
  mutate(Indiv_ID=row_number()) %>% 
  slice(rep(1:n(), each = 12)) %>% #duplicate each row 12 times 
  group_by(Indiv_ID) %>% 
  
  #weekly size 
  mutate(week=row_number(), #12 weeks for each indiv
         week_next=lead(week, order_by = week), 
         elapsed_weeks= week_next - week, #interval length in weeks
         size = if_else(germinated==1, if_else(week==1, beta, #size in week 1 is the min size (beta)
                        weibull(t = week, #otherwise use the weibull formula to calculate size
                                alpha = 40,
                                beta = beta,
                                k = k, 
                                delta = 0.5)),
                        NA) # didn't germinate
                        ) %>%
  ungroup() %>%
  mutate(size.scale = as.numeric(scale(size, scale = diff(range(size, na.rm = TRUE))))) %>% # could also try traditional scaling, but this gives range of -1 to +1
  group_by(Indiv_ID) %>%
  
  #weekly survival
  mutate(weekly.surv.prob = weekly_surv(pick(size.scale, y1surv.prob, week, established), weeks, estab_weeks), 
         alive = ifelse(established, rbinom(n(), size = 1, prob = weekly.surv.prob), established), #stay alive this week?
         alive = cummin(alive),  # if dead in a previous week, stay dead
         y1surv = min(alive)) %>%
  
  # update other traits based on death week
  mutate(size=ifelse(alive==1, size, NA), # if plant is dead, make size NA
         flowered = ifelse(sum(alive)>=8,  #only flower if alive for at least 8 weeks
           rbinom(n(), size = 1, prob = flowering.prob), NA),
         fruitPerPlant = ifelse(flowered & sum(alive)>=10, fruitPerPlant, NA)
         ) %>%
  ungroup()
```

## Extract Genotype Scores

``` r
# Extract genotype scores for each individual
genotype_scores <- pop1@gv %>%
    as_tibble() %>%
    rownames_to_column(var = "Indiv_ID") %>%
    mutate(Indiv_ID = as.integer(Indiv_ID),
         k=k100/100,
         k = ifelse(k <=0, min(k[k>0], na.rm = TRUE), k) # no negative k
) %>%
    select(Indiv_ID, everything())

# Join genotype scores with phenotypic data
pheno_weekly_merged <- pheno_weekly %>%
    left_join(genotype_scores, by = "Indiv_ID", suffix = c("_pheno", "_geno")) %>% 
  mutate(across(matches("beta|fruit"), \(x) ifelse(x < 0, min(x[x > 0], na.rm = TRUE), x))) #deal with initial size or fruit # less than 0 

pheno_weekly_merged %>% 
  select(Indiv_ID, week, size, week_next, germination.logit_geno, germination.logit_pheno)
```

```
## # A tibble: 12,000 × 6
##    Indiv_ID  week  size week_next germination.logit_geno germination.logit_pheno
##       <int> <int> <dbl>     <int>                  <dbl>                   <dbl>
##  1        1     1  8.49         2                 -0.794                  -0.520
##  2        1     2  8.78         3                 -0.794                  -0.520
##  3        1     3  8.84         4                 -0.794                  -0.520
##  4        1     4  8.90         5                 -0.794                  -0.520
##  5        1     5  8.94         6                 -0.794                  -0.520
##  6        1     6  8.99         7                 -0.794                  -0.520
##  7        1     7  9.02         8                 -0.794                  -0.520
##  8        1     8  9.06         9                 -0.794                  -0.520
##  9        1     9  9.09        10                 -0.794                  -0.520
## 10        1    10  9.13        11                 -0.794                  -0.520
## # ℹ 11,990 more rows
```

``` r
write_csv(pheno_weekly_merged, "../output/PopSim_Phenos-Genos_v4.csv")
```

