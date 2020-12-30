
<style>
div.sourceCode {
    overflow: hidden;
}
</style>

<br> <br>

> This document serves as a description of contents in `FSI_code.zip`
> and a guide to reproducing results described in “Over half of western
> United States’ most abundant tree species in decline” (Stanke et al.,
> 2020; `NCOMMS-20-20430`). We describe all system requirements and
> installation guidelines required to install the publicly available
> `rFIA` R package, which contains most of the code necessary to
> reproduce our results (all remaining code available in
> `FSI_code.zip`). We provide a demonstration of our software
> (implemented in `rFIA`) as well as specific instructions for software
> use and replication of results. Contents of `FIA` and `results`
> directories within `FSI_code.zip` are empty initially. All required
> data will be downloaded and saved upon running `.R` scripts.

# System Requirements

We have implemented methods to compute the Forest Stability Index
(described in `NCOMMS-20-20430`) using publicly available USDA Forest
Inventory and Analysis database in the open source R package, `rFIA`.
Methods in the `rFIA` R package and other code used to produce the
results outlined in `NCOMMS-20-20430` have been tested on Windows 10,
Ubuntu 18.04, and R versions 3.5.1 and 4.0.0. For more on the rFIA R
package, check out our [website](https://rfia.netlify.app/).

**The data used to produce results outlined in `NCOMMS-20-20430` may be
too large to be processed on a standard desktop computer (RAM
limitations may prevent loading/processing in R). We recommend using a
machine with at least 64GB of RAM to reproduce our results (all analysis
presented in `NCOMMS-20-20430` were completed on a machine with 384GB of
RAM).** We recognize that reviewers may not have access to these
computational resources. In such event, we have provided an additional
alternative workflow that uses a subset of the full dataset for
demonstration and review purposes (see lines 51-54 of `estimate.R`).

<br>

# Installation

Users can install the released version of `rFIA` from CRAN with:

``` r
install.packages('rFIA')
```

or the development version from GitHub with:

``` r
devtools::install_github('hunter-stanke/rFIA')
```

<br>

# Software demonstration

As an example, we will implement the Forest Stability Index using USDA
Forest Inventory and Analyis (FIA) data collected between 2001-2019 in
the state of Colorado, USA. All code shown below can be run in R
independent of the contents of `FSI_code.zip`. Processing time will vary
with network speed as data must be downloaded from an online repository,
however should not exceed 5 minutes on a normal desktop computer. All
results will be returned as `data.frame` objects in R.

We first download the FIA data for Colorado using `rFIA`:

``` r
library(rFIA)
library(dplyr)

## Dowloading and loading into R
co <- getFIA('CO')
```

We can then use the `fsi` function to estimate the Forest Stability
Index for all forest in the state using the ‘simple moving average’
estimator. As done in `NCOMMS-20-20430`, we allow slopes and intercepts
of maximum size density curves to vary by forest type using the
`scaleBy` argument:

``` r
## FSI for all forest
fsi(co, 
    method = 'sma',
    scaleBy = FORTYPCD)
```

    ## Modeling maximum size-density curve(s)...
    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 1242
    ##    Unobserved stochastic nodes: 1289
    ##    Total graph size: 13774
    ## 
    ## Initializing model

    ## # A tibble: 7 x 20
    ##    YEAR      FSI PERC_FSI FSI_STATUS FSI_INT PERC_FSI_INT PREV_RD CURR_RD
    ##   <int>    <dbl>    <dbl> <chr>        <dbl>        <dbl>   <dbl>   <dbl>
    ## 1  2012 -0.00194   -0.728 Decline    2.31e-5      0.00813   0.266   0.247
    ## 2  2013 -0.00283   -1.05  Decline    1.37e-5      0.00466   0.269   0.241
    ## 3  2014 -0.00267   -1.04  Decline    8.58e-6      0.00301   0.258   0.231
    ## 4  2015 -0.00281   -1.06  Decline    6.71e-6      0.00228   0.265   0.237
    ## 5  2016 -0.00288   -1.09  Decline    5.53e-6      0.00189   0.263   0.234
    ## 6  2017 -0.00273   -1.06  Decline    4.43e-6      0.00155   0.258   0.231
    ## 7  2018 -0.00277   -1.06  Decline    3.94e-6      0.00136   0.262   0.234
    ## # ... with 12 more variables: TPA_RATE <dbl>, BA_RATE <dbl>, REMPER <dbl>,
    ## #   FSI_VAR <dbl>, PERC_FSI_VAR <dbl>, PREV_RD_VAR <dbl>, CURR_RD_VAR <dbl>,
    ## #   TPA_RATE_VAR <dbl>, BA_RATE_VAR <dbl>, REMPER_VAR <dbl>, nPlots <dbl>,
    ## #   N <int>

where `FSI` is the estimated Forest Stability Index for all forestland
in CO by `YEAR`. `SI_STATUS` indicates the status of the population:
decline if `FSI` is negative and confidence interval (`FSI_INT`)
excludes zero; expand if `FSI` is positive and confidence interval
excludes zero; stable otherwise. Estimates can be returned at the
plot-level by specifying `byPlot = TRUE`:

``` r
## FSI plot-level
fsi(co, 
    byPlot = TRUE,
    scaleBy = FORTYPCD) %>%
   ## Forested plots only
   filter(PLOT_STATUS_CD == 1) %>%
   filter(!is.na(PREV_PLT_CN))
```

    ## Modeling maximum size-density curve(s)...
    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 1242
    ##    Unobserved stochastic nodes: 1289
    ##    Total graph size: 13774
    ## 
    ## Initializing model

    ## # A tibble: 2,661 x 15
    ##     YEAR  PLT_CN PREV_PLT_CN PLOT_STATUS_CD PLOT_STATUS pltID REMPER      FSI
    ##    <int>   <dbl>       <dbl>          <int> <chr>       <chr>  <dbl>    <dbl>
    ##  1  2012 4.04e13     3.66e12              1 Forest      1_8_~   10.1  2.06e-4
    ##  2  2012 4.04e13     3.66e12              1 Forest      1_8_~    9.9  6.25e-4
    ##  3  2012 4.04e13     3.66e12              1 Forest      1_8_~   10.1  5.58e-3
    ##  4  2012 4.04e13     3.64e12              1 Forest      1_8_~    9.8 -1.41e-3
    ##  5  2012 4.04e13     3.64e12              1 Forest      1_8_~    9.8  1.43e-2
    ##  6  2012 4.04e13     3.64e12              1 Forest      1_8_~   10    1.53e-3
    ##  7  2012 4.04e13     3.64e12              1 Forest      1_8_~   10    3.95e-4
    ##  8  2012 4.04e13     3.64e12              1 Forest      1_8_~    9.9  3.61e-3
    ##  9  2012 4.04e13     3.64e12              1 Forest      1_8_~   10    8.93e-3
    ## 10  2012 4.04e13     3.65e12              1 Forest      1_8_~    9.5 -4.71e-3
    ## # ... with 2,651 more rows, and 7 more variables: PERC_FSI <dbl>,
    ## #   PREV_RD <dbl>, CURR_RD <dbl>, PREV_TPA <dbl>, CURR_TPA <dbl>,
    ## #   PREV_BAA <dbl>, CURR_BAA <dbl>

We can group estimates by species by specifying `bySpecies = TRUE`:

``` r
## FSI by species
fsi(co, 
    method = 'sma', 
    bySpecies = TRUE,
    scaleBy = FORTYPCD)
```

    ## Modeling maximum size-density curve(s)...
    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 1242
    ##    Unobserved stochastic nodes: 1289
    ##    Total graph size: 13774
    ## 
    ## Initializing model

    ## # A tibble: 158 x 23
    ##     YEAR  SPCD COMMON_NAME SCIENTIFIC_NAME      FSI PERC_FSI FSI_STATUS FSI_INT
    ##    <int> <int> <chr>       <chr>              <dbl>    <dbl> <chr>        <dbl>
    ##  1  2012    15 white fir   Abies concolor   6.70e-4   0.676  Expand     2.38e-4
    ##  2  2012    18 corkbark f~ Abies lasiocar~  1.06e-3   1.94   Expand     5.20e-4
    ##  3  2012    19 subalpine ~ Abies lasiocar~ -9.97e-4  -0.689  Decline    1.15e-4
    ##  4  2012    65 Utah junip~ Juniperus oste~  2.74e-3   1.17   Expand     1.83e-4
    ##  5  2012    66 Rocky Moun~ Juniperus scop~ -9.90e-3  -6.70   Decline    3.84e-4
    ##  6  2012    69 oneseed ju~ Juniperus mono~  1.50e-2   8.44   Expand     1.51e-3
    ##  7  2012    93 Engelmann ~ Picea engelman~  1.90e-4   0.0880 Expand     9.63e-5
    ##  8  2012    96 blue spruce Picea pungens    4.12e-3   3.13   Expand     3.38e-4
    ##  9  2012   102 Rocky Moun~ Pinus aristata   4.34e-4   0.480  Expand     7.27e-5
    ## 10  2012   106 common or ~ Pinus edulis    -8.49e-4  -0.828  Decline    3.16e-5
    ## # ... with 148 more rows, and 15 more variables: PERC_FSI_INT <dbl>,
    ## #   PREV_RD <dbl>, CURR_RD <dbl>, TPA_RATE <dbl>, BA_RATE <dbl>, REMPER <dbl>,
    ## #   FSI_VAR <dbl>, PERC_FSI_VAR <dbl>, PREV_RD_VAR <dbl>, CURR_RD_VAR <dbl>,
    ## #   TPA_RATE_VAR <dbl>, BA_RATE_VAR <dbl>, REMPER_VAR <dbl>, nPlots <dbl>,
    ## #   N <int>

or by any other grouping variable contained in the `co` database using
the `grpBy` argument:

``` r
## Grouping by ownership class
fsi(co, 
    method = 'sma',
    grpBy = OWNGRPCD,
    scaleBy = FORTYPCD)
```

    ## Modeling maximum size-density curve(s)...
    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 1242
    ##    Unobserved stochastic nodes: 1289
    ##    Total graph size: 13774
    ## 
    ## Initializing model

    ## # A tibble: 28 x 21
    ##     YEAR OWNGRPCD      FSI PERC_FSI FSI_STATUS FSI_INT PERC_FSI_INT PREV_RD
    ##    <int>    <int>    <dbl>    <dbl> <chr>        <dbl>        <dbl>   <dbl>
    ##  1  2012       10 -3.19e-3   -1.04  Decline    5.33e-5      0.0159    0.308
    ##  2  2012       20 -3.84e-4   -0.161 Decline    3.33e-5      0.0135    0.238
    ##  3  2012       30  3.29e-4    0.128 Expand     1.14e-4      0.0453    0.257
    ##  4  2012       40 -1.23e-3   -0.570 Decline    4.30e-5      0.0187    0.216
    ##  5  2013       10 -4.46e-3   -1.43  Decline    3.08e-5      0.00881   0.311
    ##  6  2013       20 -1.18e-3   -0.497 Decline    2.14e-5      0.00844   0.238
    ##  7  2013       30  4.01e-4    0.184 Expand     5.28e-5      0.0247    0.219
    ##  8  2013       40 -1.26e-3   -0.561 Decline    2.47e-5      0.0104    0.225
    ##  9  2014       10 -4.20e-3   -1.42  Decline    1.85e-5      0.00554   0.295
    ## 10  2014       20 -1.29e-3   -0.547 Decline    1.46e-5      0.00577   0.236
    ## # ... with 18 more rows, and 13 more variables: CURR_RD <dbl>, TPA_RATE <dbl>,
    ## #   BA_RATE <dbl>, REMPER <dbl>, FSI_VAR <dbl>, PERC_FSI_VAR <dbl>,
    ## #   PREV_RD_VAR <dbl>, CURR_RD_VAR <dbl>, TPA_RATE_VAR <dbl>,
    ## #   BA_RATE_VAR <dbl>, REMPER_VAR <dbl>, nPlots <dbl>, N <int>

# Instructions for use

## *Applying the FSI to our study region*

To apply the `FSI` to our study region (i.e., western US), we simply
expand our population of interest to include Washington, Oregon,
California, Idaho, Montana, Utah, Nevada, Arizona, New Mexico, and
Colorado:

``` r
## Download and load data for all states
db <- getFIA(c('WA', 'OR', 'CA', 'NV', 'AZ',
               'NM', 'CO', 'UT', 'MT', 'ID'))

## Estimate the FSI by species across the full
## study region (range-wide indices)
rangeWide <- fsi(db, 
                 method = 'sma',
                 bySpecies = TRUE,
                 scaleBy = FORTYPCD)
```

These data are large (\~8GB), and download speeds will depend on a users
network connection. Once downloaded, processing should not exceed 10
minutes.

*Reproducing of our results*

To reproduce the results outlined in `NCOMMS-20-20430`, modify working
directory in `estimate.R`, `climate.R`, and `model.R` to the location
where `FSI_code` is unzipped. Then simply source the `.R` files listed
at the root of the `FSI_code` directory with the following (must be run
in this sequence):

``` r
## Download data and estimate the FSI by 
## populations of interest
source('./FSI_code/estimate.R')


## Run disturbance model
source('./FSI_code/disturbMod.R')
```

The above will (1) download all FIA data required for the analysis, and
save the data in the `FIA` directory. (2) Produce estimates of the FSI
for populations of interest (i.e., species). Results will be stored in
the `results` directory, with plot-level estimates found in
`results/plt`, ecoregion-scale estimates found in `results/spatial` and
range-wide estimates found in `results/regionWide`. Model results
(estimated coefficients and associated credible intervals) will be
stored in `results/model`.
