Code availability and usage for NCOMMS-20-20430
================
Hunter Stanke
12/30/2020

<style>
div.sourceCode {
    overflow: hidden;
}
</style>

<br> <br>

> This document serves as a description of contents in `FSI_code.zip`
> and a guide to reproducing results described in `NCOMMS-20-20430`. We
> describe all system requirements and installation guidelines required
> to install the publicly available `rFIA` R package, which contains
> most of the code necessary to reproduce our results (all remaining
> code available in `FSI_code.zip`). We provide a demonstration of our
> software (implemented in `rFIA`) as well as specific instructions for
> software use and replication of results. Contents of `FIA`, `PRISM`,
> and `results` directories within `FSI_code.zip` are empty initially.
> All required data will be downloaded and saved upon running `.R`
> scripts.

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
    ## 1  2012 -0.00197   -0.728 Decline    2.35e-5      0.00814   0.271   0.251
    ## 2  2013 -0.00288   -1.05  Decline    1.40e-5      0.00466   0.274   0.246
    ## 3  2014 -0.00272   -1.03  Decline    8.74e-6      0.00301   0.263   0.235
    ## 4  2015 -0.00286   -1.06  Decline    6.83e-6      0.00228   0.270   0.242
    ## 5  2016 -0.00293   -1.09  Decline    5.63e-6      0.00189   0.268   0.239
    ## 6  2017 -0.00278   -1.06  Decline    4.52e-6      0.00155   0.263   0.235
    ## 7  2018 -0.00282   -1.06  Decline    4.01e-6      0.00136   0.267   0.239
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
    ##  1  2012 4.04e13     3.66e12              1 Forest      1_8_~   10.1  2.08e-4
    ##  2  2012 4.04e13     3.66e12              1 Forest      1_8_~    9.9  6.29e-4
    ##  3  2012 4.04e13     3.66e12              1 Forest      1_8_~   10.1  5.62e-3
    ##  4  2012 4.04e13     3.64e12              1 Forest      1_8_~    9.8 -1.43e-3
    ##  5  2012 4.04e13     3.64e12              1 Forest      1_8_~    9.8  1.44e-2
    ##  6  2012 4.04e13     3.64e12              1 Forest      1_8_~   10    1.54e-3
    ##  7  2012 4.04e13     3.64e12              1 Forest      1_8_~   10    4.02e-4
    ##  8  2012 4.04e13     3.64e12              1 Forest      1_8_~    9.9  3.64e-3
    ##  9  2012 4.04e13     3.64e12              1 Forest      1_8_~   10    9.01e-3
    ## 10  2012 4.04e13     3.65e12              1 Forest      1_8_~    9.5 -4.74e-3
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
    ##  1  2012    15 white fir   Abies concolor   6.69e-4   0.678  Expand     2.38e-4
    ##  2  2012    18 corkbark f~ Abies lasiocar~  1.05e-3   1.90   Expand     5.21e-4
    ##  3  2012    19 subalpine ~ Abies lasiocar~ -9.76e-4  -0.675  Decline    1.14e-4
    ##  4  2012    65 Utah junip~ Juniperus oste~  2.74e-3   1.18   Expand     1.81e-4
    ##  5  2012    66 Rocky Moun~ Juniperus scop~ -9.79e-3  -6.69   Decline    3.80e-4
    ##  6  2012    69 oneseed ju~ Juniperus mono~  1.48e-2   8.43   Expand     1.49e-3
    ##  7  2012    93 Engelmann ~ Picea engelman~  1.86e-4   0.0866 Expand     9.52e-5
    ##  8  2012    96 blue spruce Picea pungens    4.10e-3   3.11   Expand     3.37e-4
    ##  9  2012   102 Rocky Moun~ Pinus aristata   4.37e-4   0.487  Expand     7.25e-5
    ## 10  2012   106 common or ~ Pinus edulis    -8.47e-4  -0.827  Decline    3.17e-5
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
    grpBy = OWNCD,
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

    ## # A tibble: 70 x 21
    ##     YEAR OWNCD       FSI PERC_FSI FSI_STATUS   FSI_INT PERC_FSI_INT PREV_RD
    ##    <int> <int>     <dbl>    <dbl> <chr>          <dbl>        <dbl>   <dbl>
    ##  1  2012    11  -3.18e-3   -1.03  Decline      5.88e-5       0.0176   0.307
    ##  2  2012    12 NaN        NaN     Stable     NaN           NaN      NaN    
    ##  3  2012    21   1.59e-3    0.377 Expand       3.76e-4       0.0965   0.421
    ##  4  2012    22  -5.08e-4   -0.229 Decline      3.64e-5       0.0157   0.221
    ##  5  2012    23 NaN        NaN     Stable     NaN           NaN      NaN    
    ##  6  2012    24  -1.48e-3   -0.670 Decline      2.38e-4       0.0937   0.222
    ##  7  2012    25 NaN        NaN     Stable     NaN           NaN      NaN    
    ##  8  2012    31  -8.67e-4   -0.302 Decline      1.24e-4       0.0410   0.287
    ##  9  2012    32   6.28e-3    5.93  Expand       8.82e-4       0.818    0.106
    ## 10  2012    46  -1.25e-3   -0.581 Decline      4.33e-5       0.0188   0.215
    ## # ... with 60 more rows, and 13 more variables: CURR_RD <dbl>, TPA_RATE <dbl>,
    ## #   BA_RATE <dbl>, REMPER <dbl>, FSI_VAR <dbl>, PERC_FSI_VAR <dbl>,
    ## #   PREV_RD_VAR <dbl>, CURR_RD_VAR <dbl>, TPA_RATE_VAR <dbl>,
    ## #   BA_RATE_VAR <dbl>, REMPER_VAR <dbl>, nPlots <dbl>, N <int>

# Instructions for use

## *Applying the FSI to our study region *

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

*Reproduction of our results*

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
