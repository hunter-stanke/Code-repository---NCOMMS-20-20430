##=====================================================
##=====================================================
##
## This script contains all code necessary to produce
## estimates of the Forest Stability Index used in the
## analysis outlined in NCOMMS-20-20430. Specifically,
## the code below (1) downloads the necessary Forest
## Inventory and Analysis and (2) uses those data to
## estimate species-level FSI across the entire study
## region, within ecoregion subsections & divisions,
## and at individual FIA plots. All estimates are
## saved and used in subsequent scripts.
##
## Estimation of the FSI is carried out using the 'fsi'
## function implemented in the 'rFIA' package. The source
## code for the 'fsi' function can be found in 'fsi.R'
## and is identical to that in the rFIA package. Therefore
## there is no reason to source the 'fsi.R' file if 'rFIA'
## is loaded, it is provided purely for convenience.
##
## Last modified: 14 September 2020 - Hunter Stanke
##
##====================================================
##====================================================





##===============================================================================
##  Set up your working directory (location where FSI_code is unzipped) ---------
##===============================================================================

library(rFIA)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

## Set working directory
setwd('/home/hunter/gcb/natureComm/revision')

## How many cores do you want to use?
## Default to one
cores <- 1


##===============================================================================
##  Download/ read FIA data -----------------------------------------------------
##===============================================================================
## Download data for all states
getFIA(c('WA', 'OR', 'CA', 'NV', 'AZ',
         'NM', 'CO', 'UT', 'MT', 'ID'),
       dir = './FSI_code/FIA/',
       load = FALSE)

## If the above fails due to memory limitations,
## use a subset of the study region for review
## purposes
#getFIA(c('CO'),
#       dir = './FSI_code/FIA/',
#       load = FALSE)


## Read it in - very large
im <- readFIA('./FSI_code/FIA/', nCores = cores)
im$PLOT$PLT_CN <- im$PLOT$CN


## A most recent subset to use for running through the estimators
imMR <- clipFIA(im)


## Get rid of whitespace - make codes for
## ecoregion divisions
imMR$PLOT$ECOSUBCD <- str_trim(imMR$PLOT$ECOSUBCD)
imMR$PLOT$ECODIV <- str_sub(imMR$PLOT$ECOSUBCD, 1, -4)
## One is sticky
imMR$PLOT <- imMR$PLOT %>%
  mutate(ECODIV = case_when(ECODIV == 'M332' ~ 'M33',
                            TRUE ~ ECODIV))



#===============================================================================
#########  Compute disturbance and treatment indices ---------------------------
#===============================================================================

## Don't want to drop NAs here
imMR$COND <- imMR$COND %>%
  mutate(DSTRBCD1 = replace_na(DSTRBCD1, 0),
         DSTRBCD2 = replace_na(DSTRBCD2, 0),
         DSTRBCD3 = replace_na(DSTRBCD3, 0),
         DSTRBYR1 = replace_na(DSTRBYR1, 0),
         DSTRBYR2 = replace_na(DSTRBYR2, 0),
         DSTRBYR3 = replace_na(DSTRBYR3, 0),
         TRTCD1 = replace_na(TRTCD1, 0),
         TRTCD2 = replace_na(TRTCD2, 0),
         TRTCD3 = replace_na(TRTCD3, 0),
         TRTYR1 = replace_na(TRTYR1, 0),
         TRTYR2 = replace_na(TRTYR2, 0),
         TRTYR3 = replace_na(TRTYR3, 0))

## Proportion of each plot subject to disturbance
disturb <- area(imMR, byPlot = TRUE, nCores = cores,
                grpBy = c(DSTRBCD1, DSTRBYR1, DSTRBCD2, DSTRBYR2, DSTRBCD3, DSTRBYR3))
fwrite(disturb, './FSI_code/results/plt/disturb.csv')


## Proportion of each plot subject to treatment
treat <- area(imMR, byPlot = TRUE, nCores = cores,
              grpBy = c(TRTCD1, TRTYR1, TRTCD2, TRTYR2, TRTCD3, TRTYR3))
fwrite(treat, './FSI_code/results/plt/treat.csv')


## Binary variable for treatment - these are the plots where
## harvest/ artificial regen occured on a plot
treat <- treat %>%
  as.data.frame() %>%
  filter(!c(TRTCD1 %in% c(0, 40))) %>%
  mutate(treat = 1) %>%
  distinct(pltID, treat)

## Join on the PLOT table by pltID
imMR$PLOT <- imMR$PLOT %>%
  mutate(pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_')) %>%
  left_join(treat, by = 'pltID') %>%
  mutate(treat = replace_na(treat, 0))


##===============================================================================
##  Total abundance by species --------------------------------------------------
##===============================================================================

### Get total abundance by species
baaFull <- tpa(im,
               nCores = cores,
               method = 'sma',
               bySpecies = TRUE,
               totals = TRUE)

data.table::fwrite(baaFull, './FSI_code/results/regionWide/baaFull.csv')


##===============================================================================
##  Tree size percentiles by species --------------------------------------------
##===============================================================================

## Estimating size class percentiles
spPercentiles <- imMR$TREE %>%
  left_join(select(imMR$COND, PLT_CN, CONDID, SITECLCD), by = c('PLT_CN', 'CONDID')) %>%
  filter(!is.na(DIA)) %>%
  filter(STATUSCD == 1) %>%
  ## DIA is measured to tenth of an inch
  group_by(SPCD, SITECLCD) %>%
  mutate(spSC = percent_rank(DIA),
         spSC = case_when(spSC >= 1 ~ .9999,
                          TRUE ~ spSC)) %>%
  ungroup() %>%
  select(CN, spSC)
## Panics in a mutate, sorry
spPercentiles$spSC <- makeClasses(spPercentiles$spSC,
                                        lower = 0, interval = .1, numLabs = TRUE) + .05

## Join on the rankings
imMR$TREE <- imMR$TREE %>%
  left_join(spPercentiles, by = 'CN')


##===============================================================================
##  Compute region-wide indices -------------------------------------------------
##===============================================================================

## Species range-wide estimates of FSI
full <- fsi(imMR,
            method = 'sma',
            nCores = cores,
            bySpecies = TRUE,
            areaDomain = treat == 0,
            scaleBy = c(FORTYPCD),
            returnBetas = TRUE)
fwrite(full$results, './FSI_code/results/regionWide/rangeWide.csv')
fwrite(full$betas, './FSI_code/results/regionWide/betas.csv')

## Species range-wide estimates of FSI by size class (4 inch)
fullSC <- fsi(imMR,
              method = 'sma',
              nCores = cores,
              bySpecies = TRUE,
              areaDomain = treat == 0,
              grpBy = spSC,
              scaleBy = c(FORTYPCD),
              betas = full$betas)
fwrite(fullSC, './FSI_code/results/regionWide/rangeWideSC.csv')



#===============================================================================
#########  Compute ecoregion-level indices -------------------------------------
#===============================================================================

## Subsection ---------------

## FSI by species, by ecoregion
ecosub <- rFIA:::fsi(imMR,
                     nCores = cores,
                     method = 'sma',
                     bySpecies = TRUE,
                     grpBy = ECOSUBCD,
                     areaDomain = treat == 0,
                     scaleBy = c(FORTYPCD),
                     betas = full$betas)
data.table::fwrite(ecosub, './FSI_code/results/spatial/ecosub.csv')


## Division ---------------
## FSI by species, by division
ecodiv <- rFIA:::fsi(imMR,
                     nCores = cores,
                     method = 'sma',
                     bySpecies = TRUE,
                     grpBy = ECODIV,
                     areaDomain = treat == 0,
                     scaleBy = c(FORTYPCD),
                     betas = full$betas)
data.table::fwrite(ecodiv, './FSI_code/results/spatial/ecodiv.csv')




#===============================================================================
#########  Compute plot-level indices ------------------------------------------
#===============================================================================

## FSI by species by plot
speciesPlt <- fsi(im,
                  nCores = cores,
                  byPlot = TRUE,
                  bySpecies = TRUE,
                  scaleBy = c(FORTYPCD),
                  betas = full$betas)

## Remove treatment plots
speciesPlt <- filter(speciesPlt, !c(pltID %in% treat$pltID))

data.table::fwrite(speciesPlt, './FSI_code/results/plt/speciesPlt.csv')



