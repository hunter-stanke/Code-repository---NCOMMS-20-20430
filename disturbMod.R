##=====================================================
##=====================================================
##
## This script contains all code necessary to construct
## the model of species-level response (FSI) to forest
## disturbance severity and and probability
## as outlined in NCOMMS-20-20430.
##
## Last modified: 30 July 2020 - Hunter Stanke
##
##====================================================
##====================================================



##===============================================================================
##  Set up your working directory (location where FSI_code is unzipped) ---------
##===============================================================================

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(R2jags)

## Set working directory
setwd('/home/hunter/Dropbox/gcb/natureComm/revision')


##===============================================================================
##  Read in pre-processed results (see 'estimate.R') ----------------------------
##===============================================================================

## Species plot-level FSI
species <- fread('./FSI_code/results/plt/speciesPlt.csv') %>%
  as.data.frame()

## Species abundance range-wide
speciesBAA <- fread('./FSI_code/results/regionWide/baaFull.csv') %>%
  as.data.frame()

## Plot-level disturbance codes
disturb <- fread('./FSI_code/results/plt/disturb.csv') %>%
  as.data.frame()



##===============================================================================
##  Prep data for model (build design matrix) -----------------------------------
##===============================================================================

## Id the top 8
topSpecies <- speciesBAA %>%
  filter(YEAR == 2018) %>%
  ## Removing gambel oak and one-seed juniper
  ## They are classifed as shrub-tree growth habits
  ## Doesnt really get at what we are interested in
  filter(SPCD %in% c(814,69) == FALSE) %>%
  top_n(8, TPA)

## Cleaning up names for common pinyon
species <- species %>%
  mutate(COMMON_NAME = case_when(
    COMMON_NAME == 'common or two-needle pinyon' ~ 'Common pinyon',
    TRUE ~ COMMON_NAME))

## Binary variable for disturbance
disturb <- disturb %>%
  filter(PLOT_STATUS_CD == 1) %>%
  select(-c(YEAR,DSTRBYR1, DSTRBYR2, DSTRBYR3, PLOT_STATUS_CD, nCond)) %>%
  pivot_longer(names_to = 'cutThis', values_to = 'disturb', cols = DSTRBCD1:DSTRBCD3 ) %>%
  filter(disturb > 0) %>%
  select(-c(cutThis)) %>%
  mutate(disturb = case_when(disturb %in% 10:19 ~ 'INSECT_RATE',
                             disturb %in% 20:29 ~ 'DISEASE_RATE',
                             disturb %in% 30:39 ~ 'FIRE_RATE')) %>%
  group_by(PLT_CN, disturb) %>%
  ## Will be converted to binary so value doesnt matter
  summarize(PERC_AREA = sum(PERC_AREA, na.rm = TRUE)) %>%
  mutate(PERC_AREA = if_else(PERC_AREA > 0, 1, 0)) %>%
  pivot_wider(id_cols = PLT_CN, values_from = PERC_AREA, names_from = disturb)


## Design matrix
si <- species %>%
  ## Removing non-forested plots (no BAA at previous or current measurement)
  ## PLOT_STATUS_CD != 1 indicates non-forested condition at both measurements
  filter(!is.na(FSI)) %>%
  filter(PLOT_STATUS_CD == 1) %>%
  mutate(cut = if_else(CURR_TPA == 0 & PREV_BAA == 0, 1, 0 )) %>%
  filter(cut < 1) %>%
  ## Top 8 species
  filter(SPCD %in% topSpecies$SPCD) %>%
  select(PLT_CN, PREV_PLT_CN, pltID, YEAR, REMPER, FSI, PREV_RD, COMMON_NAME, SPCD) %>%
  ## Disturbance data
  left_join(disturb, by = 'PLT_CN') %>%
  ## We want disturbance severity to be listed as 0 when no disturbance
  ## occurred, rather than NA
  mutate(FIRE_RATE = replace_na(FIRE_RATE, 0),
         INSECT_RATE = replace_na(INSECT_RATE, 0),
         DISEASE_RATE = replace_na(DISEASE_RATE, 0),
         PREV_RD = replace_na(PREV_RD, 0)) %>%
  ## Using a lognormal to model previous relative density
  ## Just log the variable here to make life easier
  ## in JAGS
  mutate(PREV_RD = log(PREV_RD + .0001))



##===========================================================================
## Prep data for model   ----------------------------------------------------
##===========================================================================

spPrep <- si %>%
  ungroup() %>%
  arrange(SPCD) %>%
  ## Need numeric group-level assignments
  mutate(g = as.numeric(as.factor(SPCD)))



## Matrix of observation-level covariates.
## t_rate, t_init, and the interaction
X = spPrep %>%
  ungroup() %>%
  select(FIRE_RATE, INSECT_RATE, DISEASE_RATE) %>%
  as.matrix()

## Set up in a list
data <- list(N = nrow(X), ## number of obs
             K = ncol(X), ## Number of individual-level predictors: disturbance rates
             J = length(unique(spPrep$g)), # Number of groups - species
             X = X, # individual level predictors: (N by K matrix)
             y = spPrep$FSI,
             z = spPrep$PREV_RD,
             remper = round(spPrep$REMPER),
             g = spPrep$g) # Numeric ID for species


##===========================================================================
## Build and run the model --------------------------------------------------
##===========================================================================

# Parameters to estimate
params <- c('alpha', 'beta', 'psi', 'sens', 'rd')

# MCMC settings
ni <- 5000
nc <- 4

# Start Gibbs sampling
jags_mod <- jags(data,
                 parameters.to.save=params,
                 model.file="./FSI_code/disturbSensitivityMod.jag",
                 n.chains=nc,
                 n.iter=ni)


##===========================================================================
## Process MCMC samples  ----------------------------------------------------
##===========================================================================

## Convert to mcmc list
jags_mcmc <- as.mcmc(jags_mod)

chains <- jags_mcmc
## Convert to data.frame
for (i in 1:length(jags_mcmc)){
  chains[[i]] <- as.data.frame(jags_mcmc[[i]])
  names(chains)[i] <- i
}
class(chains) <- 'list'


## Make it tidy
post <- bind_rows(chains) %>%
  pivot_longer(cols = everything(), names_to = 'var', values_to = 'estimate') %>%
  filter(str_detect(var, 'alpha|beta|psi|sens|rd')) %>%
  mutate(g = readr::parse_number(str_split(var, ',', simplify = TRUE)[,1]),
         varnum = readr::parse_number(str_split(var, ',', simplify = TRUE)[,2]),
         g = as.numeric(g),
         varnum = as.numeric(varnum),
         term = case_when(str_detect(var, 'beta') ~ 'beta',
                          str_detect(var, 'psi') ~ 'psi',
                          str_detect(var, 'sens') ~ 'sens',
                          str_detect(var, 'alpha') ~ 'alpha',
                          str_detect(var, 'rd') ~ 'rd'),
         disturb = case_when(varnum == 1 ~ 'fire',
                             varnum == 2 ~ 'bugs',
                             varnum == 3 ~ 'disease')) %>%
  left_join(distinct(spPrep, g, COMMON_NAME), by = 'g') %>%
  select(COMMON_NAME, term, disturb, estimate)


## Summarize
post_sum <- post %>%
  group_by(COMMON_NAME, term, disturb) %>%
  summarise(m = median(estimate),
            upper = quantile(estimate, probs = .975),
            lower = quantile(estimate, probs = .025),
            sig = if_else(upper * lower > 0, 1, 0),
            pd = if_else(m > 0,
                         length(estimate[estimate > 0]) / length(estimate),
                         length(estimate[estimate < 0]) / length(estimate)))


## Save results
fwrite(post, './FSI_code/results/mod/species_samples.csv')
fwrite(post_sum, './FSI_code/results/mod/species_summary.csv')


#
# ##===========================================================================
# ## Process MCMC samples  ----------------------------------------------------
# ##===========================================================================
#
# ## Convert to mcmc list
# jags_mcmc <- as.mcmc(jags_mod)
#
# chains <- jags_mcmc
# ## Convert to data.frame
# for (i in 1:length(jags_mcmc)){
#   chains[[i]] <- as.data.frame(jags_mcmc[[i]])
#   names(chains)[i] <- i
#
#   ## Add a column for each beta (severity) times frequency
#   for (j in 1:data$J){
#     ## Species-level
#     out <- chains[[i]][,str_detect(names(chains[[i]]), paste0(j,','))]
#
#     beta <- out[, str_detect(names(out), 'beta')]
#     psi <- out[, str_detect(names(out), 'psi')]
#
#
#     ## Sensitivity is severity times frequency
#     sens <- beta * psi
#     names(sens) <- str_replace(names(beta), 'beta', 'sens')
#
#     ## add it on
#     chains[[i]] <- bind_cols(chains[[i]], sens)
#   }
#
#   ## add it on
#   chains[[i]] <- bind_cols(chains[[i]], global_sens)
#
# }
# class(chains) <- 'list'
#
#
# ## Make it tidy
# post <- bind_rows(chains) %>%
#   pivot_longer(cols = everything(), names_to = 'var', values_to = 'estimate') %>%
#   filter(str_detect(var, 'global', negate = TRUE)) %>%
#   filter(str_detect(var, 'beta|psi|sens')) %>%
#   mutate(g = readr::parse_number(str_split(var, ',', simplify = TRUE)[,1]),
#          varnum = readr::parse_number(str_split(var, ',', simplify = TRUE)[,2]),
#          g = as.numeric(g),
#          varnum = as.numeric(varnum),
#          term = case_when(str_detect(var, 'beta') ~ 'beta',
#                           str_detect(var, 'psi') ~ 'psi',
#                           str_detect(var, 'sens') ~ 'sens'),
#          disturb = case_when(varnum == 1 ~ 'fire',
#                           varnum == 2 ~ 'bugs',
#                           varnum == 3 ~ 'disease')) %>%
#   left_join(distinct(spPrep, g, COMMON_NAME), by = 'g') %>%
#   select(COMMON_NAME, term, disturb, estimate)
#
# ## Make it tidy
# global_post <- bind_rows(chains) %>%
#   pivot_longer(cols = everything(), names_to = 'var', values_to = 'estimate') %>%
#   filter(str_detect(var, 'global')) %>%
#   filter(str_detect(var, 'beta|psi|sens')) %>%
#   mutate(varnum = readr::parse_number(str_split(var, ',', simplify = TRUE)[,1]),
#          varnum = as.numeric(varnum),
#          term = case_when(str_detect(var, 'beta') ~ 'beta',
#                           str_detect(var, 'psi') ~ 'psi',
#                           str_detect(var, 'sens') ~ 'sens'),
#          disturb = case_when(varnum == 1 ~ 'fire',
#                              varnum == 2 ~ 'bugs',
#                              varnum == 3 ~ 'disease')) %>%
#   select(term, disturb, estimate)
#
#
# ## Summarize
# post_sum <- post %>%
#   group_by(COMMON_NAME, term, disturb) %>%
#   summarise(m = median(estimate),
#             upper = quantile(estimate, probs = .975),
#             lower = quantile(estimate, probs = .025),
#             sig = if_else(upper * lower > 0, 1, 0),
#             pd = if_else(m > 0,
#                          length(estimate[estimate > 0]) / length(estimate),
#                          length(estimate[estimate < 0]) / length(estimate)))
#
# ## Summarize
# global_post_sum <- global_post %>%
#   group_by(term, disturb) %>%
#   summarise(m = median(estimate),
#             upper = quantile(estimate, probs = .975),
#             lower = quantile(estimate, probs = .025),
#             sig = if_else(upper * lower > 0, 1, 0),
#             pd = if_else(m > 0,
#                          length(estimate[estimate > 0]) / length(estimate),
#                          length(estimate[estimate < 0]) / length(estimate)))
#
#
# ## Save results
# fwrite(post, './FSI_code/results/mod/species_samples.csv')
# fwrite(post_sum, './FSI_code/results/mod/species_summary.csv')
# fwrite(global_post, './FSI_code/results/mod/global_samples.csv')
# fwrite(global_post_sum, './FSI_code/results/mod/global_summary.csv')
