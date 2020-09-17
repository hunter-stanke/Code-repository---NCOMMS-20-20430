##=====================================================
##=====================================================
##
## This script contains the source code for the 'fsi'
## function in the 'rFIA' R package. Reviewers should
## not source this file, but rather load the rFIA package
## to implement the functions listed herein (internal data
## dependencies are only met within rFIA namespace).
##
## This document is provided purely for reference in the
## review of NCOMMS-20-20430.
##
## Last modified: 28 May 2020 - Hunter Stanke
##
##====================================================
##====================================================


fsiStarter <- function(x,
                       db,
                       grpBy_quo = NULL,
                       scaleBy_quo = NULL,
                       polys = NULL,
                       returnSpatial = FALSE,
                       bySpecies = FALSE,
                       bySizeClass = FALSE,
                       landType = 'forest',
                       treeType = 'live',
                       method = 'sma',
                       lambda = .5,
                       treeDomain = NULL,
                       areaDomain = NULL,
                       totals = FALSE,
                       byPlot = FALSE,
                       useSeries = FALSE,
                       nCores = 1,
                       remote,
                       mr){
  
  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  
  if (remote){
    ## Store the original parameters here
    params <- db
    
    ## Read in one state at a time
    db <- readFIA(dir = db$dir, common = db$common,
                  tables = reqTables, states = x, ## x is the vector of state names
                  nCores = nCores)
    
    ## If a clip was specified, run it now
    if ('mostRecent' %in% names(params)){
      db <- clipFIA(db, mostRecent = params$mostRecent,
                    mask = params$mask, matchEval = params$matchEval,
                    evalid = params$evalid, designCD = params$designCD,
                    nCores = nCores)
    }
    
  } else {
    ## Really only want the required tables
    db <- db[names(db) %in% reqTables]
  }
  
  
  
  ## Need a plotCN, and a new ID
  db$PLOT <- db$PLOT %>% mutate(PLT_CN = CN,
                                pltID = paste(UNITCD, STATECD, COUNTYCD, PLOT, sep = '_'))
  db$TREE$TRE_CN <- db$TREE$CN
  ##  don't have to change original code
  #grpBy_quo <- enquo(grpBy)
  
  # Probably cheating, but it works
  if (quo_name(grpBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[10,] %>% # Just the first row
        left_join(select(db$COND, PLT_CN, names(db$COND)[names(db$COND) %in% names(db$PLOT) == FALSE]), by = 'PLT_CN') %>%
        inner_join(select(db$TREE, PLT_CN, names(db$TREE)[names(db$TREE) %in% c(names(db$PLOT), names(db$COND)) == FALSE]), by = 'PLT_CN') %>%
        select(!!grpBy_quo)
    )
    
    # If column doesnt exist, just returns 0, not a dataframe
    if (is.null(nrow(d_quo))){
      grpName <- quo_name(grpBy_quo)
      stop(paste('Columns', grpName, 'not found in PLOT, TREE, or COND tables. Did you accidentally quote the variables names? e.g. use grpBy = ECOSUBCD (correct) instead of grpBy = "ECOSUBCD". ', collapse = ', '))
    } else {
      # Convert to character
      grpBy <- names(d_quo)
    }
  } else {
    grpBy <- NULL
  }
  
  # Probably cheating, but it works
  if (quo_name(scaleBy_quo) != 'NULL'){
    ## Have to join tables to run select with this object type
    plt_quo <- filter(db$PLOT, !is.na(PLT_CN))
    ## We want a unique error message here to tell us when columns are not present in data
    d_quo <- tryCatch(
      error = function(cnd) {
        return(0)
      },
      plt_quo[10,] %>% # Just the first row
        left_join(select(db$COND, PLT_CN, names(db$COND)[names(db$COND) %in% names(db$PLOT) == FALSE]), by = 'PLT_CN') %>%
        select(!!scaleBy_quo)
    )
    
    # If column doesnt exist, just returns 0, not a dataframe
    if (is.null(nrow(d_quo))){
      grpName <- quo_name(grpBy_quo)
      stop(paste('Columns', grpName, 'not found in PLOT or COND tables. Did you accidentally quote the variables names? e.g. use scaleBy = ECOSUBCD (correct) instead of scaleBy = "ECOSUBCD". ', collapse = ', '))
    } else {
      # Convert to character
      scaleBy <- names(d_quo)
    }
  } else {
    scaleBy <- NULL
  }
  
  
  reqTables <- c('PLOT', 'TREE', 'COND', 'POP_PLOT_STRATUM_ASSGN', 'POP_ESTN_UNIT', 'POP_EVAL',
                 'POP_STRATUM', 'POP_EVAL_TYP', 'POP_EVAL_GRP')
  
  if (!is.null(polys) & first(class(polys)) %in% c('sf', 'SpatialPolygons', 'SpatialPolygonsDataFrame') == FALSE){
    stop('polys must be spatial polygons object of class sp or sf. ')
  }
  if (landType %in% c('timber', 'forest') == FALSE){
    stop('landType must be one of: "forest" or "timber".')
  }
  if (any(reqTables %in% names(db) == FALSE)){
    missT <- reqTables[reqTables %in% names(db) == FALSE]
    stop(paste('Tables', paste (as.character(missT), collapse = ', '), 'not found in object db.'))
  }
  if (str_to_upper(method) %in% c('TI', 'SMA', 'LMA', 'EMA', 'ANNUAL') == FALSE) {
    warning(paste('Method', method, 'unknown. Defaulting to Temporally Indifferent (TI).'))
  }
  
  # I like a unique ID for a plot through time
  if (byPlot) {grpBy <- c('pltID', grpBy)}
  # Save original grpBy for pretty return with spatial objects
  grpByOrig <- grpBy
  
  
  ### DEAL WITH TEXAS
  if (any(db$POP_EVAL$STATECD %in% 48)){
    ## Will require manual updates, fix your shit texas
    txIDS <- db$POP_EVAL %>%
      filter(STATECD %in% 48) %>%
      filter(END_INVYR < 2017) %>%
      filter(END_INVYR > 2006) %>%
      ## Removing any inventory that references east or west, sorry
      filter(str_detect(str_to_upper(EVAL_DESCR), 'EAST', negate = TRUE) &
               str_detect(str_to_upper(EVAL_DESCR), 'WEST', negate = TRUE))
    db$POP_EVAL <- bind_rows(filter(db$POP_EVAL, !(STATECD %in% 48)), txIDS)
  }
  
  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # # Convert polygons to an sf object
    # polys <- polys %>%
    #   as('sf')%>%
    #   mutate_if(is.factor,
    #             as.character)
    # ## A unique ID
    # polys$polyID <- 1:nrow(polys)
    #
    # # Add shapefile names to grpBy
    grpBy = c(grpBy, 'polyID')
    
    ## Make plot data spatial, projected same as polygon layer
    pltSF <- select(db$PLOT, c('LON', 'LAT', pltID)) %>%
      filter(!is.na(LAT) & !is.na(LON)) %>%
      distinct(pltID, .keep_all = TRUE)
    coordinates(pltSF) <- ~LON+LAT
    proj4string(pltSF) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    pltSF <- as(pltSF, 'sf') %>%
      st_transform(crs = st_crs(polys))
    
    ## Split up polys
    polyList <- split(polys, as.factor(polys$polyID))
    suppressWarnings({suppressMessages({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        cl <- makeCluster(nCores)
        clusterEvalQ(cl, {
          library(dplyr)
          library(stringr)
          library(rFIA)
          library(sf)
        })
        out <- parLapply(cl, X = names(polyList), fun = areal_par, pltSF, polyList)
        #stopCluster(cl) # Keep the cluster active for the next run
      } else { # Unix systems
        out <- mclapply(names(polyList), FUN = areal_par, pltSF, polyList, mc.cores = nCores)
      }
    })})
    pltSF <- bind_rows(out)
    
    # A warning
    if (length(unique(pltSF$pltID)) < 1){
      stop('No plots in db overlap with polys.')
    }
    ## Add polygon names to PLOT
    db$PLOT <- db$PLOT %>%
      left_join(select(pltSF, polyID, pltID), by = 'pltID')
    
    
    ## TO RETURN SPATIAL PLOTS
  }
  if (byPlot & returnSpatial){
    grpBy <- c(grpBy, 'LON', 'LAT')
  } # END AREAL
  
  ## Build domain indicator function which is 1 if observation meets criteria, and 0 otherwise
  # Land type domain indicator
  if (tolower(landType) == 'forest'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
  } else if (tolower(landType) == 'timber'){
    db$COND$landD <- ifelse(db$COND$COND_STATUS_CD == 1 & db$COND$SITECLCD %in% c(1, 2, 3, 4, 5, 6) & db$COND$RESERVCD == 0, 1, 0)
    # Tree Type domain indicator
    if (tolower(treeType) == 'live'){
      db$TREE$typeD <- ifelse(db$TREE$STATUSCD == 1, 1, 0)
    } else if (tolower(treeType) == 'gs'){
      db$TREE$typeD <- ifelse(db$TREE$DIA >= 5 & db$TREE$STATUSCD == 1, 1, 0)
    }
  }
  
  # update spatial domain indicator
  if(!is.null(polys)){
    db$PLOT$sp <- ifelse(db$PLOT$pltID %in% pltSF$pltID, 1, 0)
  } else {
    db$PLOT$sp <- 1
  }
  
  # User defined domain indicator for area (ex. specific forest type)
  pcEval <- left_join(db$PLOT, select(db$COND, -c('STATECD', 'UNITCD', 'COUNTYCD', 'INVYR', 'PLOT')), by = 'PLT_CN')
  #areaDomain <- substitute(areaDomain)
  pcEval$aD <- rlang::eval_tidy(areaDomain, pcEval) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(pcEval$aD)) pcEval$aD[is.na(pcEval$aD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(pcEval$aD)) pcEval$aD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  pcEval$aD <- as.numeric(pcEval$aD)
  db$COND <- left_join(db$COND, select(pcEval, c('PLT_CN', 'CONDID', 'aD')), by = c('PLT_CN', 'CONDID')) %>%
    mutate(aD_c = aD)
  aD_p <- pcEval %>%
    group_by(PLT_CN) %>%
    summarize(aD_p = as.numeric(any(aD > 0)))
  db$PLOT <- left_join(db$PLOT, aD_p, by = 'PLT_CN')
  rm(pcEval)
  
  # Same as above for tree (ex. trees > 20 ft tall)
  #treeDomain <- substitute(treeDomain)
  tD <- rlang::eval_tidy(treeDomain, db$TREE) ## LOGICAL, THIS IS THE DOMAIN INDICATOR
  if(!is.null(tD)) tD[is.na(tD)] <- 0 # Make NAs 0s. Causes bugs otherwise
  if(is.null(tD)) tD <- 1 # IF NULL IS GIVEN, THEN ALL VALUES TRUE
  db$TREE$tD <- as.numeric(tD)
  
  
  
  ## We only want inventory/ population info from t2 plots, but we need the plot tree cond data
  ## for t1 and t2
  remPlts <- db$PLOT %>%
    select(PLT_CN, PREV_PLT_CN, DESIGNCD, REMPER, PLOT_STATUS_CD) %>%
    ## Has to have a remeasurement, be in the current sample, and of the national design
    filter(!is.na(REMPER) & !is.na(PREV_PLT_CN) & PLOT_STATUS_CD != 3 & DESIGNCD == 1) %>%
    left_join(select(db$PLOT, PLT_CN, DESIGNCD, PLOT_STATUS_CD), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('', '.prev')) %>%
    ## past emasurement must be in the previous sample and of national design
    filter(PLOT_STATUS_CD.prev != 3 & DESIGNCD.prev == 1)
  
  ### Snag the EVALIDs that are needed
  db$POP_EVAL<- db$POP_EVAL %>%
    select('CN', 'END_INVYR', 'EVALID', 'ESTN_METHOD') %>%
    inner_join(select(db$POP_EVAL_TYP, c('EVAL_CN', 'EVAL_TYP')), by = c('CN' = 'EVAL_CN')) %>%
    filter(EVAL_TYP == 'EXPVOL') %>%
    filter(!is.na(END_INVYR) & !is.na(EVALID) & END_INVYR >= 2003) %>%
    distinct(END_INVYR, EVALID, .keep_all = TRUE)
  
  ## If a most-recent subset, make sure that we don't get two reporting years in
  ## western states
  if (mr) {
    db$POP_EVAL <- db$POP_EVAL %>%
      group_by(EVAL_TYP) %>%
      filter(END_INVYR == max(END_INVYR, na.rm = TRUE))
  }
  
  ## Make an annual panel ID, associated with an INVYR
  
  ### The population tables
  pops <- select(db$POP_EVAL, c('EVALID', 'ESTN_METHOD', 'CN', 'END_INVYR', 'EVAL_TYP')) %>%
    rename(EVAL_CN = CN) %>%
    left_join(select(db$POP_ESTN_UNIT, c('CN', 'EVAL_CN', 'AREA_USED', 'P1PNTCNT_EU')), by = c('EVAL_CN')) %>%
    rename(ESTN_UNIT_CN = CN) %>%
    left_join(select(db$POP_STRATUM, c('ESTN_UNIT_CN', 'EXPNS', 'P2POINTCNT', 'CN', 'P1POINTCNT', 'ADJ_FACTOR_SUBP', 'ADJ_FACTOR_MICR', "ADJ_FACTOR_MACR")), by = c('ESTN_UNIT_CN')) %>%
    rename(STRATUM_CN = CN) %>%
    left_join(select(db$POP_PLOT_STRATUM_ASSGN, c('STRATUM_CN', 'PLT_CN', 'INVYR', 'STATECD')), by = 'STRATUM_CN') %>%
    ## ONLY REMEASURED PLOTS MEETING CRITERIA ABOVE
    filter(PLT_CN %in% remPlts$PLT_CN) %>%
    ungroup() %>%
    mutate_if(is.factor,
              as.character)
  
  ### Which estimator to use?
  if (str_to_upper(method) %in% c('ANNUAL')){
    ## Want to use the year where plots are measured, no repeats
    ## Breaking this up into pre and post reporting becuase
    ## Estimation units get weird on us otherwise
    popOrig <- pops
    pops <- pops %>%
      group_by(STATECD) %>%
      filter(END_INVYR == INVYR) %>%
      ungroup()
    
    prePops <- popOrig %>%
      group_by(STATECD) %>%
      filter(INVYR < min(END_INVYR, na.rm = TRUE)) %>%
      distinct(PLT_CN, .keep_all = TRUE) %>%
      ungroup()
    
    pops <- bind_rows(pops, prePops) %>%
      mutate(YEAR = INVYR)
    
  } else {     # Otherwise temporally indifferent
    pops <- rename(pops, YEAR = END_INVYR)
  }
  
  ## P2POINTCNT column is NOT consistent for annnual estimates, plots
  ## within individual strata and est units are related to different INVYRs
  p2_INVYR <- pops %>%
    group_by(ESTN_UNIT_CN, STRATUM_CN, INVYR) %>%
    summarize(P2POINTCNT_INVYR = length(unique(PLT_CN)))
  ## Want a count of p2 points / eu, gets screwed up with grouping below
  p2eu_INVYR <- p2_INVYR %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, INVYR, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN, INVYR) %>%
    summarize(p2eu_INVYR = sum(P2POINTCNT_INVYR, na.rm = TRUE))
  p2eu <- pops %>%
    distinct(ESTN_UNIT_CN, STRATUM_CN, .keep_all = TRUE) %>%
    group_by(ESTN_UNIT_CN) %>%
    summarize(p2eu = sum(P2POINTCNT, na.rm = TRUE))
  
  ## Rejoin
  pops <- pops %>%
    left_join(p2_INVYR, by = c('ESTN_UNIT_CN', 'STRATUM_CN', 'INVYR')) %>%
    left_join(p2eu_INVYR, by = c('ESTN_UNIT_CN', 'INVYR')) %>%
    left_join(p2eu, by = 'ESTN_UNIT_CN')
  
  
  ## Recode a few of the estimation methods to make things easier below
  pops$ESTN_METHOD = recode(.x = pops$ESTN_METHOD,
                            `Post-Stratification` = 'strat',
                            `Stratified random sampling` = 'strat',
                            `Double sampling for stratification` = 'double',
                            `Simple random sampling` = 'simple',
                            `Subsampling units of unequal size` = 'simple')
  
  
  ## Add species to groups
  if (bySpecies) {
    db$TREE <- db$TREE %>%
      left_join(select(intData$REF_SPECIES_2018, c('SPCD','COMMON_NAME', 'GENUS', 'SPECIES')), by = 'SPCD') %>%
      mutate(SCIENTIFIC_NAME = paste(GENUS, SPECIES, sep = ' ')) %>%
      ungroup() %>%
      mutate_if(is.factor,
                as.character)
    grpBy <- c(grpBy, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
    grpByOrig <- c(grpByOrig, 'SPCD', 'COMMON_NAME', 'SCIENTIFIC_NAME')
  }
  
  ## Break into size classes
  if (bySizeClass){
    
    grpBy <- c(grpBy, 'sizeClass')
    grpByOrig <- c(grpByOrig, 'sizeClass')
    db$TREE$sizeClass <- makeClasses(db$TREE$DIA, interval = 2, numLabs = TRUE)
    db$TREE <- db$TREE[!is.na(db$TREE$sizeClass),]
  }
  
  
  
  ## Only the necessary plots for EVAL of interest
  remPltList <- unique(c(remPlts$PLT_CN, remPlts$PREV_PLT_CN))
  db$PLOT <- filter(db$PLOT, PLT_CN %in% remPltList)
  db$COND <- filter(db$COND, PLT_CN %in% remPltList)
  db$TREE <- filter(db$TREE, PLT_CN %in% remPltList)
  
  ## Tree basal area per acre
  db$TREE <- db$TREE %>%
    mutate(BAA = basalArea(DIA) * TPA_UNADJ)
  
  ## Narrow up the tables to the necessary variables, reduces memory load
  ## sent out the cores
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% c(grpBy, scaleBy)]
  grpC <- names(db$COND)[names(db$COND) %in% c(grpBy, scaleBy) & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% c(grpBy, scaleBy) & names(db$TREE) %in% c(grpP, grpC) == FALSE]
  
  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  db$PLOT <- select(db$PLOT, c('PLT_CN', pltID, 'REMPER', 'DESIGNCD', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR',
                               'MEASYEAR', 'MEASMON', 'MEASDAY', 'PLOT_STATUS_CD', PREV_PLT_CN, grpP, 'aD_p', 'sp'))
  db$COND <- select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD',
                               DSTRBCD1, DSTRBCD2, DSTRBCD3, TRTCD1, TRTCD2, TRTCD3))
  db$TREE <- select(db$TREE, c('PLT_CN', 'TRE_CN', 'CONDID', 'DIA', 'TPA_UNADJ', 'BAA', 'SUBP', 'TREE', grpT, 'tD', 'typeD',
                               PREVCOND, PREV_TRE_CN, STATUSCD, SPCD))
  
  
  ## Merging state and county codes
  plts <- split(db$PLOT, as.factor(paste(db$PLOT$COUNTYCD, db$PLOT$STATECD, sep = '_')))
  
  ## Summarize to plot-level
  suppressWarnings({
    ## Compute estimates in parallel -- Clusters in windows, forking otherwise
    if (Sys.info()['sysname'] == 'Windows'){
      cl <- makeCluster(nCores)
      clusterEvalQ(cl, {
        library(dplyr)
        library(stringr)
        library(rFIA)
        library(tidyr)
      })
      out <- parLapply(cl, X = names(plts), fun = fsiHelper1, plts, db, grpBy, scaleBy, byPlot)
      #stopCluster(cl) # Keep the cluster active for the next run
    } else { # Unix systems
      out <- mclapply(names(plts), FUN = fsiHelper1, plts, db, grpBy, scaleBy, byPlot, mc.cores = nCores)
    }
  })
  
  
  ## back to dataframes
  out <- unlist(out, recursive = FALSE)
  t <- bind_rows(out[names(out) == 't'])
  t1 <- bind_rows(out[names(out) == 't1'])
  a <- bind_rows(out[names(out) == 'a'])
  
  
  out <- list(t = t, t1 = t1, a = a, grpBy = grpBy, scaleBy = scaleBy, grpByOrig = grpByOrig, pops = pops)
  
}




fsi <- function(db,
                grpBy = NULL,
                polys = NULL,
                returnSpatial = FALSE,
                bySpecies = FALSE,
                bySizeClass = FALSE,
                landType = 'forest',
                treeType = 'live',
                method = 'sma',
                lambda = .5,
                treeDomain = NULL,
                areaDomain = NULL,
                totals = TRUE,
                variance = TRUE,
                byPlot = FALSE,
                useSeries = FALSE,
                scaleBy = NULL,
                betas = NULL,
                returnBetas = FALSE,
                nCores = 1) {
  
  ##  don't have to change original code
  grpBy_quo <- rlang::enquo(grpBy)
  scaleBy_quo <- rlang::enquo(scaleBy)
  areaDomain <- rlang::enquo(areaDomain)
  treeDomain <- rlang::enquo(treeDomain)
  
  ### Is DB remote?
  remote <- ifelse(class(db) == 'Remote.FIA.Database', 1, 0)
  if (remote){
    
    iter <- db$states
    
    ## In memory
  } else {
    ## Some warnings
    if (class(db) != "FIA.Database"){
      stop('db must be of class "FIA.Database". Use readFIA() to load your FIA data.')
    }
    
    ## an iterator for remote
    iter <- 1
    
  }
  
  ## Check for a most recent subset
  if (remote){
    if ('mostRecent' %in% names(db)){
      mr = db$mostRecent # logical
    } else {
      mr = FALSE
    }
    ## In-memory
  } else {
    if ('mostRecent' %in% names(db)){
      mr = TRUE
    } else {
      mr = FALSE
    }
  }
  
  ### AREAL SUMMARY PREP
  if(!is.null(polys)) {
    # Convert polygons to an sf object
    polys <- polys %>%
      as('sf') %>%
      mutate_if(is.factor,
                as.character)
    ## A unique ID
    polys$polyID <- 1:nrow(polys)
  }
  
  
  ## Run the main portion
  out <- lapply(X = iter, FUN = fsiStarter, db,
                grpBy_quo = grpBy_quo, scaleBy_quo, polys, returnSpatial,
                bySpecies, bySizeClass,
                landType, treeType, method,
                lambda, treeDomain, areaDomain,
                totals, byPlot, useSeries,
                nCores, remote, mr)
  
  ## Bring the results back
  out <- unlist(out, recursive = FALSE)
  t <- bind_rows(out[names(out) == 't'])
  t1 <- bind_rows(out[names(out) == 't1'])
  a <- bind_rows(out[names(out) == 'a'])
  grpBy <- out[names(out) == 'grpBy'][[1]]
  scaleBy <- out[names(out) == 'scaleBy'][[1]]
  grpByOrig <- out[names(out) == 'grpByOrig'][[1]]
  
  
  ## Prep the data for modeling the size-density curve
  scaleSyms <- syms(scaleBy)
  grpSyms <- syms(grpBy)
  
  ## Get groups prepped to fit model
  grpRates <- select(t1, PLT_CN, !!!scaleSyms, BA2, TPA2) %>%
    ungroup() %>%
    filter(TPA2 > 0) %>%
    tidyr::drop_na(!!!scaleSyms) %>%
    ## Stand-level variables here
    mutate(t = log(TPA2),
           b = log(BA2)) %>%
    select(t, b, PLT_CN, !!!scaleSyms)
  
  
  if (!is.null(scaleBy)){
    ## group IDS
    grpRates <- mutate(grpRates, grps = as.factor(paste(!!!scaleSyms)))
    t <- mutate(t, grps = as.factor(paste(!!!scaleSyms)))
    
    ## Remove groups with less than 10 obs
    nGrps <- grpRates %>%
      group_by(grps) %>%
      summarise(n = n())
    grpRates <- grpRates %>%
      left_join(nGrps, by = 'grps')
    
  } else {
    
    grpRates$grps <- 1
    t$grps = 1
    
  }
  
  
  if (is.null(betas)){
    ## Prep data for model
    prep <- grpRates %>%
      ungroup() %>%
      arrange(grps) %>%
      ## Need numeric group-level assignments
      mutate(grp_index = as.numeric(as.factor(grps)))
    
    
    
    
    ## If more than one group use a mixed model
    if (length(unique(grpRates$grps)) > 1){
      
      modFile <- system.file("extdata", "qrLMM.jag", package = "rFIA")
      # Parameters to estimate
      params <- c('fe_alpha', 'fe_beta', 'alpha', 'beta')
      ## Set up in a list
      data <- list(I = nrow(prep), ## number of obs
                   J = length(unique(prep$grp_index)), # Number of groups
                   y = prep$t, ## log scale TPA
                   x = prep$b, ## log scale BA
                   p = .99, ## Percentile for qr
                   grp_index = prep$grp_index) # Numeric ID for groups
      
    } else {
      
      modFile <- system.file("extdata", "qrLM.jag", package = "rFIA")
      # Parameters to estimate
      params <- c('alpha', 'beta')
      ## Set up in a list
      data <- list(I = nrow(prep), ## number of obs
                   y = prep$t, ## log scale TPA
                   x = prep$b, ## log scale BA
                   p = .99) ## Percentile for qr
      
    }
    
    # MCMC settings
    ni <- 1000
    nc <- 3
    
    print('Modeling maximum size-density curve(s)...')
    
    # Start Gibbs sampling
    jags_mod_start <- R2jags::jags(data,
                                   parameters.to.save=params,
                                   model.file=modFile,
                                   n.chains=nc,
                                   n.iter=ni)
    jags_mod <- R2jags::autojags(jags_mod_start, n.iter = 1000)
    
    
    ## Convert to mcmc list
    jags_mcmc <- coda::as.mcmc(jags_mod)
    
    chains <- jags_mcmc
    ## Convert to data.frame
    for (i in 1:length(jags_mcmc)){
      chains[[i]] <- as.data.frame(jags_mcmc[[i]])
      names(chains)[i] <- i
    }
    class(chains) <- 'list'
    
    
    if (length(unique(grpRates$grps)) > 1){
      ## Make it tidy
      betas <- bind_rows(chains) %>%
        pivot_longer(cols = everything(), names_to = 'var', values_to = 'estimate') %>%
        filter(str_detect(var, 'fe_beta|fe_alpha|deviance', negate = TRUE)) %>%
        mutate(grp_index = unlist(regmatches(var, gregexpr("\\[.+?\\]", var))),
               grp_index = as.numeric(str_sub(grp_index, 2, -2)),
               term = case_when(str_detect(var, 'alpha') ~ 'int',
                                TRUE ~ 'rate')) %>%
        left_join(distinct(prep, grp_index, grps, n), by = 'grp_index') %>%
        select(grps, term, estimate, n) %>%
        mutate(estimate = case_when(term == 'int' ~ exp(estimate),
                                    TRUE ~ estimate)) %>%
        group_by(grps, term) %>%
        summarise(mean = mean(estimate),
                  upper = quantile(estimate, probs = .975),
                  lower = quantile(estimate, probs = .025),
                  n = first(n)) %>%
        pivot_wider(id_cols = c(grps, n), names_from = term, values_from = mean:lower) %>%
        rename(int = mean_int,
               rate = mean_rate)
      
      ## Fixed effect for missing groups
      post_fe <- bind_rows(chains) %>%
        pivot_longer(cols = everything(), names_to = 'var', values_to = 'estimate') %>%
        filter(str_detect(var, 'fe_beta|fe_alpha')) %>%
        mutate(term = case_when(str_detect(var, 'alpha') ~ 'fe_int',
                                TRUE ~ 'fe_rate')) %>%
        select(term, estimate) %>%
        mutate(estimate = case_when(term == 'fe_int' ~ exp(estimate),
                                    TRUE ~ estimate)) %>%
        group_by(term) %>%
        summarise(mean = mean(estimate),
                  upper = quantile(estimate, probs = .975),
                  lower = quantile(estimate, probs = .025))%>%
        pivot_wider(names_from = term, values_from = mean:lower) %>%
        rename(fe_int = mean_fe_int,
               fe_rate = mean_fe_rate)
      
      ## Adding fixed effect info
      betas <- betas %>%
        mutate(fe_int = post_fe$fe_int,
               upper_fe_int = post_fe$upper_fe_int,
               lower_fe_int = post_fe$lower_fe_int,
               fe_rate = post_fe$fe_rate,
               upper_fe_rate = post_fe$upper_fe_rate,
               lower_fe_rate = post_fe$lower_fe_rate)
      
    } else {
      ## Make it tidy
      betas <- bind_rows(chains) %>%
        pivot_longer(cols = everything(), names_to = 'var', values_to = 'estimate') %>%
        filter(str_detect(var, 'deviance', negate = TRUE)) %>%
        mutate(term = case_when(str_detect(var, 'alpha') ~ 'int',
                                TRUE ~ 'rate')) %>%
        select(term, estimate) %>%
        mutate(estimate = case_when(term == 'int' ~ exp(estimate),
                                    TRUE ~ estimate)) %>%
        group_by(term) %>%
        summarise(mean = mean(estimate),
                  upper = quantile(estimate, probs = .975),
                  lower = quantile(estimate, probs = .025)) %>%
        pivot_wider(names_from = term, values_from = mean:lower) %>%
        rename(int = mean_int,
               rate = mean_rate)
      betas$n <- nrow(grpRates)
      betas$grps = 1
      
    }
  }
  
  
  ## If groups are missing, assume the fixed effects
  if ('fe_int' %in% names(betas)) {
    t <- t %>%
      left_join(select(betas, c(grps, int, fe_int, fe_rate, rate)), by = 'grps') %>%
      mutate(int = case_when(!is.na(int) ~ fe_int,
                             TRUE ~ int),
             rate = case_when(!is.na(rate) ~ fe_rate,
                              TRUE ~ rate)) %>%
      mutate(ba = BAA / TPA_UNADJ,
             tmax = int * (ba^rate),
             rd = TPA_UNADJ / tmax)
  } else {
    ## Add the betas onto t
    t <- t %>%
      left_join(select(betas, c(grps, int, rate)), by = 'grps') %>%
      mutate(ba = BAA / TPA_UNADJ,
             tmax = int * (ba^rate),
             rd = TPA_UNADJ / tmax)
  }
  
  
  
  
  
  if (byPlot){
    
    ## back to dataframes
    tOut <- t
    
    tOut <- tOut %>%
      ## Summing within scaleBy
      group_by(.dots = unique(c(grpBy[!c(grpBy %in% 'YEAR')], scaleBy)), YEAR, PLT_CN, PLOT_STATUS_CD, PREV_PLT_CN,
               REMPER) %>%
      summarize(PREV_RD = -sum(rd[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
                CURR_RD = sum(rd[ONEORTWO == 2] * tDI[ONEORTWO == 2], na.rm = TRUE),
                PREV_TPA = -sum(TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
                PREV_BAA = -sum(BAA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
                CURR_TPA = sum(TPA_UNADJ[ONEORTWO == 2] * tDI[ONEORTWO == 2], na.rm = TRUE),
                CURR_BAA = sum(BAA[ONEORTWO == 2] * tDI[ONEORTWO == 2], na.rm = TRUE)) %>%
      ## Summing across scaleBy
      group_by(.dots = grpBy[!c(grpBy %in% 'YEAR')], YEAR, PLT_CN, PLOT_STATUS_CD, PREV_PLT_CN,
               REMPER) %>%
      summarize(PREV_RD = mean(PREV_RD, na.rm = TRUE),
                CURR_RD = mean(CURR_RD, na.rm = TRUE),
                FSI = (CURR_RD - PREV_RD) / first(REMPER),
                PERC_FSI = FSI / PREV_RD * 100,
                PREV_TPA = sum(PREV_TPA, na.rm = TRUE),
                PREV_BAA = sum(PREV_BAA, na.rm = TRUE),
                CURR_TPA = sum(CURR_TPA, na.rm = TRUE),
                CURR_BAA = sum(CURR_BAA, na.rm = TRUE))
    
    ## If we want to use multiple remeasurements to estimate change,
    ## handle that here
    if (useSeries) {
      ## Get a unique ID for each remeasurement in the series
      nMeas <- tOut %>%
        distinct(pltID, PLT_CN, YEAR, REMPER) %>%
        group_by(pltID) %>%
        mutate(n = length(unique(PLT_CN)),
               series = min_rank(YEAR)) %>%
        ungroup() %>%
        select(pltID, PLT_CN, REMPER, n, series)
      
      ## Only if more than one remeasurement available
      if (any(nMeas$n > 1)){
        
        ## Now we loop over the unique values of n
        ## Basically have to chunk up the data each time
        ## in order to get intermediate estimates
        nRems <- unique(nMeas$n)
        remsList <- list()
        for (i in 1:length(nRems)){
          ## Temporal weights for each plot
          wgts <- nMeas %>%
            filter(series <= nRems[i] & n >= nRems[i]) %>%
            group_by(pltID) %>%
            ## Total remeasurement interval and weights for
            ## individual remeasurements
            mutate(fullRemp = sum(REMPER, na.rm = TRUE),
                   wgt = REMPER / fullRemp) %>%
            ungroup() %>%
            select(PLT_CN, n, series, wgt, fullRemp)
          
          dat <- tOut %>%
            left_join(wgts, by = c('PLT_CN')) %>%
            filter(series <= nRems[i] & n >= nRems[i]) %>%
            group_by(.dots = grpBy[grpBy %in% c('YEAR', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD') == FALSE]) %>%
            mutate(PLOT_STATUS_CD = case_when(any(PLOT_STATUS_CD == 1) ~ as.double(1),
                                              TRUE ~ as.double(PLOT_STATUS_CD))) %>%
            summarize(FSI = sum(FSI*wgt, na.rm = TRUE),
                      PLT_CN = PLT_CN[which.max(series)],
                      CURR_RD = CURR_RD[which.max(series)],
                      PREV_RD = PREV_RD[which.min(series)],
                      PREV_TPA = PREV_TPA[which.min(series)],
                      PREV_BAA = PREV_BAA[which.min(series)],
                      CURR_TPA = CURR_TPA[which.max(series)],
                      CURR_BAA = CURR_BAA[which.max(series)],
                      REMPER = first(fullRemp),
                      PLOT_STATUS_CD = first(PLOT_STATUS_CD)) %>%
            ungroup()# %>%
          # select(-c(pltID))
          remsList[[i]] <- dat
        }
        ## Bring it all back together
        dat <- bind_rows(remsList)
        
        ## Update columns in tEst
        tOut <- tOut %>%
          ungroup() %>%
          select(-c(PREV_RD:CURR_BAA, REMPER, PLOT_STATUS_CD)) %>%
          left_join(dat, by = c('PLT_CN', grpBy[!c(grpBy %in% c('YEAR', 'INVYR', 'MEASYEAR', 'PLOT_STATUS_CD'))])) %>%
          mutate(PERC_FSI = FSI / PREV_RD * 100)
      }
    }
    
    ## Make it spatial
    if (returnSpatial){
      tOut <- tOut %>%
        filter(!is.na(LAT) & !is.na(LON)) %>%
        st_as_sf(coords = c('LON', 'LAT'),
                 crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
      grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
      
    }
    
    
    tOut <- select(tOut, YEAR, PLT_CN, any_of('PREV_PLT_CN'), PLOT_STATUS_CD, grpBy[grpBy != 'YEAR'],
                   REMPER, FSI, PERC_FSI, PREV_RD, CURR_RD, PREV_TPA, CURR_TPA, PREV_BAA, CURR_BAA)
    
    
    
    ## Population estimation
  } else {
    ## back to dataframes
    pops <- bind_rows(out[names(out) == 'pops'])
    
    ## Adding YEAR to groups
    grpBy <- c('YEAR', grpBy)
    popState <- split(pops, as.factor(pops$STATECD))
    
    
    suppressWarnings({
      ## Compute estimates in parallel -- Clusters in windows, forking otherwise
      if (Sys.info()['sysname'] == 'Windows'){
        cl <- makeCluster(nCores)
        clusterEvalQ(cl, {
          library(dplyr)
          library(stringr)
          library(rFIA)
          library(tidyr)
        })
        out <- parLapply(cl, X = names(popState), fun = fsiHelper2, popState, t, a, grpBy, scaleBy, method, useSeries)
        stopCluster(cl)
      } else { # Unix systems
        out <- mclapply(names(popState), FUN = fsiHelper2, popState, t, a, grpBy, scaleBy, method, useSeries, mc.cores = nCores)
      }
    })
    ## back to dataframes
    out <- unlist(out, recursive = FALSE)
    tEst <- bind_rows(out[names(out) == 'tEst'])
    tEst <- ungroup(tEst)
    
    ##### ----------------- MOVING AVERAGES
    if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA')){
      ### ---- SIMPLE MOVING AVERAGE
      if (str_to_upper(method) == 'SMA'){
        ## Assuming a uniform weighting scheme
        wgts <- pops %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(wgt = 1 / length(unique(INVYR)))
        
        tEst <- left_join(tEst, wgts, by = 'ESTN_UNIT_CN')
        
        #### ----- Linear MOVING AVERAGE
      } else if (str_to_upper(method) == 'LMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))
        
        ## Want a number of INVYRs per EU
        neu <- wgts %>%
          group_by(ESTN_UNIT_CN) %>%
          summarize(n = sum(rank, na.rm = TRUE))
        
        ## Rejoining and computing wgts
        wgts <- wgts %>%
          left_join(neu, by = 'ESTN_UNIT_CN') %>%
          mutate(wgt = rank / n) %>%
          ungroup() %>%
          select(ESTN_UNIT_CN, INVYR, wgt)
        
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        
        #### ----- EXPONENTIAL MOVING AVERAGE
      } else if (str_to_upper(method) == 'EMA'){
        wgts <- pops %>%
          distinct(YEAR, ESTN_UNIT_CN, INVYR, .keep_all = TRUE) %>%
          arrange(YEAR, ESTN_UNIT_CN, INVYR) %>%
          group_by(as.factor(YEAR), as.factor(ESTN_UNIT_CN)) %>%
          mutate(rank = min_rank(INVYR))
        
        
        if (length(lambda) < 2){
          ## Want sum of weighitng functions
          neu <- wgts %>%
            mutate(l = lambda) %>%
            group_by(ESTN_UNIT_CN) %>%
            summarize(l = 1 - first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))
          
          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = 'ESTN_UNIT_CN') %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(ESTN_UNIT_CN, INVYR, wgt)
        } else {
          grpBy <- c('lambda', grpBy)
          ## Duplicate weights for each level of lambda
          yrWgts <- list()
          for (i in 1:length(unique(lambda))) {
            yrWgts[[i]] <- mutate(wgts, lambda = lambda[i])
          }
          wgts <- bind_rows(yrWgts)
          ## Want sum of weighitng functions
          neu <- wgts %>%
            group_by(lambda, ESTN_UNIT_CN) %>%
            summarize(l = 1 - first(lambda),
                      sumwgt = sum(l*(1-l)^(1-rank), na.rm = TRUE))
          
          ## Rejoining and computing wgts
          wgts <- wgts %>%
            left_join(neu, by = c('lambda', 'ESTN_UNIT_CN')) %>%
            mutate(wgt = l*(1-l)^(1-rank) / sumwgt) %>%
            ungroup() %>%
            select(lambda, ESTN_UNIT_CN, INVYR, wgt)
        }
        
        tEst <- left_join(tEst, wgts, by = c('ESTN_UNIT_CN', 'INVYR'))
        
      }
      
      ### Applying the weights
      tEst <- tEst %>%
        mutate_at(vars(ctEst:faEst), ~(.*wgt)) %>%
        mutate_at(vars(ctVar:cvEst_psi), ~(.*(wgt^2))) %>%
        group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
        summarize_at(vars(ctEst:plotIn_t), sum, na.rm = TRUE)
    }
    
    suppressMessages({suppressWarnings({
      ## If a clip was specified, handle the reporting years
      if (mr){
        ## If a most recent subset, ignore differences in reporting years across states
        ## instead combine most recent information from each state
        # ID mr years by group
        maxyearsT <- tEst %>%
          select(grpBy) %>%
          group_by(.dots = grpBy[!c(grpBy %in% 'YEAR')]) %>%
          summarise(YEAR = max(YEAR, na.rm = TRUE))
        
        # Combine estimates
        tEst <- tEst %>%
          ungroup() %>%
          select(-c(YEAR)) %>%
          left_join(maxyearsT, by = grpBy[!c(grpBy %in% 'YEAR')])
      }
    })})
    
    
    ##---------------------  TOTALS and RATIOS
    # Tree
    tTotal <- tEst %>%
      group_by(.dots = grpBy) %>%
      summarize_all(sum,na.rm = TRUE)
    
    
    
    ##---------------------  TOTALS and RATIOS
    suppressWarnings({
      tOut <- tTotal %>%
        group_by(.dots = grpBy) %>%
        summarize_all(sum,na.rm = TRUE) %>%
        mutate(TPA_RATE = ctEst / ptEst,
               BA_RATE = cbEst / pbEst,
               FSI = siEst / faEst,
               PERC_FSI = siEst / ra1Est,
               PREV_RD = ra1Est / faEst,
               CURR_RD = ra2Est / faEst,
               
               ## Ratio variance
               ctVar = (1/ptEst^2) * (ctVar + (TPA_RATE^2 * ptVar) - (2 * TPA_RATE * cvEst_ct)),
               cbVar = (1/pbEst^2) * (cbVar + (BA_RATE^2 * pbVar) - (2 * BA_RATE * cvEst_cb)),
               psiVar = (1/ra1Est^2) * (siVar + (PERC_FSI^2 * ra1Var) - (2 * PERC_FSI * cvEst_psi)),
               siVar = (1/faEst^2) * (siVar + (FSI^2 * faVar) - (2 * FSI * cvEst_si)),
               ra1Var = (1/faEst^2) * (ra1Var + (PREV_RD^2 * faVar) - (2 * PREV_RD * cvEst_ra1)),
               ra2Var = (1/faEst^2) * (ra2Var + (CURR_RD^2 * faVar) - (2 * CURR_RD * cvEst_ra2)),
               
               ## Make it a percent
               PERC_FSI = PERC_FSI * 100,
               psiVar = psiVar * (100^2),
               
               ## RATIO variance
               FSI_VAR = siVar,
               PERC_FSI_SE = sqrt(psiVar) / abs(PERC_FSI) * 100,
               PERC_FSI_VAR = psiVar,
               TPA_RATE_VAR = ctVar,
               BA_RATE_VAR = cbVar,
               PREV_RD_VAR = ra1Var,
               CURR_RD_VAR = ra2Var,
               
               nPlots = plotIn_t,
               N = nh,
               FSI_INT = qt(.975, df=N-1) * (sqrt(siVar)/sqrt(N)),
               PERC_FSI_INT = qt(.975, df=N-1) * (sqrt(psiVar)/sqrt(N))) %>%
        mutate(FSI_STATUS = case_when(
          FSI < 0 & FSI + FSI_INT < 0 ~ 'Decline',
          FSI < 0 & FSI + FSI_INT > 0 ~ 'Stable',
          FSI > 0 & FSI - FSI_INT > 0  ~ 'Expand',
          TRUE ~ 'Stable'
        ))
    })
    
    
    if (totals) {
      tOut <- tOut %>%
        select(grpBy, FSI, PERC_FSI, FSI_STATUS,
               FSI_INT, PERC_FSI_INT,
               PREV_RD, CURR_RD, TPA_RATE, BA_RATE,
               FSI_VAR, PERC_FSI_VAR, PREV_RD_VAR, CURR_RD_VAR,
               TPA_RATE_VAR, BA_RATE_VAR,
               nPlots, N)
      
    } else {
      tOut <- tOut %>%
        select(grpBy, FSI, PERC_FSI, FSI_STATUS,
               FSI_INT, PERC_FSI_INT,
               PREV_RD, CURR_RD, TPA_RATE, BA_RATE,
               FSI_VAR, PERC_FSI_VAR, PREV_RD_VAR, CURR_RD_VAR,
               TPA_RATE_VAR, BA_RATE_VAR,
               nPlots, N)
    }
    
    # Snag the names
    tNames <- names(tOut)[names(tOut) %in% grpBy == FALSE]
    
  }
  
  ## Pretty output
  tOut <- tOut %>%
    ungroup() %>%
    mutate_if(is.factor, as.character) %>%
    drop_na(grpBy) %>%
    arrange(YEAR) %>%
    as_tibble()
  
  # Return a spatial object
  if (!is.null(polys) & byPlot == FALSE) {
    ## NO IMPLICIT NA
    nospGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME') == FALSE])
    nospSym <- syms(nospGrp)
    tOut <- complete(tOut, !!!nospSym)
    ## If species, we don't want unique combos of variables related to same species
    ## but we do want NAs in polys where species are present
    if (length(nospGrp) < length(grpBy)){
      spGrp <- unique(grpBy[grpBy %in% c('SPCD', 'SYMBOL', 'COMMON_NAME', 'SCIENTIFIC_NAME')])
      spSym <- syms(spGrp)
      tOut <- complete(tOut, nesting(!!!nospSym))
    }
    
    suppressMessages({suppressWarnings({
      tOut <- left_join(tOut, polys, by = 'polyID') %>%
        select(c('YEAR', grpByOrig, tNames, names(polys))) %>%
        filter(!is.na(polyID))})})
    
    ## Makes it horrible to work with as a dataframe
    if (returnSpatial == FALSE) tOut <- select(tOut, -c(geometry))
  } else if (!is.null(polys) & byPlot){
    polys <- as.data.frame(polys)
    tOut <- left_join(tOut, select(polys, -c(geometry)), by = 'polyID')
  }
  
  ## For spatial plots
  if (returnSpatial & byPlot) grpBy <- grpBy[grpBy %in% c('LAT', 'LON') == FALSE]
  
  ## Above converts to tibble
  if (returnSpatial) tOut <- st_sf(tOut)
  # ## remove any duplicates in byPlot (artifact of END_INYR loop)
  if (byPlot) tOut <- unique(tOut)
  
  
  if (returnBetas) {
    tOut <- list(results = tOut, betas = betas)
  }
  
  return(tOut)
  
}





fsiHelper1 <- function(x, plts, db, grpBy, scaleBy, byPlot){
  
  ## Does not modify outside environment, just need scaleBy in here as well
  if (is.null(grpBy)){
    aGrps <- NULL
    grpBy <- scaleBy
  } else {
    aGrps <- grpBy
    grpBy <- unique(c(grpBy, scaleBy))
  }
  
  ## Selecting the plots for one county
  db$PLOT <- plts[[x]]
  
  
  ## Disturbance or treatment ever happen on the plot remeasurement? If so
  ## we want to cut it before we model the max size-density curve
  disturb <- select(db$PLOT, PLT_CN, pltID) %>%
    left_join(select(db$COND, PLT_CN, DSTRBCD1, TRTCD1), by = 'PLT_CN') %>%
    mutate(DSTRBCD1 = replace_na(DSTRBCD1, 0),
           TRTCD1 = replace_na(TRTCD1, 0)) %>%
    filter(DSTRBCD1 > 0) %>%
    ## Natural regen is ok
    filter(TRTCD1 > 0 & !c(TRTCD1 %in% 40)) %>%
    ## These plots were disturbed or treated
    distinct(PLT_CN, pltID)
  
  
  ## Which grpByNames are in which table? Helps us subset below
  grpP <- names(db$PLOT)[names(db$PLOT) %in% c(grpBy, scaleBy)]
  grpC <- names(db$COND)[names(db$COND) %in% c(grpBy, scaleBy) & names(db$COND) %in% grpP == FALSE]
  grpT <- names(db$TREE)[names(db$TREE) %in% c(grpBy, scaleBy) & names(db$TREE) %in% c(grpP, grpC) == FALSE]
  
  ## Making a treeID
  db$TREE$treID <- paste(db$TREE$SUBP, db$TREE$TREE, sep = '_')
  
  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  data <- db$PLOT %>%
    filter(DESIGNCD == 1 & PLOT_STATUS_CD != 3 & !is.na(REMPER) & !is.na(PREV_PLT_CN)) %>%
    left_join(db$COND, by = c('PLT_CN')) %>%
    left_join(db$TREE, by = c('PLT_CN', 'CONDID')) %>%
    left_join(select(db$PLOT, c('PLT_CN', 'sp', 'aD_p', 'DESIGNCD', 'PLOT_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN'), suffix = c('2', '1')) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDID', 'landD', 'aD_c', 'COND_STATUS_CD')), by = c('PREV_PLT_CN' = 'PLT_CN', 'PREVCOND' = 'CONDID'), suffix = c('2', '1')) %>%
    left_join(select(db$TREE, c('TRE_CN', all_of(grpT), treID, 'typeD', 'tD', 'TPA_UNADJ', 'BAA', 'DIA', 'STATUSCD', SPCD)), by = c('PREV_TRE_CN' = 'TRE_CN'), suffix = c('2', '1')) %>%
    filter(DESIGNCD1 == 1 & PLOT_STATUS_CD1 != 3) %>%
    mutate_if(is.factor,
              as.character)
  
  ## Comprehensive indicator function -- w/ growth accounting
  data$tDI2 <- data$landD2 * data$aD_p2 * data$aD_c2 * data$tD2 * data$typeD2 * data$sp2 *
    if_else(data$STATUSCD2 == 1, 1, 0)
  
  data$tDI1 <- data$landD1 * data$aD_p1 * data$aD_c1 * data$tD1 * data$typeD1 * data$sp1 *
    if_else(data$STATUSCD1 == 1, 1, 0)
  
  ## Comprehensive indicator function -- w/ growth accounting
  data$pDI2 <- data$landD2 * data$aD_p2 * data$aD_c2 * data$typeD2 * data$sp2 *
    if_else(data$STATUSCD2 == 1, 1, 0)
  
  data$pDI1 <- data$landD1 * data$aD_p1 * data$aD_c1 * data$typeD1 * data$sp1 *
    if_else(data$STATUSCD1 == 1, 1, 0)
  
  ## Save a copy for area calculations
  ### Only joining tables necessary to produce plot level estimates, adjusted for non-response
  aData <- select(db$PLOT, c('PLT_CN', 'PREV_PLT_CN', 'pltID', 'DESIGNCD', 'REMPER', 'STATECD', 'MACRO_BREAKPOINT_DIA', 'INVYR', 'MEASYEAR', 'MEASMON', 'MEASDAY', 'PLOT_STATUS_CD', all_of(grpP), 'aD_p', 'sp')) %>%
    filter(DESIGNCD == 1 & PLOT_STATUS_CD != 3) %>%
    left_join(select(db$COND, c('PLT_CN', 'CONDPROP_UNADJ', 'PROP_BASIS', 'COND_STATUS_CD', 'CONDID', grpC, 'aD_c', 'landD')), by = c('PLT_CN')) %>%
    mutate(aDI = landD * aD_p * aD_c * sp)
  
  
  ## PREVIOUS and CURRENT attributes
  data <- data %>%
    mutate(TPA_UNADJ1 = TPA_UNADJ1,
           TPA_UNADJ2 = TPA_UNADJ2,
           BAA1 = BAA1,
           BAA2 = BAA2,
           MORT = case_when(
             STATUSCD1 == 1 & STATUSCD2 == 2 ~ 1,
             STATUSCD1 == 1 & STATUSCD2 == 3 ~ 1,
             TRUE ~ 0),
           SURV = case_when(
             STATUSCD1 == 1 & STATUSCD2 == 1 ~ 1,
             TRUE ~ 0)
    )
  
  
  ## Just what we need
  data <- data %>%
    select(PLT_CN, PREV_PLT_CN, pltID, TRE_CN, SUBP, CONDID, TREE, CONDPROP_UNADJ,
           MEASYEAR, MACRO_BREAKPOINT_DIA, PROP_BASIS, grpP[grpP != 'PLOT_STATUS_CD'], grpC,
           REMPER, PLOT_STATUS_CD1, PLOT_STATUS_CD2,
           treID1, treID2,
           one_of(str_c(grpT,1),str_c(grpT,2)),
           tDI1, tDI2, pDI1, pDI2, STATUSCD1, STATUSCD2,
           DIA1, DIA2, BAA1, BAA2, TPA_UNADJ1, TPA_UNADJ2, SPCD1, SPCD2) %>%
    mutate(BAA1 = -(BAA1),
           TPA_UNADJ1 = -(TPA_UNADJ1)) %>%
    ## Rearrange previous values as observations
    pivot_longer(cols = -c(PLT_CN:REMPER),
                 names_to = c(".value", 'ONEORTWO'),
                 names_sep = -1) %>%
    mutate(PLOT_BASIS = case_when(
      ## When DIA is na, adjustment is NA
      is.na(DIA) ~ NA_character_,
      ## When DIA is less than 5", use microplot value
      DIA < 5 ~ 'MICR',
      ## When DIA is greater than 5", use subplot value
      DIA >= 5 & is.na(MACRO_BREAKPOINT_DIA) ~ 'SUBP',
      DIA >= 5 & DIA < MACRO_BREAKPOINT_DIA ~ 'SUBP',
      DIA >= MACRO_BREAKPOINT_DIA ~ 'MACR'))
  
  ## No zeros
  data <- data %>%
    mutate(TPA_UNADJ = replace_na(TPA_UNADJ, replace = 0),
           BAA = replace_na(BAA, replace = 0))
  
  ## Total trees for the size-density scaling
  t1 <- data %>%
    ## No disturbance/treatment plots
    filter(!c(pltID %in% disturb$pltID)) %>%
    distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
    filter(STATUSCD == 1) %>%
    filter(pDI == 1) %>%
    group_by(.dots = scaleBy, PLT_CN) %>%
    summarize(REMPER = first(REMPER),
              BAA1 = sum(-BAA[ONEORTWO == 1], na.rm = TRUE),
              TPA1 = sum(-TPA_UNADJ[ONEORTWO == 1], na.rm = TRUE),
              BAA2 = sum(BAA[ONEORTWO == 2], na.rm = TRUE),
              TPA2 = sum(TPA_UNADJ[ONEORTWO == 2], na.rm = TRUE),
              #times1 = round(TPA_UNADJ[ONEORTWO == 1]),
              skew1 = skewness(rep(DIA[ONEORTWO == 1], round(-TPA_UNADJ[ONEORTWO == 1]))),
              skew2 = skewness(rep(DIA[ONEORTWO == 2], round(TPA_UNADJ[ONEORTWO == 2])))
    ) %>%
    ## Mean BA
    mutate(BA1 = if_else(TPA1 != 0, BAA1 / TPA1, 0),
           BA2 = if_else(TPA2 != 0, BAA2 / TPA2, 0)) %>%
    ## Remove plots with high skewness
    filter(skew2 >= -1 & skew2 <= 1)
  
  
  if (byPlot){
    grpBy <- c('YEAR', grpBy)
    
    t <- data %>%
      mutate(YEAR = MEASYEAR) %>%
      distinct(PLT_CN, SUBP, TREE, ONEORTWO, .keep_all = TRUE) %>%
      select(all_of(grpBy), PLT_CN, PREV_PLT_CN, PLOT_STATUS_CD, REMPER, SUBP, TREE,
             ONEORTWO, tDI, TPA_UNADJ, BAA)
    # # Compute estimates at plot level
    # group_by(.dots = grpBy, PLT_CN, PREV_PLT_CN) %>%
    # summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
    #           REMPER = first(REMPER),
    #           PREV_BAA = sum(-BAA[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1& STATUSCD == 1], na.rm = TRUE),
    #           PREV_TPA = sum(-TPA_UNADJ[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
    #           CHNG_TPA = sum(TPA_UNADJ[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
    #           ## Sum here to avoid issues w/ zeros
    #           CHNG_BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
    #           PLOT_STATUS_CD = if_else(any(PLOT_STATUS_CD == 1), 1, 2)) %>%
    # ## Replace any NAs with zeros
    # mutate(PREV_BAA = replace_na(PREV_BAA, 0),
    #        PREV_TPA = replace_na(PREV_TPA, 0),
    #        CHNG_BAA = replace_na(CHNG_BAA, 0),
    #        CHNG_TPA = replace_na(CHNG_TPA, 0),
    #        ## T2 attributes
    #        CURR_BAA = PREV_BAA + CHNG_BAA,
    #        CURR_TPA = PREV_TPA + CHNG_TPA) %>%
    # ## Change in average tree BA and QMD
    # mutate(PREV_BA = if_else(PREV_TPA != 0, PREV_BAA / PREV_TPA, 0),
    #        CURR_BA = if_else(CURR_TPA != 0, CURR_BAA / CURR_TPA, 0),
    #        CHNG_BA = CURR_BA - PREV_BA,
    #        PREV_QMD = sqrt(PREV_BA / 0.005454154),
    #        CURR_QMD = sqrt(CURR_BA / 0.005454154),
    #        CHNG_QMD = CURR_QMD - PREV_QMD) %>%
    # ## All CHNG becomes an annual rate
    # mutate(CHNG_BAA = CHNG_BAA / REMPER,
    #        CHNG_TPA = CHNG_TPA / REMPER,
    #        CHNG_BA = CHNG_BA / REMPER,
    #        CHNG_QMD = CHNG_QMD / REMPER) %>%
    # #left_join(stems, by = c('PLT_CN', grpBy[!(grpBy %in% c('pltID', 'PLOT_STATUS_CD', 'YEAR'))])) %>%
    # ### Then divide by number of unique stems for an average
    # #mutate(CHNG_BAA = CHNG_BAA / n) %>%
    # ungroup() %>%
    # select(PLT_CN, PREV_PLT_CN, PLOT_STATUS_CD, REMPER, grpBy,
    #        PREV_TPA, PREV_BAA, PREV_BA, PREV_QMD,
    #        CHNG_TPA, CHNG_BAA, CHNG_BA, CHNG_QMD,
    #        CURR_TPA, CURR_BAA, CURR_BA, CURR_QMD,
    #        nLive)
    
    a = NULL
    
  } else {
    
    ### Plot-level estimates
    if (length(aGrps[aGrps %in% names(aData)]) < 1) {
      aGrps = NULL
    }  else {
      aGrps <- aGrps[aGrps %in% names(aData)]
    }
    
    aSyms <- syms(aGrps)
    
    
    ### Plot-level estimates
    a <- aData %>%
      ## date column
      mutate(date = paste(MEASYEAR, MEASMON, MEASDAY, sep = '-'),
             date = as.Date(date, "%Y-%m-%d")) %>%
      ## Will be lots of trees here, so CONDPROP listed multiple times
      ## Adding PROP_BASIS so we can handle adjustment factors at strata level
      distinct(PLT_CN, pltID, date, PROP_BASIS, CONDID, !!!aGrps, .keep_all = TRUE) %>%
      group_by(PLT_CN, pltID, date, PROP_BASIS, CONDID, .dots = aGrps) %>%
      summarize(CONDPROP_UNADJ = first(CONDPROP_UNADJ * aDI)) %>%
      mutate(CONDPROP_UNADJ = replace_na(CONDPROP_UNADJ, 0)) %>%
      ## Average forested area between min and max date
      group_by(pltID, PROP_BASIS, .dots = aGrps) %>%
      summarize(minDate = min(date, na.rm = TRUE),
                maxDate = max(date, na.rm = TRUE),
                amin = sum(CONDPROP_UNADJ[date == minDate], na.rm = TRUE),
                amax = sum(CONDPROP_UNADJ[date == maxDate], na.rm = TRUE),
                fa = (amin + amax) / 2) %>%
      left_join(select(ungroup(db$PLOT), PLT_CN, pltID), by = c('pltID'))
    
    ### Compute total TREES in domain of interest
    t <- data %>%
      distinct(PLT_CN, TRE_CN, ONEORTWO, .keep_all = TRUE) %>%
      select(all_of(grpBy), PLT_CN, pltID, PLOT_STATUS_CD, PLOT_BASIS, MEASYEAR, REMPER, TRE_CN,
             ONEORTWO, STATUSCD, tDI, TPA_UNADJ, BAA)
    # # Compute estimates at plot level
    # group_by(PLT_CN, PLOT_BASIS, .dots = grpBy) %>%
    # summarize(nLive = length(which(tDI[ONEORTWO == 1] > 0)), ## Number of live trees in domain of interest at previous measurement
    #           REMPER = first(REMPER),
    #           MEASYEAR = first(MEASYEAR),
    #           PREV_TPA = sum(-TPA_UNADJ[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
    #           PREV_BAA = sum(-BAA[ONEORTWO == 1 & STATUSCD == 1] * tDI[ONEORTWO == 1 & STATUSCD == 1], na.rm = TRUE),
    #           CHNG_TPA = sum(TPA_UNADJ[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
    #           CHNG_BAA = sum(BAA[STATUSCD == 1] * tDI[STATUSCD == 1], na.rm = TRUE),
    #           CURR_TPA = PREV_TPA + CHNG_TPA,
    #           CURR_BAA = PREV_BAA + CHNG_BAA,
    #           plotIn = if_else(sum(tDI, na.rm = TRUE) > 0, 1, 0))
  }
  
  pltOut <- list(t = t, a = a, t1 = t1)
  return(pltOut)
  
}






fsiHelper2 <- function(x, popState, t, a, grpBy, scaleBy, method, useSeries){
  
  ## DOES NOT MODIFY OUTSIDE ENVIRONMENT
  if (str_to_upper(method) %in% c("SMA", 'EMA', 'LMA', 'ANNUAL')) {
    grpBy <- c(grpBy, 'INVYR')
    popState[[x]]$P2POINTCNT <- popState[[x]]$P2POINTCNT_INVYR
    popState[[x]]$p2eu <- popState[[x]]$p2eu_INVYR
  }
  
  ######## ------------------ TREE ESTIMATES + CV
  aAdj <- a %>%
    ## Rejoin with population tables
    right_join(select(popState[[x]], -c(STATECD)), by = 'PLT_CN') %>%
    mutate(
      ## AREA
      aAdj = case_when(
        ## When NA, stay NA
        is.na(PROP_BASIS) ~ NA_real_,
        ## If the proportion was measured for a macroplot,
        ## use the macroplot value
        PROP_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
        ## Otherwise, use the subpplot value
        PROP_BASIS == 'SUBP' ~ ADJ_FACTOR_SUBP),
      fa = fa * aAdj) %>%
    ungroup()
  
  ## Sometimes less specific area groups
  aGrps <- unique(grpBy[grpBy %in% names(aAdj)])
  
  
  ## Strata level estimates
  tEst <- t %>%
    ungroup() %>%
    ## Converting to average tree size
    mutate(BA = if_else(ONEORTWO == 1, -BAA / TPA_UNADJ, BAA / TPA_UNADJ)) %>%
    ## Summing within scaleBy
    group_by(.dots = unique(c(grpBy[!c(grpBy %in% c('YEAR', 'INVYR'))], scaleBy)), PLT_CN, pltID, MEASYEAR, REMPER, PLOT_BASIS) %>%
    summarize(PREV_RD = -sum(rd[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
              CURR_RD = sum(rd[ONEORTWO == 2] * tDI[ONEORTWO == 2], na.rm = TRUE),
              PREV_TPA = -sum(TPA_UNADJ[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
              PREV_BA = -sum(BA[ONEORTWO == 1] * tDI[ONEORTWO == 1], na.rm = TRUE),
              CHNG_TPA = sum(TPA_UNADJ * tDI, na.rm = TRUE) / first(REMPER),
              CHNG_BA = sum(BA * tDI, na.rm = TRUE) / first(REMPER),
              plotIn_t = if_else(sum(tDI, na.rm = TRUE) > 0, 1, 0)) %>%
    ## Summing across scaleBy
    group_by(PLT_CN, pltID, MEASYEAR, PLOT_BASIS,
             REMPER, .dots = grpBy[!c(grpBy %in% c('YEAR', 'INVYR'))]) %>%
    summarize(CHNG_TPA = sum(CHNG_TPA, na.rm = TRUE),
              CHNG_BA = sum(CHNG_BA, na.rm = TRUE),
              PREV_TPA = sum(PREV_TPA, na.rm = TRUE),
              PREV_BA = sum(PREV_BA, na.rm = TRUE),
              PREV_RD = mean(PREV_RD, na.rm = TRUE),
              CURR_RD = mean(CURR_RD, na.rm = TRUE),
              plotIn_t = ifelse(sum(plotIn_t, na.rm = TRUE) >  0, 1,0)) %>%
    mutate(FSI = (CURR_RD - PREV_RD) / REMPER) %>%
    ungroup()
  
  ## If we want to use multiple remeasurements to estimate change,
  ## handle that here
  if (useSeries) {
    ## Get a unique ID for each remeasurement in the series
    nMeas <- t %>%
      distinct(pltID, PLT_CN, MEASYEAR, REMPER) %>%
      group_by(pltID) %>%
      mutate(n = length(unique(PLT_CN)),
             series = min_rank(MEASYEAR)) %>%
      ungroup() %>%
      select(pltID, PLT_CN, REMPER, n, series)
    
    ## Only if more than one remeasurement available
    if (any(nMeas$n > 1)){
      
      ## Now we loop over the unique values of n
      ## Basically have to chunk up the data each time
      ## in order to get intermediate estimates
      nRems <- unique(nMeas$n)
      remsList <- list()
      for (i in 1:length(nRems)){
        ## Temporal weights for each plot
        wgts <- nMeas %>%
          filter(series <= nRems[i] & n >= nRems[i]) %>%
          group_by(pltID) %>%
          ## Total remeasurement interval and weights for
          ## individual remeasurements
          mutate(fullRemp = sum(REMPER, na.rm = TRUE),
                 wgt = REMPER / fullRemp) %>%
          ungroup() %>%
          select(PLT_CN, n, series, wgt)
        
        dat <- tEst %>%
          left_join(wgts, by = c('PLT_CN')) %>%
          filter(series <= nRems[i] & n >= nRems[i]) %>%
          group_by(pltID, PLOT_BASIS, .dots = grpBy[!c(grpBy %in% c('YEAR', 'INVYR'))]) %>%
          summarize(FSI = sum(FSI*wgt, na.rm = TRUE),
                    CHNG_TPA = sum(CHNG_TPA*wgt, na.rm = TRUE),
                    CHNG_BA = sum(CHNG_BA*wgt, na.rm = TRUE),
                    PLT_CN = PLT_CN[which.max(series)],
                    CURR_RD = CURR_RD[which.max(series)],
                    PREV_RD = PREV_RD[which.min(series)],
                    PREV_TPA = PREV_TPA[which.min(series)],
                    PREV_BA = PREV_BA[which.min(series)],
                    plotIn_t = if_else(any(plotIn_t > 0), 1, 0)) %>%
          ungroup() %>%
          select(-c(pltID))
        remsList[[i]] <- dat
      }
      ## Bring it all back together
      dat <- bind_rows(remsList)
      
      ## Update columns in tEst
      tEst <- tEst %>%
        select(-c(CHNG_TPA:FSI)) %>%
        left_join(dat, by = c('PLT_CN', 'PLOT_BASIS', grpBy[!c(grpBy %in% c('YEAR', 'INVYR'))]))
    }
  }
  
  ## Now go to strata level and onward
  tEst <- tEst %>%
    ## Rejoin with population tables
    right_join(select(ungroup(popState[[x]]), -c(STATECD)), by = 'PLT_CN') %>%
    ungroup() %>%
    ## Need forest area to adjust SI indices
    #left_join(select(ungroup(aAdj), PLT_CN, aGrps, fa, aAdj), by = c('PLT_CN', aGrps)) %>%
    #Add adjustment factors
    mutate(tAdj = case_when(
      ## When NA, stay NA
      is.na(PLOT_BASIS) ~ NA_real_,
      ## If the proportion was measured for a macroplot,
      ## use the macroplot value
      PLOT_BASIS == 'MACR' ~ as.numeric(ADJ_FACTOR_MACR),
      ## Otherwise, use the subpplot value
      PLOT_BASIS == 'SUBP' ~ as.numeric(ADJ_FACTOR_SUBP),
      PLOT_BASIS == 'MICR' ~ as.numeric(ADJ_FACTOR_MICR)),
      CHNG_TPA = CHNG_TPA * tAdj,
      CHNG_BA = CHNG_BA * tAdj,
      CURR_RD = CURR_RD * tAdj,
      PREV_TPA = PREV_TPA * tAdj,
      PREV_BA = PREV_BA * tAdj,
      PREV_RD = PREV_RD * tAdj,
      FSI = FSI * tAdj) %>%
    ## Extra step for variance issues - summing micro, subp, and macr components
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, PLT_CN, .dots = grpBy) %>%
    summarize(CHNG_TPA = sum(CHNG_TPA, na.rm = TRUE),
              CHNG_BA = sum(CHNG_BA, na.rm = TRUE),
              PREV_TPA = sum(PREV_TPA, na.rm = TRUE),
              PREV_BA = sum(PREV_BA, na.rm = TRUE),
              PREV_RD = sum(PREV_RD, na.rm = TRUE),
              CURR_RD = sum(CURR_RD, na.rm = TRUE),
              FSI = sum(FSI, na.rm = TRUE),
              plotIn_t = ifelse(sum(plotIn_t, na.rm = TRUE) >  0, 1,0),
              nh = first(P2POINTCNT),
              p2eu = first(p2eu),
              a = first(AREA_USED),
              w = first(P1POINTCNT) / first(P1PNTCNT_EU)) %>%
    ## Add on area
    left_join(select(aAdj, ESTN_UNIT_CN, STRATUM_CN, PLT_CN, aGrps, fa),
              by = c('ESTN_UNIT_CN', 'STRATUM_CN', 'PLT_CN', aGrps)) %>%
    ## FSI is area adjusted
    mutate(si = FSI * fa,
           ra1 = PREV_RD * fa,
           ra2 = CURR_RD * fa) %>%
    ## Replace NAN w/ zeros
    mutate(si = replace_na(si, 0),
           ra1 = replace_na(ra1, 0),
           ra2 = replace_na(ra2, 0)) %>%
    ungroup() %>%
    ## Strata-level
    group_by(ESTN_UNIT_CN, ESTN_METHOD, STRATUM_CN, .dots = grpBy) %>%
    summarize(r_t = length(unique(PLT_CN)) / first(nh),
              ctStrat = mean(CHNG_TPA * r_t, na.rm = TRUE),
              cbStrat = mean(CHNG_BA * r_t, na.rm = TRUE),
              ptStrat = mean(PREV_TPA * r_t, na.rm = TRUE),
              pbStrat = mean(PREV_BA * r_t, na.rm = TRUE),
              siStrat = mean(si * r_t, na.rm = TRUE),
              ra1Strat = mean(ra1 * r_t, na.rm = TRUE),
              ra2Strat = mean(ra2 * r_t, na.rm = TRUE),
              faStrat = mean(fa * r_t, na.rm = TRUE),
              plotIn_t = sum(plotIn_t, na.rm = TRUE),
              n = n(),
              ## We don't want a vector of these values, since they are repeated
              nh = first(nh),
              a = first(a),
              w = first(w),
              p2eu = first(p2eu),
              ndif = nh - n,
              # ## Strata level variances
              ctv = stratVar(ESTN_METHOD, CHNG_TPA, ctStrat, ndif, a, nh),
              cbv = stratVar(ESTN_METHOD, CHNG_BA, cbStrat, ndif, a, nh),
              ptv = stratVar(ESTN_METHOD, PREV_TPA, ptStrat, ndif, a, nh),
              pbv = stratVar(ESTN_METHOD, PREV_BA, pbStrat, ndif, a, nh),
              siv = stratVar(ESTN_METHOD, si, siStrat, ndif, a, nh),
              ra1v = stratVar(ESTN_METHOD, ra1, ra1Strat, ndif, a, nh),
              ra2v = stratVar(ESTN_METHOD, ra2, ra2Strat, ndif, a, nh),
              fav = stratVar(ESTN_METHOD, fa, faStrat, ndif, a, nh),
              
              # Strata level covariances
              cvStrat_ct = stratVar(ESTN_METHOD, CHNG_TPA, ctStrat, ndif, a, nh, PREV_TPA, ptStrat),
              cvStrat_cb = stratVar(ESTN_METHOD, CHNG_BA, cbStrat, ndif, a, nh, PREV_BA, pbStrat),
              cvStrat_si = stratVar(ESTN_METHOD, si, siStrat, ndif, a, nh, fa, faStrat),
              cvStrat_ra1 = stratVar(ESTN_METHOD, ra1, ra1Strat, ndif, a, nh, fa, faStrat),
              cvStrat_ra2 = stratVar(ESTN_METHOD, ra2, ra2Strat, ndif, a, nh, fa, faStrat),
              cvStrat_psi = stratVar(ESTN_METHOD, si, siStrat, ndif, a, nh, ra1, ra1Strat)) %>%
    
    ## Estimation unit
    group_by(ESTN_UNIT_CN, .dots = grpBy) %>%
    summarize(ctEst = unitMean(ESTN_METHOD, a, nh, w, ctStrat),
              cbEst = unitMean(ESTN_METHOD, a, nh, w, cbStrat),
              ptEst = unitMean(ESTN_METHOD, a, nh, w, ptStrat),
              pbEst = unitMean(ESTN_METHOD, a, nh, w, pbStrat),
              siEst = unitMean(ESTN_METHOD, a, nh, w, siStrat),
              ra1Est = unitMean(ESTN_METHOD, a, nh, w, ra1Strat),
              ra2Est = unitMean(ESTN_METHOD, a, nh, w, ra2Strat),
              faEst = unitMean(ESTN_METHOD, a, nh, w, faStrat),
              
              nh = first(nh),
              # Estimation of unit variance
              ctVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ctv, ctStrat, ctEst),
              cbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, cbv, cbStrat, cbEst),
              ptVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ptv, ptStrat, ptEst),
              pbVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, pbv, pbStrat, pbEst),
              siVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, siv, siStrat, siEst),
              ra1Var = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ra1v, ra1Strat, ra1Est),
              ra2Var = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, ra2v, ra2Strat, ra2Est),
              faVar = unitVarNew(method = 'var', ESTN_METHOD, a, nh, first(p2eu), w, fav, faStrat, faEst),
              
              ## Covariances
              cvEst_ct = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ct, ctStrat, ctEst, ptStrat, ptEst),
              cvEst_cb = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_cb, cbStrat, cbEst, pbStrat, pbEst),
              cvEst_si = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_si, siStrat, siEst, faStrat, faEst),
              cvEst_ra1 = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ra1, ra1Strat, ra1Est, faStrat, faEst),
              cvEst_ra2 = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_ra2, ra2Strat, ra2Est, faStrat, faEst),
              cvEst_psi = unitVarNew(method = 'cov', ESTN_METHOD, a, nh, first(p2eu), w, cvStrat_psi, siStrat, siEst, ra1Strat, ra1Est),
              
              plotIn_t = sum(plotIn_t, na.rm = TRUE)) %>%
    ungroup()
  
  out <- list(tEst = tEst, aEst = NULL)
  
  return(out)
}
