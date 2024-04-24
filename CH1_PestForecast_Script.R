### This is the code for Chapter 1 Model Runs
### Includes: (1) Model run function, (2) Model output save function

### Load libraries:
library(tidyverse)
library(dplyr)
library(readr)
library(rjags)
library(coda)
library(MCMCvis)
library(ecoforecastR)

### Set wd - for SCC runs:
setwd("/projectnb/dietzelab/malmborg/")

### Load data:
# conditionscores <- read.csv("2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv")
# score_sds <- read.csv("2022_08_31_DATAGRAB/2022_12_7_sample_score_stddev_5k.csv")
# dmr_data <- read.csv("CHAPTER_1/DATA/2023_12_DMR_DATA_TCG.csv")
# load("Chapter_1/2024_02_JAGS_models/dmls.RData")
# load("Chapter_1/2024_02_JAGS_models/dpls.RData")
conditionscores <- read.csv("Ch1_PestForecast/data/condition_scores.csv")
score_sds <- read.csv("Ch1_PestForecast/data/score_std_devs.csv")
dmr_data <- read.csv("Ch1_PestForecast/data/DMR_data.csv")
load("Ch1_PestForecast/data/2024_04_magls.RData")
load("Ch1_PestForecast/data/dpls.RData")

### JAGS models:
# spongy_disturb_a <- read_file("2024_04_06_Ch1_JAGS_MODEL_ALPHA_VERSION.txt")
# spongy_disturb_b <- read_file("2024_04_06_Ch1_JAGS_MODEL_BETA_VERSION.txt")
spongy_disturb_a <- read_file("Ch1_PestForecast/JAGS_models/2024_04_06_Ch1_JAGS_MODEL_ALPHA_VERSION.txt")
spongy_disturb_b <- read_file("Ch1_PestForecast/JAGS_models/2024_04_06_Ch1_JAGS_MODEL_BETA_VERSION.txt")
#spongy_disturb_joint <- read_file("")

### Function for model runs:
##' @param scores forest condition score data >> .csv
##' @param distyr disturbance year >> numeric
##' @param dmr disturbance prob/mag/recov dataset from calculator >> .csv
##' @param stan_devs forest condition score standard deviations >> .csv
##' @param dmls disturbance magnitude covariate list >> .RData list
##' @param dpls disturbance probability covariate list >> .RData list
##' @param modelrun_a which model is being run for dist prob >> numeric
##' @param modelrun_b which model is being run for dist mag >> numeric
##' @param model spongy_disturb version to use >> either alpha _a version or beta _b version
##' @param covs 1 = alpha, 2 = beta, 3 = joint >> numeric
##' @param vars variables for JAGS model >> vector of strings 
##' @param iters number of iterations for JAGS run >> numeric
##' @param thin thin used on JAGS run >> numeric
##' @param diters number of iterations for DIC sampler >> numeric

spongy_jags <- function(scores, distyr, dmr, stan_devs, 
                        dmls, dpls, modelrun_a, modelrun_b, model, 
                        covs, vars, iters, thin, diters){
  cs <- scores %>%
    # Drop unwanted columns
    dplyr::select(dplyr::starts_with("X")) %>%
    # Select just June averages for best disturbance visual without seasonality
    dplyr::select(dplyr::contains(".06.")) %>%
    # rename with date, without extra characters
    dplyr::rename_with(~ str_replace_all(., c("X|_score_mean|_cs_mean" = "", 
                                              "\\." = "-")))
  # collect previous timesteps for filling in x[s]s
  # isolate data we want
  prevtime <- scores %>%
    # Drop unwanted columns
    dplyr::select(dplyr::starts_with("X")) %>%
    # Select just June averages for best disturbance visual without seasonality
    dplyr::select(dplyr::contains(".05.")) %>%
    # rename with date, without extra characters
    dplyr::rename_with(~ str_replace_all(., c("X|_score_mean|_cs_mean" = "", 
                                              "\\." = "-")))
  # disturbance year
  distyear = distyr
  # condition score disturbance year onset data
  cs_model <- cs[grep(paste0("^", as.character(distyear)), names(cs))]
  
  # disturbance data for sample
  cs <- cs[dmr$X,]
  # get disturbance years, previous timestep, and condition score value at dist year for each:
  dists <- vector()
  cs_dists <- vector()
  ptmestp <- vector()
  ptime <- vector()
  for (i in 1:nrow(cs)){
    # get disturbance year colname
    if (dmr$dpy1[i] > 0) {
      dists[i] <- colnames(cs[grep(paste0("^", distyear, sep = ""), names(cs))])
      ptmestp[i] <- colnames(prevtime[grep(paste0("^", distyear, sep = ""), names(prevtime))])
    } else {
      dists[i] <- colnames(cs[grep(paste0("^",distyear+1, sep = ""), names(cs))])
      ptmestp[i] <- colnames(prevtime[grep(paste0("^", distyear+1, sep = ""), names(prevtime))])
    }
    # get condition score associated with disturbance year
    cs_dists[i] <- cs[i, dists[i]]
    # get previous time step score for correct disturbance year
    ptime[i] <- prevtime[i, ptmestp[i]]
  }
  
  ## for inputting into model, x initial condition
  xic <- ptime
  # make ptime NAs = 0
  xic[which(is.na(xic))] <- 0
  
  
  ### standard deviations for precisions
  sds <- stan_devs %>%
    # Drop unwanted columns
    dplyr::select(dplyr::starts_with("X")) %>%
    # Select just June averages for best disturbance visual without seasonality
    dplyr::select(dplyr::contains(".06.")) %>%
    # get rid of the other stuff
    dplyr::rename_with(~ str_replace_all(., c("X|_score_stddev" = "",
                                              "\\." = "-")))
  # previous year sds for initial conditions
  prevsds <- stan_devs %>%
    # Drop unwanted columns
    dplyr::select(dplyr::starts_with("X")) %>%
    # Select just June averages for best disturbance visual without seasonality
    dplyr::select(dplyr::contains(".05.")) %>%
    # get rid of the other stuff
    dplyr::rename_with(~ str_replace_all(., c("X|_score_stddev" = "",
                                              "\\." = "-")))
  
  sds <- sds[dmr$X,]
  prevsds <- prevsds[dmr$X,]
  
  ### convert to precisions: 
  # use disturbance years to get precisions at dist year for each
  cs_sds <- vector()
  for (i in 1:nrow(sds)){
    # get condition score sds associated with disturbance year
    cs_sds[i] <- sds[i, dists[i]]
  }
  # fill in NA values before precisions calculation
  # find missing June SD sites
  miss_SD <- which(is.na(cs_sds))
  # find missing dist years for SD NA's
  miss_dist <- dists[which(is.na(cs_sds))]
  # take 75%ile for SDs across sites
  sd_75s <- c(quantile(cs_sds[which(dists == unique(dists)[1])], 0.75, na.rm = T),
              quantile(cs_sds[which(dists == unique(dists)[2])], 0.75, na.rm = T))
  # fill in NAs based on disturbance year
  for (i in 1:length(miss_SD)){
    if(miss_dist[i] == unique(dists)[1]){
      cs_sds[miss_SD[i]] <- sd_75s[1]
    } else {
      cs_sds[miss_SD[i]] <- sd_75s[2]
    }
  }
  
  # do the same for prev year
  # use disturbance years to get precisions at dist year for each
  prev_sds <- vector()
  for (i in 1:nrow(prevsds)){
    # get condition score associated with disturbance year
    prev_sds[i] <- prevsds[i, ptmestp[i]]
  }
  # fill in NA values before precisions calculation
  # find missing previous time step SD sites
  miss_pSD <- which(is.na(prev_sds))
  # find missing dist years for SD NA's
  miss_prev <- ptmestp[which(is.na(prev_sds))]
  # take 75%ile for SDs across sites
  prevsd_75s <- c(quantile(prev_sds[which(ptmestp == unique(ptmestp)[1])], 0.75, na.rm = T),
                  quantile(prev_sds[which(ptmestp == unique(ptmestp)[2])], 0.75, na.rm = T))
  # fill in NAs based on previous year
  for (i in 1:length(miss_pSD)){
    if(miss_prev[i] == unique(ptmestp)[1]){
      prev_sds[miss_pSD[i]] <- prevsd_75s[1]
    } else {
      prev_sds[miss_pSD[i]] <- prevsd_75s[2]
    }
  }
  
  # convert SDs to precisions
  cs_precs <- vector()
  for (i in 1:length(cs_sds)){
    cs_precs[i] <- 1/(cs_sds[i]^2)
  }
  
  prev_precs <- vector()
  for (i in 1:length(prev_sds)){
    prev_precs[i] <- 1/(prev_sds[i]^2)
  }
  
  # tau initial condition for model data object
  tic <- prev_precs
  
  
  ### Setting an informative R prior
  # get full time series from score data
  cs_all <- scores %>%
    # Drop unwanted columns
    dplyr::select(dplyr::starts_with("X")) %>%
    # rename with date, without extra characters
    dplyr::rename_with(~ str_replace_all(., c("X|_score_mean|_cs_mean" = "", 
                                              "\\." = "-")))
  #remove missing value columns
  cs_all <- cs_all[dmr$X,]
  cs_model <- cs_model[dmr$X,]
  # number of sites
  nsites = nrow(cs_all)
  # months of time series using for cors (5 years)
  mos <- 25
  # correlations compute
  cors = apply(cs_all, 1, function(ts){cor(ts[1:(mos-1)], ts[2:mos], use="pairwise.complete.obs")})
  # mean R:
  R_mean <- mean(cors, na.rm = T)
  # var R:
  Rvar <- var(cors, na.rm = T)
  R_prec <- 1/(sd(cors, na.rm = T)^2)
  
  
  ### initial state of model parameters:
  ### beta version for disturbance magnitude
  # dmls matrix list 
  covls <- dmls
  # test data object
  dat <- list(y = cs_dists,
              dist = as.numeric(cs_dists < -1),
              covls = covls)
  
  # make empty list to fill in magnitude covariates
  beta.init <- list()
  # make a loop to fill in inits for each model run
  for (i in 1:length(covls)){
    # object for each covariate column for multivariate lm
    v <- c()
    for (j in 2:ncol(covls[[i]])){
      # make multivariate lm call
      v[j-1] <- paste0('covls[[i]][,',j,']')
    }
    # make lm call
    lm_formula <- as.formula(paste("y ~ ",
                                   paste(v, collapse='+')))
    # run lm
    beta.init[[i]] <- lm(lm_formula, data = dat)
  }
  
  
  ### alpha version for disturbance probability
  # test object of covariates
  covls <- dpls
  # make test data object
  dat <- list(y = cs_dists,
              dist = as.numeric(cs_dists < -1),
              covls = covls)
  
  alpha.init <- list()
  # make a loop to fill in inits for each model run
  for (i in 1:length(covls)){
    # object for each covariate column for multivariate lm
    v <- c()
    for (j in 2:ncol(covls[[i]])){
      # make multivariate lm call
      v[j-1] <- paste0('covls[[i]][,',j,']')
    }
    # make lm call
    glm_formula <- as.formula(paste("dist ~ ",
                                    paste(v, collapse='+')))
    # run lm
    alpha.init[[i]] <- glm(glm_formula, data = dat,
                           family = binomial(link = "logit"))
    
  }
  
  # SET INITIAL CONDITIONS:
  # model.run = numeric, which model are you running?
  # param = character, either "alpha" or "beta" for probs/mags param
  initer <- function(model.run, param){
    # make the param object to call the right list (alpha, beta...)
    param.obj = paste0(param, ".init")
    # intercept term
    int <- coef(get(param.obj)[[model.run]])[1]
    # params for the rest of the params
    param.inits <- coef(get(param.obj)[[model.run]])[-1]
    # combine them
    param.init <- c(int, param.inits)
    return(param.init)
  }
  
  # set model number for dist prob
  modelrun_a = modelrun_a
  # set model number for dist mag
  modelrun_b = modelrun_b
  # choose model for JAGS feed
  model = model
  #set either alpha or beta for filling in data lists and inits
  covs = covs
  
  ## inits for either alpha run or beta run:
  if (covs == 1){
    # make inits object for model input 
    init <- list(R = R_mean,
                 beta0 = initer(1, "beta")[1],
                 alpha = initer(modelrun_a, "alpha"))
    # prob covariates
    dpalpha <- dpls[[modelrun_a]]
    # missing data in covariates
    missing_a <- as.numeric(rownames(dpalpha[!complete.cases(dpalpha),]))
    
    # set data object
    data = list(y = cs_model[-missing_a], ns = nsites[-missing_a],
                x_ic = xic[-missing_a], tau_ic = tic[-missing_a],
                tau_obs = cs_precs[-missing_a],
                a = dpalpha[-missing_a,],
                a0 = rep(0,ncol(dpalpha)),
                Va = solve(diag(rep(1,ncol(dpalpha)))),
                rmean = R_mean, rprec = R_prec)
    
  } else if (covs == 2) { # if model = b
    # make inits object for model input
    init <- list(R = R_mean,
                 beta = initer(modelrun_b, "beta"),
                 alpha0 = initer(1, "alpha")[1])
    # mag covariates
    dmbeta <- dmls[[modelrun_b]]
    # missing data in covariates
    missing_b <- as.numeric(rownames(dmbeta[!complete.cases(dmbeta),]))
    
    # set data object
    data = list(y = cs_model[-missing_b], ns = nsites[-missing_b],
                x_ic = xic[-missing_b], tau_ic = tic[-missing_b],
                tau_obs = cs_precs[-missing_b],
                b = dmbeta[-missing_b,],
                b0 = rep(0,ncol(dmbeta)),
                Vb = solve(diag(rep(1,ncol(dmbeta)))),
                rmean = R_mean, rprec = R_prec)
    
  } else if (covs == 3) {
    init <- list(R = R_mean,
                 alpha = initer(modelrun_a, "alpha"),
                 beta = initer(modelrun_b, "beta"))
    # prob and mag covariates
    dpalpha <- dpls[[modelrun_a]]
    dmbeta <- dmls[[modelrun_b]]
    # missing data in covariates
    missing_a <- as.numeric(rownames(dpalpha[!complete.cases(dpalpha),]))
    missing_b <- as.numeric(rownames(dmbeta[!complete.cases(dmbeta),]))
    missing <- union(missing_a, missing_b)
    
    # set data object
    data = list(y = cs_model[-missing], ns = nsites[-missing],
                x_ic = xic[-missing], tau_ic = tic[-missing],
                tau_obs = cs_precs[-missing],
                a = dpalpha[-missing,],
                b = dmbeta[-missing,],
                a0 = rep(0, ncol(dpalpha)),
                Va = solve(diag(rep(1,ncol(dpalpha)))),
                b0 = rep(0,ncol(dmbeta)),
                Vb = solve(diag(rep(1,ncol(dmbeta)))),
                rmean = R_mean, rprec = R_prec)
  }
  
  ### RUN THE MODEL
  j.pests <- jags.model(file = textConnection(model),
                        data = data,
                        inits = init,
                        n.chains = 3)
  
  # get outputs
  jpout<-coda.samples(j.pests,
                      variable.names = vars,
                      n.iter = iters,
                      thin = thin)
  
  # run DIC
  DIC <- dic.samples(j.pests, n.iter = diters)
  sum <- sum(DIC$deviance, DIC$penalty)
  
  ### Make output list
  # track metadata
  metadata <- tibble::lst(modelrun_a, modelrun_b, model, covs, data)
  # model selection
  dic <- list(DIC, sum)
  # raw jags output
  # jags <- jpout
  # model output
  out <- as.matrix(jpout)
  # combine output
  output <- tibble::lst(metadata, dic, jpout, out)
  
  return(output)
}


### Function for saving model output:
##' @param jagsmodel output from spongy_jags call
model_save <- function(jagsmodel, vartype){
  # choose file path
  filepath_outputs <- "Ch1_PestForecast/model_outputs/"
  filepath_runs <- "Ch1_PestForecast/model_runs/"
  filepath_var <- c("beta_", "alpha_", "joint_")
  # date
  date <- as.character(Sys.Date())
  # make file name
  filename_outputs <- paste0(filepath_outputs, 
                             date, 
                             "_modelrun_", filepath_var[vartype],
                             "a_", as.character(jagsmodel$metadata$modelrun_a),
                             "b_", as.character(jagsmodel$metadata$modelrun_b),
                             "_output",".csv")
  filename_runs <- paste0(filepath_runs,
                          date,
                          "_modelrun_", filepath_var[vartype],
                          "a_", as.character(jagsmodel$metadata$modelrun_a),
                          "b_", as.character(jagsmodel$metadata$modelrun_b),
                          "_data",".RData")
  
  # save outputs to folder
  write.csv(jagsmodel$out, file = filename_outputs)
  # save model selection and metadata to folder
  model_info <- jagsmodel[c('jpout','dic','metadata')]
  save(model_info, file = filename_runs)
}


### Model Runs ###

# # variables want in outputs from jags model
# vars <- c("beta", "alpha0", "R", "pa0")
# iters = 200000
# thin = 10
# diters = 20000
# 
# # 2024-04-16
# beta_model_1 <- spongy_jags(condition_scores, 2016, dmr_data, score_sds, dmls, dpls,
#                             1, spongy_disturb_b, 2, vars, iters, thin, diters)
# model_save(beta_model_run)