# This script is the JAGS model for pest disturbance mag and prob forecast
# This is the version that is being used for the SCC - DISTURBANCE PROBABILITY ("Alpha") VERSIONS

### THE LIBRARIES:
library(tidyverse)
library(dplyr)
library(rjags)
library(coda)
library(MCMCvis)
library(ecoforecastR)

### THE DATA:
# condition scores inputs
cfile <- "data/condition_scores.csv" #1995-2020, original grab
scores <- read.csv(cfile)
# isolate data we want
cs <- scores %>%
  # Drop unwanted columns
  dplyr::select(dplyr::starts_with("X")) %>%
  # Select just June averages for best disturbance visual without seasonality
  dplyr::select(dplyr::contains(".06.")) %>%
  # rename with date, without extra characters
  dplyr::rename_with(~ str_replace_all(., c("X|_score_mean|_cs_mean" = "", 
                                            "\\." = "-")))

# collect previous timesteps for filling in x[s]'s
# isolate data we want
prevtime <- scores %>%
  # Drop unwanted columns
  dplyr::select(dplyr::starts_with("X")) %>%
  # Select just June averages for best disturbance visual without seasonality
  dplyr::select(dplyr::contains(".05.")) %>%
  # rename with date, without extra characters
  dplyr::rename_with(~ str_replace_all(., c("X|_score_mean|_cs_mean" = "", 
                                            "\\." = "-")))

# number of sites
nsites = nrow(cs)

# disturbance year data
dmr_file <- "data/DMR_data.csv"
dmr <- read.csv(dmr_file)

# choose disturbance onset year
distyear <- 2016

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
sdfile <- "data/score_std_devs.csv"
stan_devs <- read.csv(sdfile)
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
# months of time series using for cors
mos <- 25
# correlations compute
cors = apply(cs_all, 1, function(ts){cor(ts[1:(mos-1)], ts[2:mos], use="pairwise.complete.obs")})
# mean R:
R_mean <- mean(cors, na.rm = T)
# var R:
Rvar <- var(cors, na.rm = T)
R_prec <- 1/(sd(cors, na.rm = T)^2)

## Loading calculated anomaly versions of covariate data
# disturbance magnitude
load("data/dmls.RData")
# disturbance probability
load("data/dpls.RData")


### THE MODEL:
#use the single time step version of the model:
spongy_disturb <- "model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model:
  y[s] ~ dnorm(mu[s], tau_obs[s])

  #### Process Model:
  muN[s] ~ dnorm(mun[s], pan)
  muD[s] ~ dnorm(mu0[s], pa0) 
  logit(D[s]) <- inprod(alpha[], a[s,])
  mu[s] <- D[s] * muD[s] + (1-D[s]) * muN[s]
  mu0[s] <- beta0
  mun[s] <- R * x[s]

  x[s] ~ dnorm(x_ic[s], tau_ic[s])          
  
}#end loop over sites

  #### Priors
  R ~ dnorm(rmean, rprec)              ##rho paramter (recovery rate)
  pa0 ~ dgamma(1,1)                    ##precision of disturbed state
  pan ~ dgamma(1,1)                    ##precision of undisturbed state
  beta0 ~ dnorm(0, 0.0001)            ##param for dist prob intercept

  ## covariate matrix: ADD WHEN READY FOR MODEL RUNS - 4/4/2024
  ##beta ~ dmnorm(b0, Vb)   ## for disturbance magnitude
  alpha ~ dmnorm(a0, Va)  ## for disturbance probability

}
"

### initial state of model parameters:
### beta version for disturbance magnitude
# test the dmls matrix list 
test <- dmls
# test data object
data <- list(y = cs_dists,
             dist = as.numeric(cs_dists < -1),
             test = test)

# make empty list to fill in magnitude covariates
beta.init <- list()
# make a loop to fill in inits for each model run
for (i in 1:length(test)){
  # object for each covariate column for multivariate lm
  v <- c()
  for (j in 2:ncol(test[[i]])){
    # make multivariate lm call
    v[j-1] <- paste0('test[[i]][,',j,']')
  }
  # make lm call
  lm_formula <- as.formula(paste("y ~ ",
                                 paste(v, collapse='+')))
  # run lm
  beta.init[[i]] <- lm(lm_formula, data = data)
}


### alpha version for disturbance probability
# test object of covariates
test <- dpls
# make test data object
data <- list(y = cs_dists,
             dist = as.numeric(cs_dists < -1),
             test = test)

alpha.init <- list()
# make a loop to fill in inits for each model run
for (i in 1:length(test)){
  # object for each covariate column for multivariate lm
  v <- c()
  for (j in 2:ncol(test[[i]])){
    # make multivariate lm call
    v[j-1] <- paste0('test[[i]][,',j,']')
  }
  # make lm call
  glm_formula <- as.formula(paste("dist ~ ",
                                  paste(v, collapse='+')))
  # run lm
  alpha.init[[i]] <- glm(glm_formula, data = data,
                         family = binomial(link = "logit"))
  
}

### SET INITIAL CONDITIONS:
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

# which model to run (1 = first model in list)
modelrun = 1

# make inits object for model input: 3/29 - without alphas
init <- list(R = R_mean,
             beta0 = initer(modelrun, "beta")[1],
             alpha = initer(modelrun, "alpha"))
# disturbance magnitude and probability covariates going in the model
# mag covariates
#dmbeta <- dmls[[modelrun]]
# prob covariates
dpalpha <- dpls[[modelrun]]

# condition score disturbance year onset data
cs_model <- cs[grep(paste0("^", distyear, sep = ""), names(cs))]

### MODEL INPUTS
# data object for model runs
data = list(y = cs_model, ns = nsites,
            x_ic = xic, tau_ic = tic,
            tau_obs = cs_precs,
            a = dmalpha,
            a0 = rep(0,4),
            Va = solve(diag(rep(1,4))),
            rmean = R_mean, rprec = R_prec)



### RUN THE MODEL
j.pests <- jags.model (file = textConnection(spongy_disturb),
                       data = data,
                       inits = init,
                       n.chains = 3)


# get outputs
jpout<-coda.samples(j.pests,
                    variable.names = c("beta", "alpha",
                                       "R", 
                                       "pa0"),
                    n.iter = 200000,
                    thin=10)

# run DIC
DIC <- dic.samples(j.pests, n,iter = 10000)
sum <- sum(DIC$deviance, DIC$penalty)