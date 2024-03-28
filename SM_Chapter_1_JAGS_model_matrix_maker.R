# CHAPTER 1
# this is a script for making matrix of models for disturbance magnitude and 
# disturbance probability parameterizations

### LIBRARIES:
library(tidyverse)
library(dplyr)


### LOAD DATA:
# the environmental variables
# for disturbance magnitude parameters
varfile_m <- "CHAPTER_1/DATA/MV_2023_12_DATA_distmag.csv"
magvars <- read.csv(varfile_m)[,-1]
# for disturbance probability parameters
varfile_p <- "CHAPTER_1/DATA/MV_2023_12_DATA_distprob.csv"
probvars <- read.csv(varfile_p)[,-1]
# load data for dmr to remove missing sites from environmental data
dmr_file <- "CHAPTER_1/DATA/2023_12_DMR_DATA_TCG.csv"
dmr <- read.csv(dmr_file)
# environmental data without missing sites
magvars <- magvars[dmr$X,]
probvars <- probvars[dmr$X,]

### Load model list data frame:
modfile <- "CHAPTER_1/2024_02_JAGS_models/2024_02_20_Dist_Mag_Models.csv"
#modfile <- "CHAPTER_1/2024_02_JAGS_models/2024_02_22_Dist_Prob_Models.csv"
# disturbance mag/prob model file
modeldf <- read.csv(modfile, header = F)


### Making list of data for models:
# empty list
covls <- list()
# mag (magvars) or prob (probvars)
covs <- magvars
#covs <- probvars
# loop over models
for (i in 1:nrow(modeldf)){
  # grab column numbers from model table
  cols <- as.numeric(modeldf[i,])  #add a row of ones to each, make i into i+1
  # remove NAs
  cols <- cols[!is.na(cols)]
  # make column of 1s for matrix multiplication step in JAGS model (intercept term)
  int <- rep(1, nrow(covs))
  # make dataframe with columns from covariates
  df <- cbind(int, covs[,c(cols)])
  # add to the list
  covls[[i]] <- df  ## have to change this line to make sure its mags or prob
}

# for saving result
magls <- covls
#probls <- covls


## informative R prior ("complicated version" 3/13)
# time series 1 object
ts1 <- matrix(data = NA, nrow = nrow(cs), ncol = yrs)
# time series 2 object
ts2 <- matrix(data = NA, nrow = nrow(cs), ncol = yrs)
# make cors vector
cors <- vector()
# make vectors of "before" and "after"
for (i in 1:nrow(cs)){
  # pick out disturbance date for each site
  td <- which(colnames(cs) == dists[i])
  # make vector of each row values:
  v <- cs[i, (td-yrs):(td-1)]
  v2 <- cs[i, (td-yrs+1):td]
  # fill in matrices of values for t-1 (ts1) and t (ts2)
  for (j in 1:ncol(v)){
    ts1[i,j] <- v[,j]
    ts2[i,j] <- v2[,j]
  }
  #correlation between two time series at each site:
  cors[i] <- cor(ts1[i,], ts2[i,], use = "complete.obs")
}
## ^^ re-did this the "less complicated way" in JAGS model script ^^



### MATRIX MULTIPLICATION MODEL ### -----

## priors and data inputs list: ----
# tau_obs[s] =  precisions from GEE condition score sd's
# b0 = matrix of betas (disturbance magnitude covariate intercepts)
# Vb = matrix of beta precisions <- uninformative <- use solve(diag())
# a0 = matrix of alphas (disturbance prob covariate intercepts)
# Va = matrix of alpha precisions <- uninformative 
# rmean = informative R prior from arima modeling
# rprec = informative R prior from arima modeling
# x_ic[s] = previous time step condition score
# tau_ic[s] = previous time step precision
# covs <- matrices of covariates

### 3/28/2024 Initial Conditions Draft loops
### initial state of model parameters:
# beta intercept initial condition
beta.init = list(lm(y ~ v, data = data),
                 lm(y ~ va, data = data))

# 3/28/2014 ---
# test the dmls matrix list 
test <- dmls
# make empty list to fill in magnitude covariates
beta.init <- list()
# make a loop to fill in inits for each model run
for (i in 1:length(test)){
  # make list for lms
  init.ls <- list()
  for (j in 2:ncol(test[[i]])){
    # run lm for each covariate in that set
    init.ls[[j-1]] <- lm(y ~ test[[i]][,j], data = data)
  }
  beta.init[[i]] <- init.ls 
  rm(init.ls)
}

#formula(covm, form)

# alpha intercept initial condition
# get disturbance binary data
dist = as.numeric(data$y < -1)
# glm analysis with binomial logit 
# alpha.init = glm(dist ~ data$z, family = binomial(link="logit"))
# add to init list
init<-list(R = R_mean,
           beta0 = coef(beta.init[[1]])[1],
           beta = c(coef(beta.init[[1]])[-1],
                    coef(beta.init[[2]])[-1]),
           alpha0 = coef(alpha.init)[1]#,
           #alpha = coef(alpha.init)[-1]
)



### 3/11/2024 DRAFT:
spdist_mm <- "model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model:
  y[s] ~ dnorm(mu[s], tau_obs[s]) 
  
  #### Process Model:
  muN[s] ~ dnorm(mun, pan)            
  muD[s] ~ dnorm(mu0[s], pa0)                     ##step 1: process model on mu0 (MAG) - b = covariates
  
  logit(D[s]) <- alpha0 + inprod(alpha[], a[s,])  ##step 2: adding process model here (PROB) - a = covariates
    
  mu[s] <- D[s] * muD[s] + (1-D[s]) * muN[s]
  mu0[s] <- inprod(beta[], b[s,])
  mun <- R * x[s]                                 ##step 3: dealing with modeling R (Chap 2 - RECOV) 

  x[s] ~ dnorm(x_ic[s], tau_ic[s])
  
}#end loop over sites
  
  #### Priors
  tau_obs ~ dgamma(t_obs, a_obs)      ##observation error (data model)
  ##tau_add ~ dgamma(a_add ,t_add)    ##process error (process model)
  R ~ dnorm(rmean, rprec)             ##rho paramter (recovery rate)  ## going to put an informative prior here
  ##beta0 ~ dnorm(-5,1)                 ##param for calculating mean of disturbed state
  ##alpha0 ~ dnorm(0,0.001)
  pa0 ~ dgamma(1,1)                   ##precision of disturbed state
  pan ~ dgamma(1,1)                   ##precision of undisturbed state
  
  ## covariate matrix:
  beta ~ dmnorm(b0, Vb)   ## for disturbance magnitude
  alpha ~ dmnorm(a0, Va)  ## for disturbance probability
  
}
"
