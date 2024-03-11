# CHAPTER 1
# this is a script for making matrix of models for disturbance magnitude and 
# disturbance probability parameterizations

### LIBRARIES:

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
#modfile <- "CHAPTER_1/2024_02_JAGS_models/2024_02_20_Dist_Mag_Models.csv"
modfile <- "CHAPTER_1/2024_02_JAGS_models/2024_02_22_Dist_Prob_Models.csv"
# disturbance mag/prob model file
modeldf <- read.csv(modfile, header = F)


### Making list of data for models:
# empty list
covls <- list()
# mag (magvars) or prob (probvars)
#covs <- magvars
covs <- probvars
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
probls <- covls



### MATRIX MULTIPLICATION MODEL ### -----
spdist_mm <- "model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model:
  y[s] ~ dnorm(mu[s], tau_obs[s])     ## fill in tau_obs[s] with precisions from GEE sd's
  
  #### Process Model:
  muN[s] <- R * x[s]                  ##step 3: dealing with modeling R (Chap 2 - RECOV)
     ##x[s] ~ dnorm(mu[s], tau_add)
  muD[s] ~ dnorm(mu0[s], pa0)         ##step 1: process model on mu0 (MAG)
    
     ## D[s] ~ dbern(p)                     ##step 2: adding process model here (PROB)
  logit(D[s]) <- alpha[1] + alpha[2]*z[s]
  ##alpha[1] ~ dnorm(0.0, 0.0001)
    
  mu[s] <- D[s] * muD[s] + (1-D[s]) * muN[s]
  mu0[s] <- beta0 + inprod(beta[], x[s,])

  ## x[s]~dnorm(x_ic, tau_ic)
  
}#end loop over sites
  
  #### Priors
  tau_obs ~ dgamma(t_obs, a_obs)     ##observation error (data model)
  ##tau_add ~ dgamma(a_add ,t_add)    ##process error (process model)
  R ~ dnorm(rmean, rprec)            ##rho paramter (recovery rate)  ## going to put an informative prior here
  ## p ~ dunif(0,1)  #disturbance probability
  beta0 ~ dnorm(-5,1) #param for calculating mean of disturbed state
  pa0 ~ dgamma(1,1) #precision of disturbed state
  
  ## covariate matrix:
  beta ~ dmnorm(b0, Vb)   ## for disturbance magnitude
  alpha ~ dmnorm(a0, Va)  ## for disturbance probability
  
}
"

## priors and data list: ----
# tau_obs[s] =  precisions from GEE condition score sd's
# b0 = matrix of betas (disturbance magnitude covariate intercepts)
# Vb = matrix of beta precisions <- need to ask Mike about how to make this
# a0 = matrix of alphas (disturbance prob covariate intercepts)
# Va = matrix of alpha precisions <- need to ask Mike 
# rmean = informative R prior from arima modeling
# rprec = informative R prior from arima modeling
