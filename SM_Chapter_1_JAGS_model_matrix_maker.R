# CHAPTER 1
# this is a script for making matrix of models for disturbance magnitude and 
# disturbance probability parameterizations

### LIBRARIES:
library(tidyverse)
library(dplyr)


### LOAD DATA:
# # the environmental variables
# # for disturbance magnitude parameters
# varfile_m <- "CHAPTER_1/DATA/MV_2023_12_DATA_distmag.csv"
# magvars <- read.csv(varfile_m)[,-1]
# # for disturbance probability parameters
# varfile_p <- "CHAPTER_1/DATA/MV_2023_12_DATA_distprob.csv"
# probvars <- read.csv(varfile_p)[,-1]
# # load data for dmr to remove missing sites from environmental data
# dmr_file <- "CHAPTER_1/DATA/2023_12_DMR_DATA_TCG.csv"
# dmr <- read.csv(dmr_file)
# # environmental data without missing sites
# magvars <- magvars[dmr$X,]
# probvars <- probvars[dmr$X,]

# Loaded from 4/29/2024 making new anomaly data with new SMAP data
#magvars <- read.csv("CHAPTER_1/2024_JAGS_models/2024_04_magvars.csv")
probvars <- read.csv("CHAPTER_1/2024_JAGS_models/2024_04_probvars.csv")

### Load model list data frame:
#modfile <- "CHAPTER_1/2024_JAGS_models/2024_02_20_Dist_Mag_Models.csv"
#modfile <- "CHAPTER_1/2024_JAGS_models/2024_02_22_Dist_Prob_Models.csv"
#modfile <- "CHAPTER_1/2024_JAGS_models/2024_04_04_Dist_Mag_Models.csv"
modfile <- "CHAPTER_1/2024_JAGS_models/2024_04_17_Dist_Prob_Models.csv"
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
#dmls <- covls
dpls <- covls

# save
save(dmls, file = "CHAPTER_1/2024_JAGS_models/2024_04_dmls.RData")
save(dpls, file = "CHAPTER_1/2024_JAGS_models/2024_04_dpls.RData")

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


### MAKING ANOMALY COVARIATE DATA
# # covariate data - for parameterizations
# # the environmental variables
# # for disturbance magnitude parameters
# varfile_m <- "CHAPTER_1/DATA/MV_2023_12_DATA_distmag.csv"
# magvars <- read.csv(varfile_m)[,-1]
# # for disturbance probability parameters
# varfile_p <- "CHAPTER_1/DATA/MV_2023_12_DATA_distprob.csv"
# probvars <- read.csv(varfile_p)[,-1]
# # environmental data without missing sites
# magvars <- magvars[dmr$X,]
# probvars <- probvars[dmr$X,]
# 
# # load disturbance magnitude and disturbance probability covariates
# # disturbance magnitude
# load("Chapter_1/2024_02_JAGS_models/magls.RData")
# # disturbance probability
# load("Chapter_1/2024_02_JAGS_models/probls.RData")
# 
# # making covariate data with anomalies rather than raw
# # function:
# anomfx<-function(x){
#   # get mean values for each column
#   means <- apply(x, 2, mean, na.rm=T)
#   # make anomaly matrix
#   anom <- matrix(NA,nrow=nrow(x), ncol=ncol(x))
#   # for each row, fill in covariate anomaly
#   for (i in 1:nrow(x)){
#     for (j in 1:ncol(x)){
#       anom[i,j] <- x[i,j] - means[j] # fill in anomalies for beta[] terms
#     }
#     anom[i,1] <- 1 # make first column 1s for beta0 term
#   }
#   return(anom)
# }
# 
# # converting covariate lists to anomaly
# # new empty list to populate with anomaly versions
# anomls <- list()
# # list being converted CHOOSE magls for disturbance magnitude/probls for disturbance probability
# covls <- probls
# #covls <- probls
# # loop for conversion
# for (i in 1:length(covls)){
#   # run cov members through the anomaly machine
#   anomls[[i]] <- anomfx(covls[[i]])
# }
# # disturbance magnitude saving anomaly version
# #dmls <- anomls
# # disturbance probability saving anomaly version
# dpls <- anomls


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

### beta version
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


### alpha version
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

# figure out how to call them when you need them:
# model.run = 1
# p.aram = "beta"
# # make the param object to call the right list (alpha, beta...)
# param.obj = paste0(p.aram, ".init")
# # make an empty vector to grab param inits
# param.init <- vector()
# for (i in 1:length(get(param.obj)[[model.run]])){
#   param.init[i] <- coef(coeftest[[i]])[-1]
# }

# beta.init = list(lm(y ~ v + va, data = data))
#                  #lm(y ~ va, data = data))
# # alpha intercept initial condition
# # get disturbance binary data
# dist = as.numeric(data$y < -1)
# # glm analysis with binomial logit 
# alpha.init = glm(dist ~ data$z, family = binomial(link="logit"))
# # add to init list
# init<-list(R = R_mean,
#            beta0 = coef(beta.init[[1]])[1], # ask Mike about beta 0 init...
#            beta = coef(beta.init[[1]])[-1],
#                     #coef(beta.init[[2]])[-1]),
#            alpha0 = coef(alpha.init)[1]#,
#            #alpha = coef(alpha.init)[-1]
#            )


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


# from the JAGS script (template):
# add to init list
init.test<-list(R = R_mean,
           #beta0 = coef(beta.init[[1]])[1],
           #beta = c(coef(beta.init[[1]])[-1],
                    #coef(beta.init[[2]])[-1]),
           #alpha0 = coef(alpha.init)[1]#,
           #alpha = coef(alpha.init)[-1]
           beta = initer(1, "beta"),
           alpha = initer(1, "alpha")
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
  
  logit(D[s]) <- alpha0 ##+ inprod(alpha[], a[s,])  ##step 2: adding process model here (PROB) - a = covariates
  ## is it : logit(D[s]) <- inprod(alpha[], a[s,])
    
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
  alpha0 ~ dnorm(0,0.001)
  pa0 ~ dgamma(1,1)                   ##precision of disturbed state
  pan ~ dgamma(1,1)                   ##precision of undisturbed state
  
  ## covariate matrix:
  beta ~ dmnorm(b0, Vb)   ## for disturbance magnitude
  #alpha ~ dmnorm(a0, Va)  ## for disturbance probability
  
}
"


# ## output saving test
# savetest <- function(test){
#   filepath <- paste0("CHAPTER_1/", "test", ".RData")
#   save(test, file = filepath)
# }
