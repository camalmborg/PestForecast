#This script is the JAGS model for pest disturbance mag and prob forecast

### THE LIBRARIES:
library(tidyverse)
library(dplyr)
library(rjags)
library(coda)
#library(MCMCvis)
#library(ecoforecastR)


### THE DATA:
# condition scores inputs
#cfile <- "2023_03_08_DATAGRAB/2023_03_08_5000_sites_sample_score_mean.csv" #2005-2022
cfile <- "2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv" #1995-2020, original grab
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
# # number of sites
# nsites = nrow(cs)

# disturbance year data
dmr_file <- "CHAPTER_1/DATA/2023_12_DMR_DATA_TCG.csv"
dmr <- read.csv(dmr_file)

# choose disturbance onset year
distyear <- 2016

# disturbance data for sample
cs <- cs[dmr$X,]
# get disturbance years and condition score value at dist year for each:
dists <- vector()
cs_dists <- vector()
for (i in 1:nrow(cs)){
  # get disturbance year colname
  if (dmr$dpy1[i] > 0) {
    dists[i] <- colnames(cs[grep(paste0("^", distyear, sep = ""), names(cs))])
  } else {
    dists[i] <- colnames(cs[grep(paste0("^",distyear+1, sep = ""), names(cs))])
  }
  # get condition score associated with disturbance year
  cs_dists[i] <- cs[i, dists[i]]
}

# # standard deviations for precisions
# sdfile <- "2022_08_31_DATAGRAB/2022_12_7_sample_score_stddev_5k.csv"
# stan_devs <- read.csv(sdfile)
# sds <- stan_devs %>%
#   # Drop unwanted columns
#   dplyr::select(dplyr::starts_with("X")) %>%
#   # Select just June averages for best disturbance visual without seasonality
#   dplyr::select(dplyr::contains(".06.")) %>%
#   # 
#   dplyr::rename_with(~ str_replace_all(., c("X|_score_stddevs" = "", 
#                                             "\\." = "-")))
# # convert to precisions -- DOESN'T WORK RIGHT YET 2/19/24
# prec_convert<-function(x){
#   prec<-matrix(NA,nrow=length(x),ncol=1)
#   for (i in 1:length(x)){
#     prec[i,]<-1/x[i]
#     if(is.na(prec[i,])){prec[i,]<-mean(prec,na.rm=T)}
#   }
#   return(prec)
# }
# # get precisions from std devs
# precs<-as.numeric(prec_convert(sds))


# covariate data - for parameterizations
# the environmental variables
# for disturbance magnitude parameters
#varfile_m <- "CHAPTER_1/DATA/MV_2023_12_DATA_distmag.csv"
#magvars <- read.csv(varfile_m)[,-1]
# for disturbance probability parameters
#varfile_p <- "CHAPTER_1/DATA/MV_2023_12_DATA_distprob.csv"
#probvars <- read.csv(varfile_p)[,-1]


### THE MODEL:
#use the single time step version of the model:
spongy_disturb <- "model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model:
  y[s] ~ dnorm(mu[s], tau_obs)
  
  #### Process Model:
  muN[s] <- R * x[s]                  ##step 3: dealing with modeling R (Chap 2 - RECOV)
  #x[s] ~ dnorm(mu[s], tau_add)
  muD[s] ~ dnorm(mu0[s], pa0)         ##step 1: process model on mu0 (MAG)
    
  D[s] ~ dbern(p)                     ##step 2: adding process model here (PROB)
  #logit(D[s]) <- alpha[1] + alpha[2]*z[s]
  #alpha[1] ~ dnorm(0.0, 0.0001)
    
  mu[s] <- D[s] * muD[s] + (1-D[s]) * muN[s]
  mu0[s] <- beta0   ## + beta[1] * variables ##

  x[s]~dnorm(x_ic, tau_ic)
  
}#end loop over sites
  
  #### Priors
  tau_obs ~ dgamma(t_obs, a_obs)     ##observation error (data model)
  #tau_add ~ dgamma(a_add ,t_add)    ##process error (process model)
  R ~ dnorm(rmean, rprec)            ##rho paramter (recovery rate)
  p ~ dunif(0,1)  #disturbance probability
  beta0 ~ dnorm(-5,1) #param for calculating mean of disturbed state
  pa0 ~ dgamma(1,1) #precision of disturbed state
  
  # beta[1]   ## COVARIATES WILL BE ADDED HERE
  
}
"

### SELECT SITES:
# random selection of sites for testing
smpl <- sample(nrow(cs), 25)
# make sample
cs_samp <- cs[smpl,]
# number of sites of sample
nsites = nrow(cs_samp)
# get dist years for cs_samp group
dist_samp <- dists[as.numeric(rownames(cs_samp))]
# make the single timestep data for each site
cs_samp_dist <- cs_dists[as.numeric(rownames(cs_samp))]

### initial state of model parameters:
init<-list()
nchain <- 3
for(j in 1:nchain){
  samp<- sample(!is.na(cs_samp),length(cs_samp),replace=TRUE)
  init[[j]]<-list(tau_add=1/var(diff(samp)),tau_obs=1/var(samp))
}


### MODEL INPUTS
# data and parameters for sites model:   #for full sample ns = nrow(scores)
data = list(y = cs_samp_dist, ns = nsites,    
              x_ic = 0, tau_ic = 0.1,
              a_obs = 0.1, t_obs = 0.1,
              #a_add = 0.1, t_add = 0.1,
              rmean = 0, rprec = 0.00001)

### RUN THE MODEL
j.pests <- jags.model (file = textConnection(spongy_disturb),
                       data = data,
                       inits = init,
                       n.chains = 3)


# jpout<-coda.samples(j.pests,
#                     variable.names = c("beta0", "tau_obs",
#                                        "D", "pa0"),
#                     n.iter = 50000,
#                     thin=2)


