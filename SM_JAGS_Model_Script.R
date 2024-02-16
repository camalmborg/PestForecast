#This script is the JAGS model for pest disturbance mag and prob forecast

### THE LIBRARIES:
library(tidyverse)
library(dplyr)
library(rjags)
library(coda)
#library(ecoforecastR)


### THE DATA:
# condition scores inputs
#cfile <- "2023_03_08_DATAGRAB/2023_03_08_5000_sites_sample_score_mean.csv" #2005-2022
cfile <- "2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv" #1995-2020, original grab
scores <- read.csv(cfile)
# isolate data we want:
cs <- scores %>%
  # Drop unwanted columns
  dplyr::select(dplyr::starts_with("X")) %>%
  # Select just June averages for best disturbance visual without seasonality
  dplyr::select(dplyr::contains(".06.")) %>%
  # 
  dplyr::rename_with(~ str_replace_all(., c("X|_score_mean|_cs_mean" = "", 
                                       "\\." = "-")))

# covariate data - for parameterizations
#varfile <- ""
#vars <- read.csv(varfile)


### SELECT SITES:
# random selection of sites for testing
smpl <- sample(nrow(cs), 5)
# make sample
cs_samp <- cs[smpl,]


### initial state of model parameters:
init<-list()
nchain <- 3
for(j in 1:nchain){
  samp<- sample(!is.na(cs_samp),length(cs_samp),replace=TRUE)
  init[[j]]<-list(tau_add=1/var(diff(samp)),tau_obs=1/var(samp))
}


### THE MODEL:
#use the single time step version of the model:
spongy_disturb <- "model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model
  y[s] ~ dnorm(x[s,t],tau_obs)
  
  #### Process Model
  muN[s,t]<-R*x[s,t-1] ##step 3: dealing with modeling R
  x[s,t] ~ dnorm(mu[s,t],tau_add)
  muD[s,t] ~ dnorm(mu0[s,t],pa0) ##step 1: process model on mu0
    
  D[s,t] ~ dbern(p) ##step 2: adding process model here
    
  mu[s,t] <- D[s,t]*muD[s,t] + (1-D[s,t])*muN[s,t]
  mu0[s,t] <- beta0 ## + beta[1] * variables

  x[s,1]~dnorm(x_ic,tau_ic)
  
}#end loop over sites
  
  #### Priors
  tau_obs ~ dgamma(t_obs,a_obs)
  tau_add ~ dgamma(a_add,t_add)
  R ~ dnorm(rmean,rprec)  #rho term
  p ~ dunif(0,1)  #disturbance probability
  beta0 ~ dnorm(-5,1) #param for calculating mean of disturbed state
  pa0 ~ dgamma(1,1) #precision of disturbed state
  # beta[1] 
  
}
"

#MODEL INPUTS
#data and parameters for sites model:
# data = list(y=cs_samp, n=NT, ns=nsites,
#               x_ic=0, tau_ic=0.1,
#               a_obs=0.1,t_obs=0.1,
#               a_add=0.1,t_add=0.1,
#               rmean=0,rprec=0.00001)

# j.pests <- jags.model (file = textConnection(spongy_disturb),
#                        data = data,
#                        inits = init,
#                        n.chains = 3)


# jpout<-coda.samples(j.pests, 
#                     variable.names = c("beta0", "beta[1]", "beta[2]",
#                                        "beta[3]", "beta[4]",
#                                        "tau_add","tau_obs", 
#                                        "pa0"),
#                     n.iter = 50000,
#                     thin=2)


