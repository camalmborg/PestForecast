#This script is the code for the JAGS model with only loop over sites not time
#to 

#THE LIBRARIES
library(rjags)
library(coda)
#library(ecoforecastR)
#library(tidyverse)
#library(dplyr)


#load full dataset (5000 sites):
#model uses condition scores, not TCG raw data
cond.scores.mo<-read.csv("2020_07_10_sample_score_mean_MONTHLY.csv")
condscores<-cond.scores.mo[,2:131]

#load hatch-feed temp,vpd,precip dataset:
#hf16<-read.csv("hf16_dataset_03_2022.csv")
hf16<-read.csv("hatchfeed16_daymet_tpv_mags_data.csv")
hf16<-hf16[,-c(3)]
hfnoX<-as.matrix(hf16[,2:ncol(hf16)])

##SELECT SITES:
#2016 dist mag sites:
cs16<-condscores[hf16$X,]
#random selection of sites for testing:
#smpl<-sample(nrow(cs16),100)
condscores.samp<-cs16[smpl,]

#number of sites, timesteps:
nsites = nrow(condscores.samp)
NT = ncol(condscores.samp)


#the model loops only over sites with to parameterize mu0 @ 2016 disturb timestep 
#timestep 2016 june >>> t = 107
#THE MODEL:
spongy_disturb <- "model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model
  y[s] ~ dnorm(x[s],tau_obs)
  
  #### Process Model
  muN[s,107]<-R*x[s]
  x[s] ~ dnorm(mu[s],tau_add)
  muD[s] ~ dnorm(mu0[s],pa0) 
  D[s] ~ dbern(p)
  mu[s] <- D[s]*muD[s] + (1-D[s])*muN[s]
  mu0[s] <- beta0
  
  ##mu0[s,t] <- beta0 + beta[1]*pcp[s,1] + beta[2]*pcp[s,2] + beta[3]*pcp[s,3] + beta[4]*pcp[s,4]
  ##beta[1]*vpd[s,1]
  ##beta[2]*vpd[s,2]
  ##beta[1]*pcp[s,1]
  ##beta[2]*pcp[s,2]
  ##beta[3]*pcp[s,3]
  ##beta[4]*pcp[s,4]
  
  x[s,1]~dnorm(x_ic,tau_ic)
  
}#end loop over sites
  
  #### Priors
  tau_obs ~ dgamma(t_obs,a_obs)
  tau_add ~ dgamma(a_add,t_add)
  R ~ dnorm(rmean,rprec)  #rho term
  p ~ dunif(0,1)  #disturbance probability
  beta0 ~ dnorm(-5,1) #param for calculating mean of disturbed state
  # beta[1] ~ dnorm(0,0.0001)
  # beta[2] ~ dnorm(0,0.0001)
  # beta[3] ~ dnorm(0,0.0001)
  # beta[4] ~ dnorm(0,0.0001)
  pa0 ~ dgamma(1,1) #precision of disturbed state
  
}
"
