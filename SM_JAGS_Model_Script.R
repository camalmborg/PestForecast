#This script updates the initial runs of the JAGS model for pest forecast
#Initially the model ran on 50 test sites, not the full 5000 (4997) sites
#New script updates the initial dataset build and includes the updated model

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
hf16<-read.csv("hf16_dataset_03_2022.csv")
hfnoX<-as.matrix(hf16[,2:ncol(hf16)])

##SELECT SITES:
#2016 dist mag sites:
cs16<-condscores[hf16$X,]
#random selection of sites for testing:
smpl<-sample(nrow(cs16),100)
condscores.samp<-cs16[smpl,]

#number of sites, timesteps:
nsites = nrow(condscores.samp)
NT = ncol(condscores.samp)

#vpd and precip feeding window data:
vpd<-hfnoX[smpl,152:153]
pcp<-hfnoX[smpl,100:101]
#make anomaly datasets:
anomfx<-function(x){
  means<-apply(x,2,mean)
  anom<-matrix(NA,nrow=nrow(x),ncol=ncol(x))
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      anom[i,j]<-x[i,j]-means[j]
    }
  }
  return(anom)
}
#anomaly output:
vpdanom<-anomfx(vpd)
pcpanom<-anomfx(pcp)


#initial state of model parameters
init<-list()
nchain <- 3
for(j in 1:nchain){
  samp<- sample(!is.na(condscores.samp),length(condscores.samp),replace=TRUE)
  init[[j]]<-list(tau_add=1/var(diff(samp)),tau_obs=1/var(samp))
}


#THE MODEL:
spongy_disturb <- "model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model
  for(t in 1:n){
    y[s,t] ~ dnorm(x[s,t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    muN[s,t]<-R*x[s,t-1] ##step 3: dealing with modeling R
    x[s,t] ~ dnorm(mu[s,t],tau_add)
    muD[s,t] ~ dnorm(mu0[s,t],pa0) ##step 1: process model on mu0
    D[s,t] ~ dbern(p) ##step 2: adding process model here
    mu[s,t] <- D[s,t]*muD[s,t] + (1-D[s,t])*muN[s,t]
    mu0[s,t] <- beta0
    ##beta[1]*vpd[s,1]
    ##beta[2]*vpd[s,2]
    ##beta[3]*pcp[s,1]
    ##beta[4]*pcp[s,2]
  }
  
  x[s,1]~dnorm(x_ic,tau_ic)
  
}#end loop over sites
  
  #### Priors
  tau_obs ~ dgamma(t_obs,a_obs)
  tau_add ~ dgamma(a_add,t_add)
  R ~ dnorm(rmean,rprec)  #rho term
  p ~ dunif(0,1)  #disturbance probability
  beta0 ~ dnorm(-5,1) #param for calculating mean of disturbed state
  #beta[1] ~ dnorm(0,0.0001)
  #beta[2] ~ dnorm(0,0.0001)
  #beta[3] ~ dnorm(0,0.0001)
  #beta[4] ~ dnorm(0,0.0001)
  pa0 ~ dgamma(1,1) #precision of disturbed state
  
}
"

#MODEL INPUTS
#data and parameters for sites model:
data = list(y=condscores.samp, n=NT, ns=nsites,
              x_ic=0, tau_ic=0.1,
              a_obs=0.1,t_obs=0.1,
              a_add=0.1,t_add=0.1,
              rmean=0,rprec=0.00001)#,
              #pcp=pcpanom, vpd=vpdanom)

j.pests <- jags.model (file = textConnection(spongy_disturb),
                       data = data,
                       inits = init,
                       n.chains = 3)

for (i in 1:10){
  jpout<-coda.samples(j.pests, 
                      variable.names = c("beta0",
                                         "tau_add","tau_obs", 
                                         "R", "p", 
                                         "pa0"),
                      n.iter = 5000)
  save(jpout, file=paste0("Model_1_mu0_jpout_", as.character(i),".RData"))
}

# jpout <-coda.samples(j.pests.1, 
#                      variable.names = c("beta0","x","y",
#                                         "tau_add","tau_obs", 
#                                         "R", "p", "D", 
#                                         "mu0", "pa0"),
#                      n.iter = 100000)

#plot(jpout)
burnin=5000
jpoutburn <- window(jpout, start=burnin)

jpouthin<-window(jpoutburn,thin=10)


out<-as.matrix(jpouthin)

#models for 509:
#out.fullenv<-jpout
#out.threevar<-jpout
#out.justprecip<-jpout
#out.anthroenv<-
#out.mu0model<-jpout

#DIC calculations:
#DIC.fullenv<-dic.samples(j.pests, n.iter=10000)
#DIC.threevar<-dic.samples(j.pests, n.iter=10000)
#DIC.justprecip<-dic.samples(j.pests, n.iter=10000)
DIC.mu0model<-dic.samples(j.pests, n.iter=10000)


