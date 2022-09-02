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

#make initial conditions for the x[s], mean junes 5 yrs prior:
#june sequence:
jseq<-seq(2,length(condscores),by=5)
junes<-condscores[,jseq[1:22]] #just up to the 2016 disturb
javgs<-as.numeric(apply(junes[,(length(junes)-5):(length(junes)-1)],1,mean,na.rm=T))

#load st devs for initial conditions x[s]:
stdevs<-read.csv("2022_08_31_DATAGRAB/2022_08_31_5k_score_stddev - 2022_08_31_5k_score_stddev.csv")
sds<-stdevs[,2:23] #just up to 2016
sdavgs<-apply(sds[,(length(sds)-5):(length(sds)-1)],1,mean,na.rm=T)
#convert to precisions:
prec_convert<-function(x){
  prec<-matrix(NA,nrow=length(x),ncol=1)
  for (i in 1:length(x)){
     prec[i,]<-1/x[i]
     if(is.na(prec[i,])){prec[i,]<-mean(prec,na.rm=T)}
  }
  return(prec)
}
precs<-as.numeric(prec_convert(sdavgs))


#load hatch-feed temp,vpd,precip dataset:
#hf16<-read.csv("hf16_dataset_03_2022.csv")
# hf16<-read.csv("hatchfeed16_daymet_tpv_mags_data.csv")
# hf16<-hf16[,-c(3)]
# hfnoX<-as.matrix(hf16[,2:ncol(hf16)])

##SELECT SITES:
#2016 dist mag sites:
# cs16<-condscores[hf16$X,]
# #random selection of sites for testing:
# #smpl<-sample(nrow(cs16),100)
# condscores.samp<-cs16[smpl,]

#number of sites, timesteps:
nsites = nrow(condscores)#.samp)
NT = ncol(condscores)#.samp)

# # #initial state of model parameters
# init<-list()
# nchain <- 3
# for(j in 1:nchain){
#   samp<- sample(!is.na(condscores),length(condscores),replace=TRUE)
#   init[[j]]<-list(tau_add=1/var(diff(samp)),tau_obs=1/var(samp))
# }

# #initial state of model parameters
init<-list()
nchain <- 3
for(j in 1:nchain){
  samp<- sample(!is.na(condscores),length(condscores),replace=TRUE)
  init[[j]]<-list(tau_obs=1/var(samp))
}


#the model loops only over sites with to parameterize mu0 @ 2016 disturb timestep 
#timestep 2016 june >>> t = 107
cond16<-condscores[,107]
prec16<-as.numeric(prec_convert(sds[,22]))
for(i in length(prec16)){
  if(is.na(prec16[i])){prec16[i]==mean(prec16,na.rm=T)}
}
#THE MODEL:
spongy_disturb <- "model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model
  y[s] ~ dnorm(mu[s],tau_obs)
  
  #### Process Model
  muN[s]<-R*x[s]
 ##x[s] ~ dnorm(mu[s],tau_add)
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
  
 x[s]~dnorm(x_ic[s],tau_ic[s])
  
}#end loop over sites
 
  ##xic~dnorm(x_ic,tau_ic)
  
  #### Priors
  tau_obs ~ dgamma(t_obs,a_obs)
##tau_add ~ dgamma(t_add,a_add)
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

#MODEL INPUTS
#data and parameters for sites model:
data = list(y=cond16, ns=nsites,
            x_ic=javgs, tau_ic=precs,
            a_obs=0.1,t_obs=mean(prec16),
            rmean=0,rprec=0.00001)#,#,
#pcp=pcpanom)#, vpd=vpdanom)
#a_add=0.1,t_add=0.1,

j.pests <- jags.model (file = textConnection(spongy_disturb),
                       data = data,
                       inits = init,
                       n.chains = 3)

jpout<-coda.samples(j.pests,
                    variable.names = c("beta0",
                                       "tau_obs",
                                       "pa0"),
                    n.iter = 100000)

out<-as.matrix(jpout)

