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

#vpd and precip  data:
vpd<-hf16[smpl,126:127] #vpd 2014-2015 hatch
pcp<-hf16[smpl,c(76,100:102)] #precip 2014-2015 feed and 2016 hatch/feed
#make anomaly datasets:
anomfx<-function(x){
  means<-apply(x,2,mean, na.rm=T)
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


# #initial state of model parameters
# init<-list()
# nchain <- 3
# for(j in 1:nchain){
#   samp<- sample(!is.na(condscores.samp),length(condscores.samp),replace=TRUE)
#   init[[j]]<-list(tau_add=1/var(diff(samp)),tau_obs=1/var(samp))
# }


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
    mu0[s,t] <- beta0 + beta[1]*pcp[s,1] + beta[2]*pcp[s,2] + beta[3]*pcp[s,3] + beta[4]*pcp[s,4]
    ##beta[1]*vpd[s,1]
    ##beta[2]*vpd[s,2]
    ##beta[1]*pcp[s,1]
    ##beta[2]*pcp[s,2]
    ##beta[3]*pcp[s,3]
    ##beta[4]*pcp[s,4]
  }
  
  x[s,1]~dnorm(x_ic,tau_ic)
  
}#end loop over sites
  
  #### Priors
  tau_obs ~ dgamma(t_obs,a_obs)
  tau_add ~ dgamma(a_add,t_add)
  R ~ dnorm(rmean,rprec)  #rho term
  p ~ dunif(0,1)  #disturbance probability
  beta0 ~ dnorm(-5,1) #param for calculating mean of disturbed state
  beta[1] ~ dnorm(0,0.0001)
  beta[2] ~ dnorm(0,0.0001)
  beta[3] ~ dnorm(0,0.0001)
  beta[4] ~ dnorm(0,0.0001)
  pa0 ~ dgamma(1,1) #precision of disturbed state
  
}
"

#MODEL INPUTS
#data and parameters for sites model:
data = list(y=condscores.samp, n=NT, ns=nsites,
              x_ic=0, tau_ic=0.1,
              a_obs=0.1,t_obs=0.1,
              a_add=0.1,t_add=0.1,
              rmean=0,rprec=0.00001,#,
              pcp=pcpanom)#, vpd=vpdanom)

j.pests <- jags.model (file = textConnection(spongy_disturb),
                       data = data,
                       inits = init,
                       n.chains = 3)

# jpout<-coda.samples(j.pests,
#                     variable.names = c("beta0", "beta"
#                                        "tau_add","tau_obs",
#                                        "pa0","x","D"),
#                     n.iter = 5000,
#                     thin=5)
# save(jpout, file=paste0("Model_1_mu0_jpout_esarun1_", as.character(i),".RData"))


# jpout<-coda.samples(j.pests, 
#                     variable.names = c("beta0", "beta[1]", "beta[2]",
#                                        "beta[3]", "beta[4]",
#                                        "tau_add","tau_obs", 
#                                        "pa0"),
#                     n.iter = 50000,
#                     thin=2)


for (i in 1:20){
  jpout<-coda.samples(j.pests,
                      variable.names = c("beta0", "beta[1]", "beta[2]",
                                         "beta[3]", "beta[4]",
                                         "tau_add","tau_obs",
                                         "pa0", "x", "D"),
                      n.iter = 20000,
                      thin=5)
  save(jpout, file=paste0("Model_2_pcp_jpout_esarun2_", as.character(i),".RData"))
}


###THE REST OF THIS CODE IS FOR 509 STUFF:
#getting jpout objects into matrix for plotting----

#getting proper names for each site from looped-over-sites mcmc output:
##' @param w mcmc object containing matrix outputs
##' @param pre prefix (variable name) for the matrix variable to be extracted
##' @param numeric boolean, whether to coerce class to numeric
parse.MatrixNames <- function(w, pre = "x", numeric = FALSE) {
  w <- sub(pre, "", w)
  w <- sub("[", "", w, fixed = TRUE)
  w <- sub("]", "", w, fixed = TRUE)
  w <- matrix(unlist(strsplit(w, ",")), nrow = length(w), byrow = TRUE)
  if (numeric) {
    class(w) <- "numeric"
  }
  colnames(w) <- c("row", "col")
  return(as.data.frame(w))
}

#load a sample mcmcobject for the model we want to use:
#model1<-"Model_1_mu0_jpout_run2_"
#model2<-"Model_2_precip5k_jpout_"
#model3<-"Model_3_fullenv5k_jpout_"
model1<-"Model_1_mu0_jpout_esarun1_"
model2<-"Model_2_pcp_jpout_esarun2_"
i=1
load(file = paste0(model2,as.character(i),".RData"))
#make empty matrix with sample mcmc object's # of columns (1 per tracked param)
out<-matrix(NA,ncol=ncol(jpout[[1]]))
#grab each mcmc object of that model and convert to matrix:
for (i in 20:19){
  load(file = paste0(model2,as.character(i),".RData"))
  jpmx<-as.matrix(jpout)
  out<-rbind(out,jpmx)
  rm(jpmx)
}
#out1<-out[-1,] #removes the NA row
out2<-out[-1,]
#out3<-out[-1,]
rm(jpout)
rm(out)

#separate parameters into x,D for plotting, and other params for summary:
param1<-c("^beta","tau_add","tau_obs", "pa0")
param2<-c("^beta","tau_add","tau_obs", "pa0")
#param3<-c("^beta", "beta[3]","beta[4]","tau_add","tau_obs", "pa0")

sel.1<-grepl(paste(param1, collapse = "|"), colnames(out1))
sel.2<-grepl(paste(param2, collapse = "|"), colnames(out2))
#sel.3<-grepl(paste(param3, collapse = "|"), colnames(out3))
params1<-out1[,sel.1]
params2<-out2[,sel.2]
#params3<-out3[,sel.3]

summ1<-as.data.frame.matrix(summary(params1))
summ2<-as.data.frame.matrix(summary(params2))
#summ3<-as.data.frame.matrix(summary(params3))



x.cols.1 <- grep("^x",colnames(out1))
ci.x.1 <- apply(out1[,x.cols.1],2,quantile,c(0.025,0.5,0.975))

x.cols.2 <- grep("^x",colnames(out2))
ci.x.2 <- apply(out2[,x.cols.2],2,quantile,c(0.025,0.5,0.975))
# 
# x.cols.3 <- grep("^x",colnames(out3))
# ci.x.3 <- apply(out3[,x.cols.3],2,quantile,c(0.025,0.5,0.975))

ci.x.names = parse.MatrixNames(colnames(ci.x.1),numeric=TRUE)
ci.d.names = parse.MatrixNames(colnames(ci.d.1))

d.cols.1 <- grep("^D",colnames(out1))
ci.d.1 <- apply(out1[,d.cols.1],2,quantile,c(0.25,0.5,0.975))
# 
d.cols.2 <- grep("^D",colnames(out2))
ci.d.2 <- apply(out2[,d.cols.2],2,quantile,c(0.25,0.5,0.975))
# 
# d.cols.3 <- grep("^D",colnames(out3))
# ci.d.3 <- apply(out3[,d.cols.3],2,quantile,c(0.25,0.5,0.975))

###FOR VISUALIZATIONS:-----
#load("modeldata.RData")
condscores.samp<-data$y

i=84
sitei = which(ci.x.names$row == i)
time=1:130
NT=length(time)
tiff(paste0("ESA_bothmodels_site2_",as.character(i),".tiff"), units="in", width=10, height=5, res=300)
plot(ci.x.1[2,sitei],type='l',ylim=c(-12,5),
     ylab="Forest Condition Score", 
     col="red",
     xlab="Month",
     cex=1)
lines(ci.x.2[2,sitei], col="blue")
#lines(ci.x.3[2,sitei], col="dark green")
ecoforecastR::ciEnvelope(time,ci.x.1[1,sitei],ci.x.1[3,sitei],col=ecoforecastR::col.alpha("indianred1",0.30))
ecoforecastR::ciEnvelope(time,ci.x.2[1,sitei],ci.x.2[3,sitei],col=ecoforecastR::col.alpha("lightblue1",0.30))
#ecoforecastR::ciEnvelope(time,ci.x.3[1,sitei],ci.x.3[3,sitei],col=ecoforecastR::col.alpha("greenyellow",0.30))
points(time,condscores.samp[i,],pch=16,cex=0.5,col="navyblue")
dev.off()


tiff(paste0("ESA_model1_site",as.character(i),".tiff"), units="in", width=8, height=3, res=300)
plot(ci.x.1[2,sitei],type='l',ylim=range(condscores.samp,na.rm=TRUE),
     ylab="Forest Condition Score",
     col="red",
     xlab="Month",
     cex=1)
ecoforecastR::ciEnvelope(time,ci.x.1[1,sitei],ci.x.1[3,sitei],col=ecoforecastR::col.alpha("indianred1",0.60))
points(time,condscores.samp[i,],pch=16,cex=0.5,col="navyblue")
dev.off()
# 
#
tiff(paste0("ESA_model2_site",as.character(i),".tiff"), units="in", width=8, height=3, res=300)
plot(ci.x.2[2,sitei],type='l',ylim=range(condscores.samp,na.rm=TRUE),
     ylab="Forest Condition Score",
     col="blue",
     xlab="Month",
     cex=1)
ecoforecastR::ciEnvelope(time,ci.x.2[1,sitei],ci.x.1[3,sitei],col=ecoforecastR::col.alpha("lightblue1",0.60))
points(time,condscores.samp[i,],pch=16,cex=0.5,col="navyblue")
dev.off()
# 
# tiff("509_Model_3_examp_site34.tiff", units="in", width=8, height=3, res=300)
# plot(ci.x.3[2,sitei],type='l',ylim=range(condscores.samp,na.rm=TRUE),
#      ylab="Forest Condition Score", 
#      col="green",
#      xlab="Month",
#      cex=1)
# ecoforecastR::ciEnvelope(time,ci.x.3[1,sitei],ci.x.1[3,sitei],col=ecoforecastR::col.alpha("greenyellow",0.60))
# points(time,condscores.samp[i,],pch=16,cex=0.5,col="navyblue")
# dev.off()



#DIC calculations:
#DIC.fullenv<-dic.samples(j.pests, n.iter=10000)
#DIC.threevar<-dic.samples(j.pests, n.iter=10000)
#DIC.justprecip<-dic.samples(j.pests, n.iter=10000)
#DIC.mu0model<-dic.samples(j.pests, n.iter=10000)

###plots of individual sites testing:-----
#plot(1:NT,condscores.samp[100,],type="l")


i=2
sitei = which(ci.d.names$row == paste0("D",i))
#sitei = sitei[2:NT]
#plot for disturbance prob all models:
tiff(paste0("Dist_Prob_allmodels_ESARUN1_", i,".tiff"),
     units="in", width=10, height=5, res=300)
#plot:
plot(time[2:NT],ci.d.1[2,sitei],type='l',ylim=range(ci.d.1,na.rm=TRUE),
     ylab="Disturbance Probability",
     col="red",
     xlab="Month",
     cex=1)
lines(time[2:NT],ci.d.2[2,sitei], col="blue")
ecoforecastR::ciEnvelope(time[2:NT],ci.d.1[1,sitei],ci.d.1[3,sitei],col=ecoforecastR::col.alpha("indianred1",0.30))
ecoforecastR::ciEnvelope(time[2:NT],ci.d.2[1,sitei],ci.d.2[3,sitei],col=ecoforecastR::col.alpha("lightblue1",0.30))
dev.off()


