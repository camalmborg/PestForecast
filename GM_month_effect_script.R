#Month effect script
#load TCG mean data:
tcgmeans<-"2020_07_10_500_sample_tcg_mean_MONTHLY.csv"
tcgtable<-read.csv(tcgmeans)
tcg<-as.matrix(tcgtable[,2:131])

#monthly steady states:
#maymeans<-apply(tcg[,seq(1,ncol(tcg),by=5)],1,mean,na.rm=T)
steadymonths<-matrix(NA,nrow=2500,ncol=5)
months<-5
for (i in 1:months){
  steadymonths[,i]<-apply(tcg[,seq(i,ncol(tcg[,90:109]),by=5)],1,mean,na.rm=T)
}
#disturbance window steady states:
disturbmonths<-matrix(NA,nrow=2500,ncol=5)
for (i in 1:months){
  disturbmonths[,i]<-apply(tcg[,seq(i,ncol(tcg[,110:119]),by=5)],1,mean,na.rm=T)
}
#recovery window:
recov<-tcg[,120:130]



#a little visual-------------
seq<-seq(1,2500,by=250)

plot(steadymonths[1,],type="l",ylim=c(0.09,0.32))
for (i in seq){
  lines(steadymonths[i,])
  lines(disturbmonths[i,],col="red")
}

plot(disturbmonths[1,],type="l",ylim=c(0.08,0.35))
for (i in seq){
  lines(disturbmonths[i,])
}

plot(recov[1,],type="l",ylim=c(0.08,0.35))
for (i in seq){
  lines(recov[i,])
}


###Model Structure draft------
#----disturbance prob model over sites------
pests.m = "
model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model
  for(t in 1:n){
    y[s,t] ~ dnorm(x[s,t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    muN[s,t]<-R*x[s,t-1] + month[m[t]]
    x[s,t] ~ dnorm(mu[s,t],tau_add)
    muD[s,t] ~ dnorm(mu0,pa0)
    D[s,t] ~ dbern(p)
    mu[s,t] <- D[s,t]*(muD[s,t] + month[m[t]]) + (1-D[s,t])*muN[s,t]  
  }
  
  x[s,1]~dnorm(x_ic,tau_ic)
  
}#end loop over sites

  ####Month effect
  for(t in 1:5){
  month[t] ~ dnorm(0,tau_month)
  }
  
  #### Priors
  tau_obs ~ dgamma(t_obs,a_obs)
  tau_add ~ dgamma(a_add,t_add)
  tau_month ~dgamma(1,0.1)
  R ~ dnorm(rmean,rprec)  #rho term
  p ~ dunif(0,1)  #disturbance probability
  mu0 ~ dnorm(-5,1) #mean of disturbed state
  pa0 ~ dgamma(1,1) #precision of disturbed state
}
"

data.m = list(y=obs, n=NT, ns=length(geoID),
              x_ic=0, tau_ic=0.1,
              a_obs=0.1,t_obs=0.1,
              a_add=0.1,t_add=0.1,
              rmean=0,rprec=0.00001,
              m=rep(1:5,length=NT))
