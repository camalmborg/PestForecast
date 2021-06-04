library(ecoforecastR)
library(rjags)

#load canopy condition score sample data--single site
#GMdat<-read.csv("GM_48_mcs.csv")

#load data for all test sites:
time=1:130
NT=length(time)
geoID=1:50
gm.data<-list()
for(i in geoID){
  file<-paste0("TEST_monthly_data/GM_",geoID[i],"_mcs.csv")
  gm.data[[i]]<-read.csv(file)
}

#load condition score means for each site:
obs<-list()
for (i in geoID){
  obs[[i]]<-gm.data[[i]]$score_mean
}

#### model with disturbance probability added:
pests = "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    muN[t]<-R*x[t-1]
    x[t] ~ dnorm(mu[t],tau_add)
    muD[t] ~ dnorm(mu0,pa0)
    D[t] ~ dbern(p)
    mu[t] <- D[t]*muD[t] + (1-D[t])*muN[t]
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(t_obs,a_obs)
  tau_add ~ dgamma(a_add,t_add)
  R ~ dnorm(rmean,rprec)
  p ~ dunif(0,1)
  mu0 ~ dnorm(-5,1) #mean of disturbed state
  pa0 ~ dgamma(1,1) #precision of disturbed state
}
"

#data and parameters for each site:
data<-list()
j.pests<-list()
j.pests.out<-list()
for (i in geoID){
data[[i]] = list(y=obs[[i]], n=NT,
            x_ic=0, tau_ic=0.1,
            a_obs=0.1,t_obs=0.1,
            a_add=0.1,t_add=0.1,
            rmean=0,rprec=0.00001)
}

#initial state of model parameters
nchain <- 3
for (i in geoID){
  samp<- sample(!is.na(obs[[i]]),length(obs[[i]]),replace=TRUE)
  for(j in 1:nchain){
    init<-list(tau_add=1/var(diff(samp)),tau_obs=1/var(samp))
    }

  j.pests[[i]] <- jags.model (file = textConnection(pests),
                      data = data[[i]],
                      inits = init,
                      n.chains = 3)
  j.pests.out[[i]] <- coda.samples (model = j.pests[[i]],
                            variable.names = c("x","tau_add","tau_obs", "R", "p", "D", "mu0", "pa0"),
                            n.iter = 5000)
}
#plot(j.pests.out)

out<-list()
x.cols<-list()
d.cols<-list()
ci.x<-list()
ci.d<-list()
for (i in geoID){
  out[[i]] <- as.matrix(j.pests.out[[i]])
  x.cols[[i]] <- grep("^x",colnames(out[[i]]))
  ci.x[[i]] <- apply(out[[i]][,x.cols[[i]]],2,quantile,c(0.025,0.5,0.975))

  d.cols[[i]] <- grep("^D",colnames(out[[i]]))
  ci.d[[i]] <- apply(out[[i]][,d.cols[[i]]],2,quantile,c(0.025,0.5,0.975))
}

#i=10
plot(time,ci.x[[i]][2,],type='n',ylim=range(obs[[i]],na.rm=TRUE),ylab="Forest Condition",main=i)
ecoforecastR::ciEnvelope(time,ci.x[[i]][1,],ci.x[[i]][3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,obs[[i]],pch="+",cex=0.5)
ecoforecastR::ciEnvelope(time[-1],ci.d[[i]][1,],ci.d[[i]][3,],col=ecoforecastR::col.alpha("hot pink",0.75))


#distributions of parameters:
rs<-matrix(NA,nrow=15000,ncol=50)
for (i in geoID){
  rs[,i]<-out[[i]][,"R"]
}
hist(rs)

ps<-matrix(NA,nrow=15000,ncol=50)
for (i in geoID){
  ps[,i]<-out[[i]][,"p"]
}
hist(ps)

mu0s<-matrix(NA,nrow=15000,ncol=50)
for (i in geoID){
  mu0s[,i]<-out[[i]][,"mu0"]
}
hist(mu0s)

pa0s<-matrix(NA,nrow=15000,ncol=50)
for (i in geoID){
  pa0s[,i]<-out[[i]][,"pa0"]
}
hist(pa0s)



