library(ecoforecastR)
library(rjags)
library(coda)

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

#load condition score means for each site into matrix:
obs<-matrix(NA,nrow=length(geoID),ncol=130)
for (i in geoID){
  obs[i,]<-gm.data[[i]]$score_mean
}

#----model with disturbance probability added:------
# pests.RW = "
# model{
#   
#   #### Data Model
#   for(t in 1:n){
#     y[t] ~ dnorm(x[t],tau_obs)
#   }
#   
#   #### Process Model
#   for(t in 2:n){
#     muN[t]<-R*x[t-1]
#     x[t] ~ dnorm(mu[t],tau_add)
#     muD[t] ~ dnorm(mu0,pa0)
#     D[t] ~ dbern(p)
#     mu[t] <- D[t]*muD[t] + (1-D[t])*muN[t]
#   }
#   
#   #### Priors
#   x[1] ~ dnorm(x_ic,tau_ic)
#   tau_obs ~ dgamma(t_obs,a_obs)
#   tau_add ~ dgamma(a_add,t_add)
#   R ~ dnorm(rmean,rprec)  #rho growth term
#   p ~ dunif(0,1)  #disturbance probability
#   mu0 ~ dnorm(-5,1) #mean of disturbed state
#   pa0 ~ dgamma(1,1) #precision of disturbed state
# }
# "

# #data for pests model:
# data.RW = list(y=obs, n=NT,
#                x_ic=0, tau_ic=0.1,
#                a_obs=0.1,t_obs=0.1,
#                a_add=0.1,t_add=0.1,
#                rmean=0,rprec=0.00001)
# 
# init<-list()
# #j.pests<-list()
# #j.pests.out<-list()
# nchain <- 3
# for (i in geoID){
#   for(j in 1:nchain){
#     samp<- sample(!is.na(obs[[i]]),length(obs[[i]]),replace=TRUE)
#     init[[j]]<-list(tau_add=1/var(diff(samp)),tau_obs=1/var(samp))
#   }
# 
# j.pests[[i]] <- jags.model (file = textConnection(pests),
#                           data = data,
#                           inits = init,
#                           n.chains = 3)
# j.pests.out[[i]] <- coda.samples (model = j.pests[[i]],
#                                 variable.names = c("x","tau_add","tau_obs", "R", "p", "D", "mu0", "pa0"),
#                                 n.iter = 5000)
# }

#----disturbance prob model over sites------
pests = "
model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model
  for(t in 1:n){
    y[s,t] ~ dnorm(x[s,t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    muN[s,t]<-R*x[s,t-1]
    x[s,t] ~ dnorm(mu[s,t],tau_add)
    muD[s,t] ~ dnorm(mu0,pa0)  #adding a process model on mu0[s,t] 
    D[s,t] ~ dbern(p) ##can add a process model here for the natural population cycle
    mu[s,t] <- D[s,t]*muD[s,t] + (1-D[s,t])*muN[s,t] 
  
  x[s,1]~dnorm(x_ic,tau_ic)
  
}#end loop over sites
  
  #### Priors
  tau_obs ~ dgamma(t_obs,a_obs)
  tau_add ~ dgamma(a_add,t_add)
  R ~ dnorm(rmean,rprec)  #rho term
  p ~ dunif(0,1)  #disturbance probability
  mu0 ~ dnorm(-5,1) #mean of disturbed state
  pa0 ~ dgamma(1,1) #precision of disturbed state
}
"

#data and parameters for sites model:
data.s = list(y=obs, n=NT, ns=length(geoID),
            x_ic=0, tau_ic=0.1,
            a_obs=0.1,t_obs=0.1,
            a_add=0.1,t_add=0.1,
            rmean=0,rprec=0.00001)

#j.pests.RE<-list()
#j.pests.out.RE<-list()

#initial state of model parameters
init<-list()
nchain <- 3
for(j in 1:nchain){
  samp<- sample(!is.na(obs),length(obs),replace=TRUE)
  init[[j]]<-list(tau_add=1/var(diff(samp)),tau_obs=1/var(samp))
}

j.pests <- jags.model (file = textConnection(pests.m),
                      data = data.m,
                      inits = init,
                      n.chains = 3)
j.pests.out <- coda.samples (model = j.pests,
                            variable.names = c("mu0","p","pa0","tau_add"),    #c("x","tau_add","tau_obs", "R", "p", "D", "mu0", "pa0"),
                            n.iter = 10000,
                            thin=10)
#j.pests.out.thin<-window(j.pests.out,thin=10)

# j.pests.out.2500<-coda.samples(model = j.pests,
#                                variable.names = c("x","tau_add","tau_obs", "R", "p", "D", "mu0", "pa0"),
#                                n.iter=10000)

###-----model output-----
out.pests <- as.matrix(j.pests.out)
out.pests.thin <-as.matrix(j.pests.out.thin)

###------diagnostic----------------
#Brooks Gelman Rubin test:
BGR<-gelman.diag(j.pests.out.thin)

#autocorrelation:
acfplot(j.pests.out)

#effective size:
#EffS<-effectiveSize(j.pests.out)
#cumuplot(j.pests.out,probs=c(0.025,0.25,0.5,0.75,0.975))


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

###splitting output-----still doesn't work
out = list(params=NULL,predict=NULL)
mfit = as.matrix(j.pests.out,chains=TRUE)
pred.cols = union(grep("x",colnames(mfit),fixed=TRUE),
                  grep("D",colnames(mfit),fixed=TRUE))
chain.col = which(colnames(mfit)=="CHAIN")
out$predict = mat2mcmc.list(mfit[,c(chain.col,pred.cols)])
out$params   = mat2mcmc.list(mfit[,-pred.cols])
#if(dic) out$DIC <- dic.samples(mc3, 2000)
#return(out)

#checking convergence
#paramconv <- coda.samples (model = j.pests,
                            #variable.names = c("mu0","p","pa0","tau_add"),
                            #n.iter = 10000)
#plot(paramconv)

#-----visualizations:------

x.cols <- grep("^x",colnames(out.pests))
ci.x <- apply(out.pests[,x.cols],2,quantile,c(0.025,0.5,0.975))
ci.x.names = parse.MatrixNames(colnames(ci.x),numeric=TRUE)


d.cols <- grep("^D",colnames(out.pests))
ci.d <- apply(out.pests[,d.cols],2,quantile,c(0.25,0.5,0.975))

i=48
sitei = which(ci.x.names$row == i)

#tiff("timeseriesexamp.tiff", units="in", width=8, height=3, res=300)
plot(ci.x[2,sitei],type='l',ylim=range(obs,na.rm=TRUE),
     ylab="Forest Condition Score", 
     col="black",
     xlab="Month",
     cex=1)
ecoforecastR::ciEnvelope(time,ci.x[1,sitei],ci.x[3,sitei],col=ecoforecastR::col.alpha("lightBlue",0.60))
points(time,obs[i,],pch="+",cex=0.5,col="navyblue")
#dev.off()

tiff("timeseriesexamp2.tiff", units="in", width=8, height=3, res=300)
plot(time,obs[i,],type='l',ylim=range(obs,na.rm=TRUE),
     ylab="Forest Condition Score", 
     col="black",
     xlab="Month",
     cex=1)
dev.off()


#ecoforecastR::ciEnvelope(time,ci.d[i,sitei][1,],ci.d[i,sitei][3,],col=ecoforecastR::col.alpha("hot pink",0.75))

# out<-list()
# x.cols<-list()
# d.cols<-list()
# ci.x<-list()
# ci.d<-list()
# for (i in geoID){
#   out[[i]] <- as.matrix(j.pests.out[[i]])
#   x.cols[[i]] <- grep("^x",colnames(out[[i]]))
#   ci.x[[i]] <- apply(out[[i]][,x.cols[[i]]],2,quantile,c(0.025,0.5,0.975))
# 
#   d.cols[[i]] <- grep("^D",colnames(out[[i]]))
#   ci.d[[i]] <- apply(out[[i]][,d.cols[[i]]],2,quantile,c(0.025,0.5,0.975))
# }
# 
# #i=10
# plot(time,ci.x[[i]][2,],type='n',ylim=range(obs[[i]],na.rm=TRUE),ylab="Forest Condition",main=i)
# ecoforecastR::ciEnvelope(time,ci.x[[i]][1,],ci.x[[i]][3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
# points(time,obs[[i]],pch="+",cex=0.5)
# ecoforecastR::ciEnvelope(time[-1],ci.d[[i]][1,],ci.d[[i]][3,],col=ecoforecastR::col.alpha("hot pink",0.75))
# 

#distributions of parameters:
rs<-matrix(NA,nrow=15000,ncol=50)
for (i in geoID){
  rs[,i]<-out.pests[,"R"]
}
hist(rs)

ps<-matrix(NA,nrow=15000,ncol=50)
for (i in geoID){
  ps[,i]<-out.pests[,"p"]
}
hist(ps)

mu0s<-matrix(NA,nrow=15000,ncol=50)
for (i in geoID){
  mu0s[,i]<-out.pests[,"mu0"]
}
hist(mu0s)

pa0s<-matrix(NA,nrow=15000,ncol=50)
for (i in geoID){
  pa0s[,i]<-out.pests[,"pa0"]
}
hist(pa0s)


#looking at tau values @ example site 48:
hist(1/sqrt(out.pests[,'tau_add']))
hist(1/sqrt(out.pests[,'tau_obs']))
plot(out.pests[,'tau_add'],out.pests[,'tau_obs'])

