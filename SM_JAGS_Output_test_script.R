#load libraries
library(rjags)
library(MCMCvis)

#testing the splitting params within loop method:

#select model:
model1<-"Model_1_mu0_jpout_run2_"
#model2<-"Model_2_precip5k_jpout_"
#model3<-"Model_3_fullenv5k_jpout_"

#select parameters:
param1<-c("beta0","tau_add","tau_obs", "pa0")
#param2<-c("beta[3]","beta[4]","tau_add","tau_obs", "pa0")
#param3<-c("beta[1]","beta[2]","beta[3]","beta[4]","tau_add","tau_obs", "pa0")

#MCMC VIS VERSION
#the following code has two loops, one for saving the parameters for the model
#and the other for saving x's and D's for certain sites

###LOOP FOR THE MODELPARAMS----
outptemp<-matrix(NA,ncol=length(param1)) #1, 2, or 3
#run loop for collecting parameters:
for (i in 20:5){
  #load jpout:
  load(file = paste0(model1,as.character(i),".RData")) #model 1, 2, or 3
  
  #split output:
  jpps<-MCMCchains(jpout, params=param1) #choose right params
  
  #grab sites for x and D:
  outptemp<-rbind(outptemp,jpps)
  
  #remove big stuff
  rm(jpps)
  rm(jpout)
}

outp<-outptemp[-1,] #removes the NA row
rm(outptemp)
#summaries of model parameters:
summ1<-as.data.frame.matrix(summary(outp))
#summ2<-as.data.frame.matrix(summary(outp))
#summ3<-as.data.frame.matrix(summary(outp))

###x and D loop----
#choose site:
siten=22
#empty matrices:
outxtemp<-matrix(NA,ncol=NT)
outdtemp<-matrix(NA,ncol=NT-1)

for (i in 20:5){
  #load jpout:
  load(file = paste0(model1,as.character(i),".RData"))
  
  #split output:
  #jpps<-MCMCchains(jpout, params=param1)
  jpds<-MCMCchains(jpout, params="D")
  jpxs<-MCMCchains(jpout, params = "x")
  
  #grab sites for x and D:
  xparam<-which(colnames(jpxs) %in% paste0("x[",siten,",",1:130,"]"))
  dparam<-which(colnames(jpds) %in% paste0("D[,",1:130,"]"))
  xsite<-jpxs[,xparam]
  dsite<-jpds[,dparam]
  #outptemp<-rbind(outptemp,jpps)
  outxtemp<-rbind(outxtemp,xsite)
  outdtemp<-rbind(outdtemp,dsite)
  
  #remove big stuff
  #rm(jpps)
  rm(jpxs)
  rm(jpds)
  rm(xsite,dsite)
  rm(jpout)
}

#outp<-outptemp[-1,] #removes the NA row
outx<-outxtemp[-1,]
outd<-outdtemp[-1,]
#rm(outptemp)
rm(outxtemp)
rm(outdtemp)


##VISUALIZATIONS-----
#if you need the data again:
#condscores.samp<-data$y

#confidence intervals:
ci.x.1 <- apply(outx,2,quantile,c(0.025,0.5,0.975))
ci.d.1 <- apply(outd,2,quantile,c(0.025,0.5,0.975))
# ci.x.2 <- apply(outx,2,quantile,c(0.025,0.5,0.975))
# ci.d.2 <- apply(outd,2,quantile,c(0.025,0.5,0.975))
# ci.x.3 <- apply(outx,2,quantile,c(0.025,0.5,0.975))
# ci.d.3 <- apply(outd,2,quantile,c(0.025,0.5,0.975))

#details for quick file naming:
siten<-22
mod<-"Mod1"
modelname<-"Model 1"

#time for plots:
time=1:130
NT=length(time)

#save single model plot details:
tiff(paste0("509", mod, "_examp_site", siten,".tiff"), 
     units="in", width=12, height=5, res=300)
plot(ci.x.1[2,],type='l',
     ylim=c(min(ci.x.1,na.rm=TRUE),max(condscores.samp,na.rm=TRUE)),
     main=paste0(modelname),
     ylab="Forest Condition Score",
     col="red",
     xlab="Month",
     cex=1)
ecoforecastR::ciEnvelope(time,ci.x.1[1,],ci.x.1[3,],col=ecoforecastR::col.alpha("indianred1",0.60))
points(time,condscores.samp[i,],pch=16,cex=0.5,col="navyblue")
dev.off()

#save single model plot details: for D's
tiff(paste0("509", mod, "_examp_site", siten,"distprob.tiff"), 
     units="in", width=10, height=5, res=300)
plot(time[2:NT],ci.d.1[2,],type='l',ylim=range(ci.d.1,na.rm=TRUE),
     main=paste0(modelname," :Disturbance Probability"),
     ylab="Disturbance Probability",
     col="red",
     xlab="Month",
     cex=1)
ecoforecastR::ciEnvelope(time[2:NT],ci.d.1[1,],ci.d.1[3,],col=ecoforecastR::col.alpha("indianred1",0.60))
dev.off()

#save all models plot details:
# tiff(paste0("509_allmodels_", siten, ".tiff"), 
#      units="in", width=10, height=5, res=300)
# plot:
# plot(time,ci.x.1[2,],type='l',ylim=c(-12,5),
#      ylab="Forest Condition Score", 
#      col="red",
#      xlab="Month",
#      cex=1)
# lines(time,ci.x.2[2,], col="blue")
# lines(time,ci.x.3[2,], col="dark green")
# ecoforecastR::ciEnvelope(time,ci.x.1[1,],ci.x.1[3,],col=ecoforecastR::col.alpha("indianred1",0.30))
# ecoforecastR::ciEnvelope(time,ci.x.2[1,sitei],ci.x.2[3,sitei],col=ecoforecastR::col.alpha("lightblue1",0.30))
# ecoforecastR::ciEnvelope(time,ci.x.3[1,sitei],ci.x.3[3,sitei],col=ecoforecastR::col.alpha("greenyellow",0.30))
# points(time,condscores.samp[siten,],pch=16,cex=0.5,col="navyblue")
# dev.off()

