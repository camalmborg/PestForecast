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

#CODA SPLIT VERSION-----
# i=1
# load(file = paste0(model1,as.character(i),".RData"))
# #make empty matrix with sample mcmc object's # of columns (1 per tracked param)
# #out<-matrix(NA,ncol=ncol(jpout[[1]]))
# 
# #outp<-matrix(NA,ncol=4)
# outxs<-matrix(NA,ncol=(ncol(jpout[[1]])-4)/2)
# #grab each mcmc object of that model and convert to matrix:
# for (i in 20:19){
#   #load jpout:
#   load(file = paste0(model1,as.character(i),".RData"))
#   
#   #split output:
#   codaSplit <- function(jpout,pattern){
#     out = list()
#     mfit = as.matrix(jpout,chains=TRUE)
#     pat.cols = grep(pattern,colnames(mfit),fixed=TRUE)
#     chain.col = which(colnames(mfit)=="CHAIN")
#     out[[1]] = mat2mcmc.list(mfit[,c(chain.col,pat.cols)])
#     out[[2]]   = mat2mcmc.list(mfit[,-pat.cols])
#     return(out)
#   }
#   
#   mat2mcmc.list <- function(w) {
#     temp <- list()
#     chain.col <- which(colnames(w) == "CHAIN")
#     for (i in unique(w[, "CHAIN"])) {
#       temp[[i]] <- coda:::as.mcmc(w[w[, "CHAIN"] == i, -chain.col])
#     }
#     return(as.mcmc.list(temp))
#   } ##to use: out <- codaSplit(jags.out,"like")
#   
#   #jpmx<-as.matrix(jpout)
#   #out<-rbind(out,jpmx)
#   #rm(jpmx)
#   #outr<-out[-1]
#   jpxs<-codaSplit(jpout,"x")
#   xs<-as.matrix(jpxs)
#   outx<-rbind(outxs,xs)
#   rm(jpxs)
# }
# outx1<-outx[-1,] #removes the NA row
# #out2<-out[-1,]
# #out3<-out[-1,]
# rm(jpout)
# rm(outp)


#MCMC VIS VERSION-----
outptemp<-matrix(NA,ncol=length(param1))
outxtemp<-matrix(NA,ncol=NT)
outdtemp<-matrix(NA,ncol=NT-1)

for (i in 20:5){
  #load jpout:
  load(file = paste0(model1,as.character(i),".RData"))
  
  #split output:
  jpps<-MCMCchains(jpout, params=param1)
  jpds<-MCMCchains(jpout, params="D")
  jpxs<-MCMCchains(jpout, params = "x")
  
  #grab sites for x and D:
  xparam<-which(colnames(jpxs) %in% paste0("x[22,",1:130,"]"))
  dparam<-which(colnames(jpds) %in% paste0("D[22,",1:130,"]"))
  xsite<-jpxs[,xparam]
  dsite<-jpds[,dparam]
  outptemp<-rbind(outptemp,jpps)
  outxtemp<-rbind(outxtemp,xsite)
  outdtemp<-rbind(outdtemp,dsite)
  
  #remove big stuff
  rm(jpps)
  rm(jpxs)
  rm(jpds)
  rm(xsite,dsite)
  rm(jpout)
}

outp<-outptemp[-1,] #removes the NA row
outx<-outxtemp[-1,]
outd<-outdtemp[-1,]
rm(outptemp)
rm(outxtemp)
rm(outdtemp)

#summaries of model parameters:
summ1<-as.data.frame.matrix(summary(outp))
#summ2<-as.data.frame.matrix(summary(outp))
#summ3<-as.data.frame.matrix(summary(outp))

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
     units="in", width=10, height=5, res=300)
plot(ci.x.1[2,],type='l',ylim=range(condscores.samp,na.rm=TRUE),
     main=paste0(modelname),
     ylab="Forest Condition Score",
     col="red",
     xlab="Month",
     cex=1)
ecoforecastR::ciEnvelope(time,ci.x.1[1,],ci.x.1[3,],col=ecoforecastR::col.alpha("indianred1",0.60))
points(time,condscores.samp[i,],pch=16,cex=0.5,col="navyblue")
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

