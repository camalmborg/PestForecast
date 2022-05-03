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


