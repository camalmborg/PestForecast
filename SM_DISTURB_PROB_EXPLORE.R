#code for disturbance probability exploratory analyses

###Load libraries:
library(pROC) #for AUC analyses
library(mgcv)

### Load condition score data:
cfile<-"2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv"
condscores<-read.csv(cfile)

spongy_ROC <- function(dmvars,dmrdat,yr,coln){
  #make empty matrix:
  rocs <- matrix(NA,nrow=ncol(dmvars[[1]]),ncol=length(dmvars))
  
  #grab just disturbance probability columns:
  dists<-dmrdat[,grep("^dp",colnames(dmrdat))]

  #loop over all members of dmvars list:
  for (i in 1:length(dmvars)){
    #first extract the list you want:
    dmvariable <- as.data.frame(dmvars[[i]][as.numeric(dmrdat$sitenum),])
    #grab column number:
    cn <- as.matrix(as.numeric(dmrdat$colnum))
    #grab exploratory response  of choice:
    y <- as.matrix(as.numeric(dists[,yr]))
    x <- as.data.frame(cbind(y,cn,dmvariable))
    vardat <- x[x$cn==coln,]
    
    #loop for filling in R2 table:  
    for (j in 1:ncol(dmvariable)){
      var.gam<-gam(vardat[,1]~s(vardat[,j+2]), data=vardat, family="binomial")
      var.roc<-roc(vardat[,1],var.gam$fitted.values)
      rocs[j,i] <-var.roc$auc
    }
  }
  #return table of AUCs
  return(rocs)

}


############ ARCHIVE:
# spongy_dp<-function(cs,){
#   junes<-cs[,grep("[:.:]06",colnames(cs))]
#   
# }

#junetest<-condscores[,grep("[:.:]06",colnames(condscores))]
#yeartest<-testfx[,grep("^dp",colnames(testfx))]
