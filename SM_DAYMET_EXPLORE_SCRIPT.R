# This code is for the daymet exploratory analyses

# necessary libraries:
library(mgcv)

dm_explore<-function(dmvars,dmrdat,dmr,){
  #loop over all members of dmvars list:
  for (i in 1:length(dmvars)){
    #first extract the list you want:
    dmvariable <- dmvars[[1]][as.numeric(dmrdat$sitenum),]
    #grab column number:
    cn <- as.matrix(as.numeric(dmr$colnum))
    #grab exploratory response  of choice:
    y <- as.matrix(as.numeric(dmr$mags))
    vardat<-as.data.frame(cbind(y,cn,dmvariable))
    
    #make empty table:
    r2s<-matrix(NA,nrow=nrow(dmvariable),ncol=ncol(dmvariable))
      
      for (j in 1:nrow(vardat)){
        for (k in 3:nrow(vardat)){
          var.gam <- gam(vardat[j,1]~s(vardat[j,k]),data=vardat)
        }
      }
    #extract R2:
    summ<-summary(var.gam)
    r2 <- summ$r.sq
    print(r2)
  }
}