# This code is for the daymet exploratory analyses

# necessary libraries:
library(mgcv)

#function for univariate analyses:
dm_explore<-function(dmvars,dmrdat,dmr,coln){
  #make empty matrix:
  r2s <- matrix(NA,nrow=ncol(dmvars[[1]]),ncol=length(dmvars))
  
  #loop over all members of dmvars list:
  for (i in 1:length(dmvars)){
    #first extract the list you want:
    dmvariable <- as.data.frame(dmvars[[i]][as.numeric(dmrdat$sitenum),])
    #grab column number:
    cn <- as.matrix(as.numeric(dmrdat$colnum))
    #grab exploratory response  of choice:
    y <- as.matrix(as.numeric(dmrdat[,dmr]))
    x <- as.data.frame(cbind(y,cn,dmvariable))
    vardat <- x[x$cn==coln,]
  
    #loop for filling in R2 table:  
      for (j in 1:ncol(dmvariable)){
        var.gam <- gam(vardat[,1]~s(vardat[,j+2]),data=vardat)
        
        #extract r2:
        summ <- summary(var.gam)
        r2s[j,i] <-summ$r.sq
      }
   }
  return(r2s)
}

testing <- dm_explore(dmvars,testfx,"mags",22)


##---------------------------------------------------------------------####
###Multivariate Analyses:



