# This new script contains all the code for running both the univariate and
# multivariate analyses for Chapter 1 - disturbance probabilities and magnitudes
# This code will run for Daymet (section 1), SMAP, VIIRS, and DEM data.

### load libraries:
library(mgcv)
library(pROC)
library(combinat)

##### UNIVARIATE ANALYSES SECTION ------------------------------------------------

### Function for Univariate Analyses AS LIST  (Daymet data):----------------------
# dmvars = data in list - Daymet data (each variable is list member)
# dmrdat = disturbance magnitude and probability data - object from DISTMAGRECOV 
#          calculation script
# dmr = which from disturbance probability, magnitude, or recovery to use for 
#       analysis - "dpy1"/"dpy2" for 2016/2017 distprob, "mags" for magnitude
# coln = column number of disturbance - 22 for 2016, 23 for 2017

SM_l_explore<-function(dmvars,dmrdat,dmr,coln){
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
  # return(vardat)
  return(r2s)
}

### Function for univariate analysis NOT in list mode:----------------------------
# var = variable of interest (DEM, VIIRS, SMAP etc.)
# dmrdat = disturbance magnitude and probability data - object from DISTMAGRECOV 
#          calculation script
# dmr = which from disturbance probability, magnitude, or recovery to use for 
#       analysis - "dpy1"/"dpy2" for 2016/2017 distprob, "mags" for magnitude
# coln = column number of disturbance - 22 for 2016, 23 for 2017

SM_var_explore<-function(var,dmrdat,dmr,coln){
  #make empty matrix:
  r2s <- matrix(NA,nrow=1,ncol=ncol(var))
  
  #loop over all members of dmvars list:
  for (i in 1:ncol(var)){
    #first extract the list you want:
    varvariables <- as.data.frame(var[as.numeric(dmrdat$sitenum),])
    #grab column number:
    cn <- as.matrix(as.numeric(dmrdat$colnum))
    #grab exploratory response  of choice:
    y <- as.matrix(as.numeric(dmrdat[,dmr]))
    x <- as.data.frame(cbind(y,cn,varvariables))
    vardat <- x[x$cn==coln,]
    
    #run the GAM:
    var.gam <- gam(vardat[,1]~s(vardat[,i+2]),data=vardat)
    
    #extract r2:
    summ <- summary(var.gam)
    r2s[,i] <-summ$r.sq
  }
  # return(vardat)
  return(r2s)
}

