# This new script contains all the code for running both the univariate and
# multivariate analyses for Chapter 1 - disturbance probabilities and magnitudes
# This code will run for Daymet (section 1), SMAP, VIIRS, and DEM data.

### load libraries:
library(mgcv)
library(pROC)
library(combinat)

##### UNIVARIATE ANALYSES SECTION ------------------------------------------------
### Function for Univariate mags Analyses AS LIST  (Daymet data):----------------------
# dmvars = data in list - Daymet data (each variable is list member)
# dmrdat = disturbance magnitude and probability data - object from DISTMAGRECOV 
#          calculation script
# dmr = which from disturbance probability, magnitude, or recovery to use for 
#       analysis - "dpy1"/"dpy2" for 2016/2017 distprob, "mags" for magnitude
# coln = column number of disturbance - 22 for 2016, 23 for 2017

spongy_l_explore<-function(dmvars,dmrdat,dmr,coln){
  #make empty matrix:
  r2s <- matrix(NA,nrow=ncol(dmvars[[1]]),ncol=length(dmvars))
  aics <- matrix(NA,nrow=ncol(dmvars[[1]]),ncol=length(dmvars))
  
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
      aics[j,i] <- var.gam$aic
    }
  }
  # return(vardat)
  #return(r2s)
  models <- as.data.frame(cbind(r2s,aics))
  return(models)
}

### Function for univariate mags analysis NOT in list mode:----------------------------
# var = variable of interest (DEM, VIIRS, SMAP etc.)
# dmrdat = disturbance magnitude and probability data - object from DISTMAGRECOV 
#          calculation script
# dmr = which from disturbance probability, magnitude, or recovery to use for 
#       analysis - "dpy1"/"dpy2" for 2016/2017 distprob, "mags" for magnitude
# coln = column number of disturbance - 22 for 2016, 23 for 2017

spongy_var_explore<-function(var,dmrdat,dmr,coln){
  #make empty matrix:
  r2s <- matrix(NA,nrow=1,ncol=ncol(var))
  aics <- vector()
  
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
    
    #extract r2 and aic:
    summ <- summary(var.gam)
    r2s[,i] <- summ$r.sq
    aics[i] <- var.gam$aic
  }
  # return(vardat)
  r2s <- t(r2s)
  models <- as.data.frame(cbind(r2s, aics))
  return(models)
}

### Functions for disturbance probability univariate analyses --------------------
### List input version (Daymet):---------------------------------------------------
# var = variable of interest (DEM, VIIRS, SMAP etc.)
# dmrdat = disturbance magnitude and probability data - object from DISTMAGRECOV 
#          calculation script
# yr = year disturbance takes place - 1/2 for 2016/2017 respectively
# coln = column number of disturbance - 22 for 2016, 23 for 2017

spongy_lv_ROC <- function(dmvars,dmrdat,yr,coln){
  #make empty matrix:
  rocs <- matrix(NA,nrow=ncol(dmvars[[1]]),ncol=length(dmvars))
  aics <- matrix(NA,nrow=ncol(dmvars[[1]]),ncol=length(dmvars))
  
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
      aics[j,i] <-var.gam$aic
    }
  }
  #return table of AUCs
  #return(rocs)
  models <- as.data.frame(cbind(rocs,aics))
  return(models)
  
}

### Matrix input version: ---------------------------------------------------------
# var = variable of interest (DEM, VIIRS, SMAP etc.)
# dmrdat = disturbance magnitude and probability data - object from DISTMAGRECOV 
#          calculation script
# yr = year disturbance takes place - 1/2 for 2016/2017 respectively
# coln = column number of disturbance - 22 for 2016, 23 for 2017
spongy_var_ROC <- function(var,dmrdat,yr,coln,nvar){
  #make empty matrix:
  rocs <- matrix(NA,nrow=ncol(var),ncol=1)
  aics <- vector()
  
  #grab just disturbance probability columns:
  dists<-dmrdat[,grep("^dp",colnames(dmrdat))]
  
  #first get data:
  vvar <- as.data.frame(var[as.numeric(dmrdat$sitenum),])
  #grab column number:
  cn <- as.matrix(as.numeric(dmrdat$colnum))
  #grab exploratory response  of choice:
  y <- as.matrix(as.numeric(dists[,yr]))
  x <- as.data.frame(cbind(y,cn,vvar))
  vardat <- x[x$cn==coln,]
  
  #loop over all members of dmvars:
  for (i in 1:ncol(var)){  #this is hard coded until further notice
    #loop for filling in R2 table:  
    var.gam<-gam(vardat[,1]~s(vardat[,i+2]), data=vardat, family="binomial")
    var.roc<-roc(var.gam$y,var.gam$fitted.values)
    #var.roc <- roc(var.gam[,1], var.gam$fitted.values)
    rocs[i,] <-var.roc$auc
    aics[i] <- var.gam$aic
  }
  #return table of AUCs
  #return(rocs)
  models <- as.data.frame(cbind(rocs,aics))
  return(models)
}

# NOTE: the roc's only work if there are sufficient non-NA values for the var.gam 
# object. I have changed the line in the code to reflect this, but I am not sure
# how to identify which points it is excluding in the analyses, merely know that it
# excludes some portion and the amounts need to match to get the ROC values.

### MULTIVARIATE ANALYSES SECTION ------------------------------------------------

### Function for Multivariate disturbance magnitude analyses:---------------------
# data = explanatory variables of interest matrix/dataframe
# nvars = number of variables in each model run (numeric)
# dmrdat = disturbance magnitude and probability data - object from DISTMAGRECOV 
#          calculation script

spongy_multi_var <- function(data, nvars, dmrdat) {
  #make combinations of nvars variables:
  combi <- t(combn(ncol(data),nvars))
  
  #make empty matrix:
  r2s <- matrix(nrow=nrow(combi), ncol=2) #need to know dims
  r2s[,1] <- 1:nrow(combi)
  
  aics <- vector()
  
  #vardattest <- list()
  #vardat loop:
  for (i in 1:nrow(combi)){
    
    #first get data:
    vvar <- as.data.frame(data[as.numeric(dmrdat$sitenum),])
    
    #make gam variables data frame:
    vardat <- as.data.frame(cbind(dmrdat$mags, vvar[,c(combi[i,])]))
    
    #make gam explantory variables list
    ex_vars <- c()
    for (j in 1:nvars){
      ex_vars[j] <- paste0('s(vardat[,', j+1, '])')
    }
    
    #make a single string:
    gam_formula <- as.formula(paste("vardat[,1] ~ ",
                                    paste(ex_vars, collapse='+')))
    
    #run gam with those data:
    mv_gam <- gam(gam_formula)
    summ <- summary(mv_gam)
    r2s[i,2] <- summ$r.sq
    aics[i] <- mv_gam$aic
  }
  
  #delta <- min(aics)
  #delAIC <- aics-delta
  models <- as.data.frame(cbind(r2s, combi, aics))
  return(models)
}


### Function for Multivariate disturbance probability analyses:-------------------
# data = explanatory variables of interest matrix/dataframe
# nvars = number of variables in each model run (numeric)
# dmrdat = disturbance magnitude and probability data - object from DISTMAGRECOV 
#          calculation script
# yr = disturbance year (1 = 2016; 2 = 2017)
spongy_multi_ROC <- function(data,nvars,dmrdat,yr){
  #make combinations of nvars variables:
  combi <- t(combn(ncol(data),nvars))
  
  #make empty matrix:
  rocs <- matrix(nrow=nrow(combi), ncol=2) #need to know dims
  rocs[,1] <- 1:nrow(combi)
  
  aics <- vector()
  
  #grab data without missing values:
  #first get data:
  vvar <- as.data.frame(data[as.numeric(dmrdat$sitenum),])
  
  #grab just disturbance probability columns:
  dists<-dmrdat[,grep("^dp",colnames(dmrdat))]
  
  #vardat loop:
  for (i in 1:nrow(combi)){
    #make gam variables data frame:
    vardat <- as.data.frame(cbind(dists[,yr], vvar[,c(combi[i,])]))
    
    #make gam explantory variables list
    ex_vars <- c()
    for (j in 1:nvars){
      ex_vars[j] <- paste0('s(vardat[,', j+1, '])')
    }
    
    #make a single string:
    gam_formula <- as.formula(paste("vardat[,1] ~ ",
                                    paste(ex_vars, collapse='+')))
    
    #run gam with those data:
    mv_gam <- gam(gam_formula, data=vardat, family="binomial")
    #mv_roc <- roc(vardat[,1], mv_gam$fitted.values)
    mv_roc<-roc(mv_gam$y,mv_gam$fitted.values)
    rocs[i,2] <- mv_roc$auc
    aics[i] <- mv_gam$aic
  }
  
  #delta <- min(aics)
  #delAIC <- aics-delta
  models <- as.data.frame(cbind(rocs, combi, aics))
  return(models)
  
}


##### RUNNING THE MODELS--------------------------------------------------------

#UNIVARIATE---------------------------------------------------------------------
#Disturbance Magnitude:
#distmag_univar_16 <- spongy_var_explore(MV_DATA, testfx, "mags", 22)
#distmag_univar_17 <- spongy_var_explore(MV_DATA, textfx, "mags", 23)

# Daymet:
distmag_daymet_monthly_univar_16 <- spongy_l_explore(dmvars_mo, testfx, "mags", 22)
distmag_daymet_monthly_univar_17 <- spongy_l_explore(dmvars_mo, testfx, "mags", 23)
distmag_daymet_monthly_cs_univar_16 <- spongy_l_explore(dmvars_mo, testfx2, "mags", 22)
distmag_daymet_monthly_cs_univar_17 <- spongy_l_explore(dmvars_mo, testfx2, "mags", 23)
# SMAP:
distmag_SMAP_univar_16 <- spongy_var_explore(SMAPdat, testfx, "mags", 22)
distmag_SMAP_univar_17 <- spongy_var_explore(SMAPdat, testfx, "mags", 23)
distmag_SMAP_cs_univar_16 <- spongy_var_explore(SMAPdat, testfx2, "mags", 22)
distmag_SMAP_cs_univar_17 <- spongy_var_explore(SMAPdat, testfx2, "mags", 23)
# VIIRS:
distmag_viirs_univar_16 <- spongy_var_explore(yr_rads, testfx, "mags", 22)
distmag_viirs_univar_17 <- spongy_var_explore(yr_rads, testfx, "mags", 23)
distmag_viirs_cs_univar_16 <- spongy_var_explore(yr_rads, testfx2, "mags", 22)
distmag_viirs_cs_univar_17 <- spongy_var_explore(yr_rads, testfx2, "mags", 23)
#DEM:
distmag_DEM_univar_16 <- spongy_var_explore(site_data, testfx, "mags", 22)
distmag_DEM_univar_17 <- spongy_var_explore(site_data, testfx, "mags", 23)
distmag_DEM_cs_univar_16 <- spongy_var_explore(site_data, testfx2, "mags", 22)
distmag_DEM_cs_univar_17 <- spongy_var_explore(site_data, testfx2, "mags", 23)

#write.csv(distmag_DEM_cs_univar_17, file="Analyses_July2023/distmag_DEM_cs_univar_17.csv")

#Disturbance Probability:
# Daymet:
distprob_daymet_monthly_univar_16 <- spongy_lv_ROC(dmvars_mo, testfx, 1, 22)
distprob_daymet_monthly_univar_17 <- spongy_lv_ROC(dmvars_mo, testfx, 2, 23)
distprob_daymet_monthly_cs_univar_16 <- spongy_lv_ROC(dmvars_mo, testfx2, 1, 22)
distprob_daymet_monthly_cs_univar_17 <- spongy_lv_ROC(dmvars_mo, testfx2, 2, 23)
# SMAP:
distprob_SMAP_univar_16 <- spongy_var_ROC(SMAPdat, testfx, 1, 22)
distprob_SMAP_univar_17 <- spongy_var_ROC(SMAPdat, testfx, 2, 23)
distprob_SMAP_cs_univar_16 <- spongy_var_ROC(SMAPdat, testfx2, 1, 22)
distprob_SMAP_cs_univar_17 <- spongy_var_ROC(SMAPdat, testfx2, 2, 23)
# VIIRS:
distprob_viirs_univar_16 <- spongy_var_ROC(yr_rads, testfx, 1, 22)
distprob_viirs_univar_17 <- spongy_var_ROC(yr_rads, testfx, 2, 23)
distprob_viirs_cs_univar_16 <- spongy_var_ROC(yr_rads, testfx2, 1, 22)
distprob_viirs_cs_univar_17 <- spongy_var_ROC(yr_rads, testfx2, 2, 23)
#DEM:
distprob_DEM_univar_16 <- spongy_var_ROC(site_data, testfx, 1, 22)
distprob_DEM_univar_17 <- spongy_var_ROC(site_data, testfx, 2, 23)
distprob_DEM_cs_univar_16 <- spongy_var_ROC(site_data, testfx2, 1, 22)
distprob_DEM_cs_univar_17 <- spongy_var_ROC(site_data, testfx2, 2, 23)

#write.csv(distprob_DEM_cs_univar_17, file="Analyses_July2023/distprob_DEM_cs_univar_17.csv")

### MULIVARIATE ----------------------------------------------------------------
#2 variables:
dist_mag_2var_tcg <- spongy_multi_var(MV_DATA, 2, testfx)
dist_mag_2var_cs <- spongy_multi_var(MV_DATA, 2, testfx2)

dist_prob_2016_2var_tcg <- spongy_multi_ROC(MV_DATA, 2, testfx, 1)
dist_prob_2016_2var_cs <- spongy_multi_ROC(MV_DATA, 2, testfx2, 1)
  
dist_prob_2017_2var_tcg <- spongy_multi_ROC(MV_DATA, 2, testfx, 2)
dist_prob_2017_2var_cs <- spongy_multi_ROC(MV_DATA, 2, testfx2, 2)

#3 variables:
dist_mag_3var_tcg <- spongy_multi_var(MV_DATA, 3, testfx)
dist_mag_3var_cs <- spongy_multi_var(MV_DATA, 3, testfx2)

dist_prob_2016_3var_tcg <- spongy_multi_ROC(MV_DATA, 3, testfx, 1)
dist_prob_2016_3var_cs <- spongy_multi_ROC(MV_DATA, 3, testfx2, 1)

dist_prob_2017_3var_tcg <- spongy_multi_ROC(MV_DATA, 3, testfx, 2)
dist_prob_2017_3var_cs <- spongy_multi_ROC(MV_DATA, 3, testfx2, 2)

write.csv(dist_prob_2016_2var_cs, file="Analyses_July2023/dist_prob_2016_2var_cs.csv")
