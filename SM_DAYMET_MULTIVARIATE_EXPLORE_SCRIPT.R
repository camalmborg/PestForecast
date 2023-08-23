#this code is for multivariate analyses for spongy moth Chapter 1

#plan:
#make matrix of combinations and set a loop to run them
#multivariate gams

# necessary libraries:
library(mgcv) #for GAM analyses
library(combinat)
library(pROC) #for AUC analyses


#load data:
load("DMVARS_MO.RData")
load("DM_MAXTEMPS_MO.RData") #maxtemp_mos<- c(46,47,48,49,50,51,58,59,60,61,62,63)
#load("DM_MINTEMPS_MO.RData") #mintemp_mos<- c(52, 53)
load("MAM_spring_precip.RData") #col 6 for 2016
load("MA_spring_precip.RData") #cols 5,6 for 2015, 2016
load("AM_spring_precip.RData") #cols 4,6 for 2014, 2016
load("JJA_precip.RData") #cols 4,5 for 2014, 2015
load("JJ_precip.RData") #cols 4,5 for 2014, 2015
load("MJJ_precip.RData") #cols 4,5 for 2014, 2015
load("MJ_precip.RData") #cols 4,5 for 2014, 2015
load("SON_precip.RData") # col 4 for 2014
load("SMAP_data.RData") # col 1-4 for Apr-July 2015 soil moisture
load("site_DEM_slope_aspect_TWI_data.RData") # col 1 for DEM data
load("viirs_annual_averages_data.RData") #not included

#load mags:
dmr <- read.csv("SM_distmagrecov_data.csv")

#combining all data into on dataframe:
MV_DATA <- cbind.data.frame(dm_maxtemps[,c(46,47,48,49,50,51,58,59,60,61,62,63)],
                 MAM_spring_precip[,6],
                 MA_spring_precip[,5:6],
                 AM_spring_precip[,5:6],
                 JJA_precip[,4:5],
                 JJ_precip[,4:5],
                 MJJ_precip[,4:5],
                 MJ_precip[,4:5],
                 SON_precip[,4],
                 SMAPdat[,1:4],
                 site_data[,1])


test_data <- cbind(dm_maxtemps[,c(60,61,62,63)],
                   MAM_spring_precip[,6],
                   SON_precip[,4])

#making a combinations loop:
#combi <- t(combn(ncol(MV_DATA),2))
#vardat <- as.data.frame(cbind(dmr$mags, MV_DATA[,c(combi[1,])]))
#mv_gam <- gam(dmr$mags~s(vardat[,1])+s(vardat[,2]),data=vardat)


### DISTURBANCE MAGNITUDE MULTIVARIATE LOOP------------------------------------------------
SM_multi_var <- function(data, nvars, dmrdat) {
  #make combinations of nvars variables:
  combi <- t(combn(ncol(data),nvars))
  
  #make empty matrix:
  r2s <- matrix(nrow=nrow(combi), ncol=2) #need to know dims
  r2s[,1] <- 1:nrow(combi)
  
  aics <- vector()
  
  #vardattest <- list()
  #vardat loop:
  for (i in 1:nrow(combi)){
    #make gam variables data frame:
    vardat <- as.data.frame(cbind(dmrdat$mags, data[,c(combi[i,])]))
    
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
  
  delta <- min(aics)
  delAIC <- aics-delta
  models <- cbind(r2s, combi, delAIC)
  return(models)
}

test <- SM_multi_var(MV_DATA, 2, dmr)
test3v <- SM_multi_var(MV_DATA, 3, dmr)

test2 <- SM_multi_var(test_data, 2, dmr)

checkmodels <- c(which(test3v[,2]>=0.3))
checking <- text3v[checkmodels,]


### DISTURBANCE PROBABILITY MULTIVARIATE LOOP: ------

SM_multi_ROC <- function(data,nvars,dmrdat,yr){
  #make combinations of nvars variables:
  combi <- t(combn(ncol(data),nvars))
  
  #make empty matrix:
  rocs <- matrix(nrow=nrow(combi), ncol=2) #need to know dims
  rocs[,1] <- 1:nrow(combi)
  
  aics <- vector()
  
  #grab just disturbance probability columns:
  dists<-dmrdat[,grep("^dp",colnames(dmrdat))]
  
  #vardat loop:
  for (i in 1:nrow(combi)){
    #make gam variables data frame:
    vardat <- as.data.frame(cbind(dists[,yr], data[,c(combi[i,])]))
    
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
    mv_roc <- roc(vardat[,1], mv_gam$fitted.values)
    rocs[i,2] <- mv_roc$auc
    aics[i] <- mv_gam$aic
  }
  
  delta <- min(aics)
  delAIC <- aics-delta
  models <- cbind(rocs, combi, delAIC)
  return(models)
  
}


testing_roc_s <- SM_multi_ROC(data,2,testfx,1) #recall if yr=1, coln=22 (2016)


### LET'S LOOK AT THESE MODELS, SHALL WE? -------------------------------------

#Univariate - 7/10/2023 Analysis
dist_mag_tcg <- read.csv("Analyses_July2023/Dist_Mag_TCG_2016_Models_R2s_AICs.csv")
dist_prob_tcg <- read.csv("Analyses_July2023/Dist_Prob_TCG_2016_Models_ROCs_AICs.csv")
dist_mag_cs <- read.csv("Analyses_July2023/Dist_Mag_CS_2016_Models_R2s_AICs.csv")
dist_prob_cs <- read.csv("Analyses_July2023/Dist_Prob_CS_2016_Models_ROCs_AICs.csv")

#multivariate - 7/11/2023 Analyses


#if the model results need to be loaded you can load them here:
## Dist Mag:
#read.csv("Analyses_June2023/2023_06_28_multivariate_analyses_3vars.csv")
## Dist Prob:
#read.csv("Analyses_June2023/2023_06_28_distprob_multivariate_analyses_3var_2016.csv")

# which model:
model <- spongy_multi_3var
model <- dist_mag_tcg

### finding the top performers by DIC/AIC:
best <- which(model$delAIC == 0)

### most common variables included:
#subset data by R2/delAIC:
topmodels <- model[model$V2 > 0.27,]
topmodels <- model[model$delAIC < 100,]
modAIC <- model[order(model$delAIC),]
modr2 <- model[order(model$V2,decreasing = T),]
names(which.max(table(topmodels[,3])))
varvars <- c(sort(unique(topmodels$V3)), 
             sort(unique(topmodels$V4)), 
             sort(unique(topmodels$V5)))
sort(unique(varvars))


### Let's try to combine the AICs for all?
### Starting with Univariate:
