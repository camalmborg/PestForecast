#this code is for multivariate analyses for spongy moth Chapter 1

#plan:
#make matrix of combinations and set a loop to run them
#multivariate gams

# necessary libraries:
library(mgcv)
library(combinat)

#load data:
load("DMVARS_MO.RData")
load("DM_MAXTEMPS_MO.RData") #maxtemp_mos<- c(46,47,48,49,50,51,58,59,60,61,62,63)
load("MAM_spring_precip.RData") #col 6 for 2016
load("MA_spring_precip.RData") #cols 5,6 for 2015, 2016
load("AM_spring_precip.RData") #cols 4,6 for 2014, 2016
load("JJA_precip.RData") #cols 4,5 for 2014, 2015
load("JJ_precip.RData") #cols 4,5 for 2014, 2015
load("MJJ_precip.RData") #cols 4,5 for 2014, 2015
load("MJ_precip.RData") #cols 4,5 for 2014, 2015
load("SON_precip.RData") # col 4 for 2014

#load mags:
dmr <- read.csv("SM_distmagrecov_data.csv")

#combining all data into on dataframe:
MV_DATA <- cbind(dm_maxtemps[,c(46,47,48,49,50,51,58,59,60,61,62,63)],
                 MAM_spring_precip[,6],
                 MA_spring_precip[,5:6],
                 AM_spring_precip[,5:6],
                 JJA_precip[,4:5],
                 JJ_precip[,4:5],
                 MJJ_precip[,4:5],
                 MJ_precip[,4:5],
                 SON_precip[,4])

#making a combinations loop:
combi <- t(combn(ncol(MV_DATA),2))
vardat <- as.data.frame(cbind(dmr$mags, MV_DATA[,c(combi[1,])]))
mv_gam <- gam(dmr$mags~s(vardat[,1])+s(vardat[,2]),data=vardat)

SM_multi_var <- function(data, nvars, dmrdat) {
  #make empty matrix:
  r2s <- matrix(nrow=1, ncol=1) #need to know dims
  
  #make combinations of nvars variables:
  combi <- t(combn(ncol(data),nvars))
  
  #vardattest <- list()
  #vardat loop:
  for (i in 1:nrow(combi)){
    #make gam variables data frame:
    vardat[[i]] <- as.data.frame(cbind(dmrdat$mags, data[,c(combi[i,])]))
    
    #make gam explantory variables list
    ex_vars <- c()
    for (j in 1:nvars){
      ex_vars[j] <- paste0('+s(vardat[,', j+1, '])')
      ex_vars[1] <- gsub('\\+','', ex_vars[1])
    }
    
    #make a single string:
    ex_var <- paste(ex_vars, collapse='')
    
    #run gam with those data:
    mv_gam <- gam(vardat[,1]~as.factor(ex_var),data=vardat)
  }
}

