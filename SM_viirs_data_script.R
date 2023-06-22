### script for VIIRS download(?) and processing
#this is the script for processing the viirs data. Data are monthly composites
#of nighttime radiance from GEE

#load libraries:
library(tidyr)
library(stringr)
library(tidyverse)

#load data:
viirs_data <- read.csv("viirs_data/2023_06_21_5000sample_viirs_3.csv")

#geographic coordinates from GEE extract:
geo<-as.data.frame(viirs_data[,".geo"])

#make lat and lon columns from .geo data:
coords<-matrix(nrow=nrow(geo),ncol=2)
for (i in 1:nrow(geo)){
  #longitudes:
  lon<-str_extract(geo[i,], "\\d+\\.*\\d*")
  coords[i,1]<-as.numeric(lon)*-1
  
  #latitudes:
  extlon<-sub(lon,"",geo[i,])
  coords[i,2]<-as.numeric(str_extract(extlon, "\\d+\\.*\\d*"))
}
colnames(coords)<-c("lon","lat")

#viirs_RC <- viirs <- viirs_data[,2:217]
viirs_rads <- viirs_data[,c(grep("avg",colnames(viirs_data)))]
viirs_coms <- viirs_data[,c(grep("cvg",colnames(viirs_data)))]

#Finding columns where all values == 0:
# zeros <- list()
# mins <- vector()
# for (i in 1:ncol(viirs_rads)){
#        zeros[[i]]<-which(viirs_rads[,i]==0)
#        mins[i]<-min(viirs_rads[,i])
#    }

#all zeros 
# NOTE: for all times cf_cvg==0, avg_rad == 0, so I am making NA's for all 0's
#making annual averages but skipping 0's:
viirs_rads[viirs_rads == 0] <- NA
#names(viirs_rads) <- gsub("[^.-0-9]", "", names(viirs_rads), fixed = TRUE)

#selecting years and averaging:
vseq <- seq(1, ncol(viirs_rads), by=12)
yr_rads <- matrix(nrow=nrow(viirs_rads), ncol=length(vseq))
for (i in 1:length(vseq)){
  yr_rads[,i] <- apply(viirs_rads[,(vseq[i]):(vseq[i]+11)], 1, mean, na.rm=T)
}

###----TIME FOR SOME ANALYSES----------------------------------------

library(mgcv)

testfx <- read.csv("SM_distmagrecov_data.csv")

#function for univariate analyses:
VIIRS_explore<-function(viirs,dmrdat,dmr,coln){
  #make empty matrix:
  r2s <- matrix(NA,nrow=1,ncol=ncol(viirs))
  
  #loop over all members of dmvars list:
  for (i in 1:ncol(viirs)){
    #first extract the list you want:
    vvariables <- as.data.frame(viirs[as.numeric(dmrdat$sitenum),])
    #grab column number:
    cn <- as.matrix(as.numeric(dmrdat$colnum))
    #grab exploratory response  of choice:
    y <- as.matrix(as.numeric(dmrdat[,dmr]))
    x <- as.data.frame(cbind(y,cn,vvariables))
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

###disturbance probability analyses:
library(mgcv)
library(pROC)

viirs_ROC <- function(viirs,dmrdat,yr,coln){
  #make empty matrix:
  rocs <- matrix(NA,nrow=ncol(viirs),ncol=1)  #hard coded rn
  
  #grab just disturbance probability columns:
  dists<-dmrdat[,grep("^dp",colnames(dmrdat))]
  
  #first get just SMAP data:
  vvar <- as.data.frame(viirs[as.numeric(dmrdat$sitenum),])
  #grab column number:
  cn <- as.matrix(as.numeric(dmrdat$colnum))
  #grab exploratory response  of choice:
  y <- as.matrix(as.numeric(dists[,yr]))
  x <- as.data.frame(cbind(y,cn,vvar))
  vardat <- x[x$cn==coln,]
  
  # #remove missing values:
  # miss <- which(is.na(vardat[,3]))
  # vardat <- vardat[-miss,]
  
  #loop over all members of dmvars:
  for (i in 1:ncol(viirs)){  #this is hard coded until further notice
    #loop for filling in R2 table:  
    var.gam<-gam(vardat[,1]~s(vardat[,i+2]), data=vardat, family="binomial")
    var.roc<-roc(vardat[,1],var.gam$fitted.values)
    rocs[i,] <-var.roc$auc
  }
  #return table of AUCs
  return(rocs)
}

viirsrocs <- viirs_ROC(yr_rads,dmrdat,1,22)

demrocs <- viirs_ROC(site_data,dmrdat,1,22)

###BEFORE RUNNING: make sure coords match between mags and dem data
testing_viirs <- VIIRS_explore(yr_rads,testfx,"mags",22)


### Attempt to use opendapr package----------------------------------

#install packages:
# library(devtools)
# devtools::install_github("ptaconet/opendapr")
## Online it says "still in development" -- wasn't able to get it to install properly


