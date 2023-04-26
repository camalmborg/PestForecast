#This code is for downloading Daymet data for forecasts. Process:
#1) Download Daymet data for each site
#2) Extract variables
#3) Get monthly values (rows) for each year (column) within list of sites


#### Load Libraries:
#library(daymetr)
library(tidyr)
library(stringr)
library(tidyverse)
#library(miceadds)


#### Load condition score .csv from GEE extract:
#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"
file <- "2023_03_08_DATAGRAB/2023_03_08_5000_sites_sample_score_mean.csv"

#condition score object:
cond.scores<-read.csv(file)
#cond.samp<-cond.scores[1:3,]


#### Function for grabbing daymet data:-----
spongy_met<-function(scores,startyr,endyr,var,filenm){
  ##Section for getting sites from GEE dataset:
  #geographic coordinates from GEE extract:
  geo<-as.data.frame(scores[,".geo"])
  
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
  
  #make dataset with condition scores and coordinates:
  nsites<-1:nrow(scores)
  sites<-as.data.frame(cbind(nsites,coords))
  ######these lines were for testing: ###
  #sites<-sites[1:5,]
  #nsites<-1:nrow(sites)
  
  ##Section for downloading daymet for each site:
  dm <- list()
  for(i in nsites){
    dm[[i]] <- daymetr::download_daymet(site = sites$nsites[i],
                                        lat = sites$lat[i],
                                        lon = sites$lon[i],
                                        start = startyr,
                                        end = endyr,
                                        internal = TRUE)
  }
  
  #use any dm# to get day of year
  doy <- dm[[1]]$data$yday
  #get all years and unique years for later:
  for (i in nsites){
    metyr=dm[[i]]$data$year
    metyears=unique(metyr)
  }
  
  #make empty matrix for all daymet variable data:
  #dm.v<-matrix(data=NA, nrow=nrow(sites))
  dm.v<-list()
  
  ##Loop for extracting daymet variables of interest (var):
  for (k in 1:length(var)){
    dmvar<-list()
    for (i in nsites){
      metyr<-dm[[i]]$data$year
      metyears<-unique(metyr)
      dmvar[[i]]<-matrix(NA,length(metyears),365)
      for(j in 1:nrow(dm[[i]]$data)){
        dmvar[[i]][as.numeric(as.factor(metyr))[j], 
                  dm[[i]]$data$yday[j]]=dm[[i]]$data[[var[k]]][j]
      }
    }
  #save data:
  #save(dmvar,file=paste0(filenm[k],".RData"))
  
  ##Section for computing monthly average variable values:
  meanvar<-list()
  #monthly breaks:
  jan<-c(doy[1:31])
  feb<-c(doy[32:59])
  mar<-c(doy[60:90])
  apr<-c(doy[91:120])
  may<-c(doy[121:151])
  jun<-c(doy[152:181])
  jul<-c(doy[182:212])
  aug<-c(doy[213:243])
  sep<-c(doy[244:273])
  oct<-c(doy[274:304])
  nov<-c(doy[305:334])
  dec<-c(doy[335:365])

  ##Loop for getting mean values per month
  for (i in nsites){
    meanvar[[i]]<-matrix(NA, nrow=12, ncol=length(metyears))
    for (j in 1:length(metyears)){
      for (s in 1:12){
        if (s==1){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,jan])
        }
        if (s==2){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,feb])
        }
        if (s==3){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,mar])
        }
        if (s==4){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,apr])
        }
        if (s==5){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,may])
        }
        if (s==6){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,jun])
        }
        if (s==7){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,jul])
        }
        if (s==8){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,aug])
        }
        if (s==9){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,sep])
        }
        if (s==10){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,oct])
        }
        if (s==11){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,nov])
        }
        if (s==12){
          meanvar[[i]][s,j]<-mean(dmvar[[i]][j,dec])
        }
      }
    }
  }
  #save meanvar data:
  #save(meanvar, file=paste0(filenm[k],"_monthly_means",".Rdata"))
  
  ### Getting mean values into 1 matrix:
  #number of years:
  alltime=c(1,1:(endyr-startyr)+1)
  #months:
  mons=1:12
  ## Loop for extracting monthly values
  x.p <- matrix(data=NA, nrow=nrow(sites))
  for (p in alltime){
    x <- matrix(data = NA, nrow=length(nsites), ncol=length(mons))
    for (m in 1:length(mons)){
      for (s in nsites){
        #months become columns, rows become sites:
        x[s,m] <- meanvar[[s]][mons[m],alltime[p]]
      }
    }
    x.p <- cbind(x.p,x)
    rm(x) #remove last loop
  }
  #remove x.p NA column, remove x.p extra variable
  envar <- x.p[,2:((length(mons)*length(alltime))+1)] #remove NA column
  rm(x.p)
  
  #add to daymet variable:
  dm.v[[filenm[k]]]<-envar
  rm(envar)
  }

  #return(envar)
  return(dm.v)
}

#(as characters):
#testing, testing, is this thing on?
dmvars<-spongy_met(cond.scores,2011,2017,
                   c("tmax..deg.c.","tmin..deg.c.","prcp..mm.day.","vp..Pa."),
                   c("maxtemp","mintemp","pcp","vpd"))

#### Getting seasonal variable values for analyses:-----
## Quick function for making sequences to extract 
## monthly values based on # of years of daymet data downloaded:
dm_me <- function(x,yrs,month){
  mcols <- vector()
  for (i in 1:length(month)){
    colpyr <- c(seq(month[i],yrs*12,by=12))
    mcols <- c(mcols,colpyr)
  }
  cols <- sort(mcols) #x = month (1=jan, 2=feb, etc.)
  dmvar <- x[,cols]
  return(dmvar)
}

# spongyvars<-cbind(dm_me(dmvars$maxtemp,2,c(1,6)),
#                   dm_me(dmvars$pcp,2,c(1,6)))

### Expanding DM_ME to make seasonal daymet averages for spring, summer, fall, winter months
#spring = march/april/may (3,4,5), summer = june/july/august (6/7/8)

#if necessary:
#library(tidyverse)
#library(dplyr)

dm_me_seasonal <- function(x,yrs,month){
  mcols <- vector()
  m_id <- vector()
  y_id <- vector()
  n_id <- c(1:(length(month)*yrs))
  for (i in 1:length(month)){
    colpyr <- c(seq(month[i],yrs*12,by=12))
    mcols <- c(mcols,colpyr)
    m_id <- c(m_id, rep(month[i], times=yrs))
    y_id <- c(y_id, rep(1:yrs))
  }
  
  #makes an id daata frame with which to attach months and years to 
  #the resulting dmvar data frame in long form:
  dm_mycols <- cbind(mcols, n_id, m_id, y_id)
  
  #grab the dmvar data:
  cols <- sort(dm_mycols[,"mcols"]) #x = month (1=jan, 2=feb, etc.)
  dmvarbs <- as.matrix(x[,cols])
  
  #make final thing
  dmcols <- t(dm_mycols)
  dmcols <- dmcols[,order(dmcols["mcols",])]
  dmv <- rbind(dmcols,dmvarbs)
  
  #average the colums:
  dmvar <- matrix(NA, nrow=nrow(dmv)-4, ncol=7)
  for (i in 1:yrs){
    dmvar[,i] <- apply(dmv[5:nrow(dmv),dmv["y_id",]==i], 1, mean, na.rm=T) 
  }
  
  return(dmvar)
}

#seasons:
spring <- c(3,4,5)
summer <- c(6,7,8)
fall <- c(9,10,11)
winter <- c(12,1,2)

spongyvars <- list(dm_me_seasonal(dmvars$maxtemp, 7, spring),
                   dm_me_seasonal(dmvars$mintemp, 7, spring),
                   dm_me_seasonal(dmvars$pcp, 7, spring),
                   dm_me_seasonal(dmvars$vpd, 7, spring))


# dm_me is used to create a data frame that includes daymet variables
# of interest - arguments are:
#1) dmvars$... <- the variable you want to extract
#2) the number of years in your analysis (typically the previous 2 years)
#3) c(....) <- the months of the year you want the values (eg. 1=jan, 2=feb)
# you cbind these and it makes a dataframe with these values where columns
# are the monthly daymet variables and rows are sites

#dm_me_seasonal is used to create a matrix (or data frame) that includes
#daymet variables of interest SEASONALLY from mean monthly values
#arguments are:
#1) x <- enter dmvars from spongymet function
#2) yrs <- number of years of the analyses (2011-2017 = 7)
#3) months <- months you want to use in analyses: 3,4,5 = Mar, Apr, May
#spring = 3,4,5 #summer = 6,7,8 eg (calendar months)

#### Getting variables for the forecast:


###WHEN RUNNING FOR REAL BE SURE TO REMOVE 5 SITE SUBSET^^^-----

#choose your variable from the daymet list, add filename
#variables:
#c("tmax..deg.c.","tmin..deg.c.","prcp..mm.day.","vp..Pa.")
#filenames:
#c("maxtemp","mintemp","pcp","vpd")

#load"maxtemp_monthly_means.RData")
#check if it works....
#maxtemp <- miceadds::load.Rdata2("maxtemp_monthly_means.RData")
#mintemp <- miceadds::load.Rdata2("mintemp_monthly_means.RData")

# it works! Charlotte is a beautiful genius!



#----------##### THIS VERSION IS FOR 15-to-15 MONTHLY DATA: #-----------######
### Function for grabbing daymet data:---
spongy_met_15s<-function(scores,startyr,endyr,var,filenm){
  ##Section for getting sites from GEE dataset:
  #geographic coordinates from GEE extract:
  geo<-as.data.frame(scores[,".geo"])
  
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
  
  #make dataset with condition scores and coordinates:
  nsites<-1:nrow(scores)
  sites<-as.data.frame(cbind(nsites,coords))
  ######these lines were for testing: ###
  #sites<-sites[1:5,]
  #nsites<-1:nrow(sites)
  
  ##Section for downloading daymet for each site:
  dm <- list()
  for(i in nsites){
    dm[[i]] <- daymetr::download_daymet(site = sites$nsites[i],
                                        lat = sites$lat[i],
                                        lon = sites$lon[i],
                                        start = startyr,
                                        end = endyr,
                                        internal = TRUE)
  }
  
  #use any dm# to get day of year
  doy <- dm[[1]]$data$yday
  #get all years and unique years for later:
  for (i in nsites){
    metyr=dm[[i]]$data$year
    metyears=unique(metyr)
  }
  
  #make empty matrix for all daymet variable data:
  #dm.v<-matrix(data=NA, nrow=nrow(sites))
  dm.v<-list()
  
  ##Loop for extracting daymet variables of interest (var):
  for (k in 1:length(var)){
    dmvar<-list()
    for (i in nsites){
      metyr<-dm[[i]]$data$year
      metyears<-unique(metyr)
      dmvar[[i]]<-matrix(NA,length(metyears),365)
      for(j in 1:nrow(dm[[i]]$data)){
        dmvar[[i]][as.numeric(as.factor(metyr))[j], 
                   dm[[i]]$data$yday[j]]=dm[[i]]$data[[var[k]]][j]
      }
    }
    #save data:
    #save(dmvar,file=paste0(filenm[k],".RData"))
    
    ##Section for computing 1/15 monthly average variable values:
    meanvar<-list()
    #monthly breaks:
    jan<-c(doy[c(1:14, 349:365)])
    feb<-c(doy[15:45])
    mar<-c(doy[46:73])
    apr<-c(doy[74:104])
    may<-c(doy[105:134])
    jun<-c(doy[135:165])
    jul<-c(doy[166:195])
    aug<-c(doy[196:226])
    sep<-c(doy[227:257])
    oct<-c(doy[258:287])
    nov<-c(doy[288:318])
    dec<-c(doy[319:348])
    
    ##Loop for getting mean values per month
    for (i in nsites){
      meanvar[[i]]<-matrix(NA, nrow=12, ncol=length(metyears))
      for (j in 1:length(metyears)){
        for (s in 1:12){
          if (s==1){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,jan])
          }
          if (s==2){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,feb])
          }
          if (s==3){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,mar])
          }
          if (s==4){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,apr])
          }
          if (s==5){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,may])
          }
          if (s==6){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,jun])
          }
          if (s==7){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,jul])
          }
          if (s==8){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,aug])
          }
          if (s==9){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,sep])
          }
          if (s==10){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,oct])
          }
          if (s==11){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,nov])
          }
          if (s==12){
            meanvar[[i]][s,j]<-mean(dmvar[[i]][j,dec])
          }
        }
      }
    }
    #save meanvar data:
    #save(meanvar, file=paste0(filenm[k],"_monthly_means",".Rdata"))
    
    ### Getting mean values into 1 matrix:
    #number of years:
    alltime=c(1,1:(endyr-startyr)+1)
    #months:
    mons=1:12
    ## Loop for extracting monthly values
    x.p <- matrix(data=NA, nrow=nrow(sites))
    for (p in alltime){
      x <- matrix(data = NA, nrow=length(nsites), ncol=length(mons))
      for (m in 1:length(mons)){
        for (s in nsites){
          #months become columns, rows become sites:
          x[s,m] <- meanvar[[s]][mons[m],alltime[p]]
        }
      }
      x.p <- cbind(x.p,x)
      rm(x) #remove last loop
    }
    #remove x.p NA column, remove x.p extra variable
    envar <- x.p[,2:((length(mons)*length(alltime))+1)] #remove NA column
    rm(x.p)
    
    #add to daymet variable:
    dm.v[[filenm[k]]]<-envar
    rm(envar)
  }
  
  #return(envar)
  return(dm.v)
}

#(as characters):
#testing, testing, is this thing on?
dmvars<-spongy_met_15s(cond.scores,2012,2017,
                   c("tmax..deg.c.","tmin..deg.c.","prcp..mm.day.","vp..Pa."),
                   c("maxtemp","mintemp","pcp","vpd"))

