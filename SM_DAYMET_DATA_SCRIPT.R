#This code is for downloading Daymet data for forecasts. Process:
#1) Download Daymet data for each site
#2) Extract variables
#3) Get monthly values (rows) for each year (column) within list of sites
#4) 


#### Load Libraries:
library(daymetr)
library(tidyr)
library(stringr)
library(miceadds)


#### Load condition score .csv from GEE extract:
file<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"

#condition score object:
cond.scores<-read.csv(file)


#### Function for grabbing daymet data:
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
  sites<-sites[1:5,]
  nsites<-1:nrow(sites)
  
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
  save(meanvar, file=paste0(filenm[k],"_monthly_means",".Rdata"))
  
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
  dm.v[[k]]<-envar
  rm(envar)
  }
  
  #return(envar)
  return(dm.v)
}

#choose your variable from the daymet list, add filename 
#(as characters):
#testing, testing, is this thing on?
dmvars<-spongy_met(cond.scores,2020,2021,c("tmax..deg.c.","tmin..deg.c."),c("maxtemp","mintemp"))


#load"maxtemp_monthly_means.RData")
#check if it works....
#maxtemp <- miceadds::load.Rdata2("maxtemp_monthly_means.RData")
#mintemp <- miceadds::load.Rdata2("mintemp_monthly_means.RData")

# it works! Charlotte is a beautiful genius!