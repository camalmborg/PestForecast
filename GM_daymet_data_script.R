### Daymet Data Grab Script
library(daymetr)
library(tidyr)
library(stringr)

###for if your data has lat lon columns:----
###load data
#load condition scores:
# file<-"2020_07_10_sample_score_mean_MONTHLY.csv"
# cond.scores<-read.csv(file)

###if you have lat and lon columns:
# coordcols<-c("lat","lon")
# coords<-cond.scores.mo[,coordcols]
# nsites<-1:nrow(cond.scores.mo)
# sites<-cbind(nsites,coords)

###if you have data with just geo col:-----
file<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"
cond.scores<-read.csv(file)
geo<-as.data.frame(cond.scores[,".geo"])

# extract coordinates from .geo column:
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
nsites<-1:nrow(cond.scores)
sites<-as.data.frame(cbind(nsites,coords))
#these lines were for testing:
sites<-sites[1:5,]
nsites<-1:nrow(sites)


###loop version for putting daymet data into a list:
#start and end
startyr <-1995
endyr<-2000

#loop for downloading daymet data
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

#function for grabbing daymet data:
spongy_met<-function(var,filenm){
  dmvar<-list()
  for (i in nsites){
    metyr<-dm[[i]]$data$year
    metyears<-unique(metyr)
    dmvar[[i]]<-matrix(NA,length(metyears),365)
    for(j in 1:nrow(dm[[i]]$data)){
      dmvar[[i]][as.numeric(as.factor(metyr))[j], 
                 dm[[i]]$data$yday[j]]=dm[[i]]$data[[var]][j]
    }
  }
  #return(dmvar)
  save(dmvar,file=paste0(filenm,".RData"))
}

#choose your variable from the daymet list, add filename 
#(as characters):
spongy_met("tmax..deg.c.","maxtemp")

#filename can include full path to folder for data

#can do for variables:
#maxtemp=tmax..deg.c.
#mintemp=tmin..deg.c.
#vpd=vp..Pa.
#precip=prcp..mm.day.

#will load as object called dmvar:
load("maxtemp.RData")

#^^^this code works now!!!


###ARCHIVED CODE-----

#these lines were for testing:
#sites<-sites[1:5,]
#nsites<-1:nrow(sites)

###put daymet data into a list
#for single site to get dm:
# dm<-daymetr::download_daymet(site = sites$nsites[1],
#                              lat = sites$lat[1],
#                              lon = sites$lon[1],
#                              start = 1995,
#                              end = 2020,
#                              internal = TRUE)


###individual variable loops:
# maxtemp<-list()
# for (i in nsites){
#   metyr= dm[[i]]$data$year
#   metyears=unique(metyr)
#   maxtemp[[i]]<-matrix(NA,length(metyears),365)
#   for(j in 1:nrow(dm[[i]]$data)){
#     maxtemp[[i]][as.numeric(as.factor(metyr))[j],dm[[i]]$data$yday[j]]=dm[[i]]$data$tmax..deg.c.[j]
#   }
#   print("done!", i)
# }
# #this code works^^^
# print("it worked!")
# #saveRDS(maxtemp, file='maxtemp.rds')
# save(maxtemp, file='maxtemp.RData')
# #load('maxtemp.RData')