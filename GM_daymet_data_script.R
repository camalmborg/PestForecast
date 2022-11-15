### Daymet Data Grab Script
library(daymetr)
library(tidyr)
library(stringr)

###load data
#load annual condition scores:
file<-"2020_07_10_sample_score_mean_MONTHLY.csv"
cond.scores.mo<-read.csv(file)
coordcols<-c("lat","lon")
coords<-cond.scores.mo[,coordcols]
nsites<-1:nrow(cond.scores.mo)
sites<-cbind(nsites,coords)

#if you have data with just geo col:
file<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"
cond.scores.an<-read.csv(file)
geo<-as.data.frame(cond.scores.an[,".geo"])
#extgeo<-gsub('[typePointcoordinates\\\\"{}::()]',"",geo)
#newgeo.w<-gsub("([0-9]+)_.*", "\\1", newgeo)

#coordinates:
coords<-matrix(nrow=nrow(geo),ncol=2)
for (i in 1:nrow(geo)){
  #longitudes:
  lon<-str_extract(geo[i,], "\\d+\\.*\\d*")
  coords[i,1]<-as.numeric(lon)*-1
  
  #latitudes:
  extlon<-sub(lon,"",geo[i,])
  coords[i,2]<-as.numeric(str_extract(extlon, "\\d+\\.*\\d*"))
}



#put daymet data into a list
#for single site to get dm:
dm<-daymetr::download_daymet(site = sites$nsites[1],
                             lat = sites$lat[1],
                             lon = sites$lon[1],
                             start = 1995,
                             end = 2020,
                             internal = TRUE)

#loop version:
dm <- list()
for(i in nsites){
  dm[[i]] <- daymetr::download_daymet(site = sites$nsites[i],
                                      lat = sites$lat[i],
                                      lon = sites$lon[i],
                                      start = 1995,
                                      end = 2020,
                                      internal = TRUE)
}

#use any dm# to get day of year
doy <- dm[[1]]$data$yday

#grabbing max temp from each site's daymet data in dm list
maxtemp<-list()
for (i in nsites){
   metyr= dm[[i]]$data$year
   metyears=unique(metyr)
   maxtemp[[i]]<-matrix(NA,length(metyears),365)
   for(j in 1:nrow(dm[[i]]$data)){
     maxtemp[[i]][as.numeric(as.factor(metyr))[j],dm[[i]]$data$yday[j]]=dm[[i]]$data$tmax..deg.c.[j]
   }
print("done!", i)
}
#this code works^^^
print("it worked!")
#saveRDS(maxtemp, file='maxtemp.rds')
save(maxtemp, file='maxtemp.RData')
#load('maxtemp.RData')


#grabbing min temp from each site's daymet data in dm list
mintemp<-list()
for (i in nsites){
   metyr= dm[[i]]$data$year
   metyears=unique(metyr)
   mintemp[[i]]<-matrix(NA,length(metyears),365)
   for(j in 1:nrow(dm[[i]]$data)){
      mintemp[[i]][as.numeric(as.factor(metyr))[j],dm[[i]]$data$yday[j]]=dm[[i]]$data$tmin..deg.c.[j]
   }
   print("done!")
}
#this code works^^^
print("it worked!")
#saveRDS(mintemp, file='mintemp.rds')
save(mintemp, file='mintemp.RData')
#load('mintemp.RData')


#grabbing vpd from each site's daymet data in dm list
vpd<-list()
for (i in nsites){
   metyr= dm[[i]]$data$year
   metyears=unique(metyr)
   vpd[[i]]<-matrix(NA,length(metyears),365)
   for(j in 1:nrow(dm[[i]]$data)){
      vpd[[i]][as.numeric(as.factor(metyr))[j],dm[[i]]$data$yday[j]]=dm[[i]]$data$vp..Pa.[j]
   }
   print("done!")
}
#this code works^^^
print("it worked!")
#saveRDS(vpd, file='vpd.rds')
save(vpd, file='vpd.RData')
#load('mintemp.RData')


precip<-list()
for (i in nsites){
  metyr=dm[[i]]$data$year
  metyears=unique(metyr)
  precip[[i]]<-matrix(NA,length(metyears),365)
  for(j in 1:nrow(dm[[i]]$data)){
    precip[[i]][as.numeric(as.factor(metyr))[j],dm[[i]]$data$yday[j]]=dm[[i]]$data$prcp..mm.day.[j]
  }
  print("done!")
}
#this code works^^^
print("it worked!")
save(precip, file='precip.RData')


solr<-list()
for (i in nsites){
  metyr=dm[[i]]$data$year
  metyears=unique(metyr)
  solr[[i]]<-matrix(NA,length(metyears),365)
  for(j in 1:nrow(dm[[i]]$data)){
    solr[[i]][as.numeric(as.factor(metyr))[j],dm[[i]]$data$yday[j]]=dm[[i]]$data$srad..W.m.2.[j]
  }
  print("done!")
}
#this code works^^^
print("it worked!")
save(solr, file='solr.RData')

#automating to run through x number of variables?
#change to May-Sep grab instead of full year? 
#   -no, we may want winter temps
#   -separate seasons as needed later

