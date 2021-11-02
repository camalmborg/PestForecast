### Daymet Data Grab Script
#load data
#load annual condition scores:
#cond.scores.an<-read.csv("2020_07_10_sample_score_mean_ANNUAL.csv")
coords<-cond.scores.an[,30:31]
nsites<-1:5000
sites<-cbind(nsites,coords)

#put daymet data into a list
dm <- list()
for(i in nsites){
  dm[[i]] <- daymetr::download_daymet(site = sites$nsites[i],
                                      lat = sites$lat[i],
                                      lon = sites$lon[i],
                                      start = 1995,
                                      end = 2020,
                                      internal = TRUE)
}

#can use any dm# to get day of year
doy <- dm[[1]]$data$yday

#grabbing max temp from each site's daymet data in dm list
maxtemp<-list()
for (i in nsites){
   metyr= dm[[i]]$data$year
   metyears=unique(metyr)
   maxtemp[[i]]<-matrix(NA,length(metyears),366)
#   #maxtemp[,i]<-dm[[i]]$data$tmax..deg.c.
#   for(j in 1:nrow(dm[[i]]$data)){
#     maxtemp[[i]][as.numeric(as.factor(metyr))[j],dm[[i]]$data$yday[j]]=dm[[i]]$data$tmax..deg.c.[j]
#   }
}

#add the doy to be the first column, each of the next columns 2-9 are each site (in alphabetical order)
#maxtemp<-cbind(doy,maxtemp)
