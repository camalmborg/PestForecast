#load packages:
library(ncdf4)
library(tidyverse)

#####GET SOIL MOISTURE DATA FROM NASA GESDISC-----
#fixed url for new NLDAS location:
url = "https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_NOAH0125_M.002"
#url = "https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_NOAH0125_H.002"
nc = nc_open(url)

#establish lat long bounding box
bounding<-read.csv("2022_03_29_latlonboundingbox.csv")

#nc lat and long:
nclat <- ncvar_get(nc,"lat")
nclon <- ncvar_get(nc,"lon")

lats <- which(nclat>min(bounding[,1]) & nclat<max(bounding[,1]))
longs <- which(nclon>min(bounding[,2]) & nclon<max(bounding[2,]))

#time:
time <- ncvar_get(nc,"time")
#for data that is Jan 1, 1979-April 1, 2022
#2015-2022:
t <- time[387:520] #match start below with first # here

#test lat long:
lat <- lats[1]
lon <- longs[1]

#scrape data for soil moisture variables:
soilm <- ncvar_get(nc,"soilm0_200cm",
                      start=c(min(longs),min(lats),387),
                      count=c(length(longs),length(lats),length(t)))


# #example lat long:
# point<-cond.scores.mo[30,]
# plat<-point$lat
# plong<-point$lon
# pt<-c(plat,plong)
# 
# gridpoints = data.frame(y = rep(seq_along(lats),times=length(longs)),
#                         x=rep(seq_along(longs),each=length(lats))) %>%
#   mutate(lat = nclat[lats[y]],lon=nclon[lats[x]])
# 
# dist = (pt[1]-gridpoints$lat)^2+(pt[2]-gridpoints$lon)^2
# row = which.min(dist)
# sm = soilm[gridpoints$x[row],gridpoints$y[row],]
# 
# plot(t,sm)
load(file = "cond_scores_mo.csv")

###MAKE A MATRIX OF SOIL MOISTURE FOR EACH POINT 
soilm.sites<-matrix(data=NA, nrow=nrow(cond.scores.mo), ncol=length(t))
for (p in 1:nrow(cond.scores.mo)){
  point<-cond.scores.mo[p,] #select site
  plat<-point$lat           #grab lat and long
  plong<-point$lon
  pt<-c(plat,plong)
  
  #find which gridpoint that site is in:
  gridpoints = data.frame(y = rep(seq_along(lats),times=length(longs)),
                          x=rep(seq_along(longs),each=length(lats))) %>%
    mutate(lat = nclat[lats[y]],lon=nclon[lats[x]])
  
  dist = (pt[1]-gridpoints$lat)^2+(pt[2]-gridpoints$lon)^2
  row = which.min(dist)
  
  #grab sm data, put in matrix:
  soilm.sites[p,]<-soilm[gridpoints$x[row],gridpoints$y[row],]
  
}

#close nc file:
nc_close(nc)

#### Make dataset with soil moisture for hatching and feeding windows
#hatch = April/May avg
#feed = June/July avg

#grab years--predisturb window:
yr11<-2:5
yr12<-14:17
yr13<-26:29
yr14<-38:41
yr15<-50:53

#soilm for each:
sm11<-soilm.sites[,yr11]
sm12<-soilm.sites[,yr12]
sm13<-soilm.sites[,yr13]
sm14<-soilm.sites[,yr14]
sm15<-soilm.sites[,yr15]

#get hatch and feed avgs for each:
hfmeans<-function(x){
  hfmean<-cbind(apply(x[,1:2],1,mean),
          apply(x[,3:4],1,mean))
  return(hfmean)
}

#hatch and feed avgs and anomaly function:
sm11hf<-hfmeans(sm11)
sm12hf<-hfmeans(sm12)
sm13hf<-hfmeans(sm13)
sm14hf<-hfmeans(sm14)
sm15hf<-hfmeans(sm15)

# sm11hf<-anomfx(hfmeans(sm11))
# sm12hf<-anomfx(hfmeans(sm12))
# sm13hf<-anomfx(hfmeans(sm13))
# sm14hf<-anomfx(hfmeans(sm14))
# sm15hf<-anomfx(hfmeans(sm15))

smpredist<-cbind(sm11hf,sm12hf,sm13hf,sm14hf,sm15hf)
#smpredistanom<-cbind(sm11hf,sm12hf,sm13hf,sm14hf,sm15hf)

#save data:
write.csv(smpredist, "soil_moisture_data.csv")
#write.csv(smpredistanom,"soil_moisture_anom_data.csv")

