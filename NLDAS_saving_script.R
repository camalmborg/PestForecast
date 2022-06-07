#load packages:
library(ncdf4)
library(tidyverse)

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



