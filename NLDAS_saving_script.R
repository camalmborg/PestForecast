#load packages:
library(ncdf4)
library(tidyverse)

#fixed url for new NLDAS location:
#url = "https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_NOAH0125_M.002"
url = "https://hydro1.gesdisc.eosdis.nasa.gov/dods/NLDAS_NOAH0125_H.002"
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
#2015-2022:
t <- time[434:518]

#test lat long:
lat <- lats[1]
lon <- longs[1]

#scrape data for soil moisture variables:
soilm <- ncvar_get(nc,"soilm0_200cm",
                      start=c(min(longs),min(lats),434),
                      count=c(length(longs),length(lats),length(t)))


#example lat long:
point<-cond.scores.mo[1,]
plat<-point$lat
plong<-point$lon
pt<-c(plat,plong)

gridpoints = data.frame(y = rep(seq_along(lats),times=length(longs)),
                        x=rep(seq_along(longs),each=length(lats))) %>%
  mutate(lat = nclat[lats[y]],lon=nclon[lats[x]])

dist = (pt[1]-gridpoints$lat)^2+(pt[2]-gridpoints$lon)^2
row = which.min(dist)
sm = soilm[gridpoints$x[row],gridpoints$y[row],]

plot(t,sm)
# for (l in lats){
#   for (g in longs){
#     
#   }
# }

#close nc file:
nc_close(nc)
