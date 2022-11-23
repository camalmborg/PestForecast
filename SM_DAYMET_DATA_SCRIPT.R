#This code is for downloading Daymet data for forecasts. Process:
#1) Download Daymet data for each site
#2) Extract variables
#3) Get monthly values
#4) 


#### Load Libraries:
library(daymetr)
library(tidyr)
library(stringr)


#### Load condition score .csv from GEE extract:
file<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"

#condition score object:
cond.scores<-read.csv(file)

#geographic coordinates from GEE extract:
geo<-as.data.frame(cond.scores[,".geo"])

#make lat and lon columns from .geo data:
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

