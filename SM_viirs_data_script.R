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
#making annual averages but skipping 0's:
viirs_rads[viirs_rads == 0] <- NA
names(viirs_rads) <- gsub("[^.-0-9]", "", names(viirs_rads), fixed = TRUE)

#selecting years and averaging:
vseq <- seq(1, ncol(viirs_rads), by=12)
yr_rads <- matrix(nrow=nrow(viirs_rads), ncol=length(seq))
for (i in 1:length(seq)){
  yr_rads[,i] <- apply(viirs_rads[,(vseq[i]):(vseq[i]+11)], 1, mean)
}


### Attempt to use opendapr package----------------------------------

#install packages:
# library(devtools)
# devtools::install_github("ptaconet/opendapr")
## Online it says "still in development" -- wasn't able to get it to install properly


