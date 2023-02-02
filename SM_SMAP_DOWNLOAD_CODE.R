#SMAP download code

####load libraries:
library(tidyr)
library(stringr)
library(leafletR) #cannot download on this version of R???
library(rgdal)

####load data:
#Load condition score .csv from GEE extract:
file<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"
#condition score object:
cond.scores<-read.csv(file)


spongy_SMAP <- function(start, end,
                        site_info,
                        geoJSON_outdir, 
                        smap_outdir){
  
  geo<-as.data.frame(site_info[,".geo"])
  
  #make lat and lon columns from .geo data:
  coords<-matrix(nrow=nrow(geo),ncol=3)
  for (i in 1:nrow(geo)){
    #longitudes:
    lon<-str_extract(geo[i,], "\\d+\\.*\\d*")
    coords[i,1]<-as.numeric(lon)*-1
    
    #latitudes:
    extlon<-sub(lon,"",geo[i,])
    coords[i,2]<-as.numeric(str_extract(extlon, "\\d+\\.*\\d*"))
    coords[i,3]<-i
  }
  colnames(coords)<-c("lon","lat","site_id")
  
  # Create geoJSON file for site
  site_GeoJSON <- data.frame(coords$lon, coords$lat) %>%
    setNames(c("lon","lat")) %>% 
    leafletR::toGeoJSON(name = site_info$name, dest = geoJSON_outdir, overwrite = TRUE) %>%
    rgdal::readOGR()
  site_GeoJSON$name = site_info$name
  site_GeoJSON = site_GeoJSON[-1] %>%
    leafletR::toGeoJSON(name = site_info$name, dest = geoJSON_outdir, overwrite = TRUE)
}