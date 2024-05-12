#let's try to make a map?

### load libraries
library(tidyverse)
library(dplyr)
library(usmap)
library(terra)
library(ggplot2)


### plot maps
#plot_usmap(include = c("MA", "CT", "RI"))

#plot_usmap(include = c("MA"))


### making data for QGIS maps
# load the dp/dm/recov covariate data
# disturbance year data
dmr_file <- "CHAPTER_1/DATA/2023_12_DMR_DATA_TCG.csv"
dmr <- read.csv(dmr_file)
dmr_dm <- read.csv("CHAPTER_1/2024_JAGS_models/2024_04_magvars.csv")
dmr_dp <- read.csv("CHAPTER_1/2024_JAGS_models/2024_04_probvars.csv")
# load condition scores for coordinates
cfile <- "2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv" #1995-2020, original grab
scores <- read.csv(cfile)[dmr$X,]
geo <- as.data.frame(scores$.geo)

#make lat and lon columns from .geo:
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

# site data
sites <- cbind(id = 1:nrow(scores),
               sys = scores$system.index,
               lon = coords[,1],
               lat = coords[,2])

# make site-based covariate datasets
sites_dmr_mags <- cbind.data.frame(sites, dmr_dm)
sites_dmr_probs <- cbind.data.frame(sites, dmr_dp)
# save them
#write.csv(sites_dmr_mags, file = "Maps/Chapter_1/Data/sites_dmr_mags.csv")
#write.csv(sites_dmr_probs, file = "Maps/Chapter_1/Data/sites_dmr_probs.csv")
write.csv(sites_dmr_mags, file = "Maps/Chapter_1/Data/sites_mv_mags_centered.csv")
write.csv(sites_dmr_probs, file = "Maps/Chapter_1/Data/sites_mv_probs_centered.csv")
