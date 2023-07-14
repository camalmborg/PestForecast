#HF Field Data Script

###ARCHIVE ###load GEE data:---------
# #tcg data:
# hfplotstcg<-read.csv("HF_2022_Field_Data/HFplots_tcg_mean.csv")
# #condition scores:
# hfplotscond<-read.csv("HF_2022_Field_Data/HFplots_score_mean.csv")
# 
# #remove extra columns:
# hftcg<-hfplotstcg[,2:131]
# hfcond<-hfplotscond[,2:131]

# 
# ###load field data:
# HF.latlon<-read.csv("HF_2022_Field_Data/2022_HF_plots_latlon.csv")
# HF.plot.data<-read.csv("HF_2022_Field_Data/HF_Plot_data.csv")
# HF.seedlings<-read.csv("HF_2022_Field_Data/HF_Seedlings_long.csv")
# HF.trees<-read.csv("HF_2022_Field_Data/HF_Tree_data.csv")
# HF.understory<-read.csv("HF_2022_Field_Data/HF_Und_ground_survey.csv")
# 
# #separating mortality sites:
# library("dplyr")
# HF.condition<-HF.trees %>% group_by(plot_2, Cond) %>% summarize(count=n())
# HF.mort<-HF.condition[HF.condition$Cond=="D",]
# HF.mort$mort<-1
# 
# #merge with plot data:
# HF.plot.data<- merge(HF.plot.data, HF.mort, all = TRUE)
# HF.plot.data$mort[is.na(HF.plot.data$mort)] <- 0

#### HARVARD FOREST DATA:--------------------------------------------------
file <- "HF_2022_Field_Data/GEE_Data/2023_05_17_hfplots_sample_tcg_mean.csv"
cfile <- "HF_2022_Field_Data/GEE_Data/2023_05_17_hfplots_sample_score_mean.csv"

#tcg and condition score objects:
tcg.values<-read.csv(file)
#tcgs<-tcg.values[,c(grep("^X",colnames(tcg.values)))] #grab with "X1996 eg) for all tcg value columns
cond.scores<-read.csv(cfile)

#separate lat lons:
#make coords object:
geo<-as.data.frame(tcg.values[,".geo"])
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


#load understory data:
hf_unds <- "HF_2022_Field_Data/Und_ground_survey.csv"

#load seedling data:
hf_seedlings <- "HF_2022_Field_Data/Seedlings_long.csv"

### HARVARD FOREST SUMMER FIELD DATA 2022 CLEANING --------------------------

#load libraries:
library(dplyr)
library(tidyverse)

###load and fix data-----
#field plot data:
field_plots <- read.csv("HF_2022_Field_Data/Plot_data  - Sheet1.csv")
field_plots <- field_plots %>%
  mutate(plot = str_replace(plot, " ", "-")) %>%
  mutate(latitude = str_replace(latitude, " N", "")) %>%
  mutate(longitude = str_replace(longitude, " W", "")) %>%
  mutate_at(c('latitude','longitude'), as.numeric) %>%
  mutate(invasives = ifelse(invasives == "no",0,1)) %>%
  mutate(recent_timber_harvest = ifelse(recent_timber_harvest == "no",0,1))
field_plots$longitude <- field_plots$longitude*-1
  #mutate(invasives = str_replace(invasives, "yes", "1")) ###need to do the 0/1 replace this this and timber harvest

#making a new lat/lon file:
# hf_lat_lon <- cbind(field_plots$longitude, field_plots$latitude)
# colnames(hf_lat_lon) <- c('longitude','latitude')
# write.csv(hf_lat_lon, file="hf_lat_lon.csv", row.names=TRUE)


#individual tree data (to get plots with mortality observed):
hf_trees <- read.csv("HF_2022_Field_data/Tree_data - Sheet1.csv")
hf_trees <- hf_trees %>% 
  mutate(plot = str_replace(plot, " ", "-")) %>%
  mutate(CondBin = ifelse(Cond == "D",1,0))
hf_condition <- hf_trees %>% 
  group_by(plot, Cond) %>% 
  summarize(count=n())
hf_mort <- hf_condition[hf_condition$Cond=="D",]
hf_mort$mort<-1

#plots with oaks present:
oaks <- c("BO", "RO", "WO")
hf_spec <- hf_trees %>%
  group_by(plot, spp) %>%
  summarize(count=n())

#count number of trees, species, and oaks in each plot:
plots <- unique(hf_spec$plot)
plot <- vector()
n_trees <- vector()
n_dead <- vector()
n_spec <- vector()
n_oaks <- vector()
n_d_oaks <- vector()
for (i in plots){
  plot[i] <- i
  hfs <- hf_spec[hf_spec$plot==i,]
  n_trees[i] <- sum(hfs$count)
  n_spec[i] <- nrow(hfs)
  oak <- rbind(hfs[hfs$spp == oaks[1],],
               hfs[hfs$spp == oaks[2],],
               hfs[hfs$spp == oaks[3],])
  n_oaks[i] <- sum(oak$count)
  
  trees <- hf_trees[hf_trees$plot==i,]
  n_dead[i] <- sum(trees$CondBin)
  oak_d <- rbind(trees[trees$spp == oaks[1],],
                 trees[trees$spp == oaks[2],],
                 trees[trees$spp == oaks[3],])
  n_d_oaks[i] <- sum(oak_d$CondBin)
  rm(hfs,oak, trees, oak_d)
}
tree_data <- cbind(plot, n_trees, n_spec, n_oaks, n_dead, n_d_oaks)

###merge with plot data:
field_plots <- merge(field_plots,hf_mort, all = TRUE)
field_plots$mort[is.na(field_plots$mort)] <- 0

field_data <- merge(field_plots, tree_data, all=TRUE)
field_data <- field_data[-which(is.na(field_data$n_trees)),c(1:11,14:19)]

###oak plots:
# hf_oaks <- as.data.frame(matrix(nrow=nrow(field_plots), ncol=length(oaks)+1))
# hf_oaks[,1]<- field_plots$plot


# for (i in 1:length(oaks)){
#   for (j in 1:row(hf_spec)){
#     hf_oaks[j,i+1] <- ifelse(unique(hf_spec$plot)[j] == hf_oaks$plot[j] & hf_spec$spp[j,] == oaks[i],
#                             1,0)
#   }
# }

#merge with remote-sensing data:
#load:
hf_mags <- read.csv("HF_2022_Field_Data/HF_mags_recov_from_GEE_data.csv")
#join:
hf_data <- cbind(field_plots,hf_mags)

### GAM time:
library(mgcv)
library(pROC)

hf_gam <- gam(hf_data$mort ~ s(hf_data$mags), 
              data=hf_data,
              family = "binomial")
hf_roc<-roc(hf_data$mags,hf_gam$fitted.values)

#summ <- 

### Making some plots:
# library(lattice)
# library(latticeExtra)
# library(tactile)

plot.roc(hf_data$mort,hf_data$mags,
         percent=T)


plot(hf_data$mags, hf_data$mort)

#bin means:
bins=seq(min(hf_data$mags)-1,max(hf_data$mags)+1,length=10)
bin = findInterval(hf_data$mags,bins)
mu = tapply(hf_data$mort,bin,mean,na.rm=TRUE)
n  = tapply(hf_data$mort,bin,length)
sigma = sqrt(mu*(1-mu)/n)
x = min(hf_data$mags)-1 + cumsum(diff(bins))
points(x,mu,col=2)
points(x,mu+2*sigma,col=3,pch="-")
points(x,mu-2*sigma,col=3,pch="-")

#next stop: plot the gam
fit <- hf_gam$fitted.values
points(hf_data$mags, fit, col=4, pch="*")
