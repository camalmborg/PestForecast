#HF Field Data Script

###ARCHIVE ###load GEE data:---------
#tcg data:
hfplotstcg<-read.csv("HF_2022_Field_Data/HFplots_tcg_mean.csv")
#condition scores:
hfplotscond<-read.csv("HF_2022_Field_Data/HFplots_score_mean.csv")

#remove extra columns:
hftcg<-hfplotstcg[,2:131]
hfcond<-hfplotscond[,2:131]

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


###recovery rate slopes for field sites:

### HARVARD FOREST SUMMER FIELD DATA 2022 --------------------------

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
  mutate(plot = str_replace(plot, " ", "-"))
hf_condition <- hf_trees %>% 
  group_by(plot, Cond) %>% 
  summarize(count=n())
hf_mort <- hf_condition[hf_condition$Cond=="D",]
hf_mort$mort<-1

###merge with plot data:
field_plots <- merge(field_plots,hf_mort, all = TRUE)
field_plots$mort[is.na(field_plots$mort)] <- 0

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