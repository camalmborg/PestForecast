#HF Field Data Script

#load libraries:
library(dplyr)
library(tidyverse)

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
hf_unds <- read.csv("HF_2022_Field_Data/Und_ground_survey.csv")
hf_unds <- hf_unds %>%
  mutate(plot = str_replace(plot, " ", "-"))
hf_ground <- hf_unds[hf_unds$type == "g",]
#hf_unds <- hf_unds[hf_unds$type == "u",]

#load seedling data:
hf_seedlings <- read.csv("HF_2022_Field_Data/Seedlings_long.csv")
hf_seedlings <- hf_seedlings %>%
  mutate(plot = str_replace(plot, " ", "-"))

### HARVARD FOREST SUMMER FIELD DATA 2022 CLEANING --------------------------

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

#individual tree data (to get plots with mortality observed):
hf_trees <- read.csv("HF_2022_Field_data/Tree_data - Sheet1.csv")
#calculate Basal Area (BA) from DBH for each tree:
hf_trees$BA <- hf_trees$dbh^2 * 0.005454
hf_trees <- hf_trees %>% 
  mutate(plot = str_replace(plot, " ", "-")) %>%
  mutate(CondBin = ifelse(Cond == "D",1,0))
hf_condition <- hf_trees %>% 
  group_by(plot, Cond) %>% 
  summarize(count=n())
hf_mort <- hf_condition[hf_condition$Cond=="D",]
hf_mort$mort<-1

#epicormic sprouts present/absent - regrowth/resprout evidence
hf_epsprouts <- hf_trees %>%
  mutate(Epi.S = ifelse(Epi.S == "Yes",1,0))
hf_epsprouts <- hf_epsprouts[which(hf_epsprouts$Epi.S==1),]
hf_sprouts <- hf_epsprouts %>%
  group_by(plot) %>%
  summarize(count=n())
hf_sprouts$sprouts<-1

#merge with plot data:
field_plots <- merge(field_plots,hf_mort, all = TRUE)
field_plots$mort[is.na(field_plots$mort)] <- 0
field_plots <- merge(field_plots,hf_sprouts, all.x = TRUE)
field_plots$sprouts[is.na(field_plots$sprouts)] <- 0

### Getting tree and oak counts in each plot:
#species counts in each plot:
oaks <- c("BO", "RO", "WO", "CO")
hf_spec <- hf_trees %>%
  group_by(plot, spp) %>%
  summarize(count=n())

#count number of trees, species, and oaks in each plot:
plots <- unique(hf_spec$plot)
plot <- vector()
n_trees <- vector()
plot_BA_from_BAF <- vector()
n_dead <- vector()
n_spec <- vector()
n_oaks <- vector()
n_d_oaks <- vector()
for (i in plots){
  plot[i] <- i
  hfs <- hf_spec[hf_spec$plot==i,]
  
  #get number of trees in each plot, number of species, and number of oaks:
  n_trees[i] <- sum(hfs$count)
  plot_BA_from_BAF[i] <- n_trees[i] * 10 #Basal area per acre >> BA = # trees in * Basal Area Factor (10) 
  n_spec[i] <- nrow(hfs)
  oak <- rbind(hfs[hfs$spp == oaks[1],],
               hfs[hfs$spp == oaks[2],],
               hfs[hfs$spp == oaks[3],],
               hfs[hfs$spp == oaks[4],])
  n_oaks[i] <- sum(oak$count)
  
  #get number of dead oaks:
  trees <- hf_trees[hf_trees$plot==i,]
  n_dead[i] <- sum(trees$CondBin)
  oak_d <- rbind(trees[trees$spp == oaks[1],],
                 trees[trees$spp == oaks[2],],
                 trees[trees$spp == oaks[3],],
                 trees[trees$spp == oaks[4],])
  n_d_oaks[i] <- sum(oak_d$CondBin)
  rm(hfs,oak, trees, oak_d)
}
tree_data <- cbind(plot, n_trees, plot_BA_from_BAF, n_spec, n_oaks, n_dead, n_d_oaks)


###basal area calculations (for % DBH):
plot <- vector()
plot_dbh <- vector()
plot_BA <- vector()
oak_dbh <- vector()
oak_BA <- vector()
dead_dbh <- vector()
dead_BA <- vector()
dead_oak_dbh <- vector()
dead_oak_BA <- vector()
for (i in plots){
  plot[i] <- i
  
  #plot dbh calculation:
  hfs <- hf_trees[hf_trees$plot==i,]
  plot_dbh[i] <- sum(hfs$dbh)
  plot_BA[i] <- sum(hfs$BA)
  oak <- rbind(hfs[hfs$spp == oaks[1],],
               hfs[hfs$spp == oaks[2],],
               hfs[hfs$spp == oaks[3],],
               hfs[hfs$spp == oaks[4],])
  oak_dbh[i] <- sum(oak$dbh)
  oak_BA[i] <- sum(oak$BA)
  
  #dbh and BA of dead trees in each plot:
  d_trees <- hfs[hfs$CondBin == 1,]
  dead_dbh[i] <- sum(d_trees$dbh)
  dead_BA[i] <- sum(d_trees$BA)
  
  #number of dead trees and dead oaks:
  n_dead[i] <- sum(hfs$CondBin)
  oak_d <- rbind(hfs[hfs$spp == oaks[1],],
                 hfs[hfs$spp == oaks[2],],
                 hfs[hfs$spp == oaks[3],],
                 hfs[hfs$spp == oaks[4],])
  n_d_oaks[i] <- sum(oak_d$CondBin)
  
  #dbh and BA of dead oaks:
  dead_oak <- oak_d[oak_d$CondBin == 1,]
  dead_oak_dbh[i] <- sum(dead_oak$dbh)
  dead_oak_BA[i] <- sum(dead_oak$BA)
  
  rm(hfs, oak, d_trees, oak_d, dead_oak)
}

tree_data <- cbind(tree_data, plot_dbh, plot_BA, dead_dbh, dead_BA, dead_oak_dbh, dead_oak_BA)


###cleaning seedling data:
plot <- vector()
n_seedlings <- data.frame()
total_seed <- vector()
for (i in 1:length(plots)){
  plot[i] <- plots[i]
  for (j in 1:3){
    hfs <- hf_seedlings[hf_seedlings$plot == plots[i],]
    hfs <- hfs[hfs$size == j,]
    n_seedlings[i,j] <- sum(hfs[,4:7])
  }
  total_seed[i] <- sum(n_seedlings[i,])
  rm(hfs)
}
seed_data <- cbind(plot, total_seed, n_seedlings)

field_data <- merge(field_plots, tree_data, all=TRUE)
field_data <- field_data[-which(is.na(field_data$n_trees)),c(1:11,14:19)]
field_data <- merge(field_data, seed_data, all=TRUE)

#add ferns data to field data:
field_data <- cbind(field_data, hf_ground$f.a)
#write.csv(field_data, file="HF_2022_Field_Data/2023_07_14_hf_field_data_cleaned.csv")



#merge with remote-sensing data:
#load:
hf_mags <- read.csv("HF_2022_Field_Data/HF_mags_recov_from_GEE_data.csv")
#join:
hf_data <- cbind(field_data,hf_mags)

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

#making a new lat/lon file:
# hf_lat_lon <- cbind(field_plots$longitude, field_plots$latitude)
# colnames(hf_lat_lon) <- c('longitude','latitude')
# write.csv(hf_lat_lon, file="hf_lat_lon.csv", row.names=TRUE)

