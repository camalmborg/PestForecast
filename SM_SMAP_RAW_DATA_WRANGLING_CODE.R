#this code is for the SMAP data- raw monthly values for April-August 2015-2017 from GEE for 
#univariate analyses. The R SMAP download code is not set up for 5000 sites of
#data and so I'm going rogue to get my initial analyses done.

##EVERYTHING here is hard coded and not super useful if you don't have my environments loaded, sorry everyone

#load libraries
library(mgcv)
library(dplyr)
library(tidyverse)
library(tidyr)
library(stringr)

#load data
#file<-"2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"
#cfile <- cfile<-"2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv"
#sfile <- "SMAP_Data/SMAP_04_08_2015_2017_try3.txt"
sfile <- "SMAP_Data/2024_04_23_SMAP_04_08_2015_2017_5000sites.csv"
SMAP.GEE <- read.table(sfile, header=T, sep=",")
#load lat long data from sampling points:
sitesfile <- "2022_03_22_5000sites_lat_long_points_for_GEE_asset.csv"
sites <- read.csv(sitesfile)
#cond.scores<-read.csv(cfile)


#get rid of ] from smap data, make numeric:
SMAPvalues <- SMAP.GEE %>%
  # Drop unwanted columns
  dplyr::select(dplyr::starts_with("X")) %>%
  # remove X from dates
  dplyr::rename_with(~ str_replace_all(., c("X" = ""))) %>%
  # put in order
  dplyr::rename_with(~ str_replace_all(., c(".*_" = "")))
# remove . in date and replace with -
colnames(SMAPvalues) <- gsub("\\.", "-0", colnames(SMAPvalues))
# make into dates
colnames(SMAPvalues) <- as.Date(paste0(names(SMAPvalues), "-01", format = '%d %m %Y'))
# sort by date
SMAPvalues <- SMAPvalues[,sort(names(SMAPvalues))]

#re-make SMAP.GEE object with just lat/long/smap:
geo<-as.data.frame(SMAP.GEE[,".geo"])
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
colnames(coords)<-c("lon","lat","site_id")
# combine with SMAPvalues
SMAPdata <- as.data.frame(cbind("site_id" = coords[,"site_id"],
                                "lon" = coords[,"lon"],
                                "lat" = coords[,"lat"],
                                SMAPvalues))

# for finding missing values:
missingSMAP <- SMAPdata[!complete.cases(SMAPvalues), 1:3]
# save CSV
write.csv(missingSMAP, file = "SMAP_Data/2024_04_23_missing_SMAP_values.csv")


#####-------------------ANALYSES------------------#####################

#function for univariate analyses:
smap_explore<-function(smap,dmrdat,dmr,coln){
  #make empty matrix:
  r2s <- matrix(NA,nrow=length(smap)-4,ncol=1)  #hard coded rn
  
  #first get just SMAP data:
  smapvar <- as.data.frame(smap[5:ncol(smap)])
  #grab column number:
  cn <- as.matrix(as.numeric(dmrdat$colnum))
  #grab exploratory response  of choice:
  y <- as.matrix(as.numeric(dmrdat[,dmr]))
  x <- as.data.frame(cbind(y,cn,smapvar))
  vardat <- x[x$cn==coln,]
  
  #loop over all members of dmvars:
  for (i in 1:15){  #this is hard coded until further notice
    #loop for filling in R2 table:  
    var.gam <- gam(vardat[,1]~s(vardat[,i+2]),data=vardat)

    #extract r2:
    summ <- summary(var.gam)
    r2s[i,] <-summ$r.sq
  }
  return(r2s)
}

testing_smap_2 <- smap_explore(SMAPmags,testfx,"mags",22)


#####-----Let's make some plots---------------########
library(mgcv)
library(lattice)
library(latticeExtra)
library(tactile)

mo<-SMAPmags[,"1"]  #what month column here
smap.gam <- gam(mags~s(mo), data = SMAPmags)

smapplot<-xyplot(mags ~ mo, data = SMAPmags,
                panel = function(x, y) {
                  ci<-predict(smap.gam, se=T)
                  ci$lower<-ci$fit-qt(0.975,smap.gam$df.null)*ci$se.fit
                  ci$upper<-ci$fit+qt(0.975,smap.gam$df.null)*ci$se.fit
                  l.ci<-cbind(smap.gam$model$mo,ci$fit,ci$lower,ci$upper)
                  l<-l.ci[order(l.ci[,1]),]
                  panel.ci(l[,1],l[,2],l[,4],l[,3],
                           fill="seagreen3",alpha = 0.3)
                  panel.xyplot(x, y, pch=20,col="seagreen")
                  panel.lines(l[,1], l[,2],lty=1, col='black', lwd=1.5)
                  summ<-summary(smap.gam)
                  r2 <- summ$r.sq
                  #f <- summ$fstatistic
                  # p <- pf(f[1],f[2],f[3],lower.tail=F)
                  panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                             x=2150,y=-0.12,cex=0.75)
                  # panel.text(labels = bquote(italic(p)==.(format(p,digits = 3))),
                  #            x=1.2,y=-0.15,cex=0.75)
                },
                ylab="Disturbance Magnitude (TCG)",
                xlab="Mean Soil Moisture (SMAP)",
)
print(smapplot)



#####-----Disturbance Probabilities Analyses----------#####
library(mgcv)
library(pROC)

smap_ROC <- function(dmvars,dmrdat,yr,coln){
  #make empty matrix:
  rocs <- matrix(NA,nrow=length(smap),ncol=1)  #hard coded rn
  
  #grab just disturbance probability columns:
  dists<-dmrdat[,grep("^dp",colnames(dmrdat))]
  
  #first get just SMAP data:
  smapvar <- as.data.frame(smap)
  #grab column number:
  cn <- as.matrix(as.numeric(dmrdat$colnum))
  #grab exploratory response  of choice:
  y <- as.matrix(as.numeric(dists[,yr]))
  x <- as.data.frame(cbind(y,cn,smapvar))
  vardat <- x[x$cn==coln,]
  
  #remove missing values:
  miss <- which(is.na(vardat[,3]))
  vardat <- vardat[-miss,]
  
  #loop over all members of dmvars:
  for (i in 1:15){  #this is hard coded until further notice
    #loop for filling in R2 table:  
    var.gam<-gam(vardat[,1]~s(vardat[,i+2]), data=vardat, family="binomial")
    var.roc<-roc(vardat[,1],var.gam$fitted.values)
    rocs[i,] <-var.roc$auc
  }
  #return table of AUCs
  return(rocs)
}

SMAProcs <- smap_ROC(dmvars,dmrdat,1,22)


### ARCHIVE:
#gsub(".", "-", colnames(SMAPvalues))  
#dplyr::rename_with(~ str_replace_all(., c(".." = "-")))
# # rename with date, without extra characters
# dplyr::rename_with(~ str_replace_all(., c("X|_score_mean|_cs_mean" = "",
#                                           "\\." = "-")))
#SMAPvalues <- SMAP.GEE[,"smp."]
#SMAPvalues <- SMAP.GEE[,(grep("^X",colnames(SMAP.GEE)))]
#SMAPvalues <- gsub("\\[|\\]", "", SMAPvalues)
#SMAPvalues <- gsub(" null","NA",SMAPvalues)
#SMAPvalues <- as.numeric(SMAPvalues)

#make coords object:
# geo<-as.data.frame(cond.scores[,".geo"])
# #make lat and lon columns from .geo data:
# coords<-matrix(nrow=nrow(geo),ncol=3)
# for (i in 1:nrow(geo)){
#   #longitudes:
#   lon<-str_extract(geo[i,], "\\d+\\.*\\d*")
#   coords[i,1]<-as.numeric(lon)*-1
#   
#   #latitudes:
#   extlon<-sub(lon,"",geo[i,])
#   coords[i,2]<-as.numeric(str_extract(extlon, "\\d+\\.*\\d*"))
#   coords[i,3]<-i
# }
# # colnames(coords)<-c("lon","lat","site_id")
# coords <- cbind(c(1:5000),sites[,2],sites[,1])
# colnames(coords)<-c("site_id","lon","lat")

# #how many sites:
# nsites=nrow(sites)
# #loop for organizing SMAP data into monthly columns 2015-2017:
# SMAPorg <- matrix(data=NA, nrow=nsites, ncol=15) #15= 3 years X 5 months
# for (i in 1:15){
#   rows<-seq(i, nrow(SMAPdata), by=15)
#   SMAPorg[,i]<-SMAPdata[rows,3]
# }


# ####matching SMAP lat/lon to mags lat/lon
# #make data frame with coords and mags:
# #magsdat <- as.data.frame(cbind(coords[-missing,],mags))
# #magscoords <- magsdat[order(magsdat[,1]),]
# 
# sortcoords <- as.data.frame(coords[order(coords[,2]),]) #with full 5000
# 
# #remove duplicate values in SMAP lat/lon:
# SMAPlonlat <- as.data.frame(unique(SMAPdata[,1:2]))
# 
# #make data frame with SMAP coords and SMAP data:
# SMAPall <- as.data.frame(cbind(SMAPlonlat,SMAPorg))
# SMAPcoords <- SMAPall[order(SMAPall[,1]),]
# SMAPsmap <- cbind(sortcoords$site_id,SMAPcoords)
# sortSMAP <- SMAPsmap[order(SMAPsmap[,1]),]
#removing correct missing values:
#sortSMAP <- sortSMAP[-missing,]

# # 4/2/2024 - missing SMAP values
# missingSMAP <- SMAPsmap[which(is.na(SMAPsmap[,4])),]
# missingSMAPcoords <- missingSMAP[,2:3]
# 
# #make new data frame:
# #SMAPmags <- cbind(magsdat, sortSMAP[,4:ncol(sortSMAP)])
# #smap<-SMAPmags[,5:ncol(SMAPmags)]
# SMAPdat <- sortSMAP[,4:ncol(sortSMAP)]
# #save data:
# write.csv(SMAPdat, file = "2023_06_26_SMAP_data_sorted.csv")
# save(SMAPdat, file = "SMAP_data.Rdata")
