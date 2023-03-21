#this code is for the SMAP data- raw monthly values for April-August 2015-2017 from GEE for 
#univariate analyses. The R SMAP download code is not set up for 5000 sites of
#data and so I'm going rogue to get my initial analyses done.

#load libraries
library(mgcv)
library(tidyr)
library(stringr)

#load data
#file<-"2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"
#cfile <- cfile<-"2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv"
sfile <- "SMAP_Data/SMAP_04_08_2015_2017_try3.txt"
SMAP.GEE <- read.table(sfile, header=T, sep=",")
#load lat long data from sampling points:
sitesfile <- "2022_03_22_5000sites_lat_long_points_for_GEE_asset.csv"
sites <- read.csv(sitesfile)
#cond.scores<-read.csv(cfile)


#get rid of ] from smap data, make numeric:
SMAPvalues <- SMAP.GEE[,"smp."]
SMAPvalues <- gsub("\\[|\\]", "", SMAPvalues)
SMAPvalues <- gsub("null","NA",SMAPvalues)
SMAPvalues <- as.numeric(SMAPvalues)

#re-make SMAP.GEE object with just lat/long/smap:
SMAPdata <- as.data.frame(cbind(SMAP.GEE[,"longitude"],
                                SMAP.GEE[,"latitude"],
                                SMAPvalues))

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
# colnames(coords)<-c("lon","lat","site_id")
coords <- cbind(c(1:5000),sites[,2],sites[,1])
colnames(coords)<-c("site_id","lon","lat")

#how many sites:
nsites=nrow(sites)
#loop for organizing SMAP data into monthly columns 2015-2017:
SMAPorg <- matrix(data=NA, nrow=nsites, ncol=15) #15= 3 years X 5 months
for (i in 1:15){
  rows<-seq(i, nrow(SMAPdata), by=15)
  SMAPorg[,i]<-SMAPdata[rows,3]
}


####matching SMAP lat/lon to mags lat/lon
#make data frame with coords and mags:
#magsdat <- as.data.frame(cbind(coords[-missing,],mags))
#magscoords <- magsdat[order(magsdat[,1]),]

sortcoords <- as.data.frame(coords[order(coords[,2]),]) #with full 5000

#remove duplicate values in SMAP lat/lon:
SMAPlonlat <- as.data.frame(unique(SMAPdata[,1:2]))

#make data frame with SMAP coords and SMAP data:
SMAPall <- as.data.frame(cbind(SMAPlonlat,SMAPorg))
SMAPcoords <- SMAPall[order(SMAPall[,1]),]
SMAPsmap <- cbind(sortcoords$site_id,SMAPcoords)
sortSMAP <- SMAPsmap[order(SMAPsmap[,1]),]
#removing correct missing values:
#sortSMAP <- sortSMAP[-missing,]


#make new data frame:
SMAPmags <- cbind(magsdat, sortSMAP[,4:ncol(sortSMAP)])
#smap<-SMAPmags[,5:ncol(SMAPmags)]


#####-------------------ANALYSES------------------#####################

#function for univariate analyses:
smap_explore<-function(smap,dmrdat,dmr,coln){
  #make empty matrix:
  r2s <- matrix(NA,nrow=length(smap)-4,ncol=1)  #hard coded rn
  
  #loop over all members of dmvars list:
  for (i in 1:15){  #this is hard coded until further notice
    #first get just SMAP data:
    smapvar <- as.data.frame(smap[5:ncol(smap)])
    #grab column number:
    cn <- as.matrix(as.numeric(dmrdat$colnum))
    #grab exploratory response  of choice:
    y <- as.matrix(as.numeric(dmrdat[,dmr]))
    x <- as.data.frame(cbind(y,cn,smapvar))
    vardat <- x[x$cn==coln,]
    
    #loop for filling in R2 table:  
    var.gam <- gam(vardat[,1]~s(vardat[,i+2]),data=vardat)

    #extract r2:
    summ <- summary(var.gam)
    r2s[i,] <-summ$r.sq
  }
  return(r2s)
}

testing_smap <- smap_explore(SMAPmags,testfx,"mags",22)
