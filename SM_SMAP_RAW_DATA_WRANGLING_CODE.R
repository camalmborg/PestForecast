#this code is for the SMAP data- raw monthly values for April-August 2015-2017 from GEE for 
#univariate analyses. The R SMAP download code is not set up for 5000 sites of
#data and so I'm going rogue to get my initial analyses done.

#load libraries

#load data
file <- "SMAP_Data/SMAP_04_08_2015_2017_try3.txt"
SMAP.GEE <- read.table(file, header=T, sep=",")
#load lat long data from sampling points:
sitesfile <- "2022_03_22_5000sites_lat_long_points_for_GEE_asset.csv"
sites <- read.csv(sitesfile)

#get rid of ] from smap data, make numeric:
SMAPvalues <- SMAP.GEE[,"smp."]
SMAPvalues <- gsub("\\[|\\]", "", SMAPvalues)
SMAPvalues <- gsub("null","NA",SMAPvalues)
SMAPvalues <- as.numeric(SMAPvalues)

#re-make SMAP.GEE object with just lat/long/smap:
SMAPdata <- as.data.frame(cbind(SMAP.GEE[,"longitude"],
                                SMAP.GEE[,"latitude"],
                                SMAPvalues))

nsites=5000
SMAPorg <- matrix(data=NA, nrow=nsites, ncol=15) #15= 3 years X 5 months
for (i in 1:15){
  rows<-seq(i, nrow(SMAPdata), by=15)
  SMAPorg[,i]<-SMAPdata[rows,3]
}
