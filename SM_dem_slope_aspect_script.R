#DEM test script

#getting country code and data:   #says "use geodata package instead"
#getData('ISO3')
#getData('alt', country= 'USA', mask=T)

#load libraries:
library(lattice)
library(geodata)
#library(rgdal)  ##deprecating?
library(rgeos)
library(terra)
library(sp)
library(sf)
library(tidyr)
library(stringr)
library(tidyverse)
#library(raster)  ##deprecated
#library(rasterVis)

#url for DEM data:
#url <- "http://doi.org/10.5067/MEaSUREs/SRTM/SRTMGL1N.003"

### GET DEM DATA AND SLOPE/ASPECT CALCULATIONS USING TERRA PACKAGE:------------

#load tif file:
dem_file <- "DEM_Data/2023_06_07_NE_MERGE_DEM.tif"
NE_rast <- rast(dem_file)

#load TWI tif:
twi_file <- "DEM_Data/2023_06_19_TWI_DZ.tif"
NE_twi <- rast(twi_file)

#load sites:
file<-"2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"
cfile<-"2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv"

#tcg and condition score objects:
tcg.values<-read.csv(file)
#tcgs<-tcg.values[,c(grep("^X",colnames(tcg.values)))] #grab with "X1996 eg) for all tcg value columns
cond.scores<-read.csv(cfile)

#geographic coordinates from GEE extract:
geo<-as.data.frame(cond.scores[,".geo"])

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
colnames(coords)<-c("x","y")

#site_lat_lon <- read.csv("2022_03_22_5000sites_lat_long_points_for_GEE_asset.csv")
#coords <- site_lat_lon[,c(2,1)] #correcting order for conversion to spatialpointsdataframe obj
#set site crs:
#site_crs <- crs(NE_rast)
#convert sites to spatial:
#coordinates(sites) <- ~y+x  ##makes a SpatialPoints object
# sites <- SpatialPointsDataFrame(coords = coords, 
#                                 data = site_lat_lon,
#                                 proj4string = CRS(site_crs))



#plot(sites)

#convert sites to vector

#calculating slope and aspect from DEM:
aspect_r <- terrain(NE_rast, "aspect", unit = "radians")
slope_r <- terrain (NE_rast, "slope", unit = "radians")

aspect_d <- terrain(NE_rast, "aspect", unit="degrees")
slope_d <- terrain(NE_rast, "slope", unit="degrees")

#stack 'em:
NE_dem_data <- c(NE_rast, slope_d, aspect_d, NE_twi)

#extracting site values: 
site_data <- terra::extract(NE_dem_data, coords)
#note: if tidyverse is employed terra needs to explictly be called for extract()

#saving site data:
write.csv(site_data, file="2023_06_26_5000sample_site_DEM_slope_aspect_TWI_data.csv")
save(site_data, file="site_DEM_slope_aspect_TWI_data.RData")

### Get DEM data ARCHIVED: USED RASTER PACKAGE-------------------------------

#this will get the DEM for the whole USA:
#raster package becoming unavailable, replaced by terra
#getData("ISO3")
#getData('alt', country='USA', mask=TRUE)
#rastertest <- raster('USA1_msk_alt.grd')

# #crop to bounding box:
# bounds <- read.csv("2022_03_29_latlonboundingbox.csv")
# extents <- c(min(bounds$y),
#              max(bounds$y),
#              min(bounds$x),
#              max(bounds$x))
# box <- as(extent(extents), 'SpatialPolygons')
# crs(box) <- "+proj=longlat +datum=WGS84 +no_defs"
# #study area:
# rastbox <- crop(rastertest, box)

# #re-project raster and convert to gtiff:
# rastboxproj <- projectRaster(rastbox, 
#                              crs= "+proj=longlat +datum=WGS84 +no_defs",
#                              format='GTiff')
# 
# 
# #plot(rastertest)
# #plot(rastboxproj)

### Make spatial dataset from csv of lat/lon points for sites--

# #load csv with lat/lon coords:
# sites <- read.csv("2022_03_22_5000sites_lat_long_points_for_GEE_asset.csv")
# #convert to spatial:
# coordinates(sites) <- ~y+x
# crs(sites) <- crs(rastbox)
# 
# #plot(sites)


### Get slope and aspect data--

# #slope:
# slopes <- terrain(rastbox, opt="slope", unit="degrees")
# #aspect:
# aspects <- terrain(rastbox, opt="aspect", unit="degrees")
# 
# #make stack:
# datastack <- stack(rastbox, slopes, aspects)
# 
# #extract site values:
# DEMdata <- extract(datastack, sites)
# #save a csv and GEOTIFF of this -- 5/19


### Univariate Analyses ------------------------------------------

#load libraries:
library(mgcv)

#load testfx object if not in environment:
testfx <- read.csv("SM_distmagrecov_data.csv")

#function for univariate analyses:
DSA_explore<-function(dem,dmrdat,dmr,coln){
  #make empty matrix:
  r2s <- matrix(NA,nrow=1,ncol=ncol(dem))
  
  #loop over all members of dmvars list:
  for (i in 1:ncol(dem)){
    #first extract the list you want:
    demvariables <- as.data.frame(dem[as.numeric(dmrdat$sitenum),])
    #grab column number:
    cn <- as.matrix(as.numeric(dmrdat$colnum))
    #grab exploratory response  of choice:
    y <- as.matrix(as.numeric(dmrdat[,dmr]))
    x <- as.data.frame(cbind(y,cn,demvariables))
    vardat <- x[x$cn==coln,]
    
    #run the GAM:
    var.gam <- gam(vardat[,1]~s(vardat[,i+2]),data=vardat)
      
    #extract r2:
    summ <- summary(var.gam)
    r2s[,i] <-summ$r.sq
  }
  # return(vardat)
  return(r2s)
}

###BEFORE RUNNING: make sure coords match between mags and dem data
testing_dem <- DSA_explore(site_data,testfx,"mags",22)


### Plots -------------------------------------------------------
#----------#### MAKING FIGURES: ####------------------------
library(lattice)
library(latticeExtra)
library(tactile)
library(mgcv)

vardat <- vardat
mo <- vardat[,4]
var.gam <- var.gam <- gam(vardat[,1]~s(mo),data=vardat)

dmplot<-xyplot(y ~ mo, data = vardat,
               panel = function(x, y) {
                 ci<-predict(var.gam, se=T)
                 ci$lower<-ci$fit-qt(0.975,var.gam$df.null)*ci$se.fit
                 ci$upper<-ci$fit+qt(0.975,var.gam$df.null)*ci$se.fit
                 l.ci<-cbind(var.gam$model$mo,ci$fit,ci$lower,ci$upper)
                 l<-l.ci[order(l.ci[,1]),]
                 panel.ci(l[,1],l[,2],l[,4],l[,3],
                          fill="seagreen3",alpha = 0.3)
                 panel.xyplot(x, y, pch=20,col="seagreen")
                 panel.lines(l[,1], l[,2],lty=1, col='black', lwd=1.5)
                 summ<-summary(var.gam)
                 r2 <- summ$r.sq
                 #f <- summ$fstatistic
                 # p <- pf(f[1],f[2],f[3],lower.tail=F)
                 panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                            x=2150,y=-0.12,cex=0.75)
                 # panel.text(labels = bquote(italic(p)==.(format(p,digits = 3))),
                 #            x=1.2,y=-0.15,cex=0.75)
               },
               ylab="Disturbance Magnitude (TCG)",
               xlab="DEM Slope",
)
print(dmplot)
