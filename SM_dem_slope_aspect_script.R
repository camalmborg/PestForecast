#DEM test script

#getting country code and data:   #says "use geodata package instead"
#getData('ISO3')
#getData('alt', country= 'USA', mask=T)

#load libraries:
library(lattice)
library(geodata)
library(rgdal)
library(rgeos)
library(terra)
#library(raster)  ##deprecated
#library(rasterVis)

#url for DEM data:
#url <- "http://doi.org/10.5067/MEaSUREs/SRTM/SRTMGL1N.003"


### Get DEM data --------------------------------------------------------------

#this will get the DEM for the whole USA:
#raster package becoming unavailable, replaced by terra
rastertest <- raster('USA1_msk_alt.tif')

#crop to bounding box:
bounds <- read.csv("2022_03_29_latlonboundingbox.csv")
extents <- c(min(bounds$y),
             max(bounds$y),
             min(bounds$x),
             max(bounds$x))
box <- as(extent(extents), 'SpatialPolygons')
crs(box) <- "+proj=longlat +datum=WGS84 +no_defs"
#study area:
rastbox <- crop(rastertest, box)

#re-project raster and convert to gtiff:
rastboxproj <- projectRaster(rastbox, 
                             crs= "+proj=longlat +datum=WGS84 +no_defs",
                             format='GTiff')


#plot(rastertest)
#plot(rastboxproj)

### Make spatial dataset from csv of lat/lon points for sites------------------

#load csv with lat/lon coords:
sites <- read.csv("2022_03_22_5000sites_lat_long_points_for_GEE_asset.csv")
#convert to spatial:
coordinates(sites) <- ~y+x
crs(sites) <- crs(rastbox)

#plot(sites)


### Get slope and aspect data--------------------------------------------------

#slope:
slopes <- terrain(rastbox, opt="slope", unit="degrees")
#aspect:
aspects <- terrain(rastbox, opt="aspect", unit="degrees")

#make stack:
datastack <- stack(rastbox, slopes, aspects)

#extract site values:
DEMdata <- extract(datastack, sites)
#save a csv and GEOTIFF of this -- 5/19


### Univariate Analyses ------------------------------------------

#load libraries:
library(mgcv)

#function for univariate analyses:
DSA_explore<-function(dem,x,y,z){
  #make empty matrix:
  r2s <- matrix(NA,nrow=ncol(dem[[1]]),ncol=length(dmvars))
  
  #loop over all members of dmvars list:
  for (i in 1:length(dmvars)){
    #first extract the list you want:
    dmvariable <- as.data.frame(dmvars[[i]][as.numeric(dmrdat$sitenum),])
    #grab column number:
    cn <- as.matrix(as.numeric(dmrdat$colnum))
    #grab exploratory response  of choice:
    y <- as.matrix(as.numeric(dmrdat[,dmr]))
    x <- as.data.frame(cbind(y,cn,dmvariable))
    vardat <- x[x$cn==coln,]
    
    #loop for filling in R2 table:  
    for (j in 1:ncol(dmvariable)){
      var.gam <- gam(vardat[,1]~s(vardat[,j+2]),data=vardat)
      
      #extract r2:
      summ <- summary(var.gam)
      r2s[j,i] <-summ$r.sq
    }
  }
  # return(vardat)
  return(r2s)
}

testing <- dm_explore(spongyvars,testfx,"mags",22)