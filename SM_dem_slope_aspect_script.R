#DEM test script

#getting country code and data:   #says "use geodata package instead"
#getData('ISO3')
#getData('alt', country= 'USA', mask=T)

#load libraries:
library(lattice)
library(geodata)
library(rgdal)
library(rgeos)
library(raster)
library(rasterVis)

#url for DEM data:
#url <- "http://doi.org/10.5067/MEaSUREs/SRTM/SRTMGL1N.003"


### Get DEM data --------------------------------------------------------------

#this will get the DEM for the whole USA:
rastertest <- raster('USA1_msk_alt.gri')

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

