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

#this will get the DEM for the whole USA:
rastertest <- raster('USA1_msk_alt.gri')

#crop to bounding box:
bounds <- read.csv("2022_03_29_latlonboundingbox.csv")
extents <- c(round(min(bounds$x)),
             round(max(bounds$x)),
             round(min(bounds$y)),
             round(max(bounds$y)))
box <- as(extent(extents), 'SpatialPolygons')
crs(box) <- "+proj=longlat +datum=WGS84 +no_defs"

#project raster:
rastertestproj2 <- projectRaster(rastertest, crs= "+proj=longlat +datum=WGS84 +no_defs",
                                format='GTiff')

#crop to study area bounds:
rastbox <- crop(rastertestproj2, box)

#plot(rastertestproj)
plot(rastbox)






