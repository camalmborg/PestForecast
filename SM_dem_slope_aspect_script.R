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
#plot(rastertest)





