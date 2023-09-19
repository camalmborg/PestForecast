### Making Figures for Chapter 1 Univariate Analyses:

### loading environmental variables data:
load("DMVARS_MO.RData") ##daymet monthly variables
load("DMVARS_SEAS.RData") ##daymet seasonal variables
load("SMAP_data.RData") # col 1-4 for Apr-July 2015 soil moisture
#load("site_DEM_slope_aspect_TWI_data.RData") # col 1 for DEM data
#load("viirs_annual_averages_data.RData")

### loading magnitude data:
#dmr <- read.csv("SM_distmagrecov_data.csv")
#dmr <- read.csv("SM_distmagrecov_data_2016_tcg.csv")
#dmr <- read.csv("SM_distmagrecov_data_2017_tcg.csv")
dmrcs <- read.csv("SM_distmagrecov_data_2016_cs.csv")
#dmrcs <- read.csv("SM_distmagrecov_data_2017_cs.csv")


### Libraries for figures:
library(ggplot2)
library(tidyverse)
library(mgcv)
library(tidygam)
library(pROC)
library(viridis)

### DISTURBANCE MAGNITUDE FIGURES: univariate ------------

# data for gam:
yvar <- dmrcs$mags
xvar <- dmvars_mo[[4]][dmrcs$sitenum,43] 
#xvar <- SMAPdat[dmrcs$sitenum,4]
data <- cbind.data.frame(yvar, xvar)
# run gam:
var_gam <- gam(yvar ~ s(xvar), 
              data=data)

#fit <- var_gam$fitted.values

# plot:
ggplot(data=data, aes(x=xvar, y=yvar)) +
  geom_point(color="brown") +
  geom_smooth(method = "gam") +
  ggtitle("July 2014 VPD")

