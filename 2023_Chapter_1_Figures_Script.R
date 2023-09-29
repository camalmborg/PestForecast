### Making Figures for Chapter 1 Univariate Analyses:

### loading environmental variables data:
#load("DMVARS_MO.RData") ##daymet monthly variables
#load("DMVARS_SEAS.RData") ##daymet seasonal variables
#load("SMAP_data.RData") # col 1-4 for Apr-July 2015 soil moisture
#load("site_DEM_slope_aspect_TWI_data.RData") # col 1 for DEM data
#load("viirs_annual_averages_data.RData")
load("CHAPTER_1/2023_09_distmag_bestunimodels_variables_data.RData")


### loading magnitude data:
#dmr <- read.csv("SM_distmagrecov_data.csv")
#dmr <- read.csv("SM_distmagrecov_data_2016_tcg.csv")
#dmr <- read.csv("SM_distmagrecov_data_2017_tcg.csv")
#dmrcs <- read.csv("SM_distmagrecov_data_2016_cs.csv")
#dmrcs <- read.csv("SM_distmagrecov_data_2017_cs.csv")
dmr_tcg <- read.csv("CHAPTER_1/2023_09_DMR_DATA_CS_2016.csv")
dmr_cs <- read.csv("CHAPTER_1/2023_09_DMR_DATA_TCG_2016.csv")

### loading best models data:
best_tcg <- read.csv("CHAPTER_1/2023_09_29_DISTMAG_2016_TCG_delAICsort.csv")
best_cs <- read.csv("CHAPTER_1/2023_09_29_DISTMAG_2016_CS_delAICsort.csv")

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



### DISTURBANCE MAGNITUDE FIGURES: multivariate

# choose data type tcg or cs
models <- best_tcg
#models <- best_cs

# take top performers:
tops <- models[1:10,grep("^VARIABLE",colnames(models))]


for (i in 1:10){
  
}
