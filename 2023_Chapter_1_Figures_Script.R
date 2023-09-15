### Making Figures for Chapter 1 Univariate Analyses:

# loading data:

load("DMVARS_MO.RData") ##daymet monthly variables
load("DMVARS_SEAS.RData") ##daymet seasonal variables
#load("SMAP_data.RData") # col 1-4 for Apr-July 2015 soil moisture
#load("site_DEM_slope_aspect_TWI_data.RData") # col 1 for DEM data
#load("viirs_annual_averages_data.RData")

#loading magnitude data:
#dmr <- read.csv("SM_distmagrecov_data.csv")
dmr <- read.csv("SM_distmagrecov_data_2016_tcg.csv")
dmr <- read.csv("SM_distmagrecov_data_2017_tcg.csv")
dmrcs <- read.csv("SM_distmagrecov_data_2016_cs.csv")
dmrcs <- read.csv("SM_distmagrecov_data_2017_cs.csv")