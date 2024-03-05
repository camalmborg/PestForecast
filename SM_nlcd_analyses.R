### SCRIPT FOR NLCD ANALYSES 03/2024

### Load libraries:
library(dplyr)
library(tidyverse)
library(mgcv)

### Load DMR data:
# disturbance probability and magnitude data
# tcg version
dmr_file <- "CHAPTER_1/DATA/2023_12_DMR_DATA_TCG.csv"
dmr <- read.csv(dmr_file)
# cs version
# the distmagrecov data CS
dmrcs <- read.csv("CHAPTER_1/DATA/2023_12_DMR_DATA_CS.csv")
# site data
# tcg values
tfile <- "2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"
tcg <- read.csv(tfile)
# condition scores
cfile <- "2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv" #1995-2020, original grab
scores <- read.csv(cfile)
# NLCD data
nlcdfile <- "data_archive/2020_07_10_500_sample_nlcd.csv"
nlcd <- read.csv(nlcdfile) ## checked if sites match, they do, we're good!
# select sites without missing dmr data
nlcd <- nlcd[dmr$X,]
# convert landcover type to factor
nlcd$landcover <- as.factor(nlcd$landcover)
# select unique landcovers for reference later
nlcd_cover <- unique(nlcd$landcover)


### run some GAMs


# notes:
# convert categories to factor
# 
