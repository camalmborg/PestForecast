### SCRIPT FOR NLCD ANALYSES 03/2024

### Load libraries:
library(dplyr)
library(tidyverse)
library(mgcv)
library(pROC)


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
#nlcd$landcover <- as.factor(nlcd$landcover)
# select unique landcovers for reference later
nlcd_cover <- unique(nlcd$landcover)


### run some GAMs
# make data object
# nlcd variable
lc <- nlcd$landcover
# change 41 (deciduous forest class) to 1 for setting reference class in gam
lc <- replace(lc, lc == 41, 1)
# column number for disturbance years
cn <- as.matrix(as.numeric(dmr$colnum))
# which response variable? "mags" for dist mag, "dpy1" or "dpy2" for dist prob
yvar <- "mags"
#yvar <- "dpy1"
#yvar <- "dpy2"
y <- as.matrix(as.numeric(dmr[,yvar]))
# make data frame to convert to data
x <- as.data.frame(cbind(y, cn, lc))
colnames(x) <- c("y", "cn", "lc")
# data object for gams
#select disturbance year column number -- 22 for 2016, 23 for 2017
coln = 22
#coln = 23
dat <- x[x$cn == coln,]
# sort to put reference class first:
dat <- dat[order(dat$lc),]

# run model
# disturbance magnitude
nlcd_gam <- gam(y ~ as.factor(lc), data = dat, method = "REML")
# disturbance probability
#nlcd_gam <- gam(y ~ as.factor(lc), data = dat, family = "binomial", method = "REML")

# check model summary
summ <- summary(nlcd_gam)
r2 <- summ$r.sq
#aic <- nlcd_gam$aic

# if doing disturbance probability, ROCs:
# run ROC
nlcd_roc <- roc(dat$y, nlcd_gam$fitted.values)
roc <- nlcd_roc$auc
#aic <- nlcd_gam$aic


# notes:
# convert categories to factor
# 
