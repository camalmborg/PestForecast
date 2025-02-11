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
nlcdfile <- "CHAPTER_1/DATA/2020_07_10_sample_nlcd.csv"
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
ptc <- nlcd$percent_tree_cover
# change 41 (deciduous forest class) to 1 for setting reference class in gam
lc <- replace(lc, lc == 41, 1)
# column number for disturbance years
dmr <- dmrcs
cn <- as.matrix(as.numeric(dmr$colnum))
# which response variable? "mags" for dist mag, "dpy1" or "dpy2" for dist prob
#yvar <- "mags"
yvar <- "dpy1"
#yvar <- "dpy2"
y <- as.matrix(as.numeric(dmr[,yvar]))
# make data frame to convert to data
x <- as.data.frame(cbind(y, cn, lc, ptc))
colnames(x) <- c("y", "cn", "lc", "ptc")
# data object for gams
#select disturbance year column number -- 22 for 2016, 23 for 2017
coln = 22
#coln = 23
dat <- x[x$cn == coln,]
# sort to put reference class first:
dat <- dat[order(dat$lc),]

# run model
# disturbance magnitude
nlcd_gam <- gam(y ~ as.factor(lc), data = dat)
# check model summary
summ <- summary(nlcd_gam)
r2 <- summ$r.sq
aic <- nlcd_gam$aic
#ran 2/11/2025: r2 = 0.186744, aic = -4278.522635 for y = mags tcg
#               r2 = 0.093025, aic = 8104.702825 for y = mags CS
#ran 2/11/2025: r2 = 0.090457, aic = 2002.571447 for y = dpy1 tcg
#               r2 = 0.0707096, aic = 1771.436604 for y = dpy1 CS

# do with percent tree cover
nlcd_ptc_gam <- gam(y ~ as.factor(ptc), data = dat)
summ <- summary(nlcd_ptc_gam)
r2 <- summ$r.sq
aic <- nlcd_ptc_gam$aic
#ran 2/11/2025 r2 = 0.0984026, aic = -4084.3285187 for y = mags tcg
#              r2 = 0.0500345, aic = 8195.4280983 for y = mags cs
#ran 2/11/2025 r2 = 0.0457663, aic = 2098.2506737 for y = dpy1 tcg
#              r2 = 0.0399439, aic = 1838.1852318 for y = dpy1 cs

# disturbance probability
nlcd_gam <- gam(y ~ as.factor(lc), data = dat, family = "binomial", method = "REML")

# if doing disturbance probability, ROCs:
# run ROC
nlcd_roc <- multiclass.roc(dmr$dpy1, lc)
roc <- nlcd_roc$auc
aic <- nlcd_gam$aic
#ran 2/11/25 roc = 0.454, aic = 1938.28555
#            roc = 0.0399943, aic - 1771.436604
# notes:
# convert categories to factor
# 
