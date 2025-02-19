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
dmrtcg <- read.csv(dmr_file)
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
nlcd <- nlcd[dmrtcg$X,]
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
#dmr <- dmrtcg
dmr <- dmrcs
cn <- as.matrix(as.numeric(dmr$colnum))
# which response variable? "mags" for dist mag, "dpy1" or "dpy2" for dist prob
#yvar <- "mags"
#yvar <- "dpy1"   # for binomial rocs
yvar <- "dpy2"  # for binomial rocs
y <- as.matrix(as.numeric(dmr[,yvar]))
# make data frame to convert to data
x <- as.data.frame(cbind(y, cn, lc, ptc))
colnames(x) <- c("y", "cn", "lc", "ptc")
# data object for gams
#select disturbance year column number -- 22 for 2016, 23 for 2017
#coln = 22 #if using dpy1 (2016)
coln = 23 #if using dpy2 (2017)
dat <- x[x$cn == coln,]
# sort to put reference class first:
dat <- dat[order(dat$lc),]

# run models
# disturbance magnitude ----
nlcd_gam <- gam(y ~ as.factor(lc), data = dat)
# check model summary
summ <- summary(nlcd_gam)
r2 <- summ$r.sq
aic <- nlcd_gam$aic
#ran 2/11/2025: r2 = 0.186744, aic = -4278.522635 for y = mags tcg 2016
#               r2 = 0.093025, aic = 8104.702825 for y = mags CS 2016
#               r2 = 0.05882254795, aic = -8017.91914090441 for y = mags tcg 2017
#               r2 = 0.031488367423, aic = 14661.36496287 for y mags CS 2017


# do with percent tree cover
nlcd_ptc_gam <- gam(y ~ as.factor(ptc), data = dat)
summ <- summary(nlcd_ptc_gam)
r2 <- summ$r.sq
aic <- nlcd_ptc_gam$aic
#ran 2/11/2025 r2 = 0.0984026, aic = -4084.3285187 for y = mags tcg 2016
#              r2 = 0.0500345, aic = 8195.4280983 for y = mags CS 2016
#              r2 = 0.02555262694, aic = -7897.868254918 for y = mags tcg 2017
#              r2 = 0.01242232680, aic = 14733.4843811874 for y = mags CS 2017



# disturbance probability ----
nlcd_gam <- gam(y ~ as.factor(lc), data = dat, family = "binomial", method = "REML")

# if doing disturbance probability, ROCs:
# for land cover class (factor data)
# run ROC
nlcd_roc <- multiclass.roc(dmr$dpy2, lc)  #multiclass for factor data
roc <- nlcd_roc$auc
aic <- nlcd_gam$aic
#ran 2/11/25 roc = 0.4535, aic = 1938.28555 tcg 2016 dyp1
#            roc = 0.4535, aic = 1748.74447876079 cs 2016 dpy1
#            roc = 0.4507, aic = 4327.88200002787 tcg 2017 dpy2 tcg 2017
#            roc = 0.4516, aic = 4393.57364950584 cs 2017 dpy2

nlcd_gam <- gam(y ~ ptc, data = dat, family = "binomial", method = "REML")
# run ROC
# for percent tree cover:
nlcd_roc <- roc(dmr$dpy2, ptc)
roc <- nlcd_roc$auc
aic <- nlcd_gam$aic
#ran 2/11/25 roc = 0.5260, aic = 2030.5319533483 for tcg 2016 dpy1
#            roc = 0.5197, aic = 1803.09157065308 for cs 2016 dpy1
#            roc = 0.5632, aic = 4356.18818246671 for tcg 2017 dpy2
#            roc = 0.5691, aic = 4439.44611045792 for cs 2017 dpy2






# notes:
# convert categories to factor
# 
