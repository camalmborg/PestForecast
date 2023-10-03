### Making Figures for Chapter 1 Univariate Analyses:

### loading environmental variables data:
#load("DMVARS_MO.RData") ##daymet monthly variables
#load("DMVARS_SEAS.RData") ##daymet seasonal variables
#load("SMAP_data.RData") # col 1-4 for Apr-July 2015 soil moisture
#load("site_DEM_slope_aspect_TWI_data.RData") # col 1 for DEM data
#load("viirs_annual_averages_data.RData")
load("CHAPTER_1/2023_09_distmag_bestunimodels_variables_data.RData")
load("CHAPTER_1/2023_09_distprob_bestunimodels_variables_data.RData")

### loading magnitude data:
#dmr <- read.csv("SM_distmagrecov_data.csv")
#dmr <- read.csv("SM_distmagrecov_data_2016_tcg.csv")
#dmr <- read.csv("SM_distmagrecov_data_2017_tcg.csv")
#dmrcs <- read.csv("SM_distmagrecov_data_2016_cs.csv")
#dmrcs <- read.csv("SM_distmagrecov_data_2017_cs.csv")
dmr_tcg <- read.csv("CHAPTER_1/2023_09_DMR_DATA_TCG.csv")
dmr_cs <- read.csv("CHAPTER_1/2023_09_DMR_DATA_CS.csv")

### loading best models data:
best_tcg <- read.csv("CHAPTER_1/2023_09_29_DISTMAG_2016_TCG_delAICsort.csv")
best_cs <- read.csv("CHAPTER_1/2023_09_29_DISTMAG_2016_CS_delAICsort.csv")

### Libraries for figures:
library(ggplot2)
library(tidyverse)
library(mgcv)
library(mgcViz)
library(gratia)
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



### DISTURBANCE MAGNITUDE FIGURES: multivariate ------

# MAKING THE GAM:
# choose data type tcg or cs:
#models <- dm_bio_best_cs
models <- dp17_bio_best_tcg

# choose matching mags data:
#mags <- dmr_tcg[,c("sitenum","mags")]
#mags <- dmr_cs[,c("sitenum","mags")]
probs <- dmr_tcg[,c("sitenum", "dpy2")]
#probs <- dmr_cs[,c("sitenum", "dpy2")]

# take top performers:
tops <- models[1:25,grep("^VARIABLE",colnames(models))]

# where should the plots be saved?
plot_folder <- "Analyses_September2023/Figures_Multi/dp17_bio_best/"
plot_type <- "tcg/"
#plot_type <- "cs/"

### dist mag plots: -----
for (i in 1:nrow(tops)){
  #start the loop here >>>
  # grab model variable column numbers:
  v1 <- as.numeric(tops[i, 1])
  v2 <- as.numeric(tops[i, 2])
  v3 <- as.numeric(tops[i, 3])
  
  # grab variables:
  vars <- cbind(MV_2023_09_DATA[,c(v1,v2,v3)])
  var_nms <- as.character(c(v1,v2,v3))
  plot_nm <- paste0("2023_09_plot_", as.character(i), "_vars_",
                    paste(var_nms, collapse = "_"))

  # combine data:
  data <- cbind.data.frame(mags$mags, vars[mags$sitenum,])
  colnames(data) <- c("mags", "var1", "var2", "var3")

  # make the formula:
  ex_vars <- c()
  for (j in 1:ncol(data)-1){
    ex_vars[j] <- paste0('s(', names(data[j+1]), ')')
  }

  #make a single string:
  gam_formula <- as.formula(paste("mags ~",
                                  paste(ex_vars, collapse='+')))

  # make the model:
  dm_gam <- gam(gam_formula, data=data, method="REML")
  #gam_list[[i]] <- dm_gam

  # FIGURE SECTION:
  # tiff(paste0(plot_folder,plot_type,plot_nm,".tiff"),
  #      units="in", width=6, height=4, res=200)
  plot <- draw(dm_gam)
  ggsave(plot, file = paste0(plot_folder,plot_type,plot_nm,".tiff"),
         units="in", width=6, height=4, dpi=200)
  #dev.off()
}


### dist prob plots:-----
for (i in 1:nrow(tops)){
  #start the loop here >>>
  # grab model variable column numbers:
  v1 <- as.numeric(tops[i, 1])
  v2 <- as.numeric(tops[i, 2])
  v3 <- as.numeric(tops[i, 3])
  
  # grab variables:
  vars <- cbind(MV_2023_09_DATA_DP[,c(v1,v2,v3)])
  var_nms <- as.character(c(v1,v2,v3))
  plot_nm <- paste0("2023_09_plot_", as.character(i), "_vars_",
                    paste(var_nms, collapse = "_"))
  
  # combine data:
  data <- cbind.data.frame(probs[,2], vars[probs$sitenum,])
  colnames(data) <- c("dpy", "var1", "var2", "var3")
  
  # make the formula:
  ex_vars <- c()
  for (j in 1:ncol(data)-1){
    ex_vars[j] <- paste0('s(', names(data[j+1]), ')')
  }
  
  #make a single string:
  gam_formula <- as.formula(paste("dpy ~",
                                  paste(ex_vars, collapse='+')))
  
  # make the model:
  dm_gam <- gam(gam_formula, data=data, family = "binomial", method="REML")
  #gam_list[[i]] <- dm_gam
  
  # FIGURE SECTION:
  # tiff(paste0(plot_folder,plot_type,plot_nm,".tiff"),
  #      units="in", width=6, height=4, res=200)
  plot <- draw(dm_gam)
  ggsave(plot, file = paste0(plot_folder,plot_type,plot_nm,".tiff"),
         units="in", width=6, height=4, dpi=200)
  #dev.off()
}
