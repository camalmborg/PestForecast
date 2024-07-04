# JAGS model output script - for extracting data for maps
# and for making final figures for Chapter 1

### load libraries:
library(tidyverse)
library(dplyr)
library(readr)
library(boot)
library(MCMCvis)
library(ecoforecastR)
library(ggplot2)
library(hrbrthemes)
library(pROC)
library(knitr)
library(kableExtra)
library(webshot)


### load environmental data
# dmls object:
load("CHAPTER_1/2024_JAGS_models/2024_05_dmls.RData")
# dpls object:
load("CHAPTER_1/2024_JAGS_models/2024_05_dpls.RData")
# find sites with missing output values:
# choose dmls obj with all environmental data
find_miss <- dmls[[2]]
# find missing values
missing <- as.numeric(rownames(find_miss[!complete.cases(find_miss),]))
rm(find_miss)
# should have 4990 sites total (4997 for multivariate analyses, 4990 for JAGS)

### load TGC and CS dist prob and dist mags data
# tcg dist-mag-recov object:
dmr_tcg <- read.csv("CHAPTER_1/DATA/2023_12_DMR_DATA_TCG.csv")[-missing,]
# cs dist-mag-recov object:
dmr_cs <- read.csv("CHAPTER_1/DATA/2023_12_DMR_DATA_CS.csv")[-missing,]
# Landsat product data:
scores <- read.csv("CHAPTER_1/DATA/2022_12_7_sample_score_mean_5k.csv")[dmr_cs$X,]
# getting geographic data for mapping
# isolate geographic data from Landsat product
geo <- as.data.frame(scores$.geo)
#make lat and lon columns from .geo:
coords<-matrix(nrow=nrow(geo),ncol=3)
for (i in 1:nrow(geo)){
  #longitudes:
  lon<-str_extract(geo[i,], "\\d+\\.*\\d*")
  coords[i,1]<-as.numeric(lon)*-1
  #latitudes:
  extlon<-sub(lon,"",geo[i,])
  coords[i,2]<-as.numeric(str_extract(extlon, "\\d+\\.*\\d*"))
  coords[i,3]<-i
}
# site data
sites <- cbind(id = 1:nrow(scores),
               sys = scores$system.index,
               lon = coords[,1],
               lat = coords[,2])
colnames(coords) <- c("lon", "lat", "site")
rm(geo)
#write.csv(coords, "CHAPTER_1/DATA/site_coordinates.csv")


### ALPHA MODELS --- using output for computing disturbance predictions:
# prepare the predicted values:
# load best alpha model outputs
runpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_runs/"
outpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_outputs/"
# model data
runfile <- "A_best_alpha/2024-05-10_modelrun_alpha_a_10_b_1_data.RData"
#runfile <- "2024-06-21_modelrun_joint_a_1_b_6_data.RData"
run <- load(paste0(runpath, runfile))
# model output
outfile <- "A_best_alpha/2024-05-10_modelrun_alpha_a_10_b_1_output.csv"
#outfile <- "2024-06-21_modelrun_joint_a_1_b_6_output.csv"
out <- read.csv(paste0(outpath, outfile))

# data from model inputs for computing dist prob and dist mag pred/obs
data <- model_info$metadata$data

# covariates
env <- data$a
# predicted intercept and slope means:
preds <- apply(out, 2, mean)
# get probabilities
# which outputs are alpha covariates
ps <- grep("^alpha", names(preds))
# predicted disturbances
Ed <- inv.logit(as.matrix(env) %*% as.matrix(preds[ps]))
# observed disturbances:
# dist16 <- dmr_cs$dpy1
# dist17 <- dmr_cs$dpy2
# make empty matrix to fill with 2016 and 2017 disturbance data
obsdist <- c(rep(NA, nrow(dmr_cs)))
# add 1s for disturbances in each year:
dists <- c(which(dmr_cs$dpy1 == 1), 
           which(dmr_cs$dpy2 == 1))
obsdist[dists] <- 1
obsdist[is.na(obsdist)] <- 0

# percentage of observed disturbances:
oD <- sum(dist16 + dist17) / nrow(dmr_cs) * 100
oD16 <- sum(dist16) / nrow(dmr_cs) * 100
# percentage of predicted disturbance:
pD <- length(which(Ed > 0.95)) / nrow(dmr_cs) * 100

# save data for maps:
# add coordinates
probs_pd_obs <- cbind(coords[,'lon'], coords[,'lat'],
                     obsdist, Ed)
colnames(probs_pd_obs) <- c("lon", "lat", "obs", "pred")
write.csv(probs_pd_obs, "Maps/Chapter_1/Data/probs_pred_obs.csv")

  
### BETA MODELS --- using output for computing disturbance magnitudes:
# load model outputs
runpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_runs/"
outpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_outputs/"

runfile <- "A_best_beta/2024-05-16_modelrun_beta_a_1_b_4_data.RData"
run <- load(paste0(runpath, runfile))
# model output
outfile <- "A_best_beta/2024-05-16_modelrun_beta_a_1_b_4_output.csv"
out <- read.csv(paste0(outpath, outfile))

# data from model inputs for computing mu
data <- model_info$metadata$data

# covariates
env <- data$b
# predicted intercept and slopes means
preds <- apply(out, 2, mean)
# find beta param estimates from outputs
ps <- grep("^beta", names(preds))
# predicted magnitudes
Emu0 <- as.matrix(env) %*% as.matrix(preds[ps])
# observed minimum values from disturbance
obsdist <- dmr_cs$mins
# prior timestep value
prior <- data$x_ic
# predicted score value - magnitude = disturbance
preddist <- prior - Emu0

# save data for maps:
# add coordinates
mags_pd_obs <- cbind(coords[,'lon'], coords[,'lat'],
                     obsdist, preddist)
colnames(mags_pd_obs) <- c("lon", "lat", "obs", "pred")
write.csv(mags_pd_obs, "Maps/Chapter_1/Data/mags_pred_obs.csv")


### JOINT MODELS --- using output from best joint model:
# load model outputs
runpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_runs/"
outpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_outputs/"

runfile <- "A_best_joint/2024-06-21_modelrun_joint_a_4_b_2_data.RData"
run <- load(paste0(runpath, runfile))
# model output
outfile <- "A_best_joint/2024-06-21_modelrun_joint_a_4_b_2_output.csv"
out <- read.csv(paste0(outpath, outfile))

# data from model inputs for computing dist prob and dist mag pred/obs
data <- model_info$metadata$data

# covariates
enva <- data$a
envb <- data$b
# predicted intercept and slope means
preds <- apply(out, 2, mean)
# get probabilities
# which outputs are alpha covariates
psa <- grep("^alpha", names(preds))
psb <- grep("^beta", names(preds))
# predicted disturbances:
Ed <- inv.logit(as.matrix(enva) %*% as.matrix(preds[psa]))
# observed values:
# make empty matrix to fill with 2016 and 2017 disturbance data
obsdist_p <- c(rep(NA, nrow(dmr_cs)))
# add 1s for disturbances in each year:
dists <- c(which(dmr_cs$dpy1 == 1), 
           which(dmr_cs$dpy2 == 1))
obsdist_p[dists] <- 1
obsdist_p[is.na(obsdist_p)] <- 0
# predicted magnitudes:
Emu0 <- as.matrix(envb) %*% as.matrix(preds[psb])
# observed values:
# observed minimum values from disturbance
obsdist_m <- dmr_cs$mins
# prior timestep value
prior <- data$x_ic
# predicted score value - magnitude = disturbance
preddist_m <- prior - Emu0

# save data for maps:
# add coordinates
joint_pd_obs <- cbind(coords[,'lon'], coords[,'lat'],
                      obsdist_p, Ed,
                      obsdist_m, preddist_m)
colnames(joint_pd_obs) <- c("lon", "lat", "a_obs", "a_pred", "b_obs", "b_pred")
write.csv(joint_pd_obs, "Maps/Chapter_1/Data/joint_pred_obs.csv")


### Example time series for conceptual figure ---------------------------
just_scores <- scores %>%
  

##### ARCHIVE ##### -----------------------------------------------------
### PLOTTING AND TABLE STUFF FOR EFI CONFERENCE

# ARCHIVED
#Emu <- 
# muN[s] ~ dnorm(mun[s], pan)
# muD[s] ~ dnorm(mu0[s], pa0) 
# logit(D[s]) <- alpha0
# mu[s] <- D[s] * muD[s] + (1-D[s]) * muN[s]
# mu0[s] <- inprod(beta[], b[s,])
# mun[s] <- R * x[s]
# x[s] ~ dnorm(x_ic[s], tau_ic[s])


# # plot for pred vs obs dist mag
# plot_data = cbind.data.frame(june_2016_cs_m, Emu0)
# plot_data = cbind.data.frame(mags_2016_cs_correct, Emu0)
# #plot_lm <- lm(plot_data$Emu0 ~ plot_data$june_2016_cs_m)
# plot_lm <- lm(plot_data$mags_2016_cs ~ plot_data$Emu0)
# 
# mags_plot <- ggplot(plot_data, aes(x = plot_data[,2], y = plot_data[,1])) +
#   geom_point(color = "grey37") +
#   geom_abline(color = "red", lwd = 1)
#   #geom_smooth(method=lm , color="red", se=FALSE)
# print(mags_plot)
# sqrt(mean(plot_data$Emu0-plot_data$mags_2016_cs, na.rm=T)^2)
# 
# # plot for ROC/AUC
# plot_data_roc = cbind.data.frame(roc_fit_m, Ed)
# plot_roc <- roc(plot_data_roc$roc_fit_m ~ plot_data_roc$Ed)
# plot_auc <- plot_roc$auc
# 
# auc_plot <- ggplot(plot_data_roc, aes(y = plot_data_roc$roc_fit_m,
#                                       x = plot_data_roc$Ed)) +
#   geom_point(color = "grey37") + 
#   #geom_smooth(se = F)
#   abline()
# print(auc_plot)
# 
# # RMSE
# sqrt(mean(plot_data_roc$roc_fit_m-plot_data_roc$Ed)^2)
# 
# 
# ### Tables
# 
# #prep table
# summ = summary(best_alpha_jpout)
# means = formatC(summ$statistics[,"Mean"], 
#                 digits = 3,
#                 format = 'f')
# alpha_table = as.data.frame(matrix(NA, length(means)+1, 1))
# rownames(alpha_table) <- c(names(means)[1],"Int",
#                            names(dpls[[1]][2:4]),
#                            names(means)[6:7],
#                            "DIC")
# colnames(alpha_table) <- "Value"
# alpha_table[,1] <- c(means, formatC(best_alpha_dic[[2]][1],
#                                     digits = 3,
#                                     format = 'f'))
# kbl(alpha_table) %>%
#   kable_styling()
# #save_kable(alpha_table, "Chapter_1/2024_06_alpha_table.png")
# 
# 
# 
# # beta table
# summ = summary(best_beta_jpout)
# means = formatC(summ$statistics[,"Mean"], 
#                 digits = 3,
#                 format = 'f')
# beta_table = as.data.frame(matrix(NA, length(means)+1, 1))
# rownames(beta_table) <- c(names(means)[1:2], "Int",
#                           names(dmls[[1]][2:5]),
#                           names(means[8]), "DIC")
# colnames(beta_table) <- "Value"
# beta_table[,1] <- c(means, formatC(best_beta_dic[[2]][1],
#                                     digits = 3,
#                                     format = 'f'))
# kbl(beta_table) %>%
#   kable_styling()
