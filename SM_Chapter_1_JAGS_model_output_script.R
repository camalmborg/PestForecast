# JAGS model output script - for extracting data for maps
# and for making final figures for Chapter 1

### load libraries:
library(tidyverse)
library(dplyr)
library(readr)
library(boot)
library(rjags)
#library(MCMCvis)
#library(ecoforecastR)
library(ggplot2)
library(hrbrthemes)
library(mgcv)
library(pROC)
library(knitr)
library(kableExtra)
library(webshot2)
library(superheat)
library(RColorBrewer)


### load environmental data-----
# dmls object:
load("CHAPTER_1/2024_09_JAGS_models/2024_09_dmls.RData") #making new figures 9/19
# dpls object:
load("CHAPTER_1/2024_09_JAGS_models/2024_09_dpls.RData")
load("CHAPTER_1/2024_09_JAGS_models/2024_09_dpls_2.RData")
# find sites with missing output values:
# choose dmls obj with all environmental data
find_miss <- dmls[[4]] # one with SMAP
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


### ALPHA MODELS using output for computing disturbance predictions: -----
# prepare the predicted values:
# load best alpha model outputs
#runpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_runs/"
#outpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_outputs/"
outpath <- "CHAPTER_1/2024_09_JAGS_models/model_outputs/alpha/"
runpath <- "CHAPTER_1/2024_09_JAGS_models/model_runs/alpha/"
# model data
runfile <- "2024-09-09_modelrun_alpha_a_21_b_1_data.RData"
run <- load(paste0(runpath, runfile))
# model output
outfile <- "2024-09-09_modelrun_alpha_a_21_b_1_output.csv"
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
dist16 <- dmr_cs$dpy1
#dist17 <- dmr_cs$dpy2
# make empty matrix to fill with 2016 and 2017 disturbance data
# obsdist <- c(rep(NA, nrow(dmr_cs)))
# # add 1s for disturbances in each year:
# dists <- c(which(dmr_cs$dpy1 == 1), 
#            which(dmr_cs$dpy2 == 1))
# obsdist[dists] <- 1
# obsdist[is.na(obsdist)] <- 0

# percentage of observed disturbances:
# oD <- sum(dist16 + dist17) / nrow(dmr_cs) * 100
# oD16 <- sum(dist16) / nrow(dmr_cs) * 100
# # percentage of predicted disturbance:
# pD <- length(which(Ed > 0.95)) / nrow(dmr_cs) * 100

# save data for maps:
# add coordinates
probs_pd_obs <- cbind.data.frame(coords[,'lon'], coords[,'lat'],
                     dist16, Ed)
colnames(probs_pd_obs) <- c("lon", "lat", "obs", "pred")
head(probs_pd_obs[probs_pd_obs$obs == 1,])
#write.csv(probs_pd_obs, "Maps/Chapter_1/Data/2024_09_probs_pred_obs_test_a_4.csv")

  
### BETA MODELS --- using output for computing disturbance magnitudes:-----
# load model outputs
# runpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_runs/"
# outpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_outputs/"
runpath <- "CHAPTER_1/2024_09_JAGS_models/model_runs/beta/"
outpath <- "CHAPTER_1/2024_09_JAGS_models/model_outputs/beta/"

runfile <- "2024-09-09_modelrun_beta_a_17_b_1_data.RData"
run <- load(paste0(runpath, runfile))
# model output
outfile <- "2024-09-09_modelrun_beta_a_17_b_1_output.csv"
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
preddist <- prior + Emu0

# save data for maps:
# add coordinates
mags_pd_obs <- cbind(coords[,'lon'], coords[,'lat'],
                     obsdist, preddist)
colnames(mags_pd_obs) <- c("lon", "lat", "obs", "pred")

#mags_map_dat <- cbind(coords[,'lon'], coords[,'lat'], Emu0)
write.csv(mags_pd_obs, "Maps/Chapter_1/Data/2024_09_mags_pred_obs.csv")


### JOINT MODELS --- using output from best joint model:-----
# load model outputs
#runpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_runs/"
#outpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_outputs/"
runpath <- "CHAPTER_1/2024_09_JAGS_models/model_runs/"
outpath <- "CHAPTER_1/2024_09_JAGS_models/model_outputs/"

runfile <- "2024-09-19_modelrun_joint_a_2_b_2_data.RData"
run <- load(paste0(runpath, runfile))
# model output
outfile <- "2024-09-19_modelrun_joint_a_2_b_2_output.csv"
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
#obsdist_p <- c(rep(NA, nrow(dmr_cs)))
# add 1s for disturbances in each year:
obsdist_p <- dmr_cs$dpy1
# dists <- c(which(dmr_cs$dpy1 == 1), 
#            which(dmr_cs$dpy2 == 1))
# obsdist_p[dists] <- 1
# obsdist_p[is.na(obsdist_p)] <- 0
# predicted magnitudes:
Emu0 <- as.matrix(envb) %*% as.matrix(preds[psb])
# observed values:
# observed minimum values from disturbance
obsdist_m <- dmr_cs$mins
# prior timestep value
prior <- data$x_ic
# predicted score value - magnitude = disturbance
preddist_m <- prior + Emu0

# save data for maps:
# add coordinates
joint_pd_obs <- cbind(coords[,'lon'], coords[,'lat'],
                      obsdist_p, Ed,
                      obsdist_m, preddist_m)
colnames(joint_pd_obs) <- c("lon", "lat", "a_obs", "a_pred", "b_obs", "b_pred")
#write.csv(joint_pd_obs, "Maps/Chapter_1/Data/2024_09_joint_pred_obs.csv")
joint_mags_test <- cbind(coords[,'lon'], coords[,'lat'], Emu0)

### Predicted/Observed plots: -----
# DP ------
# get environmental data
var <- env[,-1]
nvars = ncol(var)
# get disturbance probability data
dists<-dmr_cs[,grep("^dp",colnames(dmr_cs))]
yr = 1 # for 2016, 2 for 2017
#make gam data frame:
vardat <- as.data.frame(cbind(dists[,yr], var))

#make gam explantory variables list
ex_vars <- c()
for (j in 1:nvars){
  ex_vars[j] <- paste0('s(vardat[,', j+1, '])')
}
#make a single string:
gam_formula <- as.formula(paste("vardat[,1] ~ ",
                                paste(ex_vars, collapse='+')))
#run gam with those data:
mv_gam <- gam(gam_formula, data=vardat, family="binomial")
mv_roc<-roc(mv_gam$y,mv_gam$fitted.values)
roc_fit_m <- mv_gam$fitted.values

# plot for ROC/AUC
plot_data_roc = cbind.data.frame(roc_fit_m, Ed)
plot_roc <- roc(plot_data_roc$Ed, plot_data_roc$roc_fit_m)
plot_auc <- plot_roc$auc

auc_plot <- ggplot(plot_data_roc, aes(y = plot_data_roc$roc_fit_m,
                                      x = plot_data_roc$Ed)) +
  geom_point(color = "grey37") +
  #geom_smooth(se = F)
  abline()
print(auc_plot)

# RMSE
sqrt(mean(plot_data_roc$roc_fit_m-plot_data_roc$Ed)^2)


# DM: ------
# data
#obsdist <- obsdist_m  # for joint
#preddist <- preddist_m
plot_data <- cbind.data.frame(obsdist, preddist)
pred_obs_lm <- lm(plot_data$obsdist ~ plot_data$preddist, plot_data)
summ = summary(pred_obs_lm)

#png(mags_plot, "2024_07_distmag_pred_v_obs_joint.png",
    #width = 6, height = 4, units = "in", res = 300)
mags_plot <- ggplot(plot_data, aes(x = plot_data[,2], y = plot_data[,1]),
                    xlim(-20,10), ylim(-20,10)) +
  geom_point(color = "grey50", size = 1) +
  geom_abline(color = "firebrick", lwd = 1) +
  labs(x = "Forest Condition (Predicted Score)",
       y = "Forest Condition (Observed Score)",
       title = "Disturbance Magnitude Predicted vs Observed") +
  annotate("text",
           x = 3.5, y = -12,
           #x = 6, y = -12, #for joint
           label = paste("R-squared:", round(summ$adj.r.squared, 3)),
           color = "navyblue", size = 3) +
  theme_classic()
  #geom_smooth(method=lm , color="red", se=FALSE)
mags_plot
#dev.off()
# save plot
ggsave(
  filename = "2024_07_pred_v_obs_mags_beta.png",
  plot = mags_plot,
  device = "png",
  width = 7,
  height = 4,
  units = "in",#c("in", "in"),
  dpi = 300 #"print"
)

#sqrt(mean(plot_data$preddist-plot_data$obsdist, na.rm=T)^2)


### Example time series for conceptual figure ---------------------------
# get the condition scores for time series examples
# using just june months for smoothing
just_scores <- scores %>% 
  # Drop unwanted columns
  dplyr::select(dplyr::starts_with("X")) %>%
  # Select just June averages for best disturbance visual without seasonality
  dplyr::select(dplyr::contains(".06.")) %>%
  # site column for flipping longer
  mutate(sites = as.factor(1:nrow(x = .)), .before = dplyr::everything()) %>%
  # flip
  tidyr:: pivot_longer(cols = starts_with("X")) %>%
  # rename with date, without extra characters
  dplyr::mutate(date_char = str_replace_all(string = name,
                                            c("X|_score_mean|_cs_mean" = "", 
                                            "\\." = "-"))) %>%
  # year and month info
  dplyr:: mutate(
    year = as.numeric(stringr::str_sub(string = date_char, start = 1, end = 4)),
    month = as.numeric(stringr::str_sub(string = date_char, start = 6, end = 7))) %>%
  # Get a real date column
  dplyr::mutate(date = as.Date(date_char)) %>%
  # filter time amount
  filter(date > as.Date("2004-06-01")) %>%
  # select columns
  select(sites, date, value)

# example site number
# running some samples:
eg <- which(dmr_cs$mags > 8 & dmr_cs$mags < 8.3 & dmr_cs$colnum == 22)
#site_num <- 3080   
site_num = eg[16]   
# filter for just site number
eg_site <- just_scores[just_scores$sites == site_num,]


# make a nice looking plot
time_series <- ggplot(data = eg_site, aes(x = date, y = value)) +
  # add background color for disturbance years:
  # geom_rect(data = NULL, aes(xmin = as.Date("2016-05-01"), xmax = as.Date("2020-07-01"),
  #                            ymin = -Inf, ymax = Inf), alpha = 0.1, show.legend = FALSE,
  #           fill = "rosybrown1", color = NA) +
  # add line
  geom_line() +
  # add points
  geom_point(pch = 16, size = 1.5, color = "black",
             alpha = 0.8, position = position_dodge(width = 0.2)) +
  # add dashes connecting missing bits
  geom_line(data = filter(eg_site, !is.na(value)),
            #linetype = "dashed", 
            #linewidth = 0.3, 
            show.legend = FALSE,
            position = position_dodge(width = 0.2)) +
  # axis titles
  xlab("Date") +
  ylab("Condition Score") +
  # change plot theme
  theme(
    # remove grid
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    # remove background
    panel.background = element_blank(), 
    # remove sides
    axis.line = element_line(colour = "black"),
    # axis title x and y
    axis.title = element_text(color = "black")
  )
print(time_series)
# save plot
ggsave(
  filename = "2024_07_25_forest_condition_eg_plot_no_box.png",
  plot = time_series,
  device = "png",
  width = 6,
  height = 2.5,
  units = "in",#c("in", "in"),
  dpi = 300 #"print"
)
  
##### TABLES ##### ------------------------------------------------------
### ALPHA best model PARAMS TABLE -----
# #prep table
# identify the model:
meta <- model_info$metadata
alph = 1  # for alpha model
# covariate names:
cov_names_a <- names(dpls[[alph]])
# summary - param means:
summ = summary(model_info$jpout)
means = formatC(summ$statistics[,"Mean"], 
                digits = 3,
                format = 'f')
# make alpha param table:
alpha_means <- as.data.frame(matrix(NA, length(cov_names_a) + 2, 1))
rownames(alpha_means) <- c("Dist Prob Intercept",
                           cov_names_a[-1],
                           "Dist Mag Intercept",
                           "Process Error")
colnames(alpha_means) <- "Disturbance Probability Model Parameter Means"

# compile final results:
alpha_means[,1] <- c(means[grep("^alpha", names(means))],
                     means[grep("^beta", names(means))],
                     means[grep("^p", names(means))])

# make them pretty:
kbl(alpha_means) %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = T, 
                position = "left")

### BETA best model PARAMS TABLE -----
# identify the model:
meta <- model_info$metadata
bet = 1  # for beta model
# covariate names:
cov_names_b <- names(dmls[[bet]])
# summary - param means:
summ = summary(model_info$jpout)
means = formatC(summ$statistics[,"Mean"], 
                digits = 3,
                format = 'f')
# make beta param table:
beta_means <- as.data.frame(matrix(NA, length(cov_names_b) + 2, 1))
rownames(beta_means) <- c("Dist Mag Intercept",
                           cov_names_b[-1],
                           "Dist Prob Intercept",
                           "Process Error")
colnames(beta_means) <- "Disturbance Magnitude Model Parameter Means"

# compile final results:
beta_means[,1] <- c(means[grep("^beta", names(means))],
                     means[grep("^alpha", names(means))],
                     means[grep("^p", names(means))])

# make them pretty:
kbl(beta_means) %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = T, 
                position = "left")


### JOINT best model PARAMS TABLE -----
# identify the model:
meta <- model_info$metadata
# get covariates
alph <- meta$modelrun_a
bet <- meta$modelrun_b
# covariate names:
cov_names_a <- names(dpls[[alph]])
cov_names_b <- names(dmls[[bet]]) 
# change one row name to not be duplicate for table
nm <- intersect(cov_names_a, cov_names_b)
cov_names_b[which(cov_names_b == nm[2])] <- paste0(nm[2],"_")
# model summary
summ = summary(model_info$jpout)
means = formatC(summ$statistics[,"Mean"], 
                digits = 3,
                format = 'f')
# make param table:
joint_means <- as.data.frame(matrix(NA, length(cov_names_a) + length(cov_names_b) + 1, 1))
rownames(joint_means) <- c("Dist Prob Intercept",
                          cov_names_a[-1],
                          "Dist Mag Intercept",
                          cov_names_b[-1],
                          "Process Error")
colnames(joint_means) <- "Joint Model Parameter Means"

# compile final results:
joint_means[,1] <- c(means[grep("^alpha", names(means))],
                    means[grep("^beta", names(means))],
                    means[grep("^p", names(means))])
# make table:
kbl(joint_means) %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = T, 
                position = "left")



### MULTIVAR ANALYSES tables ------
# compiling a tables of Rs and AUCs and delAICs:
# load results:
alpha_mv_results <- read.csv("CHAPTER_1/2024_01_RESULTS/best_distprob16_models_tcg.csv")
alpha_mv_results$dAIC <- min(alpha_mv_results$aics) - alpha_mv_results$aics
beta_mv_results <- read.csv("CHAPTER_1/2024_01_RESULTS/best_distmag_models_tcg.csv")
beta_mv_results$dAIC <- min(beta_mv_results$aics) - beta_mv_results$aics
# order by delAIC
alpha_mv_results <- alpha_mv_results[order(alpha_mv_results$dAIC, decreasing = TRUE),]
beta_mv_results <- beta_mv_results[order(beta_mv_results$dAIC, decreasing = TRUE),]
# covariate names:
all_covs_a <- names(dpls[[2]][-1])
all_covs_b <- names(dmls[[2]][-1])
# set up tables:
prob_mv_table <- as.data.frame(matrix(NA, 5, length(all_covs_a)+2))
mags_mv_table <- as.data.frame(matrix(NA, 5, length(all_covs_b)+2))

# fill in tables:
# ALPHA
rownames(prob_mv_table) <- as.character(c(1:5))
colnames(prob_mv_table) <- c( paste0("\u394", "AIC"),
                              "AUC",
                              all_covs_a)
prob_mv_table[,1] <- alpha_mv_results$dAIC[1:5]
prob_mv_table[,2] <- alpha_mv_results$aucs[1:5]
for (i in 1:5){
  # extract data
  vars <- alpha_mv_results$model_vars[i]
  # convert to numeric
  vn <- str_split(vars, ',') %>% unlist %>%
    str_replace('\\D*(\\d*)\\D*', '\\1') %>%
    as.numeric()
  # fill in X's for params that were used in that model
  for (j in 1:length(vn)){
    prob_mv_table[i,vn[j] + 2] <- "X"
  }
}

# BETA
rownames(mags_mv_table) <- as.character(c(1:5))
colnames(mags_mv_table) <- c("R sq", 
                             paste0("\u394", "AIC"),
                             all_covs_b)

mags_mv_table[,1] <- beta_mv_results$dAIC[1:5]
mags_mv_table[,2] <- beta_mv_results$r2s[1:5]
for (i in 1:5){
  # extract data
  vars <- beta_mv_results$model_vars[i]
  # convert to numeric
  vn <- str_split(vars, ',') %>% unlist %>%
    str_replace('\\D*(\\d*)\\D*', '\\1') %>%
    as.numeric()
  # fill in X's for params that were used in that model
  for (j in 1:length(vn)){
    mags_mv_table[i,vn[j] + 2] <- "X"
  }
}




### JAGS MODEL ANALYSES tables -----
# # load results:
#JAGS_models <- read.csv("CHAPTER_1/2024_03_JAGS_models/2024_BEST_JAGS_MODELS_RESULTS.csv")
# JAGS_models <- JAGS_models %>% mutate_if(is.numeric, round, digits = 3)
# JAGS_models$dDIC_all <- min(JAGS_models$DIC) - JAGS_models$DIC
# # sort by delDICs:
# JAGS_models <- JAGS_models[order(JAGS_models$dDIC_all, decreasing = T),]
# 
# ## heatmap figure
# # organize data into numeric-only columns:
# JAGS_params <- cbind(a1_int = JAGS_models$a1_int,
#                      JAGS_models[,grep("^a_", colnames(JAGS_models))],
#                      b1_int = JAGS_models$b1_int,
#                      JAGS_models[,grep("^b_", colnames(JAGS_models))])
# # add row names and nicer-looking column names
# rownames(JAGS_params) <- c(paste0(JAGS_models$model_type," ", JAGS_models$model_rank))
# JAGS_params <- JAGS_params %>%
#   rename('Alpha Intercept' = a1_int) %>%
#   rename('Beta Intercept' = b1_int)
# colnames(JAGS_params) <- c(gsub("\\.", replacement = " ", colnames(JAGS_params)))
# colnames(JAGS_params) <- c(gsub("a_", "", colnames(JAGS_params)))
# colnames(JAGS_params) <- c(gsub("b_", "", colnames(JAGS_params)))

# 10/7/2024 BEST JAGS MODELS 
#setwd("/projectnb/dietzelab/malmborg/")

mnames <- c(
  "2024-10-02_modelrun_joint_a_21_b_14",
  "2024-10-04_modelrun_joint_a_18_b_14",
  "2024-10-04_modelrun_joint_a_20_b_17",
  "2024-10-04_modelrun_joint_a_20_b_15",
  "2024-10-05_modelrun_joint_a_20_b_14",
  "2024-09-25_modelrun_alpha_a_10_b_1",
  "2024-09-26_modelrun_alpha_a_18_b_1",
  "2024-09-27_modelrun_alpha_a_20_b_1",
  "2024-09-27_modelrun_alpha_a_21_b_1",
  "2024-09-11_modelrun_beta_a_1_b_13",
  "2024-09-11_modelrun_beta_a_1_b_14",
  "2024-09-12_modelrun_beta_a_1_b_15",
  "2024-09-13_modelrun_beta_a_1_b_16",
  "2024-09-13_modelrun_beta_a_1_b_17")

columns <- c("R", "alpha0", "alpha[1]", "alpha[2]", "alpha[3]", "alpha[4]",
             "alpha[5]", "alpha[6]", "alpha[7]", "alpha[8]", "beta0", 
             "beta[1]", "beta[2]", "beta[3]","beta[4]",
             "beta[5]", "beta[6]", "beta[7]", "beta[8]", 
             "beta[9]", "pa0", "DIC")

best_models <- matrix(data = NA, nrow = length(mnames), ncol = length(columns))
a_vars <- matrix(data = NA, nrow = length(mnames), ncol = length(grep("^a", columns))-1)
b_vars <- matrix(data = NA, nrow = length(mnames), ncol = length(grep("^b", columns))-1)
colnames(best_models) <- columns
colnames(a_vars) <- columns[grep("^a", columns)][-1]
colnames(b_vars) <- columns[grep("^b", columns)][-1]

for (i in 1:length(mnames)){
  # load model
  model <- (paste0("CHAPTER_1/2024_09_JAGS_models/best_models/",
                   mnames[i],"_data.RData"))
  load(model)
  # remove burn in
  jpout <- model_info$jpout
  burnin <- 100000
  jburn <- window(jpout, start = burnin)
  # extract means
  summ <- summary(jburn)
  # fill in DIC
  best_models[i,"DIC"] <- model_info$dic[[2]]
  # put names in data frame
  for (j in names(summ$statistics[,"Mean"])){
    best_models[i, j] <- summ$statistics[j,"Mean"]
  }
  # get variable names
  if (model_info$metadata$covs == 1){
    for (k in 1:length(names(model_info$metadata$data$a))){
      a_vars[i, k] <- names(model_info$metadata$data$a)[k]
      } 
    } else if (model_info$metadata$covs == 2) {
      for (k in 1:length(names(model_info$metadata$data$b))){
        b_vars[i, k] <- names(model_info$metadata$data$b)[k]
        }
      } else if (model_info$metadata$covs == 3){
        for (k in 1:length(names(model_info$metadata$data$a))){
          a_vars[i, k] <- names(model_info$metadata$data$a)[k]
          }
        for (k in 1:length(names(model_info$metadata$data$b))){
          b_vars[i, k] <- names(model_info$metadata$data$b)[k]
    }
  }
}
rm(jburn, jpout, model_info, summ, burnin, i, j, k, model)

# round results to 3 digits
best_models <- as.data.frame(best_models) %>% mutate_if(is.numeric, round, digits = 3)
# calculate delDIC
best_models$delDIC <- min(best_models$DIC) - best_models$DIC
# sort
best_models <- best_models[order(best_models$delDIC, decreasing = T),]
# make new alpha and beta intercept columns
best_models[which(is.na(best_models$alpha0)),'alpha0'] <- c(best_models$`alpha[1]`[which(!is.na(best_models$`alpha[1]`))])
best_models[which(is.na(best_models$beta0)), 'beta0'] <- c(best_models$`beta[1]`[which(!is.na(best_models$`beta[1]`))])

# matching models to correct variable parameters

# put covariate names in one table
a_vars[which(a_vars[,1] == "int"),1] <- "a_int"
b_vars[which(b_vars[,1] == "int"),1] <- "b_int"
vars <- as.data.frame(cbind(a_vars, b_vars))
# make matrix to fill
model_params <- as.data.frame(matrix(data = NA, nrow = length(mnames), ncol = ncol(vars)))
colnames(model_params) <- c(as.character(vars[which(!is.na(vars["alpha[8]"]))[1], grep("^a", names(vars))]),
                            as.character(vars[which(!is.na(vars["beta[9]"]))[1], grep("^b", names(vars))]))
# fill in matrix values
for (i in 1:nrow(model_params)){
  for(j in 1:length(names(vars))){
    var <- vars[i,j]
    name <- names(vars)[j]
    model_params[i,var] <- best_models
  }
}
# make final table
#JAGS_params <- best_models[,c()]

# make values for NAs to be able to print numbers in table
JAGS_params[is.na(JAGS_params)] <- 1000
# make this dark grey:
# all NAs/1000s
JAGS_params.col <- JAGS_params > 999
# set all values that satisfy the condition to color
JAGS_params.col <- gsub("TRUE", "grey45", JAGS_params.col)
# set all values that do not satisfy the condition to "black"
JAGS_params.col <- gsub("FALSE", "black", JAGS_params.col)
# convert to matrix
JAGS_params.col <- matrix(JAGS_params.col, ncol = ncol(JAGS_params))
# add row and column names to make broken-up figures:
#rownames(JAGS_params.col) <- rownames(JAGS_params)
colnames(JAGS_params.col) <- colnames(JAGS_params)

# make heatmap:
png("2024_10_07_JAGS_params_heatmap.png",
    width = 10, height = 8, units = "in", res = 300)
superheat(as.matrix(t(JAGS_params)), 
          scale = FALSE, # the scale normalizes to mean 0 SD 1
          left.label.text.size = 4,
          bottom.label.text.size = 4,
          bottom.label.text.angle = 90,
          heat.pal = c("dodgerblue", "skyblue1", "white", "pink1","firebrick1"),
          heat.pal.values = c(0, 0.45, 0.5, 0.55, 1),
          heat.lim = c(-6,6),
          legend.num.ticks = 8,
          column.title = "Models",
          row.title = "Parameters",
          column.title.size = 4,
          row.title.size = 4,
          X.text = as.matrix(t(JAGS_params)),
          X.text.col = t(JAGS_params.col),
          heat.na.col = "grey45",
          X.text.size = 3,
          extreme.values.na = TRUE,
          left.label.size = 0.3,
          bottom.label.size = 0.15,
          legend.vspace = 0.01)
dev.off()

# three superheats for 3-part figure:
# ALPHA
png("2024_JAGS_alpha_heatmap1.png",
    width = 6.5, height = 3.5, units = "in", res = 300)
superheat(as.matrix(t(JAGS_alpha)), 
          scale = FALSE, # the scale normalizes to mean 0 SD 1
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          #bottom.label.text.angle = 90,
          heat.pal = c("dodgerblue", "skyblue1", "white", "pink1","firebrick1"),
          heat.pal.values = c(0, 0.45, 0.5, 0.55, 1),
          heat.lim = c(-6,6),
          legend.num.ticks = 8,
          column.title = "Alpha Models",
          row.title = "Parameters",
          column.title.size = 4,
          row.title.size = 4,
          X.text = as.matrix(t(JAGS_alpha)),
          X.text.col = t(JAGS_alpha.col),
          heat.na.col = "grey45",
          X.text.size = 3,
          extreme.values.na = TRUE,
          left.label.size = 0.3,
          bottom.label.size = 0.15,
          legend = FALSE)
dev.off()
# BETA
png("2024_JAGS_beta_heatmap1.png",
    width = 7, height = 4, units = "in", res = 300)
superheat(as.matrix(t(JAGS_beta)), 
          scale = FALSE, # the scale normalizes to mean 0 SD 1
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          #bottom.label.text.angle = 90,
          heat.pal = c("dodgerblue", "skyblue1", "white", "pink1","firebrick1"),
          heat.pal.values = c(0, 0.45, 0.5, 0.55, 1),
          heat.lim = c(-6,6),
          legend.num.ticks = 8,
          column.title = "Beta Models",
          row.title = "Parameters",
          column.title.size = 4,
          row.title.size = 4,
          X.text = as.matrix(t(JAGS_beta)),
          X.text.col = t(JAGS_beta.col),
          heat.na.col = "grey45",
          X.text.size = 3,
          extreme.values.na = TRUE,
          left.label.size = 0.3,
          bottom.label.size = 0.15,
          legend = FALSE)
dev.off()
# JOINT
png("2024_JAGS_joint_heatmap1.png",
    width = 8, height = 8, units = "in", res = 300)
superheat(as.matrix(t(JAGS_joint)), 
          scale = FALSE, # the scale normalizes to mean 0 SD 1
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          #bottom.label.text.angle = 90,
          heat.pal = c("dodgerblue", "skyblue1", "white", "pink1","firebrick1"),
          heat.pal.values = c(0, 0.45, 0.5, 0.55, 1),
          heat.lim = c(-6,6),
          legend.num.ticks = 8,
          column.title = "Joint Models",
          row.title = "Parameters",
          column.title.size = 4,
          row.title.size = 4,
          X.text = as.matrix(t(JAGS_joint)),
          X.text.col = t(JAGS_joint.col),
          heat.na.col = "grey45",
          X.text.size = 3,
          extreme.values.na = TRUE,
          left.label.size = 0.3,
          bottom.label.size = 0.1,
          legend.vspace = 0.01)
dev.off()

## make the table:
# # columns with params
# param_cols <- names(JAGS_models[,7:21])
# # assemble table
# best_jags <- kbl(JAGS_models)%>%
#   kable_styling(font_size = 8,
#                 full_width = F)
# # column colors
# for (i in 7:21){
#   best_jags <- column_spec(best_jags, i,
#                            color = "white",
#                            background = spec_color(JAGS_models[,i]))
# }
# 
# best_jags




#save_kable(best_jags, file = "best_jags.png")
  
  
# kbl(JAGS_models)%>%
  # # for every param column
  # reduce(
  #   which(names(JAGS_models) %in% param_cols),
  #   ~ column_spec(.x, .y,
  #                 color = "white",
  #                 background = spec_color())
  # )


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
