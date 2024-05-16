# JAGS model output script - for extracting data for maps
# and for making final figures for Chapter 1

### load libraries:
library(tidyverse)
library(dplyr)
library(readr)
library(boot)
library(MCMCvis)
library(ecoforecastR)
library(ggplot)
library(knitr)

### load model output:
runpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_runs/"
outpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_outputs/"

runfile <- "A_best_alpha/2024-05-10_modelrun_alpha_a_10_b_1_data.RData"
run <- load(paste0(runpath, runfile))

outfile <- "A_best_alpha/2024-05-10_modelrun_alpha_a_10_b_1_output.csv"
out <- read.csv(paste0(outpath, outfile))

### load environmental data:
load("CHAPTER_1/2024_JAGS_models/2024_04_dmls.RData")  # 5/16 - need to make BEST MODELS version for next set
load("CHAPTER_1/2024_JAGS_models/2024_04_dpls.RData")

### ALPHA MODELS --- using output for computing disturbance predictions:
# number of samples
niter = nrow(out)
nsite = nrow(env)
# environmental data
env <- dpls[[10]]

# get intercept and slope means
preds <- apply(out, 2, mean)

#get probabilities
Ed <- inv.logit(as.matrix(env[,1:4]) %*% as.matrix(preds[3:6]))


### BETA MODELS --- using output for computing disturbance magnitudes:
# environmental data
env <- dmls[[13]][,c(1,2,4,6,7)] ## 5/16 best performer now ID'd - need to make new dpls and dmls objects
