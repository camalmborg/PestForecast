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

### load environmental data:
# load("CHAPTER_1/2024_JAGS_models/2024_04_dmls.RData")  # 5/16 - need to make BEST MODELS version for next set
# load("CHAPTER_1/2024_JAGS_models/2024_04_dpls.RData")

### ALPHA MODELS --- using output for computing disturbance predictions:
# load model outputs
runpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_runs/"
outpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_outputs/"
# model data
runfile <- "A_best_alpha/2024-05-10_modelrun_alpha_a_10_b_1_data.RData"
run <- load(paste0(runpath, runfile))
# model output
outfile <- "A_best_alpha/2024-05-10_modelrun_alpha_a_10_b_1_output.csv"
out <- read.csv(paste0(outpath, outfile))

# data from model inputs for computing mu
data <- model_info$metadata$data
# covariates
env <- data$a

# predicted intercept and slope means
preds <- apply(out, 2, mean)

#get probabilities
Ed <- inv.logit(as.matrix(env) %*% as.matrix(preds[3:6]))


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
# environmental data
# env <- dmls[[13]][,c(1,2,4,6,7)] ## 5/16 best performer now ID'd - need to make new dpls and dmls objects
env <- data$b

# predicted intercept and slopes means
preds <- apply(out, 2, mean)

#Ex <- dnorm(nrow(env), mean = mean(data$x_ic), sd = mean(data$tau_ic))
#Emun <- dnorm(nrow(env), data$rmean, data$rprec)

Emu0 <- as.matrix(env) %*% as.matrix(preds[4:(length(preds)-1)])

#Emu <- 
# muN[s] ~ dnorm(mun[s], pan)
# muD[s] ~ dnorm(mu0[s], pa0) 
# logit(D[s]) <- alpha0
# mu[s] <- D[s] * muD[s] + (1-D[s]) * muN[s]
# mu0[s] <- inprod(beta[], b[s,])
# mun[s] <- R * x[s]
# x[s] ~ dnorm(x_ic[s], tau_ic[s])

