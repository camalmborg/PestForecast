# a scripte for viewing MCMC outputs and diagnostics

### THE LIBRARIES:
library(tidyverse)
library(dplyr)
library(rjags)
library(coda)
library(MCMCvis)
library(ecoforecastR)


### FROM SCC FILES:
# load model outputs
runpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_runs/"
outpath <- "CHAPTER_1/2024_JAGS_models/Best_Models_from_SCC/model_outputs/"
# model data
#runfile <- "A_best_alpha/2024-05-12_modelrun_alpha_a_2_b_1_data.RData"
#runfile <- "A_best_beta/2024-05-06_modelrun_beta_a_1_b_15_data.RData"
runfile <- "2024-06-25_modelrun_joint_a_4_b_6_data.RData"
run <- load(paste0(runpath, runfile))
# model output
#outfile <- "A_best_alpha/2024-05-16_modelrun_alpha_a_1_b_1_output.csv"
#outfile <- "A_best_beta/2024-05-16_modelrun_beta_a_1_b_4_output.csv"
#outfile <- "2024-06-19_modelrun_joint_a_1_b_6_output.csv"
#out <- read.csv(paste0(outpath, outfile))

# data from model inputs for computing mu
data <- model_info$metadata$data
# covariates
enva <- data$a
envb <- data$b
# DIC
dic <- model_info$dic
# JAGS output 
jpout <- model_info$jpout

### Let's do some diagnostics:
# plot trace and density plots
plot(jpout)

# GBR 
gelman.diag(jpout)

# GBR plot
gelman.plot(jpout)

# discarding burn in
# set burn in based on GBR
burnin <- 100000
# remove burn in
jburn <- window(jpout, start = burnin)
# plot burn in output
plot(jburn)

# autocorrelation plots
#acfplot(jburn)  # ask Mike about this one

# effective sample size
effectiveSize(jburn)

# summary
summary(jburn)

# correlation plots
# convert MCMC to matrix
out <- as.matrix(jburn)
# pairs plot
pairs(out)
# correlations summary
cor(out)


### Model selection step:
# DIC
# DIC <- dic.samples(j.pests, n,iter = 10000)
# sum <- sum(DIC$deviance, DIC$penalty)



### JAGS Output Figure Code ###
# For when it is time to make figures
# getting proper names for mcmc outputs
##' @param w mcmc object containing matrix outputs
##' @param pre prefix (variable name) for the matrix variable to be extracted
##' @param numeric boolean, whether to coerce class to numeric
parse.MatrixNames <- function(w, pre = "x", numeric = FALSE) {
  w <- sub(pre, "", w)
  w <- sub("[", "", w, fixed = TRUE)
  w <- sub("]", "", w, fixed = TRUE)
  w <- matrix(unlist(strsplit(w, ",")), nrow = length(w), byrow = TRUE)
  if (numeric) {
    class(w) <- "numeric"
  }
  colnames(w) <- c("row", "col")
  return(as.data.frame(w))
}

### ARCHIVE::::
# ### OUT object:
# # model name
# jagmod = j.pests
# # iterations and thinning
# niter = 100000
# jthin = 1
# # variable names object
# jagvars <- c("beta0", "alpha0",
#              "beta[1]", "alpha[1]",
#              "R",
#              "pa0")
# 
# # run from model in SM_JAGS_Model_script
# jpout<-coda.samples(jagmod,
#                     variable.names = jagvars,
#                     n.iter = niter,
#                     thin = jthin)