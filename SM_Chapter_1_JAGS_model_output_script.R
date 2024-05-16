# JAGS model output script - for extracting data for maps
# and for making final figures for Chapter 1

### load libraries:
library(tidyverse)
library(dplyr)
library(readr)
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