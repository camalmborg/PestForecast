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
#runfile <- "2024-06-21_modelrun_joint_a_1_b_6_data.RData"
run <- load(paste0(runpath, runfile))
# model output
outfile <- "A_best_alpha/2024-05-10_modelrun_alpha_a_10_b_1_output.csv"
#outfile <- "2024-06-21_modelrun_joint_a_1_b_6_output.csv"
out <- read.csv(paste0(outpath, outfile))

# data from model inputs for computing mu
data <- model_info$metadata$data
# covariates
enva <- data$a
envb <- data$b

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



### plotting:

# libraries
library(ggplot2)
library(hrbrthemes)
library(pROC)

# plot for pred vs obs dist mag
plot_data = cbind.data.frame(june_2016_cs_m, Emu0)
plot_data = cbind.data.frame(mags_2016_cs_correct, Emu0)
#plot_lm <- lm(plot_data$Emu0 ~ plot_data$june_2016_cs_m)
plot_lm <- lm(plot_data$mags_2016_cs ~ plot_data$Emu0)

mags_plot <- ggplot(plot_data, aes(x = plot_data[,2], y = plot_data[,1])) +
  geom_point(color = "grey37") +
  geom_abline(color = "red", lwd = 1)
  #geom_smooth(method=lm , color="red", se=FALSE)
print(mags_plot)
sqrt(mean(plot_data$Emu0-plot_data$mags_2016_cs, na.rm=T)^2)

# plot for ROC/AUC
plot_data_roc = cbind.data.frame(roc_fit_m, Ed)
plot_roc <- roc(plot_data_roc$roc_fit_m ~ plot_data_roc$Ed)
plot_auc <- plot_roc$auc

auc_plot <- ggplot(plot_data_roc, aes(y = plot_data_roc$roc_fit_m,
                                      x = plot_data_roc$Ed)) +
  geom_point(color = "grey37") + 
  #geom_smooth(se = F)
  abline()
print(auc_plot)

# RMSE
sqrt(mean(plot_data_roc$roc_fit_m-plot_data_roc$Ed)^2)


### Tables
# libraries
library(kableExtra)
library(knitr)
library(webshot)

#prep table
summ = summary(best_alpha_jpout)
means = formatC(summ$statistics[,"Mean"], 
                digits = 3,
                format = 'f')
alpha_table = as.data.frame(matrix(NA, length(means)+1, 1))
rownames(alpha_table) <- c(names(means)[1],"Int",
                           names(dpls[[1]][2:4]),
                           names(means)[6:7],
                           "DIC")
colnames(alpha_table) <- "Value"
alpha_table[,1] <- c(means, formatC(best_alpha_dic[[2]][1],
                                    digits = 3,
                                    format = 'f'))
kbl(alpha_table) %>%
  kable_styling()
#save_kable(alpha_table, "Chapter_1/2024_06_alpha_table.png")



# beta table
summ = summary(best_beta_jpout)
means = formatC(summ$statistics[,"Mean"], 
                digits = 3,
                format = 'f')
beta_table = as.data.frame(matrix(NA, length(means)+1, 1))
rownames(beta_table) <- c(names(means)[1:2], "Int",
                          names(dmls[[1]][2:5]),
                          names(means[8]), "DIC")
colnames(beta_table) <- "Value"
beta_table[,1] <- c(means, formatC(best_beta_dic[[2]][1],
                                    digits = 3,
                                    format = 'f'))
kbl(beta_table) %>%
  kable_styling()
