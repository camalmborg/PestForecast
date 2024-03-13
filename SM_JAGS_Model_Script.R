#This script is the JAGS model for pest disturbance mag and prob forecast

### THE LIBRARIES:
library(tidyverse)
library(dplyr)
library(rjags)
library(coda)
#library(MCMCvis)
#library(ecoforecastR)


### THE DATA:
# condition scores inputs
#cfile <- "2023_03_08_DATAGRAB/2023_03_08_5000_sites_sample_score_mean.csv" #2005-2022
cfile <- "2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv" #1995-2020, original grab
scores <- read.csv(cfile)
# isolate data we want
cs <- scores %>%
  # Drop unwanted columns
  dplyr::select(dplyr::starts_with("X")) %>%
  # Select just June averages for best disturbance visual without seasonality
  dplyr::select(dplyr::contains(".06.")) %>%
  # rename with date, without extra characters
  dplyr::rename_with(~ str_replace_all(., c("X|_score_mean|_cs_mean" = "", 
                                       "\\." = "-")))

# collect previous timesteps for filling in x[s]s
# isolate data we want
prevtime <- scores %>%
  # Drop unwanted columns
  dplyr::select(dplyr::starts_with("X")) %>%
  # Select just June averages for best disturbance visual without seasonality
  dplyr::select(dplyr::contains(".05.")) %>%
  # rename with date, without extra characters
  dplyr::rename_with(~ str_replace_all(., c("X|_score_mean|_cs_mean" = "", 
                                            "\\." = "-")))

# # number of sites
# nsites = nrow(cs)

# disturbance year data
dmr_file <- "CHAPTER_1/DATA/2023_12_DMR_DATA_TCG.csv"
dmr <- read.csv(dmr_file)

# choose disturbance onset year
distyear <- 2016

# disturbance data for sample
cs <- cs[dmr$X,]
# get disturbance years, previous timestep, and condition score value at dist year for each:
dists <- vector()
cs_dists <- vector()
ptmestp <- vector()
ptime <- vector()
for (i in 1:nrow(cs)){
  # get disturbance year colname
  if (dmr$dpy1[i] > 0) {
    dists[i] <- colnames(cs[grep(paste0("^", distyear, sep = ""), names(cs))])
    ptmestp[i] <- colnames(prevtime[grep(paste0("^", distyear, sep = ""), names(prevtime))])
  } else {
    dists[i] <- colnames(cs[grep(paste0("^",distyear+1, sep = ""), names(cs))])
    ptmestp[i] <- colnames(prevtime[grep(paste0("^", distyear+1, sep = ""), names(prevtime))])
  }
  # get condition score associated with disturbance year
  cs_dists[i] <- cs[i, dists[i]]
  # get previous time step score for correct disturbance year
  ptime[i] <- prevtime[i, ptmestp[i]]
}

### standard deviations for precisions
sdfile <- "2022_08_31_DATAGRAB/2022_12_7_sample_score_stddev_5k.csv"
stan_devs <- read.csv(sdfile)
sds <- stan_devs %>%
  # Drop unwanted columns
  dplyr::select(dplyr::starts_with("X")) %>%
  # Select just June averages for best disturbance visual without seasonality
  dplyr::select(dplyr::contains(".06.")) %>%
  # get rid of the other stuff
  dplyr::rename_with(~ str_replace_all(., c("X|_score_stddev" = "",
                                            "\\." = "-")))
# previous year sds for initial conditions
prevsds <- stan_devs %>%
  # Drop unwanted columns
  dplyr::select(dplyr::starts_with("X")) %>%
  # Select just June averages for best disturbance visual without seasonality
  dplyr::select(dplyr::contains(".05.")) %>%
  # get rid of the other stuff
  dplyr::rename_with(~ str_replace_all(., c("X|_score_stddev" = "",
                                            "\\." = "-")))

sds <- sds[dmr$X,]
prevsds <- prevsds[dmr$X,]

### convert to precisions: 
# use disturbance years to get precisions at dist year for each
cs_sds <- vector()
for (i in 1:nrow(sds)){
  # get condition score associated with disturbance year
  cs_sds[i] <- sds[i, dists[i]]
}
# fill in NA values before precisions calculation
# find missing June SD sites
miss_SD <- which(is.na(cs_sds))
# find missing dist years for SD NA's
miss_dist <- dists[which(is.na(cs_sds))]
# take 75%ile for SDs across sites
sd_75s <- c(quantile(cs_sds[which(dists == unique(dists)[1])], 0.75, na.rm = T),
            quantile(cs_sds[which(dists == unique(dists)[2])], 0.75, na.rm = T))
# fill in NAs based on disturbance year
for (i in 1:length(miss_SD)){
  if(miss_dist[i] == unique(dists)[1]){
    cs_sds[miss_SD[i]] <- sd_75s[1]
  } else {
    cs_sds[miss_SD[i]] <- sd_75s[2]
  }
}

# do the same for prev year
# use disturbance years to get precisions at dist year for each
prev_sds <- vector()
for (i in 1:nrow(prevsds)){
  # get condition score associated with disturbance year
  prev_sds[i] <- prevsds[i, ptmestp[i]]
}
# fill in NA values before precisions calculation
# find missing June SD sites
miss_pSD <- which(is.na(prev_sds))
# find missing dist years for SD NA's
miss_prev <- dists[which(is.na(prev_sds))]
# take 75%ile for SDs across sites
prevsd_75s <- c(quantile(prev_sds[which(ptmestp == unique(ptmestp)[1])], 0.75, na.rm = T),
                quantile(prev_sds[which(ptmestp == unique(ptmestp)[2])], 0.75, na.rm = T))
# fill in NAs based on disturbance year
for (i in 1:length(miss_pSD)){
  if(miss_prev[i] == unique(ptmestp)[1]){
    cs_sds[miss_prev[i]] <- prevsd_75s[1]
  } else {
    cs_sds[miss_prev[i]] <- prevsd_75s[2]
  }
}

# convert SDs to precisions
cs_precs <- vector()
for (i in 1:length(cs_sds)){
  cs_precs[i] <- 1/(cs_sds[i]^2)
}

prev_precs <- vector()
for (i in 1:length(prev_sds)){
  prev_precs[i] <- 1/(prev_sds[i]^2)
}


### Setting an informative R prior
#yrs prior for before dist time series
yrs <- 15
# time series 1 object
ts1 <- matrix(data = NA, nrow = nrow(cs), ncol = yrs)
# time series 2 object
ts2 <- matrix(data = NA, nrow = nrow(cs), ncol = yrs)
# make cors vector
cors <- vector()
# make vectors of "before" and "after"
for (i in 1:nrow(cs)){
  # pick out disturbance date for each site
  td <- which(colnames(cs) == dists[i])
  # make vector of each row values:
  v <- cs[i, (td-yrs):(td-1)]
  v2 <- cs[i, (td-yrs+1):td]
  # fill in matrices of values for t-1 (ts1) and t (ts2)
  for (j in 1:ncol(v)){
    ts1[i,j] <- v[,j]
    ts2[i,j] <- v2[,j]
  }
  #correlation between two time series at each site:
  cors[i] <- cor(ts1[i,], ts2[i,], use = "complete.obs")
}
# mean R:
Rmean <- mean(cors, na.rm = T)
# var R:
Rvar <- var(cors, na.rm = T)
Rprec <- 1/(sd(cors, na.rm = T)^2)


# # covariate data - for parameterizations
# # the environmental variables
# # for disturbance magnitude parameters
# varfile_m <- "CHAPTER_1/DATA/MV_2023_12_DATA_distmag.csv"
# magvars <- read.csv(varfile_m)[,-1]
# # for disturbance probability parameters
# varfile_p <- "CHAPTER_1/DATA/MV_2023_12_DATA_distprob.csv"
# probvars <- read.csv(varfile_p)[,-1]
# # environmental data without missing sites
# magvars <- magvars[dmr$X,]
# probvars <- probvars[dmr$X,]

# load disturbance magnitude and disturbance probability covariates
# # disturbance magnitude
# load("Chapter_1/2024_02_JAGS_modls/magsls.RData")
# # disturbance probability
# load("Chapter_1/2024_02_JAGS_models/probls.RData")

# # making covariate data with anomalies rather than raw
# # function:
# anomfx<-function(x){
#   # get mean values for each column
#   means <- apply(x, 2, mean, na.rm=T)
#   # make anomaly matrix
#   anom <- matrix(NA,nrow=nrow(x), ncol=ncol(x))
#   # for each row, fill in covariate anomaly
#   for (i in 1:nrow(x)){
#     for (j in 1:ncol(x)){
#       anom[i,j] <- x[i,j] - means[j] # fill in anomalies for beta[] terms
#     }
#     anom[i,1] <- 1 # make first column 1s for beta0 term
#   }
#   return(anom)
# }
# 
# # converting covariate lists to anomaly
# # new empty list to populate with anomaly versions
# anomls <- list()
# # list being converted CHOOSE magls for disturbance magnitude/probls for disturbance probability
# #covls <- magls
# covls <- probls
# # loop for conversion
# for (i in 1:length(covls)){
#   # run cov members through the anomaly machine
#   anomls[[i]] <- anomfx(covls[[i]])
# }
# # disturbance magnitude saving anomaly version
# #dmls <- anomls
# # disturbance probability saving anomaly version
# dpls <- anomls

# Loading anomaly versions of covariate data 
# disturbance magnitude
load("Chapter_1/2024_02_JAGS_models/dmls.RData")
# disturbance probability
load("Chapter_1/2024_02_JAGS_models/dpls.RData")


### THE MODEL:
#use the single time step version of the model:
spongy_disturb <- "model{

###Loop over individual sites
for (s in 1:ns){

  #### Data Model:
  y[s] ~ dnorm(mu[s], tau_obs[s])  ##3/1 - when adding precisions >> tau_obs[s]
  
  #### Process Model:
  muN[s] <- R * x[s]                    ##step 3: dealing with modeling R (Chap 2 - RECOV)
  #x[s] ~ dnorm(mu[s], tau_add)
  muD[s] ~ dnorm(mu0[s], pa0)           ##step 1: process model on mu0 (MAG)
    
  ##D[s] ~ dbern(p)                     ##step 2: adding process model here (PROB)
  logit(D[s]) <- alpha0 + (alpha[1]*z[s])
    
  mu[s] <- D[s] * muD[s] + (1-D[s]) * muN[s]
  mu0[s] <- beta0 + (beta[1] * v[s]) ##+ (beta[2] * va[s]) + (beta[3] * vb[s])

  ##x[s] ~ dnorm(x_ic, tau_ic) ##3/1 - this is where you would put previous timepoint
  
}#end loop over sites
  
  #### Priors
  #tau_obs ~ dgamma(t_obs, a_obs)     ##observation error (data model)
  #tau_add ~ dgamma(a_add ,t_add)     ##process error (process model)
  R ~ dnorm(rmean, rprec)             ##rho paramter (recovery rate)
  #p ~ dunif(0,1)                     ##disturbance probability
  beta0 ~ dnorm(-5,1)                 ##param for calculating mean of disturbed state
  alpha0 ~ dnorm(0, 0.0001) 
  pa0 ~ dgamma(1,1) #precision of disturbed state
  
  ## COVARIATES WILL BE ADDED HERE
  beta[1] ~ dnorm(0, 0.0001)  
  #beta[2] ~ dnorm(0, 0.0001)
  #beta[3] ~ dnorm(0, 0.0001)
  
  alpha[1] ~ dnorm(0.0, 0.0001)
  
  ## covariate matrix: ADD WHEN READY FOR MODEL RUNS - 3/6/24
  ##beta ~ dmnorm(b0, Vb)   ## for disturbance magnitude
  ##alpha ~ dmnorm(a0, Va)  ## for disturbance probability
  
}
"


### SELECT SITES:
# random selection of sites for testing (before using full sample)
smpl <- sample(nrow(cs), 50)
# make sample
cs_samp <- cs[smpl,]
# number of sites of sample
nsites = nrow(cs_samp)
# get dist years for cs_samp group
dist_samp <- dists[as.numeric(rownames(cs_samp))]
# make the single timestep data for each site
cs_samp_dist <- cs_dists[as.numeric(rownames(cs_samp))]

# make same sample of covariate data for testing individual beta[] parameters 2/27/24
# disturbance magnitude
dmbeta <- dmls[[1]][smpl,2]
# disturbance probability
dpalpha <- dpls[[1]][smpl,2]

# precisions samples
cs_prec_samp <- cs_precs[smpl]


### initial state of model parameters:
# init<-list()
# nchain <- 3
# for(j in 1:nchain){
#   samp<- sample(!is.na(cs_samp),length(cs_samp),replace=TRUE)
#   init[[j]]<-list(tau_add=1/var(diff(samp)),tau_obs=1/var(samp))
# }


### MODEL INPUTS
# data and parameters for sites model:   #for full sample ns = nrow(scores)
data = list(y = cs_samp_dist, ns = nsites,    
              #x_ic = 0, tau_ic = 0.1,
              a_obs = 0.1, t_obs = 0.1,
              #a_add = 0.1, t_add = 0.1, # for precisions from forest condition tool
              rmean = 0, rprec = 0.00001,
              v = dmbeta, z = dpalpha)  #dmbeta is sample covariate data for testing convergence with individual covs 2/27/24

# draft data object for runs with 
data = list(y = cs_samp_dist, ns = nsites,
            x_ic = x, tau_ic = tic,
            tau_obs = cs_prec_samp,
            rmean = 0, rprec = 0.00001,
            v = dmbeta, z = dpalpha,
            rmean = Rmean, rprec = Rprec)

### RUN THE MODEL
j.pests <- jags.model (file = textConnection(spongy_disturb),
                       data = data,
                       inits = init,
                       n.chains = 3)


# running on 2/27/2024 for dist mag param convergence check with covariate(s) added
jpout<-coda.samples(j.pests,
                    variable.names = c("beta0", "alpha0",
                                       "beta[1]", "alpha[1]",
                                       "tau_obs",
                                       "pa0"),
                    n.iter = 100000,
                    thin=2)

