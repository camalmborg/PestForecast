# CHAPTER 1
# this is a script for making matrix of models for disturbance magnitude and 
# disturbance probability parameterizations

### LIBRARIES:

### LOAD DATA:
# the environmental variables
# for disturbance magnitude parameters
varfile_m <- "CHAPTER_1/DATA/MV_2023_12_DATA_distmag.csv"
magvars <- read.csv(varfile_m)[,-1]
# for disturbance probability parameters
varfile_p <- "CHAPTER_1/DATA/MV_2023_12_DATA_distprob.csv"
probvars <- read.csv(varfile_p)[,-1]
# load data for dmr to remove missing sites from environmental data
dmr_file <- "CHAPTER_1/DATA/2023_12_DMR_DATA_TCG.csv"
dmr <- read.csv(dmr_file)
# environmental data without missing sites
magvars <- magvars[dmr$X,]
probvars <- probvars[dmr$X,]

### Load model list data frame:
#modfile <- "CHAPTER_1/2024_02_JAGS_models/2024_02_20_Dist_Mag_Models.csv"
modfile <- "CHAPTER_1/2024_02_JAGS_models/2024_02_22_Dist_Prob_Models.csv"
# disturbance mag/prob model file
modeldf <- read.csv(modfile, header = F)


### Making list of data for models:
# empty list
covls <- list()
# mag (magvars) or prob (probvars)
#covs <- magvars
covs <- probvars
# loop over models
for (i in 1:nrow(modeldf)){
  # grab column numbers from model table
  cols <- as.numeric(modeldf[i,])  #add a row of ones to each, make i into i+1
  # remove NAs
  cols <- cols[!is.na(cols)]
  # make column of 1s for matrix multiplication step in JAGS model (intercept term)
  int <- rep(1, nrow(covs))
  # make dataframe with columns from covariates
  df <- cbind(int, covs[,c(cols)])
  # add to the list
  covls[[i]] <- df  ## have to change this line to make sure its mags or prob
}

# for saving result
magls <- covls
probls <- covls
