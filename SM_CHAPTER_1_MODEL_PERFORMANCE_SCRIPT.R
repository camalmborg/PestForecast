# This is the script for fucking around and finding out the best performing models
# CHAPTER 1 BEST MODELS FOR JAGS HINDCASTS

# libraries:
library(dplyr)
library(tidyverse)

# loading the model results masterlist files:
distmag_uni <- read.csv("CHAPTER_1/CHAPTER 1-DISTMAG-UNI-GAM-ANALYSES-TCG-MASTERLIST.csv")
distmag_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTMAG-MULTI-GAM-ANALYSES-TCG-MASTERLIST.csv")
#distprob16_uni <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2016-UNI-AUC-ANALYSES-TCG-MASTERLIST.csv")
#distprob16_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2016-MULTI-AUC-ANALYSES-TCG-MASTERLIST.csv")
#distprob17_uni <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2017-UNI-AUC-ANALYSES-TCG-MASTERLIST.csv")
#distprob17_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2017-MULTI-AUC-ANALYSES-TCG-MASTERLIST.csv")

### Making uni and multi into something we can Rbind
# choose group you would like to use from above:
uni <- distmag_uni
multi <- distmag_multi

# univariate variable and monthyear columns combine:
uni <- uni %>%
  mutate(MONTHYEAR = str_replace(MONTHYEAR, " ", "-")) %>%
  mutate(VARIABLE1 = paste(VARIABLE1, MONTHYEAR, sep="_")) %>%
  select(-c(MONTHYEAR))

# make new dataframe with 2 new NA columsn for VARIABLE2 and VARIABLE3:
v2v3 <- data.frame(matrix(nrow=nrow(uni), ncol=2, data=NA))
colnames(v2v3) <- c("VARIABLE2","VARIABLE3")
# combine:
uni <- cbind.data.frame(uni$MODEL_ID, uni$VARIABLE1,
                           v2v3$VARIABLE2,v2v3$VARIABLE3,
                           uni$R2,uni$AIC)
colnames(uni) <- names(multi)

# Combining them into one big dataframe:
model <- rbind(uni, multi)

### Calculate delAIC:
model$delAIC <- model$AIC - min(model$AIC)

### Choosing best performers:

