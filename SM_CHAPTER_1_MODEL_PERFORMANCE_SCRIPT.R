# This is the script for fucking around and finding out the best performing models
# CHAPTER 1 BEST MODELS FOR JAGS HINDCASTS

# libraries:
library(dplyr)
library(tidyverse)

# loading the model results masterlist files:
#distmag_uni <- read.csv("CHAPTER_1/CHAPTER 1-DISTMAG-UNI-GAM-ANALYSES-TCG-MASTERLIST.csv")
#distmag_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTMAG-MULTI-GAM-ANALYSES-TCG-MASTERLIST.csv")
#distprob16_uni <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2016-UNI-AUC-ANALYSES-TCG-MASTERLIST.csv")
#distprob16_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2016-MULTI-AUC-ANALYSES-TCG-MASTERLIST.csv")
#distprob17_uni <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2017-UNI-AUC-ANALYSES-TCG-MASTERLIST.csv")
#distprob17_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2017-MULTI-AUC-ANALYSES-TCG-MASTERLIST.csv")

#distmag_uni <- read.csv("CHAPTER_1/CHAPTER 1-DISTMAG-UNI-GAM-ANALYSES-CS-MASTERLIST.csv")
#distmag_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTMAG-MULTI-GAM-ANALYSES-CS-MASTERLIST.csv")
#distprob16_uni <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2016-UNI-AUC-ANALYSES-CS-MASTERLIST.csv")
#distprob16_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2016-MULTI-AUC-ANALYSES-CS-MASTERLIST.csv")
#distprob17_uni <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2017-UNI-AUC-ANALYSES-CS-MASTERLIST.csv")
#distprob17_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB2017-MULTI-AUC-ANALYSES-CS-MASTERLIST.csv")

distmag_uni <- read.csv("CHAPTER_1/2023_September/2023_09_UNI_Dist_Mag_TCG_2016_Models_R2s_AICs.csv")
distmagcs_uni <- read.csv("CHAPTER_1/2023_September/2023_09_UNI_Dist_Mag_CS_2016_Models_R2s_AICs.csv")

### Making uni and multi into something we can Rbind
# choose group you would like to use from above:
uni <- distmag_uni
multi <- distprob16_multi

# DIST MAG: univariate variable and monthyear columns combine: -------------
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

# DIST PROB: univariate variable and monthyear columns combine: -------------
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
                        uni$AUC,uni$AIC)
colnames(uni) <- names(multi)


### Combining them into one big dataframe:
model <- rbind(uni, multi)

### Calculate delAIC: -----------
model$delAIC <- model$AIC - min(model$AIC)

### Choosing best performers:
modAIC <- model[order(model$delAIC),]
modr2 <- model[order(model$R2),]
#modauc <- model[order(model$AUC),]

#write.csv(modAIC, file="CHAPTER_1/dist_mag_cs_models_AICs_sorted.csv")
#write.csv(modr2, file="CHAPTER_1/dist_mag_cs_models_r2s_sorted.csv")
write.csv(modAIC, file="CHAPTER_1/dist_prob_2016_cs_models_AICs_sorted.csv")
write.csv(modauc, file="CHAPTER_1/dist_prob_2016_cs_models_auc_sorted.csv")
