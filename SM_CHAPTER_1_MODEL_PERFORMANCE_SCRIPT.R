# This is the script for fucking around and finding out the best performing models
# CHAPTER 1 BEST MODELS FOR JAGS HINDCASTS

# libraries:
library(dplyr)
library(tidyverse)

# loading the model results masterlist files:

# distmag_uni <- read.csv("CHAPTER_1/2023_September/2023_09_UNI_Dist_Mag_TCG_2016_Models_R2s_AICs.csv")
# distmagcs_uni <- read.csv("CHAPTER_1/2023_September/2023_09_UNI_Dist_Mag_CS_2016_Models_R2s_AICs.csv")
# distmag_multi <- read.csv("CHAPTER_1/2023_September/CHAPTER 1-DISTMAG2016-MULTI-GAM-ANALYSES-TCG-MASTERLIST.csv")
# distmagcs_multi <- read.csv("CHAPTER_1/2023_September/CHAPTER 1-DISTMAG2016- MULTI-GAM-ANALYSES-CS-MASTERLIST.csv")

# distprob16tcg_uni <-read.csv()
# distprob16cs_uni <- read.csv()
# distprob17tcg_uni <- read.csv()
# distprob17cs_uni <- read.csv()
  
distprob16cs_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB-2016-MULTI-AUC-ANALYSES-CS-MASTERLIST.csv")
distprob16tcg_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB-2016-MULTI-AUC-ANALYSES-TCG-MASTERLIST.csv")
distprob17cs_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB-2017-MULTI-AUC-ANALYSES-CS-MASTERLIST.csv")
distprob17tcg_multi <- read.csv("CHAPTER_1/CHAPTER 1-DISTPROB-2017-MULTI-AUC-ANALYSES-TCG-MASTERLIST.csv")
  
### Making uni and multi into something we can Rbind
# choose group you would like to use from above:

# uni <- distmagcs_uni
# multi <- distmagcs_multi

multi <- distprob16cs_multi


# DIST MAG: univariate variable and monthyear columns combine: -------------
# uni <- uni %>%
#   mutate(MONTHYEAR = str_replace(MONTHYEAR, " ", "-")) %>%
#   mutate(VARIABLE1 = paste(VARIABLE1, MONTHYEAR, sep="_")) %>%
#   select(-c(MONTHYEAR))

uni <- uni %>%
  unite(VARIABLE1, c("VarYear", "VarMonth", "Variable")) %>%
  rename(MODEL_ID = Model)

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
model$delAIC <- model$aics - min(model$aics)

### Choosing best performers:
modAIC <- model[order(model$delAIC),]
#modr2 <- model[order(model$R2),]
modauc <- model[order(model$AUC),]
#modauc <- model[order(model$ROC),]

write.csv(modAIC, file="CHAPTER_1/2023_09_29_DISTMAG_2016_CS_delAICsort.csv")
write.csv(modr2, file="CHAPTER_1/2023_09_29_DISTMAG_2016_CS_r2sort.csv")
#write.csv(modAIC, file="CHAPTER_1/dist_prob_2017_cs_models_AICs_sorted.csv")
#write.csv(modauc, file="CHAPTER_1/dist_prob_2017_cs_models_auc_sorted.csv")


##### BIOLOGICAL PLAUSIBILITY 10/03/2023:
distmag_vars <- c(2, 6, 7, 11, 12, 13, 14, 15, 16, 17,
                  22, 23, 24, 25, 31, 32, 33)
distprob_vars <- c( 2, 7, 8, 9, 10, 11, 12,
                    18, 19, 20, 21, 22, 23, 24, 25,
                    32, 33, 34, 35, 38, 39, 40)

dm_bio_best_tcg <- best_tcg[best_tcg$VARIABLE1 %in% distmag_vars &
                       best_tcg$VARIABLE2 %in% distmag_vars &
                       best_tcg$VARIABLE3 %in% distmag_vars,]
dm_bio_best_cs <- best_cs[best_cs$VARIABLE1 %in% distmag_vars &
                             best_cs$VARIABLE2 %in% distmag_vars &
                             best_cs$VARIABLE3 %in% distmag_vars,]

dp16_bio_best_tcg <- distprob16tcg_multi[distprob16tcg_multi$VARIABLE1 %in% distprob_vars &
                                           distprob16tcg_multi$VARIABLE2 %in% distprob_vars &
                                           distprob16tcg_multi$VARIABLE3 %in% distprob_vars,]
dp16_bio_best_cs <- distprob16cs_multi[distprob16cs_multi$VARIABLE1 %in% distprob_vars &
                                          distprob16cs_multi$VARIABLE2 %in% distprob_vars &
                                          distprob16cs_multi$VARIABLE3 %in% distprob_vars,]

dp17_bio_best_tcg <- distprob17tcg_multi[distprob17tcg_multi$VARIABLE1 %in% distprob_vars &
                                           distprob17tcg_multi$VARIABLE2 %in% distprob_vars &
                                           distprob17tcg_multi$VARIABLE3 %in% distprob_vars,]
dp17_bio_best_cs <- distprob16cs_multi[distprob17cs_multi$VARIABLE1 %in% distprob_vars &
                                         distprob17cs_multi$VARIABLE2 %in% distprob_vars &
                                         distprob17cs_multi$VARIABLE3 %in% distprob_vars,]
