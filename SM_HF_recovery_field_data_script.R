#HF Field Data Script

###load GEE data:
#tcg data:
hfplotstcg<-read.csv("HF_2022_Field_Data/HFplots_tcg_mean.csv")
#condition scores:
hfplotscond<-read.csv("HF_2022_Field_Data/HFplots_score_mean.csv")

#remove extra columns:
hftcg<-hfplotstcg[,2:131]
hfcond<-hfplotscond[,2:131]


###load field data:
HF.latlon<-read.csv("HF_2022_Field_Data/2022_HF_plots_latlon.csv")
HF.plot.data<-read.csv("HF_2022_Field_Data/HF_Plot_data.csv")
HF.seedlings<-read.csv("HF_2022_Field_Data/HF_Seedlings_long.csv")
HF.trees<-read.csv("HF_2022_Field_Data/HF_Tree_data.csv")
HF.understory<-read.csv("HF_2022_Field_Data/HF_Und_ground_survey.csv")

#separating mortality sites:
library("dplyr")
HF.condition<-HF.trees %>% group_by(plot_2, Cond) %>% summarize(count=n())
HF.mort<-HF.condition[HF.condition$Cond=="D",]
HF.mort$mort<-1

#merge with plot data:
HF.plot.data<- merge(HF.plot.data, HF.mort, all = TRUE)
HF.plot.data$mort[is.na(HF.plot.data$mort)] <- 0

