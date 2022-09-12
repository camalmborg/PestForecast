#HF Field Data Script

#load GEE data:
hfplotstcg<-read.csv("HF_2022_Field_Data/HFplots_tcg_mean.csv")
hfplotscond<-read.csv("HF_2022_Field_Data/HFplots_score_mean.csv")

hftcg<-hfplotstcg[,2:131]
hfcond<-hfplotscond[,2:131]
