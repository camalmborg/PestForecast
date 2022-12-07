#### This is the script for processing GEE data from the forest condition tool

#### Load condition score .csv from GEE extract:
#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"
#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_tcg_mean - 2022_08_31_5k_tcg_mean.csv"
file<-"2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"

#condition score object:
#cond.scores<-read.csv(file)
tcg.values<-read.csv(file)
tcgs<-tcg.values[,c(grep("^X",colnames(tcg.values)))] #grab with "X1996 eg) for all tcg value columns

#grab june values:


