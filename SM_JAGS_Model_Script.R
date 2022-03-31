#This script updates the initial runs of the JAGS model for pest forecast
#Initially the model ran on 50 test sites, not the full 5000 (4997) sites
#New script updates the initial dataset build and includes the updated model

#load full dataset (5000 sites):
#model uses condition scores, not TCG raw data
cond.scores.mo<-read.csv("2020_07_10_sample_score_mean_MONTHLY.csv")
condscores<-cond.scores.mo[,2:131]
condscores.samp<-condscores[1:250,]

#number of sites
nsites = nrow(condscores.samp)
