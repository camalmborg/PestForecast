#code for disturbance probability exploratory analyses

### Load condition score data:
cfile<-"2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv"
condscores<-read.csv(cfile)

spongy_dp<-function(cs,){
  junes<-cs[,grep("[:.:]06",colnames(cs))]
  
}

junetest<-condscores[,grep("[:.:]06",colnames(condscores))]
