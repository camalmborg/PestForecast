#### This is the script for processing GEE data from the forest condition tool

#### Load condition score .csv from GEE extract:
#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"
#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_tcg_mean - 2022_08_31_5k_tcg_mean.csv"
file<-"2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"

#condition score object:
#cond.scores<-read.csv(file)
tcg.values<-read.csv(file)
tcgs<-tcg.values[,c(grep("^X",colnames(tcg.values)))] #grab with "X1996 eg) for all tcg value columns

#### Function for computing disturbance magnitudes, probabilities, and recovery rates:
spongy_mpr<-function(tcg,distyr){
  #first get just june values:
  junes<-seq(2,130,by=5)
  tcgjune<-as.matrix(tcg[,junes])
  
  #identify steady state (mean tcg previous 3 years):
  prevyrs<-tcgjune[,grep(as.character(distyr-3),colnames(tcgjune)):
                     grep(as.character(distyr-1),colnames(tcgjune))]
  steadys<-as.matrix(apply(prevyrs[,1:3],1,mean,na.rm=T))
  
  #identify missing data columns:
  NA16<-which(is.na(tcgjune[,grep(as.character(distyr),colnames(tcgjune))]))
  NA17<-which(is.na(tcgjune[,grep(as.character(distyr+1),colnames(tcgjune))]))
  missing<-intersect(NA16,NA17)
  tcgjune<-tcgjune[-missing,]
  steadys<-steadys[-missing]
  
  
}


