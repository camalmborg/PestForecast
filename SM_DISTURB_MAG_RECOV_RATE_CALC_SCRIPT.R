#### This is the script for processing GEE data from the forest condition tool

#### Load condition score .csv from GEE extract:
#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"
#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_tcg_mean - 2022_08_31_5k_tcg_mean.csv"
file<-"2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"

#condition score object:
#cond.scores<-read.csv(file)
tcg.values<-read.csv(file)
#tcgs<-tcg.values[,c(grep("^X",colnames(tcg.values)))] #grab with "X1996 eg) for all tcg value columns

#### Function for computing disturbance magnitudes, probabilities, and recovery rates:
spongy_mpr<-function(tcg,distyr){
  #get just the tcg values from data frame:
  tcgs<-tcg[,c(grep("^X",colnames(tcg)))]
  
  #first get just june values:
  junes<-seq(2,130,by=5)
  tcgjune<-as.matrix(tcgs[,junes])
  
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
  
  ### THE CALCULATION LOOP:
  recov.rate<-vector()
  slope<-list()
  mins<-vector()
  colnum<-vector()
  recovcol<-vector()
  end<-length(tcgjune[1,])
  nsite<-length(tcgjune[,1])
  for (i in 1:nsite){
    mins[i]<-min(tcgjune[i,grep(as.character(distyr),colnames(tcgjune)):
                           grep(as.character(distyr+1),colnames(tcgjune))],na.rm=T)        #grabs min forest condition score for each site
    colnum[i]<-which(tcgjune[i,]==mins[i])           #grabs time step at which min score appears in disturbance window
    if(is.na(colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1])){
       recovcol[i]<-colnum[i] + 2
     } else {
       recovcol[i]<-colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1]
     }
     if(recovcol[i]>end){
       recov<-tcgjune[i,colnum[i]:end]
     } else {
       recov<-tcgjune[i,colnum[i]:recovcol[i]]
     }
    # #recov <- tcg[i,colnum[i]:recovcol[i]]      
     ind<-1:length(recov)                      #grabs length of recov rate   
     slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
     recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
  }
  rm(slope)
  
  #calculate disturbance magnitude:
  mags<-steadys-mins
  
  #calculate recovery time:
  recov.time<-mags/recov.rate
  
  #combine data for regressions:
  tcg.mx<-cbind(tcg[-missing,1],steadys,colnum,mins,mags,recov.rate,recov.time)
  tcg.m<-as.data.frame(tcg.mx)
  colnames(tcg.m)<-c("id","steady","colnum","mins","mags","recov.rate","recov.time")
  return(tcg.m)
}

testfx<-spongy_mpr(tcg.values,2016)
