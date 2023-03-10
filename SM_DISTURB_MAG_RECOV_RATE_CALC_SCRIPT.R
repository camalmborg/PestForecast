#### This is the script for processing GEE data from the forest condition tool

#### Load condition score .csv from GEE extract:
#cfile<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"
#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_tcg_mean - 2022_08_31_5k_tcg_mean.csv"
file<-"2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"
#file <- "2023_03_08_DATAGRAB/2023_03_08_5000_sites_sample_tcg_mean.csv"
#cfile <- "2023_03_08_DATAGRAB/2023_03_08_5000_sites_sample_score_mean.csv"
cfile<-"2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv"


#tcg and condition score objects:
tcg.values<-read.csv(file)
#tcgs<-tcg.values[,c(grep("^X",colnames(tcg.values)))] #grab with "X1996 eg) for all tcg value columns
cond.scores<-read.csv(cfile)

#### Function for computing disturbance magnitudes, probabilities, and recovery rates:
spongy_mpr<-function(tcg,cs,distyr){
  #get just the tcg values from data frame:
  tcgs<-tcg[,c(grep("^X",colnames(tcg)))]
  sitenum<-as.matrix(1:nrow(tcgs))
  
  #first get just june values:
  junes<-seq(2,ncol(tcgs),by=5) #15/15s = seq starts w 3, 1/1s = seq starts w 2
  tcgjune<-as.matrix(tcgs[,junes])

  #identify steady state (mean tcg previous 3 years):
  prevyrs<-tcgjune[,grep(as.character(distyr-3),colnames(tcgjune)):
                     grep(as.character(distyr-1),colnames(tcgjune))]
  steadys<-as.matrix(apply(prevyrs[,1:3],1,mean,na.rm=T))
  
  #identify missing data columns:
  NAy1<-which(is.na(tcgjune[,grep(as.character(distyr),colnames(tcgjune))]))
  NAy2<-which(is.na(tcgjune[,grep(as.character(distyr+1),colnames(tcgjune))]))
  missing<-intersect(NAy1,NAy2)
  if (length(missing) == 0){
    tcgjune<-tcgjune
    steadys<-steadys
  } else {
  tcgjune<-tcgjune[-missing,]
  steadys<-steadys[-missing]
  }
  
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
  
  #calculate disturbance magnitude/original tcg:
  magsdiv<-mags/steadys
  
  #calculate recovery time:
  recov.time<-mags/recov.rate
  
  #combine data into dataframe:
  if (length(missing) == 0){
    tcg.mx<-cbind(tcg[,1],sitenum,steadys,colnum,mins,mags,magsdiv,recov.rate,recov.time)
  } else{
  tcg.mx<-cbind(tcg[-missing,1],sitenum[-missing],steadys,colnum,mins,mags,magsdiv,recov.rate,recov.time)
  }
  
  tcg.m<-as.data.frame(tcg.mx)
  colnames(tcg.m)<-c("id","sitenum","steady","colnum","mins","mags","magsdiv","recov.rate","recov.time")
  #return(tcg.m)
  
  ### Condition scores for disturbance probabilties
  csj<-cs[,grep("[:.:]06",colnames(cs))]
  if (length(missing)==0){
    csj<-csj
  } else {
    csj<-csj[-missing,]
  }
  dmpr<-cbind(tcg.m,csj)
  #LOOP FOR DISTURBANCE PROBABILITY:
  distprob<-matrix(NA,nrow=nrow(dmpr),ncol=2)
  d = quantile(dmpr[,grep(as.character(distyr-5),colnames(dmpr)):
                       grep(as.character(distyr-1),colnames(dmpr))],
               c(0.01),na.rm=T) #selecting dist sites 1% quant in prev 5 yr
  #determine if disturbance (cond score < d threshold) occurs
  for (i in 1:nrow(dmpr)){
    if (dmpr[i,]$colnum == min(dmpr$colnum) & dmpr[i,grep(as.character(distyr),colnames(dmpr))] <= d){
      distprob[i,1]<- 1
    } else if (dmpr[i,]$colnum == max(dmpr$colnum) & dmpr[i,grep(as.character(distyr+1),colnames(dmpr))] <= d){
      distprob[i,2]<- 1
    } else {
      distprob[i,] <- 0
    }
    distprob[is.na(distprob)] <- 0
  }
  
  #add disturbance probabilities to dataframe:
  tcg.m$dpy1<-distprob[,1]
  tcg.m$dpy2<-distprob[,2]
  
  return(tcg.m)
}

testfx<-spongy_mpr(tcg.values,cond.scores,2016)

#once again I am a beautiful genius!!!!




###------------15s/15s version--------------------------------------------#
#### Function for computing disturbance magnitudes, probabilities, and recovery rates:
spongy_mpr_15s<-function(tcg,cs,distyr){
  #get just the tcg values from data frame:
  tcgs<-tcg[,c(grep("^X",colnames(tcg)))]
  sitenum<-as.matrix(1:nrow(tcgs))
  
  #first get just june values:
  junejul<-seq(4,ncol(tcgs),by=5) #15/15s = seq starts w 3, 1/1s = seq starts w 2
  tcgjune<-as.matrix(tcgs[,junejul])
  
  #identify steady state (mean tcg previous 3 years):
  prevyrs<-tcgjune[,grep(as.character(distyr-3),colnames(tcgjune)):
                     grep(as.character(distyr-1),colnames(tcgjune))]
  steadys<-as.matrix(apply(prevyrs[,1:3],1,mean,na.rm=T))
  
  #identify missing data columns:
  NAy1<-which(is.na(tcgjune[,grep(as.character(distyr),colnames(tcgjune))]))
  NAy2<-which(is.na(tcgjune[,grep(as.character(distyr+1),colnames(tcgjune))]))
  missing<-intersect(NAy1,NAy2)
  if (length(missing) == 0){
    tcgjune<-tcgjune
    steadys<-steadys
  } else {
    tcgjune<-tcgjune[-missing,]
    steadys<-steadys[-missing]
  }
  
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
    #disturb cols:
    colnum[i]<-which(tcgjune[i,grep(as.character(distyr),colnames(tcgjune)):
                               grep(as.character(distyr+1),colnames(tcgjune))]==mins[i])           #grabs time step at which min score appears in disturbance window
  #   if(is.na(colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1])){
  #     recovcol[i]<-colnum[i] + 2
  #   } else {
  #     recovcol[i]<-colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1]
  #   }
  #   if(recovcol[i]>end){
  #     recov<-tcgjune[i,colnum[i]:end]
  #   } else {
  #     recov<-tcgjune[i,colnum[i]:recovcol[i]]
  #   }
  #   # #recov <- tcg[i,colnum[i]:recovcol[i]]
  #   ind<-1:length(recov)                      #grabs length of recov rate
  #   slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  #   recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
  }
  rm(slope)
  
  #calculate disturbance magnitude:
  mags<-steadys-mins
  
  #calculate disturbance magnitude/original tcg:
  magsdiv<-mags/steadys
  
  #calculate recovery time:
  recov.time<-mags/recov.rate
  
  #combine data into dataframe:
  if (length(missing) == 0){
    tcg.mx<-cbind(tcg[,1],sitenum,steadys,colnum,mins,mags,magsdiv,recov.rate,recov.time)
  } else{
    tcg.mx<-cbind(tcg[-missing,1],sitenum[-missing],steadys,colnum,mins,mags,magsdiv,recov.rate,recov.time)
  }
  
  tcg.m<-as.data.frame(tcg.mx)
  colnames(tcg.m)<-c("id","sitenum","steady","colnum","mins","mags","magsdiv","recov.rate","recov.time")
  #return(tcg.m)
  
  ### Condition scores for disturbance probabilties
  csj<-cs[,grep("[:.:]06",colnames(cs))]
  if (length(missing)==0){
    csj<-csj
  } else {
    csj<-csj[-missing,]
  }
  dmpr<-cbind(tcg.m,csj)
  #LOOP FOR DISTURBANCE PROBABILITY:
  distprob<-matrix(NA,nrow=nrow(dmpr),ncol=2)
  d = quantile(dmpr[,grep(as.character(distyr-5),colnames(dmpr)):
                      grep(as.character(distyr-1),colnames(dmpr))],
               c(0.01),na.rm=T) #selecting dist sites 1% quant in prev 5 yr
  #determine if disturbance (cond score < d threshold) occurs
  for (i in 1:nrow(dmpr)){
    if (dmpr[i,]$colnum == min(dmpr$colnum) & dmpr[i,grep(as.character(distyr),colnames(dmpr))] <= d){
      distprob[i,1]<- 1
    # } else if (dmpr[i,]$colnum == max(dmpr$colnum) & dmpr[i,grep(as.character(distyr+1),colnames(dmpr))] <= d){
    #   distprob[i,2]<- 1
    } else {
      distprob[i,] <- 0
    }
    distprob[is.na(distprob)] <- 0
  }
  
  #add disturbance probabilities to dataframe:
  tcg.m$dpy1<-distprob[,1]
  tcg.m$dpy2<-distprob[,2]
  
  return(tcg.m)
}

testfx<-spongy_mpr_15s(tcg.values,cond.scores,2016)

#this is still not working^^