#### This is the script for processing GEE data from the forest condition tool

#### Load condition score .csv from GEE extract:
#DONT USE:#cfile<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv" #because these are annuals
#DONT USE:#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_tcg_mean - 2022_08_31_5k_tcg_mean.csv"
file<-"2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"
cfile<-"2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv"
#file <- "2023_03_08_DATAGRAB/2023_5000_points_sample_tcg_mean.csv"
#cfile <- "2023_03_08_DATAGRAB/2023_5000_points_sample_score_mean.csv"

#### HARVARD FOREST DATA:
#file <- "HF_2022_Field_Data/GEE_Data/2023_05_17_hfplots_sample_tcg_mean.csv"
#cfile <- "HF_2022_Field_Data/GEE_Data/2023_05_17_hfplots_sample_score_mean.csv"
#file <- "HF_2022_Field_Data/GEE_Data/2023_08_31_HF_sample_tcg_mean.csv"
#cfile <- "HF_2022_Field_Data/GEE_Data/2023_08_31_HF_sample_score_mean.csv"

#tcg and condition score objects:
tcg.values<-read.csv(file)
#tcgs<-tcg.values[,c(grep("^X",colnames(tcg.values)))] #grab with "X1996 eg) for all tcg value columns
cond.scores<-read.csv(cfile)

#### Function for computing disturbance magnitudes, probabilities, and recovery rates:-------
#tcg = tcg scores GEE data for all sites
#cs = condition score GEE data for all sites
#distyr = year of outbreak (2016)
#monthnum = 3 if tcg and condscore data include April, 2 if data does not include April
#seqnum = number of growing season months in tcg/cs data (April-Aug = 5, April-Sep = 6, etc.)
spongy_mpr<-function(tcg,cs,distyr,monthnum,seqnum){
  #get just the tcg values from data frame:
  tcgs<-tcg[,c(grep("^X",colnames(tcg)))]
  sitenum<-as.matrix(1:nrow(tcgs))
  
  #first get just june values:
  junes<-seq(monthnum,ncol(tcgs),by=seqnum)  #3 is if tcg/condscores data have april in them, 2 if data starts w/May
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
  recov<-list()
  for (i in 1:nsite){
    mins[i]<-min(tcgjune[i,grep(as.character(distyr),colnames(tcgjune)):
                           grep(as.character(distyr+1),colnames(tcgjune))],na.rm=T)        #grabs min forest condition score for each site
    colnum[i]<-which(tcgjune[i,grep(as.character(distyr),colnames(tcgjune)):
                               grep(as.character(distyr+1),colnames(tcgjune))]==mins[i])           #grabs time step at which min score appears in disturbance window
    colnum[i] <- ifelse(colnum[i]==1, (grep(as.character(distyr),colnames(tcgjune))),
                  (grep(as.character(distyr),colnames(tcgjune))+1))
    
    ## RECOVERY, ADD BACK IN FOR CHAP 2:
    if(is.na(colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1])){
       recovcol[i]<-colnum[i] + 2 #puts recovcol outside the disturbance window
     } else {
       recovcol[i]<-colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1]
     }

     if(recovcol[i]>end){
       recov<-tcgjune[i,colnum[i]:end]
     } else {
       recov<-tcgjune[i,colnum[i]:recovcol[i]]
     }
    # ##recov <- tcg[i,colnum[i]:recovcol[i]]
     ind<-1:length(recov)                      #grabs length of recov rate
     slope[[i]]<-lm(recov~ind)   #runs the lm to find slope of recov rate, saves output
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
    tcg.m<-data.frame(tcg[,1],sitenum,steadys,colnum,mins,mags,magsdiv,recov.rate,recov.time)
  } else{
  tcg.m<-data.frame(tcg[-missing,1],sitenum[-missing],steadys,colnum,mins,mags,magsdiv,recov.rate,recov.time)
  }

  #tcg.m<-as.data.frame(tcg.mx)
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
  d = quantile(dmpr[,grep(as.character(distyr-6),colnames(dmpr)):
                       grep(as.character(distyr-1),colnames(dmpr))],
               c(0.01),na.rm=T) #selecting dist sites 1% quant in prev 5 yr. 
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

testfx<-spongy_mpr(tcg.values,cond.scores, 2016, 2, 5)
testfx2 <- spongy_mpr(cond.scores, cond.scores, 2016, 2, 5)

#hf_mags_2 <- spongy_mpr(tcg.values, cond.scores, 2016)
#once again I am a beautiful genius!!!!

#### DOUBLE DISTURB CALC VERSION 8/31/2023: -------------------

#will grab mins/mags in both years:
spongy_mpr_2<-function(tcg,cs,distyr,monthnum,seqnum){
  #get just the tcg values from data frame:
  tcgs<-tcg[,c(grep("^X",colnames(tcg)))]
  sitenum<-as.matrix(1:nrow(tcgs))
  
  #first get just june values:
  junes<-seq(monthnum,ncol(tcgs),by=seqnum)  #3 is if tcg/condscores data have april in them, 2 if data starts w/May
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
  #recov.rate<-vector()
  #slope<-list()
  mins<-vector()
  colnum<-vector()
  #recovcol<-vector()
  end<-length(tcgjune[1,])
  nsite<-length(tcgjune[,1])
  #recov<-list()
  for (i in 1:nsite){
    mins[i]<-min(tcgjune[i,grep(as.character(distyr),colnames(tcgjune)):
                           grep(as.character(distyr+1),colnames(tcgjune))],na.rm=T)        #grabs min forest condition score for each site
    colnum[i]<-which(tcgjune[i,grep(as.character(distyr),colnames(tcgjune)):
                               grep(as.character(distyr+1),colnames(tcgjune))]==mins[i])           #grabs time step at which min score appears in disturbance window
    colnum[i] <- ifelse(colnum[i]==1, (grep(as.character(distyr),colnames(tcgjune))),
                        (grep(as.character(distyr),colnames(tcgjune))+1))
    
    ### RECOVERY, ADD BACK IN FOR CHAP 2:
    # if(is.na(colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1])){
    #    recovcol[i]<-colnum[i] + 2 #puts recovcol outside the disturbance window
    #  } else {
    #    recovcol[i]<-colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1]
    #  }
    # 
    #  if(recovcol[i]>end){
    #    recov<-tcgjune[i,colnum[i]:end]
    #  } else {
    #    recov<-tcgjune[i,colnum[i]:recovcol[i]]
    #  }
    # # ##recov <- tcg[i,colnum[i]:recovcol[i]]
    #  ind<-1:length(recov)                      #grabs length of recov rate
    #  slope[[i]]<-lm(recov~ind)   #runs the lm to find slope of recov rate, saves output
    #  recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
  }
  #rm(slope)
  
  #calculate disturbance magnitude:
  mags<-steadys-mins
  mags16 <- steadys-(tcgjune[,grep(as.character(distyr),colnames(tcgjune))])
  mags17 <- steadys-(tcgjune[,grep(as.character(distyr+1),colnames(tcgjune))])
  
  #calculate disturbance magnitude/original tcg:
  magsdiv<-mags/steadys
  
  # #calculate recovery time:
  # recov.time<-mags/recov.rate
  
  #combine data into dataframe:
  if (length(missing) == 0){
    tcg.m<-data.frame(tcg[,1],sitenum,steadys,colnum,mins,mags,magsdiv,mags16,mags17)#,recov.rate,recov.time)
  } else{
    tcg.m<-data.frame(tcg[-missing,1],sitenum[-missing],steadys,colnum,mins,mags,magsdiv,mags16,mags17)#,recov.rate,recov.time)
  }
  
  #tcg.m<-as.data.frame(tcg.mx)
  colnames(tcg.m)<-c("id","sitenum","steady","colnum","mins","mags","magsdiv","mags16","mags17")
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
  d = quantile(dmpr[,grep(as.character(distyr-6),colnames(dmpr)):
                      grep(as.character(distyr-1),colnames(dmpr))],
               c(0.01),na.rm=T) #selecting dist sites 1% quant in prev 5 yr. 
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

hf_fx_tcg <- spongy_mpr_2(tcg.values, cond.scores, 2016, 2, 5)
hf_fx_cs <- spongy_mpr_2(cond.scores, cond.scores, 2016, 2, 5)
