#### This is the script for processing GEE data from the forest condition tool

#### Load condition score .csv from GEE extract:
#cfile<-"2022_08_31_DATAGRAB/2022_08_31_5k_score_mean - 2022_08_31_5k_score_mean.csv"
#file<-"2022_08_31_DATAGRAB/2022_08_31_5k_tcg_mean - 2022_08_31_5k_tcg_mean.csv"
#file<-"2022_08_31_DATAGRAB/2022_12_7_sample_tcg_mean_5k.csv"
#cfile<-"2022_08_31_DATAGRAB/2022_12_7_sample_score_mean_5k.csv"

#### HARVARD FOREST DATA:
file <- "HF_2022_Field_Data/GEE_Data/2023_05_17_hfplots_sample_tcg_mean.csv"
cfile <- "HF_2022_Field_Data/GEE_Data/2023_05_17_hfplots_sample_score_mean.csv"

#tcg and condition score objects:
tcg.values<-read.csv(file)
#tcgs<-tcg.values[,c(grep("^X",colnames(tcg.values)))] #grab with "X1996 eg) for all tcg value columns
cond.scores<-read.csv(cfile)

#### Function for computing disturbance magnitudes, probabilities, and recovery rates:
spongy_mpr_cs<-function(tcg,cs,distyr){
  #get just the tcg values from data frame:
  tcgs<-tcg[,c(grep("^X",colnames(tcg)))]
  sitenum<-as.matrix(1:nrow(tcgs))
  
  #get just cs values from data frame:
  css <- cs[,c(grep("^X",colnames(cs)))]
  
  #first get just june values:
  junes<-seq(3,ncol(tcgs),by=5)  #3 is if tcg/condscores data have april in them, 2 if data starts w/May
  #tcgjune<-as.matrix(tcgs[,junes])
  csjune <- as.matrix(css[,junes])
  
  #identify steady state (mean tcg previous 3 years):
  prevyrs<-csjune[,grep(as.character(distyr-3),colnames(csjune)):
                     grep(as.character(distyr-1),colnames(csjune))]
  steadys<-as.matrix(apply(prevyrs[,1:3],1,mean,na.rm=T))
  
  #identify missing data columns:
  NAy1<-which(is.na(csjune[,grep(as.character(distyr),colnames(csjune))]))
  NAy2<-which(is.na(csjune[,grep(as.character(distyr+1),colnames(csjune))]))
  
  missing<-intersect(NAy1,NAy2) 
  
  if (length(missing) == 0){
    csjune<-csjune
    steadys<-steadys
  } else {
    csjune<-csjune[-missing,]
    steadys<-steadys[-missing]
  }
  
  ### THE CALCULATION LOOP:
  #recov.rate<-vector()
  slope<-list()
  mins<-vector()
  colnum<-vector()
  #recovcol<-vector()
  end<-length(csjune[1,])
  nsite<-length(csjune[,1])
  #recov<-list()
  for (i in 1:nsite){
    mins[i]<-min(csjune[i,grep(as.character(distyr),colnames(csjune)):
                           grep(as.character(distyr+1),colnames(csjune))],na.rm=T)        #grabs min forest condition score for each site
    colnum[i]<-which(csjune[i,grep(as.character(distyr),colnames(csjune)):
                               grep(as.character(distyr+1),colnames(csjune))]==mins[i])           #grabs time step at which min score appears in disturbance window
    colnum[i] <- ifelse(colnum[i]==1, (grep(as.character(distyr),colnames(csjune))),
                        (grep(as.character(distyr),colnames(csjune))+1))
    
    # if(is.na(colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1])){
    #   recovcol[i]<-colnum[i] + 2 #why?
    # } else {
    #   recovcol[i]<-colnum[i] + which(tcgjune[i,colnum[i]:end]>=steadys[i])[1]
    # }
    # 
    # if(recovcol[i]>end){
    #   recov<-tcgjune[i,colnum[i]:end]
    # } else {
    #   recov<-tcgjune[i,colnum[i]:recovcol[i]]
    # }
    # # ##recov <- tcg[i,colnum[i]:recovcol[i]]
    # ind<-1:length(recov)                      #grabs length of recov rate
    # slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
    # recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
  }
  rm(slope)
  
  #calculate disturbance magnitude:
  mags<-steadys-mins
  
  #calculate disturbance magnitude/original tcg:
  magsdiv<-mags/steadys
  
  #calculate recovery time:
  #recov.time<-mags/recov.rate
  
  #combine data into dataframe:
  if (length(missing) == 0){
    cs.mx<-cbind(css[,1],sitenum,steadys,colnum,mins,mags,magsdiv)#,recov.rate,recov.time)
  } else{
    cs.mx<-cbind(css[-missing,1],sitenum[-missing],steadys,colnum,mins,mags,magsdiv)#,recov.rate,recov.time)
  }
  
  cs.m<-as.data.frame(cs.mx)
  colnames(cs.m)<-c("id","sitenum","steady","colnum","mins","mags","magsdiv")#,"recov.rate","recov.time")
  #return(tcg.m)
  
  ### Condition scores for disturbance probabilties
  csj<-cs[,grep("[:.:]06",colnames(cs))]
  if (length(missing)==0){
    csj<-csj
  } else {
    csj<-csj[-missing,]
  }
  dmpr<-cbind(cs.m,csj)
  #LOOP FOR DISTURBANCE PROBABILITY:
  distprob<-matrix(NA,nrow=nrow(dmpr),ncol=2)
  d = quantile(dmpr[,grep(as.character(distyr-5),colnames(dmpr)):
                      grep(as.character(distyr-1),colnames(dmpr))],
               c(0.1),na.rm=T) #selecting dist sites 1% quant in prev 5 yr. For other data: originally 0.01 
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
  cs.m$dpy1<-distprob[,1]
  cs.m$dpy2<-distprob[,2]
  
  return(cs.m)
}

testfx_t<-spongy_mpr_cs(tcg.values,cond.scores,2016)

