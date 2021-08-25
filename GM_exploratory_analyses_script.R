####load data for 5000 sites:---------
#load monthly condition scores:
cond.scores.mo<-read.csv("2020_07_10_sample_score_mean_MONTHLY.csv")
cond.mo<-cond.scores.mo[,2:131]

#load annual condition scores:
cond.scores.an<-read.csv("2020_07_10_sample_score_mean_ANNUAL.csv")
cond.an<-cond.scores.an[,2:27]

#load monthly tcg:
tcg.month<-read.csv("2020_07_10_sample_tcg_mean_MONTHLY.csv")
tcg.mo<-tcg.month[,2:131]



####remove may:-----
nomay<-seq(1,130,by=5)
tcg.nomay<-as.matrix(tcg.mo[,-nomay])

#### steady state previous 3 years 2015------
#find the non-GM disturb condition "steady state"
steady<-apply(tcg.nomay[,73:84],1,mean,na.rm=T)
#steadyall<-apply(tcg.nomay[,1:84],1,mean,na.rm=T)


####the TCG 2016-onward MONTHLY version (no Mays):-----
#ind<-list()
recov.rate<-vector()
slope<-list()
mins<-vector()
colnum<-vector()
recovcol<-vector()
end<-length(tcg.nomay[1,])
nsite<-length(tcg.nomay[,1])
for (i in 1:nsite){
  mins[i]<-min(tcg.nomay[i,85:92],na.rm=T)          #grabs min forest condition score for each site
  colnum[i]<-which(tcg.nomay[i,]==mins[i])           #grabs time step at which min score appears in disturbance window
  if(is.na(colnum[i] + which(tcg.nomay[i,colnum[i]:end]>=steady[i])[1])){
    recovcol[i]<-colnum[i] + 2
  } else {
    recovcol[i]<-colnum[i] + which(tcg.nomay[i,colnum[i]:end]>=steady[i])[1]
  }
  if(recovcol[i]>end){
    recov<-tcg.nomay[i,colnum[i]:end]
  } else {
    recov<-tcg.nomay[i,colnum[i]:recovcol[i]]
  }
  #recov <- tcg[i,colnum[i]:recovcol[i]]      
  ind<-1:length(recov)                      #grabs length of recov rate   
  slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
}

#magnitudes:
mags<-steady-mins
recov.time<-mags/recov.rate




#### THE JUNE VERSION:-----
####just junes
junes<-seq(2,130,by=5)
tcg.june<-as.matrix(tcg.mo[,junes])

###steady state
steady<-apply(tcg.june[,19:21],1,mean,na.rm=T)

recov.rate<-vector()
slope<-list()
mins<-vector()
colnum<-vector()
recovcol<-vector()
end<-length(tcg.june[1,])
nsite<-length(tcg.june[,1])
for (i in 1:nsite){
  mins[i]<-min(tcg.june[i,21:23],na.rm=T)          #grabs min forest condition score for each site
  colnum[i]<-which(tcg.june[i,]==mins[i])           #grabs time step at which min score appears in disturbance window
  if(is.na(colnum[i] + which(tcg.june[i,colnum[i]:end]>=steady[i])[1])){
    recovcol[i]<-colnum[i] + 2
  } else {
    recovcol[i]<-colnum[i] + which(tcg.june[i,colnum[i]:end]>=steady[i])[1]
  }
  if(recovcol[i]>end){
    recov<-tcg.june[i,colnum[i]:end]
  } else {
    recov<-tcg.june[i,colnum[i]:recovcol[i]]
  }
  #recov <- tcg[i,colnum[i]:recovcol[i]]      
  ind<-1:length(recov)                      #grabs length of recov rate   
  slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
}

#magnitudes:
mags<-steady-mins
recov.time<-mags/recov.rate
