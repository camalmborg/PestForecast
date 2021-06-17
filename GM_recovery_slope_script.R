#load data for all test sites:
time=1:130
NT=length(time)
geoID=1:50
gm.data<-list()
for(i in geoID){
  file<-paste0("TEST_monthly_data/GM_",geoID[i],"_mcs.csv")
  gm.data[[i]]<-read.csv(file)
}

#load condition score means for each site into matrix:
obs<-matrix(NA,nrow=length(geoID),ncol=130)
for (i in geoID){
  obs[i,]<-gm.data[[i]]$score_mean
}


#calculating slope for each recovery rate (one site):
min<-min(obs[1,],na.rm=T)
colnum<-which(obs[1,]==min)
recov<-obs[1,colnum:130]
ind<-1:length(colnum:130)
slope<-lm(recov~ind)
#abline(slope)
recov.rate<-slope$coefficients["ind"]


#calculating slope for each recovery rate (loop):
recov.rate=matrix(NA,nrow=50,ncol=1)
colnum<-matrix(NA,nrow=50,ncol=1)
for (i in geoID){
  min<-min(obs[i,],na.rm=T)
  colnum[i,]<-which(obs[i,]==min)
  recov<-obs[i,colnum[i,]:130]
  ind<-1:length(colnum[i,]:130)
  slope<-lm(recov~ind)
  recov.rate[i,]<-slope$coefficients["ind"]
}

obs.recov<-cbind(obs,colnum,recov.rate)
