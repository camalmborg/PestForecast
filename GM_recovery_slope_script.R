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


#now let's see how to find the minimum value
min<-min(obs[1,],na.rm=T)
colnum<-which(obs[1,]==min)
recov<-obs[1,colnum:130]
ind<-1:length(colnum:130)
slope<-lm(recov~ind)
abline(slope)


