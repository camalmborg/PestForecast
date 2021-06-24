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
# min<-min(obs[1,],na.rm=T)
# colnum<-which(obs[1,]==min)
# recov<-obs[1,colnum:130]
# ind<-1:length(colnum:130)
# slope<-lm(recov~ind)
# #abline(slope)
# recov.rate<-slope$coefficients["ind"]

#calculating slope for each recovery rate (loop):
recov.rate<-matrix(NA,nrow=50,ncol=1)
slope<-list()
mins<-matrix(NA,nrow=50,ncol=1)
colnum<-matrix(NA,nrow=50,ncol=1)
for (i in geoID){
  mins[i,]<-min(obs[i,],na.rm=T)            #grabs min forest condition score for each site
  colnum[i,]<-which(obs[i,]==mins[i,])      #grabs time step at which min score appears
  recov<-obs[i,colnum[i,]:(colnum[i,]+10)]  #+10=2 year recovery period
  ind<-1:length(recov)                      #grabs length of recov rate (now 11)
  slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  recov.rate[i,]<-slope[[i]]$coefficients["ind"] #stores recovery rate
}

# ##saving the if-else chunk:
# if((colnum[i,])>120){ 
#   recov<-obs[i,colnum[i,]:130]
# } else {
#   recov<-obs[i,colnum[i,]:(colnum[i,]+10)]     #+10=2 year recovery period
# }

#the 2016-onward version:
recov.rate<-matrix(NA,nrow=50,ncol=1)
slope<-list()
mins<-matrix(NA,nrow=50,ncol=1)
colnum<-matrix(NA,nrow=50,ncol=1)
for (i in geoID){
  mins[i,]<-min(obs[i,105:130],na.rm=T)          #grabs min forest condition score for each site
  colnum[i,]<-which(obs[i,]==mins[i,])           #grabs time step at which min score appears
  recov <- obs[i,colnum[i,]:min(colnum[i,]+10,130)]  
  ind<-1:length(recov)                           #grabs length of recov rate
  slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  recov.rate[i,]<-slope[[i]]$coefficients["ind"] #stores recovery rate
}

obs.recov<-cbind(obs,colnum,mins,recov.rate)


#-----TASSLED CAP GREENNES VERSION--------------------------
#load TCG mean data:
tcgmeans<-"2021_06_23_sample_tcg_mean_GM_50_calib_points.csv"
tcgtable<-read.csv(tcgmeans)
tcg<-as.matrix(tcgtable[1:50,2:131])

#the TCG 2016-onward version:
recov.rate<-matrix(NA,nrow=50,ncol=1)
slope<-list()
mins<-matrix(NA,nrow=50,ncol=1)
colnum<-matrix(NA,nrow=50,ncol=1)
for (i in geoID){
  mins[i,]<-min(tcg[i,105:130],na.rm=T)          #grabs min forest condition score for each site
  colnum[i,]<-which(tcg[i,]==mins[i,])           #grabs time step at which min score appears
  recov <- tcg[i,colnum[i,]:min(colnum[i,]+10,130)]  
  ind<-1:length(recov)                           #grabs length of recov rate
  slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  recov.rate[i,]<-slope[[i]]$coefficients["ind"] #stores recovery rate
}

#find the non-GM disturb condition "steady state"
steady<-apply(tcg[,1:104],1,mean,na.rm=T)

#magnitudes:
mags<-steady-mins
recov.time<-mags/recov.rate

tcg.recov<-cbind(tcg,colnum,mins,steady,mags,recov.rate,recov.time)



#-----visualizations---------
plot(time,tcg[1,],type="l",ylim=(c(-0.05,0.4)))
for (i in geoID){
  lines(time,tcg[i,],col=i)
}

plot(time[105:130],tcg[1,105:130],type="l",ylim=c(-0.01,0.4))
for (i in geoID){
  lines(time[105:130],tcg[i,105:130],col=i)
}  

#sites with negative recovery rates
which(recov.rate<0)

hist(steady)
hist(mags)