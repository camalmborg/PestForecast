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


#------Forest Condition Score version:---------------

if(FALSE){
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

#find the non-GM disturb condition "steady state"
steady<-apply(obs[,90:105],1,mean,na.rm=T)

#magnitudes:
mags<-steady-mins
recov.time<-mags/recov.rate

#previous year's greenness:
endprevyr<-colnum-1
startprevyr<-endprevyr-5
prevyr<-matrix(NA,nrow=50,ncol=1)
for (i in geoID){
  prevyr[i,]<-(mean(na.omit(tcg[i,startprevyr[i,]:endprevyr[i,]])))
}

plot(prevyr,mags)
}


#-----TASSLED CAP GREENNESS VERSION--------------------------
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
  mins[i,]<-min(tcg[i,110:118],na.rm=T)          #grabs min forest condition score for each site
  colnum[i,]<-which(tcg[i,]==mins[i,])           #grabs time step at which min score appears
  recov <- tcg[i,colnum[i,]:min(colnum[i,]+15,130)]  
  ind<-1:length(recov)                           #grabs length of recov rate
  slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  recov.rate[i,]<-slope[[i]]$coefficients["ind"] #stores recovery rate
}

#find the non-GM disturb condition "steady state"
steady<-apply(tcg[,90:110],1,mean,na.rm=T) 
steadyall<-apply(tcg[,1:110],1,mean,na.rm=T)

#magnitudes:
mags<-steady-mins
recov.time<-mags/recov.rate

tcg.recov.mx<-cbind(steadyall,steady,colnum,mins,mags,recov.rate,recov.time)
tcg.recov<-as.data.frame(tcg.recov.mx)
colnames(tcg.recov)<-c("Steady State (All Years)","Steady State 2012-2015","Column Number","Minimum TCG Value","Disturbance Magnitude","Recovery Rate","Recovery Time")

if (FALSE){
#previous year's greenness analysis:
endprevyr<-colnum-1
startprevyr<-endprevyr-5
prevyr<-matrix(NA,nrow=50,ncol=1)
for (i in geoID){
  prevyr[i,]<-(mean(na.omit(tcg[i,startprevyr[i,]:endprevyr[i,]])))
}

#reduction in greenness from previous year using magnitude:
defol<-prevyr-mags

#make a nice table of it all:
recov.view.mx<-cbind(steadyall,steady,prevyr,mags,defol,mins,recov.rate,recov.time)
recov.view<-as.data.frame(recov.view.mx)
colnames(recov.view)<-c("steadyall","steady","prevyr","mags","defol","mins","recov.rate","recov.time")
}



#oh, you know we got the plots:
plot(steady,mags)
plot(steady,mins)
plot(steady,recov.rate)
plot(steady,recov.time)

plot(mags,recov.rate)
plot(mags,recov.time)
plot(mags,steady)

plot(mins,recov.time)
plot(mins,recov.rate)
plot(mins,mags)
plot(mins,steady)


hist(steady)
hist(mags)
hist(mins)
hist(recov.rate)
hist(recov.time)


defol.lm<-lm(mags~steady,tcg.recov)
recov.rate.lm<-lm(recov.rate~steady,tcg.recov)
recov.time.lm<-lm(recov.time~steady,tcg.recov)

plot(steady,mags)
abline(defol.lm)

plot(steady,recov.rate)
abline(recov.rate.lm)

plot(steady,recov.time)
abline(recov.time.lm)

#-----some visualizations----------------------------------
plot(time,tcg[1,],type="l",ylim=(c(-0.05,0.4)))
for (i in geoID){
  lines(time,tcg[i,],col=i)
}

plot(time[100:125],tcg[1,100:125],type="l",ylim=c(-0.01,0.4))
for (i in geoID){
  lines(time[100:125],tcg[i,100:125],col=i)
}  

#for recov.rate<0
plot(time,tcg[12,],type="l",ylim=c(-0.01,0.4),col=12)
lines(time,tcg[46,],col=46)

plot(time[105:130],tcg[12,105:130],type="l",ylim=c(-0.01,0.4),col=12)
lines(time[105:130],tcg[46,105:130],col=46)

#sites with negative recovery rates
which(recov.rate<0)
which(recov.time<0)

