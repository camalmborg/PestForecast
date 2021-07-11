#load data for all test sites:
time=1:130
NT=length(time)
geoID=1:2500
gm.data<-list()
for(i in geoID){
  file<-paste0("TEST_monthly_data/GM_",geoID[i],"_mcs.csv")
  gm.data[[i]]<-read.csv(file)
}

#load condition score means for each site into matrix:
obsm<-matrix(NA,nrow=length(geoID),ncol=130)
for (i in geoID){
  obs[i,]<-gm.data[[i]]$score_mean
}


#load annual condition scores:
data.an<-read.csv("2021_07_09_sample_score_mean_ANNUAL_50sites.csv")
gm.data.an<-data.an[,2:27]

#load NLCD data
nlcd.dat<-read.csv("2020_07_11_sample_nlcd_500.csv")
landcover<-nlcd.dat[,5]
treecover<-nlcd.dat[,8]

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

  
#find the non-GM disturb condition "steady state"
steady<-apply(obs[,90:105],1,mean,na.rm=T)
  

#calculating slope for each recovery rate (loop):
recov.rate<-matrix(NA,nrow=50,ncol=1)
ind<-list()
slope<-list()
mins<-vector()
colnum<-vector()
for (i in geoID){
  mins[i]<-min(obs[i],na.rm=T)            #grabs min forest condition score for each site
  colnum[i]<-which(obs[i]==mins[i])      #grabs time step at which min score appears
  #recov<-obs[i,colnum[i,]:min(obs[i,colnum[i,]+15],130)  #+10=2 year recovery period
  #ind<-1:length(recov)                      #grabs length of recov rate (now 11)
  #slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  #recov.rate[i,]<-slope[[i]]$coefficients["ind"] #stores recovery rate
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

###THE MONTHLY DATA VERSION:
#load TCG mean data:
tcgmeans<-"2020_07_10_500_sample_tcg_mean_MONTHLY.csv"
tcgtable<-read.csv(tcgmeans)
tcg<-as.matrix(tcgtable[,2:131])

#find the non-GM disturb condition "steady state"
steady<-apply(tcg[,90:110],1,mean,na.rm=T)
steadyall<-apply(tcg[,1:110],1,mean,na.rm=T)

#monthly steady states
#maymeans<-apply(tcg[,seq(1,ncol(tcg),by=5)],1,mean,na.rm=T)
steadymonths<-matrix(NA,nrow=4999,ncol=5)
months<-5
for (i in 1:months){
  steadymonths[,i]<-apply(tcg[,seq(i,ncol(tcg),by=5)],1,mean,na.rm=T)
}

#the TCG 2016-onward MONTHLY version:
#ind<-list()
recov.rate<-vector()
slope<-list()
mins<-vector()
colnum<-vector()
recovcol<-vector()
for (i in geoID){
  mins[i]<-min(tcg[i,110:118],na.rm=T)          #grabs min forest condition score for each site
  colnum[i]<-which(tcg[i,]==mins[i])           #grabs time step at which min score appears in disturbance window
  if(is.na(colnum[i] + which(tcg[i,colnum[i]:130]>=steady[i])[1])){
    recovcol[i]<-colnum[i] + 2
  } else {
    recovcol[i]<-colnum[i] + which(tcg[i,colnum[i]:130]>=steady[i])[1]
  }
  if(recovcol[i]>130){
    recov<-tcg[i,colnum[i]:130]
  } else {
    recov<- tcg[i,colnum[i]:recovcol[i]]
  }
  #recov <- tcg[i,colnum[i]:recovcol[i]]      
  ind<-1:length(recov)                      #grabs length of recov rate   
  slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
}



###ANNUAL DATA VERSION:
#tcgmeans<-"2021_07_08_tcg_mean_annual_50sites.csv"
#tcgtable<-read.csv(tcgmeans)
#tcg<-as.matrix(tcgtable[1:50,2:27])

#steady<-apply(tcg[,19:21],1,mean,na.rm=T)
#steadyall<-apply(tcg[,1:26],1,mean,na.rm=T)

#the TCG 2016-onward version:
# recov.rate<-vector()
# slope<-list()
# #ind<-list()
# mins<-vector()
# colnum<-vector()
# recovcol<-vector()
# for (i in geoID){
#   mins[i]<-min(tcg[i,22:23],na.rm=T)          #grabs min forest condition score for each site
#   colnum[i]<-which(tcg[i,]==mins[i])           #grabs time step at which min score appears in disturbance window
#   if(is.na(colnum[i] + which(tcg[i,colnum[i]:26]>=steady[i])[1])){
#     recovcol[i]<-colnum[i] + 2
#   } else {
#     recovcol[i]<-colnum[i] + which(tcg[i,colnum[i]:26]>=steady[i])[1]
#   }
#   #recovcol[i,]<-which(tcg[i,colnum[i,]:130]>steady[i])[1]
#   #recovcol[i,]<-colnum[i,]+which(tcg[i,colnum[i,]:130]>steady[i])[1]
#   recov <- tcg[i,colnum[i]:recovcol[i]]      
#   ind<-1:length(recov)                      #grabs length of recov rate   
#   slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
#   recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
# }



#magnitudes:
mags<-steady-mins
recov.time<-mags/recov.rate

tcg.recov.mx<-cbind(steadyall,steady,colnum,mins,mags,recov.rate,recov.time)
tcg.recov<-as.data.frame(tcg.recov.mx)
tcg.recov$NLCD <- factor(landcover, levels=c(41,42,43,90,21),
                              labels=c("Deciduous", "Evergreen","Mixed","Woody Wetland","Developed:open space"))
tcg.recov$percentcover<-treecover
#colnames(tcg.recov)<-c("Steady State (All Years)","Steady State 2012-2015","Column Number","Minimum TCG Value","Disturbance Magnitude","Recovery Rate","Recovery Time","NLCD Type","Percent Tree Cover")


if (FALSE){
#previous year's greenness analysis:
endprevyr<-colnum-1
startprevyr<-endprevyr-5
prevyr<-vector()
for (i in geoID){
  prevyr[i]<-(mean(na.omit(tcg[i,startprevyr[i]:endprevyr[i]])))
}

#reduction in greenness from previous year using magnitude:
defol<-prevyr-mags

#make a nice table of it all:
#recov.view.mx<-cbind(steadyall,steady,prevyr,mags,defol,mins,recov.rate,recov.time)
#recov.view<-as.data.frame(recov.view.mx)
#colnames(recov.view)<-c("steadyall","steady","prevyr","mags","defol","mins","recov.rate","recov.time")

recov.view.mx<-cbind(steadyall,steady,prevyr,mags,defol,recov.rate,recov.time)
recov.view<-as.data.frame(tcg.recov.mx)
recov.view$NLCD <- factor(landcover, levels=c(41,42,43,90,21),
                         labels=c("Deciduous", "Evergreen","Mixed","Woody Wetland","Developed:open space"))
recov.view$percentcover<-treecover
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
recov.time.lm<-lm(recov.time~mags,tcg.recov)
mins.lm<-lm(mins ~ steady, data = tcg.recov)
magrecov.lm<-lm(recov.rate ~ mags, data = tcg.recov)
magtree.lm<-lm(mags~percentcover,data=tcg.recov)
#defolrecov.lm<-lm(recov.rate~defol,recov.view)
#prevyr.lm<-lm(recov.rate~prevyr,recov.view)


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

plot(time[110:122],tcg[1,110:122],type="l",ylim=c(-0.01,0.4))
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

#plotting steady states for each site
plot(steadymonths[1,],type="l",ylim=c(0.125,0.32))
for (i in 1:50){
  lines(steadymonths[i,],col=i)
}