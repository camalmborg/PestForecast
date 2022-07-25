#script for grabbing daymet data based on insect phenology
#pre-larval "hatch" time: months april, may
#larval "feed" time: months june, july
#daymet variables: maxtemp, precip, vpd > average monthly values aggregated
#so that hatch= mean variable values for april-may, feed = mean values june-july

#MAKING THE DATASET-----
#to get data for 5 prior years 2011-2015:
prior5 <- c(17:21)
#to get data for disturbance window 2016-2017:
distwind <- c(22:23)
#to get whole time:
alltime <- c(1:26)
#number of sites and years and months:
nsite=5000
nsites = 1:5000
nyears = 1:26
ms = 1:12

x.p <- matrix(data=NA, nrow=5000)
for (p in 1:length(alltime)){
  x <- matrix(data = NA, nrow=5000, ncol=length(ms))
  for (m in 1:length(ms)){
    for (s in nsites){
      x[s,m] <- meanvar[[s]][ms[m],alltime[p]]
    }
  }
  x.p <- cbind(x.p,x)
  rm(x) #remove last loop
}

envar <- x.p[,2:((length(ms)*length(nyears))+1)] #remove NA column
rm(x.p) #remove redundant variable

#remove missing data rows:
var <- envar#[-missing,]
rm(envar) #remove redundant matrix

#var mags is the mag and recov rate data with
#the daymet data attached (eg. mean monthly precip)
#for number of years * months specified
#full set is all months (12) * all years (26) = 312 months
#varmags<-as.data.frame(cbind(tcg.recov, var))
#mags<-tcg.recov[,'mags']
#colnum<-tcg.recov[,'colnum']

#quick function for making sequences to extract monthly values:
seqfx<-function(x){
  seq(x,312,by=12) #x = month (1=jan, 2=feb, etc.)
}

#life cycle stage groups:
#april-may "pre-larval" hatching conditions
#june-july larval feeding conditions
aprilmay<-sort(c(seqfx(4),seqfx(5)))
junejuly<-sort(c(seqfx(6),seqfx(7)))

#winter months section:
dect<-c(seqfx(12))[17:23] #[17:23] takes 2011-2017 months
jant<-c(seqfx(1))[17:23]
febt<-c(seqfx(2))[17:23]
mart<-c(seqfx(3))[17:23]

#springsummer:
#prelarval:
varhatch<-var[,aprilmay]
#larval:
varfeed<-var[,junejuly]

#winter temps data:
wintermonths<-c(dect,jant,febt,mart)
varwint<-var[,wintermonths]
colnames(varwint)<-as.character(wintermonths)
#load distmags data:
distmagrecov<-read.csv("SM_distmagrecov_data.csv")
#make dataset:
varwintmags <- as.data.frame(cbind(distmagrecov$mags,distmagrecov$colnum,varwint))
colnames(varwintmags)<-c("mags", "colnum",
                    "dec2011","dec2012","dec2013","dec2014","dec2015","dec2016","dec2017",
                    "jan2011","jan2012","jan2013","jan2014","jan2015","jan2016","jan2017",
                    "feb2011","feb2012","feb2013","feb2014","feb2015","feb2016","feb2017",
                    "mar2011","mar2012","mar2013","mar2014","mar2015","mar2016","mar2017")




#grab seasonal months (2 monthly values * 26 years = 52 cols)
#grab each month of each year's season separately: varhatch or varfeed
mo1<-seq(1,ncol(varhatch),by=2)
mo2<-seq(2,ncol(varhatch),by=2)
mo3<-seq(1,ncol(varfeed),by=2)
mo4<-seq(2,ncol(varfeed),by=2)

#make containers for them:
v1<-varhatch[,mo1]
v2<-varhatch[,mo2]
v3<-varfeed[,mo3]
v4<-varfeed[,mo4]

#run loop to average each three-month season (whatever varseason is...)
vhstage<-matrix(NA, nrow=nsite, ncol=length(nyears))
vfstage<-matrix(NA, nrow=nsite, ncol=length(nyears))
#result is matrix with seasonal averages for each year:
#vhstage<-matrix(NA, nrow=nsite-(length(missing)), ncol=length(nyears))
#vfstage<-matrix(NA, nrow=nsite-(length(missing)), ncol=length(nyears))
for (i in 1:length(nyears)){
  vhall<-cbind(v1[,i],v2[,i])
  vfall<-cbind(v3[,i],v4[,i])
  vhmean<-apply(vhall,1,mean)
  vfmean<-apply(vfall,1,mean)
  vhstage[,i]<-vhmean
  vfstage[,i]<-vfmean
}

#lines for saving: t=temp, p=precip, v=vpd
vhp<-vhstage
vfp<-vfstage

#make data frames:
hatchallvar<-as.data.frame(cbind(mags,colnum,vht,vhp,vhv))
feedallvar<-as.data.frame(cbind(mags,colnum,vft,vfp,vfv))
hatchallvar<-as.data.frame(cbind(vht,vhp,vhv))
feedallvar<-as.data.frame(cbind(vft,vfp,vfv))
#write.csv(hf, "HF_Daymet_TPV_DATA.csv")
write.csv(hatchallvar,"hatch_daymet_allvar.csv")
write.csv(feedallvar,"feed_daymet_allvar.csv")

#fixing hatch/feed mess up:
hf.2<-hf
cols.2<-c(27:52,79:104,131:156)
feedallvar<-hf.2[,cols.2]
hatchallvar<-hf.2[-cols.2]

#separate 2016 and 2017 disturbance years:
hatch16<-hatchallvar[hatchallvar$colnum==22,]
hatch17<-hatchallvar[hatchallvar$colnum==23,]
feed16<-feedallvar[feedallvar$colnum==22,]
feed17<-feedallvar[feedallvar$colnum==23,]

#make big dataset with all vars for 2016 sites:
hf<-as.data.frame(cbind(vht,vft,vhp,vfp,vhv,vfv))
hf16<-hf[hf$colnum==22,]

### RUNNING THE GAMs-----
library(mgcv)

#variables to include:
mp<-hf16[,100]
mp2<-hf16[,101]
mp3<-hf16[,74]
mp4<-hf16[,75]
mv<-hf16[,126]
mv2<-hf16[,127]
mv3<-hf16[,152]
mv4<-hf16[,153]
mt<-hf16[,24]
mt2<-hf16[,50]
mt3<-hf16[,22]
mt4<-hf16[,23]
mt5<-hf16[,48]
mt6<-hf16[,49]

#data set:
# vardat = hatch16
# vardat = feed16
# vardat = hf16
vardat = varwintmags[varwintmags$colnum == 21,] #21=2015 dist, 22=2016 dist, 23=2017 dist


#var.gam <- gam(mags~s(mp)+s(mp2)+s(mp3)+s(mp4)+s(mv)+s(mv2)+s(mv3)+s(mv4), data = vardat)
var.gam <- gam(mags ~ s(dec2011), data=vardat)
summ<-summary(var.gam)
r2 <- summ$r.sq
print(r2)

