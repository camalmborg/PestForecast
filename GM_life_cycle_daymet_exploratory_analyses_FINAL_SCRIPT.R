#script for grabbing daymet data based on insect phenology
#pre-larval "hatch" time: months april, may
#larval "feed" time: months june, july
#daymet variables: maxtemp, precip, vpd > average monthly values aggregated
#so that hatch= mean variable values for april-may, feed = mean values june-july

#to get data for 5 prior years 2011-2015:
prior5 <- c(17:21)
#to get data for disturbance window 2016-2017:
distwind <- c(22:23)
#to get whole time:
alltime <- c(1:26)
#number of sites and years and months:
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
var <- envar[-missing,]
rm(envar) #remove redundant matrix

#var mags is the mag and recov rate data with
#the daymet data attached (eg. mean monthly precip)
#for number of years * months specified
#full set is all months (12) * all years (26) = 312 months
#varmags<-as.data.frame(cbind(tcg.recov, var))
mags<-tcg.recov[,'mags']
colnum<-tcg.recov[,'colnum']

#quick function for making sequences to extract monthly values:
seqfx<-function(x){
  seq(x,312,by=12) #x = month (1=jan, 2=feb, etc.)
}

#life cycle stage groups:
#april-may "pre-larval" hatching conditions
#june-july larval feeding conditions
aprilmay<-sort(c(seqfx(4),seqfx(5)))
junejuly<-sort(c(seqfx(6),seqfx(7)))

#prelarval:
varhatch<-var[,aprilmay]
#larval:
varfeed<-var[,junejuly]

#grab seasonal months (2 monthly values * 26 years = 52 cols)
#grab each month of each year's season separately: varhatch or varfeed
mo1<-seq(1,ncol(varhatch),by=2)
mo2<-seq(2,ncol(varhatch),by=2)
mo3<-seq(1,ncol(varfeed),by=2)
mo4<-seq(2,ncol(varfeed),by=2)

#make containers for them:
v1<-varhatch[,mo1]
v2<-varhatch[,mo2]
v3<-varhatch[,mo1]
v4<-varhatch[,mo2]

#run loop to average each three-month season (whatever varseason is...)
#result is matrix with seasonal averages for each year:
vhstage<-matrix(NA, nrow=nsite-(length(missing)), ncol=length(nyears))
vfstage<-matrix(NA, nrow=nsite-(length(missing)), ncol=length(nyears))
for (i in 1:length(nyears)){
  vhall<-cbind(v1[,i],v2[,i])
  vfall<-cbind(v3[,i],v4[,i])
  vhmean<-apply(vhall,1,mean)
  vfmean<-apply(vfall,1,mean)
  vhstage[,i]<-vhmean
  vfstage[,i]<-vfmean
}

#lines for saving: t=temp, p=precip, v=vpd
vhv<-vhstage
vfv<-vfstage

#make data frames:
hatchallvar<-as.data.frame(cbind(mags,colnum,vht,vhp,vhv))
feedallvar<-as.data.frame(cbind(mags,colnum,vft,vfp,vfv))
#separate 2016 and 2017 disturbance years:
hatch16<-hatchallvar[hatchallvar$colnum==22,]
hatch17<-hatchallvar[hatchallvar$colnum==23,]
feed16<-feedallvar[feedallvar$colnum==22,]
feed17<-feedallvar[feedallvar$colnum==23,]

#make big dataset with all vars for 2016 sites:
hf<-as.data.frame(cbind(mags,colnum,vht,vft,vhp,vfp,vhv,vfv))
hf16<-hf[hf$colnum==22,]

### RUNNING THE GAMs
library(mgcv)

#variables to include:
mp<-hf16[,100]
mp2<-hf16[,101]
mp3<-hf16[,102]
mp4<-hf16[,76]
mv<-hf16[,152]
mv2<-hf16[,153]
mv3<-hf16[,128]
#mv4<-hf16[,154]
mt<-hf16[,24]
mt2<-hf16[,50]
mt3<-hf16[,22]
mt4<-hf16[,23]
mt5<-hf16[,48]
mt6<-hf16[,49]

#data set:
# vardat = hatch16
# vardat = feed16
vardat = hf16

var.gam <- gam(mags~s(mp)+s(mp2)+s(mp3)+s(mp4)+s(mv)+s(mv2)+s(mv3)+s(mt)+s(mt2)+s(mt3)+s(mt4)+s(mt5)+s(mt6), data = vardat)
summ<-summary(var.gam)
r2 <- summ$r.sq
print(r2)

print('git works now?')
