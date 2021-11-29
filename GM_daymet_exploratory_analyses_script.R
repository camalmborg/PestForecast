#daymet maxtemp exploratory analyses

#load data:
#load('maxtemp.RData')

#things I might need for loops if not in environment:
#nsites=5000
#metyears=1:26
#doy <- dm[[1]]$data$yday
# for (i in nsites){
#   metyr= dm[[i]]$data$year
#   metyears=unique(metyr)
# }


###First Maxtemp Run:-----------
#remove the 366 day:
maxtempx<-list()
for (i in 1:nsites){
  maxtempx[[i]]<-maxtemp[[i]][1:26,1:365]
}
#replace old data:
maxtemp<-maxtempx
#save(maxtemp, file='maxtemp.RData')

avgyrtemp<-matrix(NA, nrow=nsites, ncol=length(metyears))
for (i in 1:nsites){
  for (j in 1:length(metyears)){
    avgyrtemp[i,j]<-mean(maxtemp[[i]][j,])
  }
}

##caluculating seasons temps:
#seasonal breaks:
winter<-c(doy[1:59],doy[336:365])
spring<-c(doy[60:151])
summer<-c(doy[152:244])
fall<-c(doy[245:335])
#make one object:
seasons<-c(winter,spring,summer,fall)

avgsmxtemp<-list()
for (i in 1:nsites){
  avgsmxtemp[[i]]<-matrix(NA, nrow=4, ncol=length(metyears))
  for (j in 1:length(metyears)){
    for (s in 1:4){
      if (s==1){
        avgsmxtemp[[i]][s,j]<-mean(maxtemp[[i]][j,winter])
      }
      if (s==2){
        avgsmxtemp[[i]][s,j]<-mean(maxtemp[[i]][j,spring])
      }
      if (s==3){
        avgsmxtemp[[i]][s,j]<-mean(maxtemp[[i]][j,summer])
      }
      if (s==4){
        avgsmxtemp[[i]][s,j]<-mean(maxtemp[[i]][j,fall])
      }
    }
  }
}


#plot a sample (maxtemp):
samp<-sample(nsites,500)
plot(avgsmxtemp[[1]][1,], tcg.june[1,], pch=20, xlim=c(-3,8),ylim=c(-0.1,0.5))
for (s in 1:nsites){
  points(avgsmxtemp[[s]][1,], tcg.june[s,], pch=20)
}



# plot(avgsmxtemp[[1]][1,], type='l')
# for (i in 1:nsites){
#   lines(avgsmxtemp[[i]][1,],type='l',col=i)
# }




##making seasons:------------
#seasonal breaks:
winter<-c(doy[1:59],doy[336:365])
spring<-c(doy[60:151])
summer<-c(doy[152:244])
fall<-c(doy[245:335])
#make one object:
#seasons<-c(winter,spring,summer,fall)


##making months:-----------
#monthly breaks:
jan<-c(doy[1:31])
feb<-c(doy[32:59])
mar<-c(doy[60:90])
apr<-c(doy[91:120])
may<-c(doy[121:151])
jun<-c(doy[152:181])
jul<-c(doy[182:212])
aug<-c(doy[213:243])
sep<-c(doy[244:273])
oct<-c(doy[274:304])
nov<-c(doy[305:334])
dec<-c(doy[335:365])

#------make useful for all variables:--------------

###-----
#loading the variable and dropping empty matrix column
#load variable:
#load('maxtemp.RData')
#load('mintemp.RData')
load('precip.RData')
#load('solr.RData')
#load('vpd.RData')

#change this depending on which variable you are using:
var=precip

#remove the 366 day:
varx<-list()
for (i in nsites){
  varx[[i]]<-var[[i]][1:26,1:365]
}
#replace old data:
var<-varx
rm(varx)
#save(maxtemp, file='maxtemp.RData')

###FOR AVERAGES, SEASONAL-----
# avgyrvar<-matrix(NA, nrow=nsites, ncol=length(metyears))
# for (i in 1:nsites){
#   for (j in 1:length(metyears)){
#     avgyrvar[i,j]<-mean(var[[i]][j,])
#   }
# }

avgsvar<-list()
for (i in 1:nsites){
  avgsvar[[i]]<-matrix(NA, nrow=4, ncol=length(metyears))
  for (j in 1:length(metyears)){
    for (s in 1:4){
      if (s==1){
        avgsvar[[i]][s,j]<-mean(var[[i]][j,winter])
      }
      if (s==2){
        avgsvar[[i]][s,j]<-mean(var[[i]][j,spring])
      }
      if (s==3){
        avgsvar[[i]][s,j]<-mean(var[[i]][j,summer])
      }
      if (s==4){
        avgsvar[[i]][s,j]<-mean(var[[i]][j,fall])
      }
    }
  }
}

###PLOTS-----
#plot a sample with tcg as x:
samp<-sample(nsites,500)
xl=c(-3,7.5) #maxtemp, winter
yl=c(min(tcg.june,na.rm=T),max(tcg.june,na.rm=T))
plot(avgsvar[[1]][1,], tcg.june[1,], pch=20, ylim=yl,xlim=xl)
for (s in samp){
  points(avgsvar[[s]][1,], tcg.june[s,], pch=20)
}

#with condition scores as x:
samp<-sample(nsites,500)
xl=c(20,30) #maxtemp, winter
yl=c(min(cond.june,na.rm=T),max(cond.june,na.rm=T))
plot(avgsvar[[1]][3,], cond.june[1,], pch=20, ylim=yl, xlim=xl)
for (s in samp){
  points(avgsvar[[s]][3,], cond.june[s,], pch=20)
}

#with disturbance mag as x:
#samp<-sample(nsites,500)
xl=c(-0.5,6) #maxtemp, winter
yl=c(0,max(mags,na.rm=T))
plot(avgsvar[[1]][1,23], mags[1], pch=20, ylim=yl, xlim=xl)
for (s in 1:nsites){
  points(avgsvar[[s]][1,23], mags[s], pch=20)
}





#-----
#
#-----
###FOR MAX VALUES, MONTHLY:-----
maxvar<-list()
for (i in nsites){
  maxvar[[i]]<-matrix(NA, nrow=12, ncol=length(metyears))
  for (j in 1:length(metyears)){
    for (s in 1:12){
      if (s==1){
        maxvar[[i]][s,j]<-max(var[[i]][j,jan])
      }
      if (s==2){
        maxvar[[i]][s,j]<-max(var[[i]][j,feb])
      }
      if (s==3){
        maxvar[[i]][s,j]<-max(var[[i]][j,mar])
      }
      if (s==4){
        maxvar[[i]][s,j]<-max(var[[i]][j,apr])
      }
      if (s==5){
        maxvar[[i]][s,j]<-max(var[[i]][j,may])
      }
      if (s==6){
        maxvar[[i]][s,j]<-max(var[[i]][j,jun])
      }
      if (s==7){
        maxvar[[i]][s,j]<-max(var[[i]][j,jul])
      }
      if (s==8){
        maxvar[[i]][s,j]<-max(var[[i]][j,aug])
      }
      if (s==9){
        maxvar[[i]][s,j]<-max(var[[i]][j,sep])
      }
      if (s==10){
        maxvar[[i]][s,j]<-max(var[[i]][j,oct])
      }
      if (s==11){
        maxvar[[i]][s,j]<-max(var[[i]][j,nov])
      }
      if (s==12){
        maxvar[[i]][s,j]<-max(var[[i]][j,dec])
      }
    }
  }
}



###FOR MIN VALUES, MONTHLY:-----
minvar<-list()
for (i in nsites){
  minvar[[i]]<-matrix(NA, nrow=12, ncol=length(metyears))
  for (j in 1:length(metyears)){
    for (s in 1:12){
      if (s==1){
        minvar[[i]][s,j]<-min(var[[i]][j,jan])
      }
      if (s==2){
        minvar[[i]][s,j]<-min(var[[i]][j,feb])
      }
      if (s==3){
        minvar[[i]][s,j]<-min(var[[i]][j,mar])
      }
      if (s==4){
        minvar[[i]][s,j]<-min(var[[i]][j,apr])
      }
      if (s==5){
        minvar[[i]][s,j]<-min(var[[i]][j,may])
      }
      if (s==6){
        minvar[[i]][s,j]<-min(var[[i]][j,jun])
      }
      if (s==7){
        minvar[[i]][s,j]<-min(var[[i]][j,jul])
      }
      if (s==8){
        minvar[[i]][s,j]<-min(var[[i]][j,aug])
      }
      if (s==9){
        minvar[[i]][s,j]<-min(var[[i]][j,sep])
      }
      if (s==10){
        minvar[[i]][s,j]<-min(var[[i]][j,oct])
      }
      if (s==11){
        minvar[[i]][s,j]<-min(var[[i]][j,nov])
      }
      if (s==12){
        minvar[[i]][s,j]<-min(var[[i]][j,dec])
      }
    }
  }
}



###FOR MEAN VALUES, MONTHLY:-----
meanvar<-list()
for (i in nsites){
  meanvar[[i]]<-matrix(NA, nrow=12, ncol=length(metyears))
  for (j in 1:length(metyears)){
    for (s in 1:12){
      if (s==1){
        meanvar[[i]][s,j]<-mean(var[[i]][j,jan])
      }
      if (s==2){
        meanvar[[i]][s,j]<-mean(var[[i]][j,feb])
      }
      if (s==3){
        meanvar[[i]][s,j]<-mean(var[[i]][j,mar])
      }
      if (s==4){
        meanvar[[i]][s,j]<-mean(var[[i]][j,apr])
      }
      if (s==5){
        meanvar[[i]][s,j]<-mean(var[[i]][j,may])
      }
      if (s==6){
        meanvar[[i]][s,j]<-mean(var[[i]][j,jun])
      }
      if (s==7){
        meanvar[[i]][s,j]<-mean(var[[i]][j,jul])
      }
      if (s==8){
        meanvar[[i]][s,j]<-mean(var[[i]][j,aug])
      }
      if (s==9){
        meanvar[[i]][s,j]<-mean(var[[i]][j,sep])
      }
      if (s==10){
        meanvar[[i]][s,j]<-mean(var[[i]][j,oct])
      }
      if (s==11){
        meanvar[[i]][s,j]<-mean(var[[i]][j,nov])
      }
      if (s==12){
        meanvar[[i]][s,j]<-mean(var[[i]][j,dec])
      }
    }
  }
}





###FOR MEDIAN VALUES, MONTHLY:-----
medvar<-list()
for (i in nsites){
  medvar[[i]]<-matrix(NA, nrow=12, ncol=length(metyears))
  for (j in 1:length(metyears)){
    for (s in 1:12){
      if (s==1){
        medvar[[i]][s,j]<-median(var[[i]][j,jan])
      }
      if (s==2){
        medvar[[i]][s,j]<-median(var[[i]][j,feb])
      }
      if (s==3){
        medvar[[i]][s,j]<-median(var[[i]][j,mar])
      }
      if (s==4){
        medvar[[i]][s,j]<-median(var[[i]][j,apr])
      }
      if (s==5){
        medvar[[i]][s,j]<-median(var[[i]][j,may])
      }
      if (s==6){
        medvar[[i]][s,j]<-median(var[[i]][j,jun])
      }
      if (s==7){
        medvar[[i]][s,j]<-median(var[[i]][j,jul])
      }
      if (s==8){
        medvar[[i]][s,j]<-median(var[[i]][j,aug])
      }
      if (s==9){
        medvar[[i]][s,j]<-median(var[[i]][j,sep])
      }
      if (s==10){
        medvar[[i]][s,j]<-median(var[[i]][j,oct])
      }
      if (s==11){
        medvar[[i]][s,j]<-median(var[[i]][j,nov])
      }
      if (s==12){
        medvar[[i]][s,j]<-median(var[[i]][j,dec])
      }
    }
  }
}


#SAVING DATA FOR EASY REUSE:-----
# save(maxvar, file='maxtempmax.RData')
# save(minvar, file='maxtempmin.RData')
# save(meanvar, file='maxtempmean.Rdata')
# save(medvar, file='maxtempmed.Rdata')
# 
# save(maxvar, file='mintempmax.RData')
# save(minvar, file='mintempmin.RData')
# save(meanvar, file='mintempmean.Rdata')
# save(medvar, file='mintempmed.Rdata')
# 
# save(maxvar, file='precipmax.RData')
# save(minvar, file='precipmin.RData')
# save(meanvar, file='precipmean.Rdata')
# save(medvar, file='precipmed.Rdata')
# 

