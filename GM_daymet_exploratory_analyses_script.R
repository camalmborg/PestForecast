#daymet maxtemp exploratory analyses

#load data:
load('maxtemp.RData')

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


##making months:
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

#load variable:
load('maxtemp.RData')
#load('mintemp.RData')
#load('precip.RData')
#load('solr.RData')
#load('vpd.RData')

#change this depending on which variable you are using:
var=maxtemp

#remove the 366 day:
varx<-list()
for (i in 1:nsites){
  varx[[i]]<-var[[i]][1:26,1:365]
}
#replace old data:
var<-varx
rm(varx)
#save(maxtemp, file='maxtemp.RData')

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

