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

plot(avgyrtemp[1,], type='l', ylim=c(8,20))
for (i in 1:nsites){
  lines(avgyrtemp[i,],type='l',col=i)
}

# for (i in 1:length(metyears)){
#   avgyrtemp[,i]<-mean(maxtemp[[1]][i,])
# }



