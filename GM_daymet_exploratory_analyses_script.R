### Daymet Exploratory Analyses Automation ###
#
# MEANVAR TO VARIABLE #####
#meanvar will load list with  mean of:
  #precip<-meanvar
  #maxtemp<-meanvar
  #mintemp<-meanvar
  #solr<-meanvar
  #vpd<-meanvar
  #rm(meanvar) #gets rid of redundant variable

#
### Extracting daymet data from list #####
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

var <- x.p[,2:((length(ms)*length(nyears))+1)] #remove NA column
rm(x.p) #remove redundant variable
