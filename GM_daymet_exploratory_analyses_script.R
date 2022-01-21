### Daymet Exploratory Analyses Automation ###
#
# MEANVAR TO VARIABLE #####
#load("")
#meanvar will load list with  mean of:
  #precip<-meanvar
  #maxtemp<-meanvar
  #mintemp<-meanvar
  #solr<-meanvar
  #vpd<-meanvar
  #rm(meanvar) #gets rid of redundant variable

#
### Extracting daymet data from list #####

#load desired variable from storage:
#meanvar<-''
#maxvar
#minvar
#medvar
#etc.

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
varmags<-as.data.frame(cbind(tcg.recov, var))
#
#
### Make testing windows for daymet data #####

#quick function for making sequences to extract monthly values:
seqfx<-function(x){
  seq(x,312,by=12) #x = month (1=jan, 2=feb, etc.)
}

#seasons:
winter<-sort(c(seqfx(12),seqfx(1),seqfx(2)))
spring<-sort(c(seqfx(3),seqfx(4),seqfx(5)))
summer<-sort(c(seqfx(6),seqfx(7),seqfx(8)))
fall<-sort(c(seqfx(9),seqfx(10),seqfx(11)))

###combo seasons
#april-august window:
sprsum<-sort(c(seqfx(4),seqfx(5),seqfx(6),seqfx(7),seqfx(8)))
#december-march:
wintspr<-sort(c(seqfx(12),seqfx(1),seqfx(2),seqfx(3)))
#still need to figure out how to make december work with these, 
#not in "order" since december is month 12, jan is month 1

###environmental variables for each season
#spring-summer season:
varsprsum<-var[,sprsum]
#winter-spring season:
varwintspr<-var[,wintspr]

###testing windows ###   #can cbind to tcg.recov for plots
#spring-summer group:
predistvar<-varsprsum[,81:105]
distwindvar<-varsprsum[,106:115]
#recovvar<-varsprsum[,116:130]

#winter-spring group:
#predistvar<-varwintspr[]
#distwindvar<-varwintspr[]
#
#
### Plotting section for testing relationships #####
