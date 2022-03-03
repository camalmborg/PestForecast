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
#1/25/2021:
#load('vpdmean.RData')
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

###INTERACTIONS:------
#pdv.temp<-predistvar
#pdv.precip<-predistvar
pdv.vpd<-predistvar
#dwv.temp<-distwindvar
#dwv.precip<-distwindvar
dwv.vpd<-distwindvar

# #winter-spring group:
# predistvar<-varwintspr[,65:84]
# distwindvar<-varwintspr[,85:92]
#
#
### Plotting section for testing relationships #####
#
#join var to tcg.recov or just mags:
# predistvarmags<-as.data.frame(cbind(tcg.recov$ID,tcg.recov$lat,tcg.recov$lon,
#                        tcg.recov$steady,tcg.recov$mins,tcg.recov$steady,
#                        tcg.recov$mags,tcg.recov$colnum,tcg.recov$recov.rate,predistvar))
predistvarmags<-as.data.frame(cbind(mags,mins,colnum,predistvar))
pdvm16<-predistvarmags[predistvarmags$colnum==22,]
pdvm17<-predistvarmags[predistvarmags$colnum==23,]
# 
# distwindvarmags<-as.data.frame(cbind(tcg.recov$ID,tcg.recov$lat,tcg.recov$lon,
#                        tcg.recov$steady,tcg.recov$mins,tcg.recov$steady,
#                        tcg.recov$mags,tcg.recov$colnum,tcg.recov$recov.rate,distwindvar))
distwindvarmags<-as.data.frame(cbind(mags,mins,colnum,distwindvar))
dwvm16<-distwindvarmags[distwindvarmags$colnum==22,]
dwvm17<-distwindvarmags[distwindvarmags$colnum==23,]

###INTERACTIONS:
pdvm<-as.data.frame(cbind(mags,mins,colnum,pdv.temp,pdv.precip,pdv.vpd))
dwvm<-as.data.frame(cbind(mags,mins,colnum,dwv.temp,dwv.precip,dwv.vpd))

pdvm16<-pdvm[pdvm$colnum==22,]
dwvm16<-dwvm[dwvm$colnum==22,]
pdvm17<-pdvm[pdvm$colnum==23,]
dwvm17<-dwvm[dwvm$colnum==23,]

### GAM PLOT #####
#load libaries:
# library(lattice)
# library(latticeExtra)
# library(tactile)
# library(mgcv)

#Univariate:
#mo<-predistvarmags[,20]  #whatever month we are using
#mo<-distwindvarmags[,]
#mo<-pdvm16[,28]
#mo<-pdvm17[,28]
#mo<-dwvm16[,6]
#mo<-dwvm17[,13]

##Interactions (Multivariate):
mt<-pdvm16[,]
#mt<-pdvm17[,]
mp<-pdvm16[,]
#mp<-pdvm17[,]
mv<-pdvm16[,]
#mv<-pdvm17[,]

#make the gam:
#vardat = predistvarmags
#vardat = distwindvarmags
vardat = pdvm16
#vardat = pdvm17
#vardat = dwvm16
#vardat = dwvm17
var.gam <- gam(mags~s(mo), data = vardat)
#plot the gam:
varplot<-xyplot(mags ~ mo, data = vardat,
                   panel = function(x, y) {
                     ci<-predict(var.gam, se=T)
                     ci$lower<-ci$fit-qt(0.975,var.gam$df.null)*ci$se.fit
                     ci$upper<-ci$fit+qt(0.975,var.gam$df.null)*ci$se.fit
                     l.ci<-cbind(var.gam$model$mo,ci$fit,ci$lower,ci$upper)
                     l<-l.ci[order(l.ci[,1]),]
                     panel.ci(l[,1],l[,2],l[,4],l[,3],
                              fill="gray",alpha = 0.3)
                     panel.xyplot(x, y, pch=20,col="gray")
                     panel.lines(l[,1], l[,2],lty=1, col='black', lwd=1.5)
                   })
print(varplot)

summ<-summary(var.gam)
r2 <- summ$r.sq
print(r2)


#
#
### plotting with colnums: -----
#data:
# mo<-predistvarmags[,2]  #whatever month we are using
# #mo<-distwindvarmags[,6]
# plot(mo,mags,col=colnum,pch=20)

### Interactions