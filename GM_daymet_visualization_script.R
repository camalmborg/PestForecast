#loading libraries:
library(ggplot2)
library(lattice)
library(latticeExtra)
library(tactile)
library(mgcv)

###PLOTS-----
# #plot a sample with tcg as x:
# samp<-sample(nsites,500)
# xl=c(1,25) #maxtemp max
# yl=c(min(tcg.june,na.rm=T),max(tcg.june,na.rm=T))
# plot(maxvar[[1]][1,], tcg.june[1,], pch=20, ylim=yl, xlim=xl)
# for (s in samp){
#   points(maxvar[[s]][1,], tcg.june[s,], pch=20)
# }

#with condition scores as x:
#samp<-sample(nsites,500)
# xl=c(5,33) #maxtemp max
# yl=c(min(cond.june,na.rm=T),max(cond.june,na.rm=T))
# plot(maxvar[[1]][3,], cond.june[1,], pch=20, ylim=yl, xlim=xl)
# for (s in nsites){
#   points(maxvar[[s]][3,], cond.june[s,], pch=20)
# }

#load desired var:
#load('maxtempmean.RData')

#with disturbance mag as x:
#samp<-sample(nsites,500)
# xl=c(22,50) 
# yl=c(0,max(mags,na.rm=T))
# plot(meanvar[[1]][5,23], mags[1], pch=20, ylim=yl, xlim=xl)
# for (s in nsites){
#   points(maxvar[[s]][5,23], mags[s], pch=20)
# }

###PLOTTING:-----------------------------------------------

#prepare data for plotting:
#code to extract each year's values
ms <- c(4:8) 
x <- matrix(data = NA, nrow=5000, ncol=5)
for (m in 1:length(ms)){ 
  for (s in nsites){
    x[s,m] <- meanvar[[s]][ms[m],17]
  }
}
x.17 <- x
#x.18 <- x
#x.19 <- x
#x.20 <- x
#x.21 <- x
#such that eg. x.17[,2]=May 2011

#Cmakes big matrix with 5 months per 5 prior years for precip data (Apr-Aug, 2011-2015)
prior5 <- c(17:21)

x.p <- matrix(data=NA, nrow=5000)
for (p in 1:length(prior5)){
  x <- matrix(data = NA, nrow=5000, ncol=5)
  for (m in 1:length(ms)){
    for (s in nsites){
      x[s,m] <- meanvar[[s]][ms[m],prior5[p]]
    }
  }
  x.p <- cbind(x.p,x)
}

priorprecip<-x.p[,2:26]
colnames(priorprecip)<-c('apr11','may11','jun11','jul11','aug11',
                         'apr12','may12','jun12','jul12','aug12',
                         'apr13','may13','jun13','jul13','aug13',
                         'apr14','may14','jun14','jul14','aug14',
                         'apr15','may15','jun15','jul15','aug15')

#doing this with 2011 data to start:
precipmags<-as.data.frame(cbind(mags, x.17))
colnames(precipmags)<- c('mags', 'apr', 'may', 'jun', 'jul', 'aug')
may.lm <- lm(mags~may, data = precipmags)
may.gam <- gam(mags~may, data = precipmags)

#Plotting the precip change over 5 years Apr-Aug:
# plot(priorprecip[1,], type ='l', ylim=c(0,max(priorprecip)), xlim=c(0,25))
# for (i in nsites){
#   lines(priorprecip[i,],)
# }

pmay17<-xyplot(mags ~ may, data = precipmags,
                     panel = function(x, y) {
                       ci<-predict(may.lm, interval="confidence")
                       upper<-ci[,3]
                       lower<-ci[,2]
                       panel.ci(x, y, upper, lower, fill="gray48",alpha = 0.4)
                       panel.xyplot(x, y, pch=16,col="black")
                       panel.abline(lm(y ~ x))
                       # summ<-summary(may.lm)
                       # r2 <- summ$adj.r.squared
                       # f <- summ$fstatistic
                       # p <- pf(f[1],f[2],f[3],lower.tail=F)
                       # panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                       #             )#x=0.13,y=0.09,cex=0.75)
                       # panel.text(labels = bquote(italic(p)==.(format(p,digits = 3))),
                      #              x=0.13,y=0.07,cex=0.75)
                      }#,
                     # ylab="Recovery Rate",
                     # xlab="TCG Steady State (2012-2015)",
)
print(pmay17)
