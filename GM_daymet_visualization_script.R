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

#prepare data for plotting:-----
#code to extract each year's values
ms <- c(4:8)
x <- matrix(data = NA, nrow=5000, ncol=5)
for (m in 1:length(ms)){
  for (s in nsites){
    x[s,m] <- meanvar[[s]][ms[m],17]
  }
}
# x.17 <- x
# #x.18 <- x
# #x.19 <- x
# #x.20 <- x
# #x.21 <- x
# #such that eg. x.17[,2]=May 2011

#Makes big matrix with 5 months per 5 prior years for precip data (Apr-Aug, 2011-2015)
prior5 <- c(17:21)

nsites.5k=1:5000
x.p <- matrix(data=NA, nrow=5000)
for (p in 1:length(prior5)){
  x <- matrix(data = NA, nrow=5000, ncol=5)
  for (m in 1:length(ms)){
    for (s in nsites.5k){
      x[s,m] <- meanvar[[s]][ms[m],prior5[p]]
    }
  }
  x.p <- cbind(x.p,x)
}

#priorprecip<-x.p[,2:26]
priorvpd<-x.p[,2:26]
# colnames(priorprecip)<-c('mags','apr11','may11','jun11','jul11','aug11',
#                          'apr12','may12','jun12','jul12','aug12',
#                          'apr13','may13','jun13','jul13','aug13',
#                          'apr14','may14','jun14','jul14','aug14',
#                          'apr15','may15','jun15','jul15','aug15')

#doing this with 2011 data to start:
#precipmags<-as.data.frame(cbind(mags, priorprecip))
vpdmags<-as.data.frame(cbind(magsdat$mags,priorvpd))
colnames(vpdmags)<-c('mags','apr11','may11','jun11','jul11','aug11',
                         'apr12','may12','jun12','jul12','aug12',
                         'apr13','may13','jun13','jul13','aug13',
                         'apr14','may14','jun14','jul14','aug14',
                         'apr15','may15','jun15','jul15','aug15')
#may.lm <- lm(mags~may, data = precipmags)
#may.gam <- gam(mags~s(may), data = precipmags)
# may.loess <- loess(mags~may, data=precipmags)
# loess.df<-cbind(may.loess$x, may.loess$fitted)
# l.df<-loess.df[order(loess.df[,1]),]

#Plotting the precip change over 5 years Apr-Aug:
plot(priorvpd[1,], type='l', ylim=c(0,max(na.omit(priorvpd))), xlim=c(0,25))
for (i in nsites){
  lines(priorvpd[i,],)
}

# mos<-1:25
# yr<-c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)
# pprecip<-apply(priorprecip,2,mean)
# #plot(mos,pprecip)
# precipmeans<-as.data.frame(cbind(yr,mos,pprecip))

#lattice plot example:-----
# pmay17<-xyplot(mags ~ may, data = precipmags,
#                      panel = function(x, y) {
#                          # ci<-predict(may.loess, se=T)
#                          # ci$lower<-ci$fit-qt(0.975,ci$df)*ci$se
#                          # ci$upper<-ci$fit+qt(0.975,ci$df)*ci$se
#                          # lloess<-cbind(may.loess$x,ci$fit,ci$lower,ci$upper)
#                          # l<-lloess[order(lloess[,1]),]
#                            ci<-predict(may.gam, se=T)
#                            ci$lower<-ci$fit-qt(0.975,may.gam$df.null)*ci$se.fit
#                            ci$upper<-ci$fit+qt(0.975,may.gam$df.null)*ci$se.fit
#                            l.ci<-cbind(may.gam$model$may,ci$fit,ci$lower,ci$upper)
#                            l<-l.ci[order(l.ci[,1]),]
#                        panel.ci(l[,1],l[,2],l[,4],l[,3],
#                                 fill="royalblue1",alpha = 0.3)
#                        panel.xyplot(x, y, pch=20,col="royalblue3")
#                        panel.lines(l[,1], l[,2],lty=1, col='black', lwd=1.5)
#                        #panel.lines(l[,1], l[,3], lty=2, col='black')
#                        #panel.lines(l[,1], l[,4], lty=2, col='black')
#                        summ<-summary(may.lm)
#                        r2 <- summ$adj.r.squared
#                        f <- summ$fstatistic
#                        p <- pf(f[1],f[2],f[3],lower.tail=F)
#                        panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
#                                     x=1.2,y=-0.1,cex=0.75)
#                        panel.text(labels = bquote(italic(p)==.(format(p,digits = 3))),
#                                     x=1.2,y=-0.15,cex=0.75)
#                       }#,
#                      # ylab="Disturbance Magnitude (TCG)",
#                      # xlab="Mean Temperature",
# )
# print(pmay17)


###FINAL PLOTS:--------------------------------------------------

#prepare inputs:
#input data is precipmags (or whatever variable-mags)
#mo will be the month we are using for the plot
#columns for reference:
#  1'mags',2'apr11',3'may11',4'jun11',5'jul11',6'aug11',
#  7'apr12',8'may12',9'jun12',10'jul12',11'aug12',
# 12'apr13',13'may13',14'jun13',15'jul13',16'aug13',
# 17 'apr14',18'may14',19'jun14',20'jul14',21'aug14',
# 22'apr15',23'may15',24'jun15',25'jul15',26'aug15')

mo<-vpdmags[,25]  #whatever column from precipmags we are using
vpd.gam <- gam(mags~s(mo), data = vpdmags)
  
vpdplot<-xyplot(mags ~ mo, data = vpdmags,
               panel = function(x, y) {
                 ci<-predict(vpd.gam, se=T)
                 ci$lower<-ci$fit-qt(0.975,vpd.gam$df.null)*ci$se.fit
                 ci$upper<-ci$fit+qt(0.975,vpd.gam$df.null)*ci$se.fit
                 l.ci<-cbind(vpd.gam$model$mo,ci$fit,ci$lower,ci$upper)
                 l<-l.ci[order(l.ci[,1]),]
                 panel.ci(l[,1],l[,2],l[,4],l[,3],
                          fill="seagreen3",alpha = 0.3)
                 panel.xyplot(x, y, pch=20,col="seagreen")
                 panel.lines(l[,1], l[,2],lty=1, col='black', lwd=1.5)
                 summ<-summary(vpd.gam)
                 r2 <- summ$r.sq
                 #f <- summ$fstatistic
                # p <- pf(f[1],f[2],f[3],lower.tail=F)
                 panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                            x=2150,y=-0.12,cex=0.75)
                 # panel.text(labels = bquote(italic(p)==.(format(p,digits = 3))),
                 #            x=1.2,y=-0.15,cex=0.75)
               },
                ylab="Disturbance Magnitude (TCG)",
                xlab="Mean VPD",
)
print(vpdplot)

#saving plots:
tiff("Plots_509/vpdplot_jul2015.tiff", units="in", width=8, height=5, res=300)
print(vpdplot)
dev.off()

# preciptrend<-xyplot(pprecip ~ mos, data = precipmeans,
#                    panel = function(x, y) {
#                     panel.xyplot(x, y, pch=16,col=yr)
#                     panel.abline(lm(y ~ x), col='royalblue4')
#                    },
#                    key=,
#                    ylab="Mean Precipitation (mm)",
#                    xlab="Month (2011-2015, Apr-Aug)",
# )
# print(preciptrend)

# tiff("preciptrend5yr_color.tiff", units="in", width=7, height=5, res=300)
# print(preciptrend)
# dev.off()
