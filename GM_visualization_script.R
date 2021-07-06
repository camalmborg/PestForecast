###VISUALIZATION CODE###

#load plotting libraries
library(ggplot2)
library(lattice)
library(latticeExtra)
library(tactile)


#------Lattice versions:-----------------------------------

#recovery rate vs steady state plot:
recovlattice<-xyplot(recov.rate ~ steady, data = tcg.recov,
                    panel = function(x, y) {
                    ci<-predict(recov.rate.lm, interval="confidence")
                      upper<-ci[,3]
                      lower<-ci[,2]
                    panel.ci(x, y, upper, lower, fill="gray48",alpha = 0.4)
                    panel.xyplot(x, y, pch=16,col="black")
                    panel.abline(lm(y ~ x))
                      rsumm<-summary(recov.rate.lm)
                      r2 <- rsumm$adj.r.squared
                      f <- rsumm$fstatistic
                      p <- pf(f[1],f[2],f[3],lower.tail=F)
                    panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                  x=0.131,y=0.0055,cex=0.75)
                    panel.text(labels = bquote(italic(p)==.(format(p,digits=3))),
                                  x=0.13,y=0.0045,cex=0.75)
                             },
                    ylab="Recovery Rate",
                    xlab="TCG Steady State (2012-2015)")

print(recovlattice)
#save plot as tiff file:
tiff("recovlattice_50sites.tiff", units="in", width=7, height=5, res=300)
print(recovlattice)
dev.off()


#disturbance magnitude vs steady state plot:
defollattice<-xyplot(mags ~ steady, data = tcg.recov,
                     panel = function(x, y) {
                       ci<-predict(defol.lm, interval="confidence")
                       upper<-ci[,3]
                       lower<-ci[,2]
                       panel.ci(x, y, upper, lower, fill="gray48",alpha = 0.4)
                       panel.xyplot(x, y, pch=16,col="black")
                       panel.abline(lm(y ~ x))
                       dsumm<-summary(defol.lm)
                       r2 <- dsumm$adj.r.squared
                       f <- dsumm$fstatistic
                       p <- pf(f[1],f[2],f[3],lower.tail=F)
                       panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                  x=0.13,y=0.18,cex=0.75)
                       panel.text(labels = bquote(italic(p)==.(format(p,digits=3))),
                                  x=0.1325,y=0.15,cex=0.75)
                     },
                     ylab="Disturbance Magnitude",
                     xlab="TCG Steady State (2012-2015)",
                     ylim=c(-0.015,0.21))

print(defollattice)
#save plot as tiff file:
tiff("defollattice_50sites.tiff", units="in", width=7, height=5, res=300)
print(defollattice)
dev.off()


#recovery time vs steady state plot:
rtimelattice<-xyplot(recov.time ~ steady, data = tcg.recov,
                     panel = function(x, y) {
                       ci<-predict(recov.time.lm, interval="confidence")
                       upper<-ci[,3]
                       lower<-ci[,2]
                       panel.ci(x, y, upper, lower, fill="gray48",alpha = 0.4)
                       panel.xyplot(x, y, pch=16,col="black")
                       panel.abline(lm(y ~ x))
                       rtsumm<-summary(recov.time.lm)
                       r2 <- rtsumm$adj.r.squared
                       f <- rtsumm$fstatistic
                       p <- pf(f[1],f[2],f[3],lower.tail=F)
                       panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                  x=0.132,y=-50,cex=0.75)
                       panel.text(labels = bquote(italic(p)==.(format(p,digits=3))),
                                  x=0.13,y=-85,cex=0.75)
                     },
                     ylab="Recovery Time",
                     xlab="TCG Steady State (2012-2015)")

print(rtimelattice)
#save plot as tiff file:
tiff("rtimelattice_50sites.tiff", units="in", width=7, height=5, res=300)
print(rtimelattice)
dev.off()


#recovery rate vs disturbance magnitude plot:
magrecovlattice<-xyplot(recov.rate ~ mags, data = tcg.recov,
                     panel = function(x, y) {
                       ci<-predict(magrecov.lm, interval="confidence")
                       upper<-ci[,3]
                       lower<-ci[,2]
                       panel.ci(x, y, upper, lower, fill="gray48",alpha = 0.4)
                       panel.xyplot(x, y, pch=16,col="black")
                       panel.abline(lm(y ~ x))
                       mrsumm<-summary(magrecov.lm)
                       r2 <- mrsumm$adj.r.squared
                       f <- mrsumm$fstatistic
                       p <- pf(f[1],f[2],f[3],lower.tail=F)
                       panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                  x=0.177,y=0.001,cex=0.75)
                       panel.text(labels = bquote(italic(p)==.(format(p,digits=3))),
                                  x=0.18,y=0,cex=0.75)
                     },
                     ylab="Recovery Rate",
                     xlab="Disturbance Magnitude")

print(magrecovlattice)

tiff("magrecovlattice_50sites.tiff", units="in", width=7, height=5, res=300)
print(magrecovlattice)
dev.off()


#minimum TCG value vs steady state plot:
minslattice<-xyplot(mins ~ steady, data = tcg.recov,
                        panel = function(x, y) {
                          ci<-predict(mins.lm, interval="confidence")
                          upper<-ci[,3]
                          lower<-ci[,2]
                          panel.ci(x, y, upper, lower, fill="gray48",alpha = 0.4)
                          panel.xyplot(x, y, pch=16,col="black")
                          panel.abline(lm(y ~ x))
                          minsumm<-summary(mins.lm)
                          r2 <- minsumm$adj.r.squared
                          f <- minsumm$fstatistic
                          p <- pf(f[1],f[2],f[3],lower.tail=F)
                          panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                     x=0.13,y=0.19,cex=0.75)
                          panel.text(labels = bquote(italic(p)==.(format(p,digits=3))),
                                     x=0.132,y=0.17,cex=0.75)
                        },
                        ylab="Minimum TCG value",
                        xlab="TCG Steady State (2012-2015)")

print(minslattice)

tiff("minslattice_50sites.tiff", units="in", width=7, height=5, res=300)
print(minslattice)
dev.off()




#------GGPLOT VERSION-------------------------------------------
if (FALSE){
recovsumm<-summary(recov.rate.lm)

recovplot<-ggplot(data = tcg.recov, aes(x=`Steady State 2012-2015`, y=`Recovery Rate`))+geom_point()
recovplot+geom_abline(intercept = recovsumm$coefficients[1],slope=recovsumm$coefficients[2])
}
