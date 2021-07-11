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
                    #ci<-predict(recov.rate.lm, interval="confidence")
                      #upper<-ci[,3]
                      #lower<-ci[,2]
                    #panel.ci(x, y, upper, lower, fill="gray48",alpha = 0.4)
                    panel.xyplot(x, y, pch=16,col="black")
                    panel.abline(lm(y ~ x))
                      rsumm<-summary(recov.rate.lm)
                      r2 <- rsumm$adj.r.squared
                      f <- rsumm$fstatistic
                      p <- pf(f[1],f[2],f[3],lower.tail=F)
                    panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                  x=0.13,y=0.06,cex=0.75)
                    panel.text(labels = bquote(italic(p)==.(format(p,digits = 3))),
                                  x=0.13,y=0.05,cex=0.75)
                             },
                    ylab="Recovery Rate",
                    xlab="TCG Steady State (2012-2015)",
                    )

print(recovlattice)
#save plot as tiff file:
tiff("recovlattice_monthlyNEW.tiff", units="in", width=7, height=5, res=300)
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
                     )

print(defollattice)
#save plot as tiff file:
tiff("defollattice_50sites.tiff", units="in", width=7, height=5, res=300)
print(defollattice)
dev.off()


#recovery time vs steady state plot:
rtimelattice<-xyplot(recov.time ~ mags, data = tcg.recov,
                     panel = function(x, y) {
                       # ci<-predict(recov.time.lm, interval="confidence")
                       # upper<-ci[,3]
                       # lower<-ci[,2]
                       # panel.ci(x, y, upper, lower, fill="gray48",alpha = 0.4)
                       panel.xyplot(x, y, pch=16,col="black")
                       panel.abline(lm(y ~ x))
                       rtsumm<-summary(recov.time.lm)
                       r2 <- rtsumm$adj.r.squared
                       f <- rtsumm$fstatistic
                       p <- pf(f[1],f[2],f[3],lower.tail=F)
                       panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                  x=0,y=200,cex=0.75)
                       panel.text(labels = bquote(italic(p)==.(format(p,digits=3))),
                                  x=0,y=170,cex=0.75)
                     },
                     ylab="Recovery Time (growing season months)",
                     xlab="Disturbance Magnitude")

print(rtimelattice)
#save plot as tiff file:
tiff("rtimelattice.tiff", units="in", width=7, height=5, res=300)
print(rtimelattice)
dev.off()


#recovery rate vs disturbance magnitude plot:
magrecovlattice<-xyplot(recov.rate ~ mags, data = tcg.recov,
                     panel = function(x, y) {
                       #ci<-predict(magrecov.lm, interval="confidence")
                       #upper<-ci[,3]
                       #lower<-ci[,2]
                       #panel.ci(x, y, upper, lower, fill="gray48",alpha = 0.4)
                       panel.xyplot(x, y, pch=16,col="black")
                       panel.abline(lm(y ~ x))
                       mrsumm<-summary(magrecov.lm)
                       r2 <- mrsumm$adj.r.squared
                       f <- mrsumm$fstatistic
                       p <- pf(f[1],f[2],f[3],lower.tail=F)
                       panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                  x=0,y=0.20,cex=0.75)
                       panel.text(labels = bquote(italic(p)==.(format(p,digits=3))),
                                  x=0,y=0.18,cex=0.75)
                     },
                     ylab="Recovery Rate",
                     xlab="Disturbance Magnitude")

print(magrecovlattice)

tiff("magrecovlattice_NEW.tiff", units="in", width=7, height=5, res=300)
print(magrecovlattice)
dev.off()

#recovery rate vs defoliation from previous year plot:
prevdefollattice<-xyplot(recov.rate ~ defol, data = recov.view,
                        panel = function(x, y) {
                          panel.xyplot(x, y, pch=16,col="black")
                        },
                        ylab="Recovery Rate",
                        xlab="Defoliation")

print(prevdefollattice)

tiff("prevdefollattice_NEW.tiff", units="in", width=7, height=5, res=300)
print(prevdefollattice)
dev.off()

#recovery rate vs previous years greenness plot:
prevyrlattice<-xyplot(recov.rate ~ prevyr, data = recov.view,
                      panel = function(x, y) {
                        panel.xyplot(x, y, pch=16,col="black")
                      },
                      ylab="Recovery Rate",
                      xlab="Previous Year Greenness")

print(prevyrlattice)

tiff("prevyrlattice_NEW.tiff", units="in", width=7, height=5, res=300)
print(prevyrlattice)
dev.off()

#magnitude vs percent tree cover
magtreelattice<-xyplot(mags ~ percentcover | NLCD, data = tcg.recov,
                      panel = function(x, y) {
                        panel.xyplot(x, y, pch=16,col="black")
                      },
                      ylab="Recovery Rate",
                      xlab="Percent Tree Cover")

print(magtreelattice)

tiff("magtreelattice_NEW.tiff", units="in", width=7, height=5, res=300)
print(magtreelattice)
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



mags_hist_lattice <- histogram(~ mags, data = tcg.recov,
                               xlab="Disturbance Magnitude",
                               breaks=250,
                               col="gray48")
print(mags_hist_lattice)

tiff("mags_hist_lattice.tiff", units="in", width=7, height=4, res=300)
print(mags_hist_lattice)
dev.off()


recov.rate_hist_lattice <- histogram(~ recov.rate, data = tcg.recov,
                               xlab="Recovery Rate",
                               breaks=250,
                               col="gray48")
print(recov.rate_hist_lattice)

tiff("recov.rate_hist_lattice.tiff", units="in", width=7, height=4, res=300)
print(recov.rate_hist_lattice)
dev.off()


steady_hist_lattice <- histogram(~ steady, data = tcg.recov,
                                     xlab="Steady State TCG Values",
                                     breaks=250,
                                     col="gray48")
print(steady_hist_lattice)

tiff("steady_hist_lattice.tiff", units="in", width=7, height=4, res=300)
print(steady_hist_lattice)
dev.off()


mins_hist_lattice <- histogram(~ mins, data = tcg.recov,
                                 xlab="Disturbance Window Minimum TCG Values",
                                 breaks=250,
                                 col="gray48")
print(mins_hist_lattice)

tiff("mins_hist_lattice.tiff", units="in", width=7, height=4, res=300)
print(mins_hist_lattice)
dev.off()



#------GGPLOT VERSION-------------------------------------------
if (FALSE){
recovsumm<-summary(recov.rate.lm)

recovplot<-ggplot(data = tcg.recov, aes(x=`Steady State 2012-2015`, y=`Recovery Rate`))+geom_point()
recovplot+geom_abline(intercept = recovsumm$coefficients[1],slope=recovsumm$coefficients[2])

ggplot(tcg.recov, aes(x=`Steady State 2012-2015`, 
                      y=`Disturbance Magnitude`,
                      color=as.factor(`NLCD Type`))) + geom_point(shape=1) +
  scale_colour_hue(l=50)+  # Use a slightly darker palette than normal
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE,    # Don't add shaded confidence region
              fullrange=TRUE) # Extend regression lines
}

#------Categorical------------------------------------------------

magrecovlattice_NLCD<-xyplot(recov.rate ~ mags | NLCD, data = tcg.recov,
                        panel = function(x, y) {
                          panel.xyplot(x, y, pch=16,col="black")
                          mr.lm<-lm(y~x)
                          mrsumm<-summary(mr.lm)
                          r2 <- mrsumm$adj.r.squared
                          f <- mrsumm$fstatistic
                          p <- pf(f[1],f[2],f[3],lower.tail=F)
                          #panel.abline(a = mr.lm$coefficients[1], 
                                       #b = mr.lm$coefficients[2])
                          #panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                     #x=0.177,y=0.015,cex=0.75)
                          #panel.text(labels = bquote(italic(p)==.(format(p,digits=3))),
                                     #x=0.176,y=0.008,cex=0.75)
                        },
                        ylab="Recovery Rate",
                        xlab="Disturbance Magnitude")

print(magrecovlattice_NLCD)

tiff("magslattice_NLCD.tiff", units="in", width=7, height=5, res=300)
print(magrecovlattice_NLCD)
dev.off()


recovtimelattice_NLCD<-xyplot(recov.time ~ mags | NLCD, data = tcg.recov,
                             panel = function(x, y) {
                               panel.xyplot(x, y, pch=16,col="black")
                               mr.lm<-lm(y~x)
                               mrsumm<-summary(mr.lm)
                               r2 <- mrsumm$adj.r.squared
                               f <- mrsumm$fstatistic
                               p <- pf(f[1],f[2],f[3],lower.tail=F)
                               #panel.abline(a = mr.lm$coefficients[1], 
                               #b = mr.lm$coefficients[2])
                               #panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                               #x=0.177,y=0.015,cex=0.75)
                               #panel.text(labels = bquote(italic(p)==.(format(p,digits=3))),
                               #x=0.176,y=0.008,cex=0.75)
                             },
                             ylab="Recovery Time",
                             xlab="Disturbance Magnitude",
                             layout=c(3,2))

print(recovtimelattice_NLCD)

tiff("recovtimelattice_NLCD.tiff", units="in", width=7, height=5, res=300)
print(recovtimelattice_NLCD)
dev.off()

