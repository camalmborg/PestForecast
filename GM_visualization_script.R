###VISUALIZATION CODE###

#load plotting libraries
library(ggplot2)
library(lattice)
library(latticeExtra)


#------Lattice version:-----------------------------------

recovlattice<-xyplot(recov.rate ~ steady, data = tcg.recov,
                    panel = function(x, y) {
                               panel.xyplot(x, y, pch=16,col="black")
                               panel.abline(lm(y ~ x))
                               rsumm<-summary(recov.rate.lm)
                               r2 <- rsumm$adj.r.squared
                               panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                                          x=0.13,y=0.0055,cex=0.75)
                             },
                    ylab="Recovery Rate",
                    xlab="TCG Steady State (2012-2015)")

print(recovlattice)

tiff("recovlattice.tiff", units="in", width=7, height=5, res=300)
print(recovlattice)
dev.off()


#------GGPLOT VERSION-------------------------------------------
if (FALSE){
recovsumm<-summary(recov.rate.lm)

recovplot<-ggplot(data = tcg.recov, aes(x=`Steady State 2012-2015`, y=`Recovery Rate`))+geom_point()
recovplot+geom_abline(intercept = recovsumm$coefficients[1],slope=recovsumm$coefficients[2])
}
