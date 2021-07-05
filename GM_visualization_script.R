###VISUALIZATION CODE###

#load plotting libraries
library(ggplot2)
library(lattice)
library(latticeExtra)

recovlattice<-xyplot(recov.rate ~ steady, data = tcg.recov,
                    panel = function(x, y) {
                               panel.xyplot(x, y)
                               panel.abline(lm(y ~ x))
                             }, 
                    ylab="Recovery Rate",
                    xlab="TCG Steady State (2012-2015)")

print(recovlattice)



#------GGPLOT VERSION--------
recovsumm<-summary(recov.rate.lm)

recovplot<-ggplot(data = tcg.recov, aes(x=`Steady State 2012-2015`, y=`Recovery Rate`))+geom_point()
recovplot+geom_abline(intercept = recovsumm$coefficients[1],slope=recovsumm$coefficients[2])

