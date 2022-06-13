####load data
#soil moisture data: predisturb window (2011-2015, hatch/feed (apr/may, june/jul))
soilmdata<-read.csv("soil_moisture_data.csv")[,2:11]

#mags data:
distmagrecov<-read.csv("SM_distmagrecov_data.csv")

####Running the GAMs
#load library
library(mgcv)

#data set:
vardat = cbind(distmagrecov$mags,soilmdata)

#choose soilm year eg. soilm[,1]==hatch soilm 2011
smh13<-vardat[,5]
smh14<-vardat[,7]
smh15<-vardat[,9]
smf13<-vardat[,6]
smf14<-vardat[,8]
smf15<-vardat[,10]


var.gam <- gam(mags~s(smh13), data = vardat)
summ<-summary(var.gam)
r2 <- summ$r.sq
print(r2)
