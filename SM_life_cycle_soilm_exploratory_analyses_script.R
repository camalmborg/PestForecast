####load data
#soil moisture data: predisturb window (2011-2015, hatch/feed (apr/may, june/jul))
soilmdata<-read.csv("soil_moisture_data.csv")[,2:11]

#mags data:
distmagrecov<-read.csv("SM_distmagrecov_data.csv")

####Running the GAMs
#load library
library(mgcv)

#data set:
magssoilm = cbind(distmagrecov$colnum,distmagrecov$mags,soilmdata)
colnames(magssoilm)<-c("colnum","mags","h11","f11","h12","f12","h13","f13",
                       "h14","f14","h15","f15")

vardat = magssoilm[magssoilm$colnum==23,]

#choose soilm year eg. soilm[,1]==hatch soilm 2011
# smh13<-vardat[,7]
# smh14<-vardat[,9]
# smh15<-vardat[,11]
# smf13<-vardat[,8]
# smf14<-vardat[,10]
# smf15<-vardat[,12]


var.gam <- gam(vardat[,2]~s(f15), data = vardat)
summ<-summary(var.gam)
r2 <- summ$r.sq
print(r2)

#thers<-matrix(NA, nrow=50, ncol=1)
#thers[8,]<-r2

#plot(vardat$h13,vardat$mags)
