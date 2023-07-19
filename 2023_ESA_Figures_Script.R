#Figures code for ESA 2023

# Harvard Forest data - Landsat time series at some of the sites:
site <- 204
#tcg <- tcg.values
tcg <- cond.scores

#get just the tcg values from data frame:
tcgs<-tcg[,c(grep("^X",colnames(tcg)))]

#first get just june values:
#junes<-seq(3,ncol(tcgs),by=5)  #3 is if tcg/condscores data have april in them, 2 if data starts w/May
#tcgjune<-as.matrix(tcgs[,junes])

y <- tcgs[site, 31:79]
x <- length(y[1,])
#samp_time_series <- tcgjune[site,]
#par(mar=c(1,1,1,1))
plot(1:x, y[1,])

