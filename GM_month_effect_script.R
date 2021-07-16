#Month effect script
#load TCG mean data:
tcgmeans<-"2020_07_10_500_sample_tcg_mean_MONTHLY.csv"
tcgtable<-read.csv(tcgmeans)
tcg<-as.matrix(tcgtable[,2:131])

#monthly steady states
#maymeans<-apply(tcg[,seq(1,ncol(tcg),by=5)],1,mean,na.rm=T)
steadymonths<-matrix(NA,nrow=2500,ncol=5)
months<-5
for (i in 1:months){
  steadymonths[,i]<-apply(tcg[,seq(i,ncol(tcg),by=5)],1,mean,na.rm=T)
}
sites<-1:2500
month<-c("May","June","July","August","September")
steadymonth<-as.data.frame(steadymonths)
#steadymonth<-cbind(sites,steadymonth)
colnames(steadymonth)=c("May","June","July","August","September")


