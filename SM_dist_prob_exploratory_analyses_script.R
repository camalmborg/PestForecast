#script for disturbance probability exploratory analyses

#load condition scores (finding where prob of dist = 0,1)
cond.scores.mo<-read.csv("2020_07_10_sample_score_mean_MONTHLY.csv")
condscores<-cond.scores.mo[,2:131]
jjunes<-c(102,107,112)
junes<-seq(2,length(condscores[1,]),by=5)
js<-condscores[,junes]

#Look at all june scores:
jconds<-condscores[,junes]
jmeans<-apply(jconds,2,mean,na.rm=T)

#previous 5 yr mean group:
prev5<-jconds[,16:20]
p5<-apply(prev5,1,mean)
bin1.2010.4<-apply(prev5[1:1000,],1,mean)
bin2.2010.4<-apply(prev5[1001:2000,],1,mean)
bin3.2010.4<-apply(prev5[2001:3000,],1,mean)
bin4.2010.4<-apply(prev5[3001:4000,],1,mean)
bin5.2010.4<-apply(prev5[4001:5000,],1,mean)

#2015 group:
minj2015<-min(jconds$X2015.06.01_score_mean[1:1000],na.rm=T)
bin1.2015<-jconds$X2015.06.01_score_mean[1:1000]
bin2.2015<-jconds$X2015.06.01_score_mean[1001:2000]
bin3.2015<-jconds$X2015.06.01_score_mean[2001:3000]
bin4.2015<-jconds$X2015.06.01_score_mean[3001:4000]
bin5.2015<-jconds$X2015.06.01_score_mean[4001:5000]

#2016 group:
minj2016<-min(jconds$X2016.06.01_score_mean[1:1000],na.rm=T)
bin1.2016<-jconds$X2016.06.01_score_mean[1:1000]
bin2.2016<-jconds$X2016.06.01_score_mean[1001:2000]
bin3.2016<-jconds$X2016.06.01_score_mean[2001:3000]
bin4.2016<-jconds$X2016.06.01_score_mean[3001:4000]
bin5.2016<-jconds$X2016.06.01_score_mean[4001:5000]


minj2017<-min(jconds$X2017.06.01_score_mean[1:1000],na.rm=T)
bin1.2017<-jconds$X2017.06.01_score_mean[1:1000]
bin2.2017<-jconds$X2017.06.01_score_mean[1001:2000]
bin3.2017<-jconds$X2017.06.01_score_mean[2001:3000]
bin4.2017<-jconds$X2017.06.01_score_mean[3001:4000]
bin5.2017<-jconds$X2017.06.01_score_mean[4001:5000]


#load tcg data for column numbers:
distmagsdata<-read.csv("SM_distmagrecov_data.csv")
dmdjs<-cbind(distmagsdata,js)

bd1<-dmdjs$mags[1:1000]
bd2<-dmdjs$mags[1001:2000]
bd3<-dmdjs$mags[2001:3000]
bd4<-dmdjs$mags[3001:4000]
bd5<-dmdjs$mags[4001:5000]

#adding 0,1 for disturbance occurrence:
#make empty matrix for 0,1 data:
distprob<-matrix(NA,nrow=nrow(condscores),ncol=1)
#disturbance threshold:
d = quantile(dmdjs[,23:27],c(0.01),na.rm=T) #selecting dist sites 1% quant in prev 5 yr
#determine if disturbance (condscore < d threshold) occurs
for (i in 1:nrow(condscores)){
  if (dmdjs[i,]$colnum == 21 & dmdjs[i,]$X2015.06.01_score_mean <= d){
    distprob[i,]<- 1
    } else if (dmdjs[i,]$colnum == 22 & dmdjs[i,]$X2016.06.01_score_mean <= d){
        distprob[i,]<- 1
    } else if (dmdjs[i,]$colnum == 23 & dmdjs[i,]$X2017.06.01_score_mean <= d){
        distprob[i,]<- 1
    } else {
        distprob[i,] <- 0
            }
}

disturbed<-dmdjs[which(distprob==1),]
csdist<-condscores[which(distprob==1),]



##### Running analyses -----
##load library
library(mgcv)

##load environmental data
dmhatch<-read.csv("hatch_daymet_allvar.csv")
dmfeed<-read.csv("feed_daymet_allvar.csv")
soilm<-read.csv("soil_moisture_data.csv")[-missing,2:11]
soilmh<-soilm[,seq(1,length(soilm),by=2)]
soilmf<-soilm[,seq(2,length(soilm),by=2)]

##create vardat
#choose the variable you want to test:
var = dmfeed[,50]
vardat = as.data.frame(cbind(distprob[-missing,], var))
colnames(vardat)<-c("distprob","var")

#plot:
plot(vardat$var,vardat$distprob)
lineseq<-seq(0,max(vardat$var,na.rm=T),by=0.1)
lines(lineseq,pnorm(lineseq,mean(vardat$var),sd(vardat$var)),type='l')

#run the glms:
var.glm <- glm(distprob~var, data = vardat, family=binomial)
summ<-summary(var.glm)
#print(summ$coefficients)

#McFaddens R2
#with(summary(var.glm), 1 - deviance/null.deviance)
