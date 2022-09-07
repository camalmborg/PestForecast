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

# #adding 0,1 for disturbance occurrence:
# #make empty matrix for 0,1 data:
# distprob<-matrix(NA,nrow=nrow(condscores),ncol=1)
# #disturbance threshold:
# d = quantile(dmdjs[,23:27],c(0.01),na.rm=T) #selecting dist sites 1% quant in prev 5 yr
# #determine if disturbance (condscore < d threshold) occurs
# for (i in 1:nrow(condscores)){
#   if (dmdjs[i,]$colnum == 21 & dmdjs[i,]$X2015.06.01_score_mean <= d){
#     distprob[i,]<- 1
#     } else if (dmdjs[i,]$colnum == 22 & dmdjs[i,]$X2016.06.01_score_mean <= d){
#         distprob[i,]<- 1
#     } else if (dmdjs[i,]$colnum == 23 & dmdjs[i,]$X2017.06.01_score_mean <= d){
#         distprob[i,]<- 1
#     } else {
#         distprob[i,] <- 0
#             }
# }

####time series version:-----
distprob<-matrix(NA,nrow=nrow(condscores),ncol=3)
#disturbance threshold:
d = quantile(dmdjs[,23:27],c(0.01),na.rm=T) #selecting dist sites 1% quant in prev 5 yr
#determine if disturbance (condscore < d threshold) occurs
for (i in 1:nrow(condscores)){
  if (dmdjs[i,]$colnum == 21 & dmdjs[i,]$X2015.06.01_score_mean <= d){
    distprob[i,1]<- 1
  } else if (dmdjs[i,]$colnum == 22 & dmdjs[i,]$X2016.06.01_score_mean <= d){
    distprob[i,2]<- 1
  } else if (dmdjs[i,]$colnum == 23 & dmdjs[i,]$X2017.06.01_score_mean <= d){
    distprob[i,3]<- 1
  } else {
    distprob[i,] <- 0
  }
  distprob[is.na(distprob)] <- 0
}
##-----
#disturbed<-dmdjs[which(distprob==1),]
#csdist<-condscores[which(distprob==1),]



##### Running analyses -----

##load environmental data
dmhatch<-read.csv("hatch_daymet_allvar.csv")[,-1]
dmfeed<-read.csv("feed_daymet_allvar.csv")[,-1]
#soilm<-read.csv("soil_moisture_data.csv")[-missing,2:11]
#soilmh<-soilm[,seq(1,length(soilm),by=2)]
#soilmf<-soilm[,seq(2,length(soilm),by=2)]

##create vardat: 
#choose the variable you want to test (FROM DAYMET SET):
vs=c(16:20,42:46,68:72) #pre-dist (2010-2014)
#vs=c(21:23,47:49,73:75) #dist (2015-2017)

#choose hatch or feed:
vars = dmhatch[,vs]
vardat = as.data.frame(cbind(distprob, vars))
#pre-d colnames:
colnames(vardat)<-c("dist2015","dist2016","dist2017",
                    "temp2010","temp2011","temp2012","temp2013","temp2014",
                    "pcp2010","pcp2011","pcp2012","pcp2013","pcp2014",
                    "vpd2010","vpd2011","vpd2012","vpd2013","vpd2014")
#dist colnames:
# colnames(vardat)<-c("dist2015","dist2016","dist2017",
#                     "temp2015","temp2016","temp2017",
#                     "pcp2015","pcp2016","pcp2017",
#                     "vpd2015","vpd2016","vpd2017")


###AUC analyses:
library(pROC)
library(mgcv)

##run the gams:
var.gam<-gam(dist2016~s(pcp2014), data=vardat, family="binomial")
roc<-roc(vardat$dist2016,var.gam$fitted.values)#, plot=T)
print(roc$auc)


## plotting:
#plot(vardat$pcp2010,vardat$dist2016)
#lines(var.gam$fitted.values)

# #plot.gam(var.gam)
# summ<-summary(var.gam)
# r2 <- summ$r.sq
# print(r2)

# #plot:
#varr<-
# plot(vardat$var2,vardat$distprob)
# lineseq<-seq(0,max(vardat$var2,na.rm=T),by=0.1)
# lines(lineseq,pnorm(lineseq,mean(vardat$var2),sd(vardat$var2)),type='l')

# #run the glms:
#var.glm <- glm(distprob~var, data = vardat, family=binomial)
# summ<-summary(var.glm)
# #print(summ$coefficients)
#plot(var.glm)

#McFaddens R2
#with(summary(var.glm), 1 - deviance/null.deviance)
