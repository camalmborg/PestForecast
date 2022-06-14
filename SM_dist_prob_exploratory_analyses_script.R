#script for disturbance probability exploratory analyses

#load condition scores (finding where prob of dist = 0,1)
cond.scores.mo<-read.csv("2020_07_10_sample_score_mean_MONTHLY.csv")
condscores<-cond.scores.mo[,2:131]
junes<-c(102,107,112)
js<-condscores[,junes]

#load tcg data for column numbers:
distmagsdata<-read.csv("SM_distmagrecov_data.csv")
dmdjs<-cbind(distmagsdata,js)

#adding 0,1 for disturbance occurrence:
#make empty matrix for 0,1 data:
distprob<-matrix(NA,nrow=nrow(condscores),ncol=1)
#disturbance threshold:
d = quantile(dmdjs[,8:10],c(0.25),na.rm=T) #selecting 75%ile of dist
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
var = dmhatch[,20]
vardat = as.data.frame(cbind(distprob[-missing,], var))
colnames(vardat)<-c("distprob","var")

plot(vardat$var,vardat$distprob)

#run the glms
var.glm <- glm(distprob~var, data = vardat, family="binomial")
summ<-summary(var.glm)
#r2 <- summ$r.sq
#print(r2)
