#This script is for wrangling the anthropogenic variables

#load datasets
#Human population density (county scale):
popdens<-read.csv("Anthro_Var_Data/County_Pop_Data_MA_CT_RI.csv")
#sample point counties:
pointcty<-read.csv("Anthro_Var_Data/point_counties.csv")
#tcg.recov data to run exploratory analyses:
tcg.recov<-read.csv("SM_distmagrecov_data.csv")

library(dplyr)
ctypop<-left_join(pointcty, popdens, by = c("STATE_NAME" = "STATE_NAME", "NAME" = "NAME"))

#2016 outbreak year values:
ctypop$X2016 <- gsub(",","",ctypop$X2016)
ctypop$X2016<-as.integer(ctypop$X2016)
pop2016<-ctypop$X2016[-missing]

#2016 outbreak year values:
ctypop$X2017 <- gsub(",","",ctypop$X2017)
ctypop$X2017<-as.integer(ctypop$X2017)
pop2017<-ctypop$X2017[-missing]


#Exploratory Analyses----
library(mgcv)
plot(pop2016,tcg.recov$mags)

#the data for the gam:
magspop<-as.data.frame(cbind(pop2016,pop2017,tcg.recov$mags,tcg.recov$recov.rate,
               tcg.recov$colnum))
colnames(magspop)<-c("pop2016","pop2017","mags","recovrate","colnum")

#only 2016:
magspop16<-magspop[magspop$colnum==22,]
#only 2017:
magspop17<-magspop[magspop$colnum==23,]


#THE GAM
agam<-gam(mags ~ pop2016, data=magspop)
#plot.gam(agam, all.terms = T)

#visualization:----
# library(lattice)
# library(latticeExtra)
# library(tactile)
#plot the gam:
gamplot<-xyplot(mags ~ pop2016, data = magspop,
                panel = function(x, y) {
                  ci<-predict(agam, se=T)
                  ci$lower<-ci$fit-qt(0.975,agam$df.null)*ci$se.fit
                  ci$upper<-ci$fit+qt(0.975,agam$df.null)*ci$se.fit
                  l.ci<-cbind(agam$model$pop2016,ci$fit,ci$lower,ci$upper)
                  l<-l.ci[order(l.ci[,1]),]
                  panel.ci(l[,1],l[,2],l[,4],l[,3],
                           fill="orange",alpha = 0.4)
                  panel.xyplot(x, y, pch=20,col="orange3")
                  panel.lines(l[,1], l[,2],lty=1, col='black', lwd=1.5)
                })
print(gamplot)

summ<-summary(agam)
r2<-summ$r.sq
print(r2)
