### TCG 2016 & 2017 DATASET VERSION ###

##### Load TCG data here#####
#monthly tcg:
tcgfile <- "2020_07_10_sample_tcg_mean_MONTHLY.csv"
tcg.month<-read.csv(tcgfile)
tcg.mo<-tcg.month[,2:131]

#just june values of tcg:
junes<-seq(2,130,by=5)
tcg.june.all<-as.matrix(tcg.mo[,junes])

#JUNE 2016:
june16<-tcg.june[,22]]
NA16<-which(is.na(june16))

#remove rows with missing values for disturbance window
NA16<-which(is.na(june16))
NA17<-which(is.na(june17))
missing<-intersect(NA16,NA17)
tcg.june<-tcg.june.all[-missing,]

june17<-tcg.june[,23
