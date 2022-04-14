### TCG 2016 & 2017 DATASET VERSION ###

##### TCG VALUES FOR JUNE 16-17 #####
#monthly tcg:
tcgfile <- "2020_07_10_sample_tcg_mean_MONTHLY.csv"
tcg.month<-read.csv(tcgfile)
tcg.mo<-tcg.month[,2:131]

#number of sites:
nsites<-5000

#just june values of tcg:
juneseq<-seq(2,130,by=5)
tcg.june.all<-as.matrix(tcg.mo[,juneseq])

#JUNE 2016:
june16<-as.matrix(tcg.june.all[1:nsites,22])
#removing NA values:
NA16<-which(is.na(june16))
tcgjune16<-june16[-NA16,]

#JUNE 2017:
june17<-as.matrix(tcg.june.all[1:nsites,23])
NA17<-which(is.na(june17))
tcgjune17<-june16[-NA17,]

#still have to calculate individual steady states for 16 and 17

#Next step: find min values for each site
##### MINIMUM TCG VALUES #####
recov.rate<-vector()
slope<-list()
mins<-vector()
colnum<-vector()
recovcol<-vector()
#end<-length(tcg.june[1,])
jyear<-tcgjune16 #or tcgjune17
colnum<-22 #or 23 for tcg17
nsite<-length(jyear)
for (i in 1:nsite){
  mins[i]<-min(jyear[i])        #grabs min tcg for each site
  #recov <- tcg[i,colnum[i]:recovcol[i]]      
  #ind<-1:length(recov)                      #grabs length of recov rate   
  #slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  #recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
}

#magnitudes:
mags<-steady-mins
mags.s<-steadyall-mins
recov.time<-mags/recov.rate
recov.time.s<-mags.s/recov.rate

#combine data for regressions:
tcg.recov.mx<-cbind(steady,colnum,mins,mags,recov.rate,recov.time)
tcg.recov<-as.data.frame(tcg.recov.mx)
#save a copy:
write.csv(tcg.recov,file="Desktop/CM_BostonU/D_Research/PestForecast/SM_distmagrecov_data.csv")
# 
# ###regressions
# #recov rate as function of magnitude of disturbance:
# mrrreg<-lm(recov.rate~mags,tcg.recov)
# #recov rate as function of magnitude and site greenness prior to disturbance:
# anotherreg<-lm(recov.rate~mags+steady,tcg.recov)