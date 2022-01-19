####load data for 5000 sites:---------
#load monthly condition scores:
cond.scores.mo<-read.csv("2020_07_10_sample_score_mean_MONTHLY.csv")
cond.mo<-cond.scores.mo[,2:131]

#load annual condition scores:
cond.scores.an<-read.csv("2020_07_10_sample_score_mean_ANNUAL.csv")
cond.an<-cond.scores.an[,2:27]

#load monthly tcg:
tcg.month<-read.csv("2020_07_10_sample_tcg_mean_MONTHLY.csv")
tcg.mo<-tcg.month[,2:131]


# ####remove may:-----
# nomay<-seq(1,130,by=5)
# tcg.nomay<-as.matrix(tcg.mo[,-nomay])
# 
# #### steady state previous 3 years 2015------
# #find the non-GM disturb condition "steady state"
# steady<-apply(tcg.nomay[,73:84],1,mean,na.rm=T)
# #steadyall<-apply(tcg.nomay[,1:84],1,mean,na.rm=T)
# 
# 
# ####the TCG 2016-onward MONTHLY version (no Mays):-----
# #ind<-list()
# recov.rate<-vector()
# slope<-list()
# mins<-vector()
# colnum<-vector()
# recovcol<-vector()
# end<-length(tcg.nomay[1,])
# nsite<-length(tcg.nomay[,1])
# for (i in 1:nsite){
#   mins[i]<-min(tcg.nomay[i,85:92],na.rm=T)          #grabs min forest condition score for each site
#   colnum[i]<-which(tcg.nomay[i,]==mins[i])           #grabs time step at which min score appears in disturbance window
#   if(is.na(colnum[i] + which(tcg.nomay[i,colnum[i]:end]>=steady[i])[1])){
#     recovcol[i]<-colnum[i] + 2
#   } else {
#     recovcol[i]<-colnum[i] + which(tcg.nomay[i,colnum[i]:end]>=steady[i])[1]
#   }
#   if(recovcol[i]>end){
#     recov<-tcg.nomay[i,colnum[i]:end]
#   } else {
#     recov<-tcg.nomay[i,colnum[i]:recovcol[i]]
#   }
#   #recov <- tcg[i,colnum[i]:recovcol[i]]      
#   ind<-1:length(recov)                      #grabs length of recov rate   
#   slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
#   recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
# }
# 
# #magnitudes:
# mags<-steady-mins
# recov.time<-mags/recov.rate




#### THE JUNE VERSION:-----
####just junes
junes<-seq(2,130,by=5)
tcg.june.all<-as.matrix(tcg.mo[,junes])
#cond.june<-as.matrix(cond.mo[,junes])

###steady state
steadys<-as.matrix(apply(tcg.june.all[,19:21],1,mean,na.rm=T))
steadyall<-apply(tcg.june.all[,1:26],1,mean,na.rm=T)
prevyr<-tcg.june.all[,21]
june16<-tcg.june.all[,22]
june17<-tcg.june.all[,23]

###remove rows with missing values for disturbance window
NA16<-which(is.na(june16))
NA17<-which(is.na(june17))
missing<-intersect(NA16,NA17)
tcg.june<-tcg.june.all[-missing,]
steady<-steadys[-missing,]

###site identification:
#testing missing cols
#tcgidsj<-cbind(tcgids,june16,june17)
#tcgjidsmiss<-tcgidsj[missing,]
idcols<-c(1,132,134:136)
tcgids<-cbind(tcg.month[,idcols],june16,june17)
tcgid<-tcgids[-missing,]

### THE LOOP: calculating recovery rates
recov.rate<-vector()
slope<-list()
mins<-vector()
colnum<-vector()
recovcol<-vector()
end<-length(tcg.june[1,])
nsite<-length(tcg.june[,1])
for (i in 1:nsite){
  mins[i]<-min(tcg.june[i,22:23],na.rm=T)          #grabs min forest condition score for each site
  colnum[i]<-which(tcg.june[i,]==mins[i])           #grabs time step at which min score appears in disturbance window
  if(is.na(colnum[i] + which(tcg.june[i,colnum[i]:end]>=steadyall[i])[1])){
    recovcol[i]<-colnum[i] + 2
  } else {
    recovcol[i]<-colnum[i] + which(tcg.june[i,colnum[i]:end]>=steadyall[i])[1]
  }
  if(recovcol[i]>end){
    recov<-tcg.june[i,colnum[i]:end]
  } else {
    recov<-tcg.june[i,colnum[i]:recovcol[i]]
  }
  #recov <- tcg[i,colnum[i]:recovcol[i]]      
  ind<-1:length(recov)                      #grabs length of recov rate   
  slope[[i]]<-lm(recov~ind)                      #runs the lm to find slope of recov rate, saves output
  recov.rate[i]<-slope[[i]]$coefficients["ind"] #stores recovery rate
}

#magnitudes:
mags<-steady-mins
#mags.s<-steadyall-mins

recov.time<-mags/recov.rate
#recov.time.s<-mags.s/recov.rate

#combine data for regressions:
tcg.recov.mx<-cbind(tcgid,steady,colnum,mins,mags,recov.rate,recov.time)
tcg.recov<-as.data.frame(tcg.recov.mx)

####Regressions#####
#recov rate as function of magnitude of disturbance:
mrrreg<-lm(recov.rate~mags,tcg.recov)
#recov rate as function of magnitude and site greenness prior to disturbance:
anotherreg<-lm(recov.rate~mags+steady,tcg.recov)

#
#
#
### Plotting exercises #####
plot(tcg.recov$mags, tcg.recov$recov.rate, 
     col=c('red', 'blue'))
# legend('bottomleft',
#        legend=c('22','23'), 
#        col=c('red','blue'))
j16<-tcg.recov[tcg.recov$colnum==22,]
j17<-tcg.recov[tcg.recov$colnum==23,]
hist(j17$mags,breaks=75,col=rgb(1,0,0,0.5))
hist(j16$mags,breaks=75,col=rgb(0,0,1,0.5), add=T)
