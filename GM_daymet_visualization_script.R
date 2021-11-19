###PLOTS-----
# #plot a sample with tcg as x:
# samp<-sample(nsites,500)
# xl=c(1,25) #maxtemp max
# yl=c(min(tcg.june,na.rm=T),max(tcg.june,na.rm=T))
# plot(maxvar[[1]][1,], tcg.june[1,], pch=20, ylim=yl, xlim=xl)
# for (s in samp){
#   points(maxvar[[s]][1,], tcg.june[s,], pch=20)
# }

#with condition scores as x:
samp<-sample(nsites,500)
xl=c(5,33) #maxtemp max
yl=c(min(cond.june,na.rm=T),max(cond.june,na.rm=T))
plot(maxvar[[1]][3,], cond.june[1,], pch=20, ylim=yl, xlim=xl)
for (s in nsites){
  points(maxvar[[s]][3,], cond.june[s,], pch=20)
}

#with disturbance mag as x:
#samp<-sample(nsites,500)
xl=c(25,37) #maxtemp max
yl=c(0,max(mags,na.rm=T))
plot(maxvar[[1]][6,23], mags[1], pch=20, ylim=yl, xlim=xl)
for (s in nsites){
  points(maxvar[[s]][6,23], mags[s], pch=20)
}

