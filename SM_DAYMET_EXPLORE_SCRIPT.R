# This code is for the daymet exploratory analyses

#load data if not in environment:
#load("DMVARS_MO.RData")

# 11/26/2024
load("2024_08_dm_grab.RData")
names <- read.table("colnamesbcimstupid.txt")
varcolnames <- as.character(names[,2])
# data frame the daymet data:
daymet <- cbind.data.frame(dmvars[[1]], dmvars[[2]], dmvars[[3]], dmvars[[4]])
colnames(daymet) <- varcolnames
write.csv(daymet, "2024_11_26_Daymet_For_Supp.csv")

# add column names:
years <- c("2014", "2015", "2016")
variable <- c("maxtemp", "mintemp", "precip", "vpd")
month <- as.character(1:12)

# 
#   for (i in 1:length(variable)){
#     for (j in 1:length(years)){
#       for (k in 1:length(month)){
#          print(paste0(years[j],"_",
#                month[k],"_",
#                variable[i]))
#       }
#     }
#   }
# printed to console, copied and pasted as a txt and loaded in - has names


# for (i in 1:length(variable)){
#   for(i in 1:length(years)){
#     print(paste0(years[j],"_",variable[i],"_"))
#   }
# }

# necessary libraries:
library(mgcv)

#function for univariate analyses:
dm_explore<-function(dmvars,dmrdat,dmr,coln){
  #make empty matrix:
  r2s <- matrix(NA,nrow=ncol(dmvars[[1]]),ncol=length(dmvars))
  
  #loop over all members of dmvars list:
  for (i in 1:length(dmvars)){
    #first extract the list you want:
    dmvariable <- as.data.frame(dmvars[[i]][as.numeric(dmrdat$sitenum),])
    #grab column number:
    cn <- as.matrix(as.numeric(dmrdat$colnum))
    #grab exploratory response  of choice:
    y <- as.matrix(as.numeric(dmrdat[,dmr]))
    x <- as.data.frame(cbind(y,cn,dmvariable))
    vardat <- x[x$cn==coln,]
  
    #loop for filling in R2 table:  
      for (j in 1:ncol(dmvariable)){
        var.gam <- gam(vardat[,1]~s(vardat[,j+2]),data=vardat)
        
        #extract r2:
        summ <- summary(var.gam)
        r2s[j,i] <-summ$r.sq
      }
  }
 # return(vardat)
  return(r2s)
}

testing <- dm_explore(dmvars_mo,testfx,"mags",22)

#season <- c("SPRING", "SUMMER", "FALL", "WINTER")
#sprmonths <- c("MarApr", "AprMay", "MayJunJul")
#summonths <- c("MayJun", "JunJul", "JulAug")
#winmonths <- c("NovDec","DecJan", "JanFeb")

#filename <- "2023_05_10_Daymet_monthly_analysis_2011_2017_r2s.csv"
#filename <- paste0("Analyses_Daymet_seasonal/2023_05_15_daymet_seasonal_analyses_", "ANNUAL", ".csv")
#write.csv(testing, file=filename)


#----------#### MAKING FIGURES: ####------------------------
library(lattice)
library(latticeExtra)
library(tactile)
library(mgcv)

vardat <- vardat
mo <- vardat[,8]
var.gam <- var.gam <- gam(vardat[,1]~s(mo),data=vardat)

dmplot<-xyplot(y ~ mo, data = vardat,
                panel = function(x, y) {
                  ci<-predict(var.gam, se=T)
                  ci$lower<-ci$fit-qt(0.975,var.gam$df.null)*ci$se.fit
                  ci$upper<-ci$fit+qt(0.975,var.gam$df.null)*ci$se.fit
                  l.ci<-cbind(var.gam$model$mo,ci$fit,ci$lower,ci$upper)
                  l<-l.ci[order(l.ci[,1]),]
                  panel.ci(l[,1],l[,2],l[,4],l[,3],
                           fill="seagreen3",alpha = 0.3)
                  panel.xyplot(x, y, pch=20,col="seagreen")
                  panel.lines(l[,1], l[,2],lty=1, col='black', lwd=1.5)
                  summ<-summary(var.gam)
                  r2 <- summ$r.sq
                  #f <- summ$fstatistic
                  # p <- pf(f[1],f[2],f[3],lower.tail=F)
                  panel.text(labels = bquote(italic(R)^2 ==.(format(r2,digits = 3))),
                             x=2150,y=-0.12,cex=0.75)
                  # panel.text(labels = bquote(italic(p)==.(format(p,digits = 3))),
                  #            x=1.2,y=-0.15,cex=0.75)
                },
                ylab="Disturbance Magnitude (TCG)",
                xlab="maxtemp feb/mar 2016",
)
print(dmplot)

#saving plots:
# tiff("Plots_509/vpdplot_jul2015.tiff", units="in", width=8, height=5, res=300)
# print(vpdplot)
# dev.off()

##---------------------------------------------------------------------####

# Making a new seasonal calculation script:

# 3 month combos version ----
seqfx<-function(x,y){
  seq(x,(y*12),by=12) #x = month (1=jan, 2=feb, etc.)
}

#grab variable of interest from daymet list:
var <- dmvars_mo[[1]]
#seasons:
spring <- c(3,4,5)
#make the season sequences:
for (i in length(sns)){
  
}
spring <- sort(c(seqfx(spring[1],8),
                 seqfx(spring[2],8),
                 seqfx(spring[3],8)))

