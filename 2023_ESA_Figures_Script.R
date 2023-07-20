#Figures code for ESA 2023

#load libraries:
library(tidyverse)
library(ggplot2)
library(dplyr)


# Harvard Forest data - Landsat time series at some of the sites:
site <- 145
#tcg <- tcg.values
tcg <- cond.scores

#get just the tcg values from data frame:
tcgs<-tcg[,c(grep("^X",colnames(tcg)))]

#first get just june values:
#junes<-seq(3,ncol(tcgs),by=5)  #3 is if tcg/condscores data have april in them, 2 if data starts w/May
#tcgjune<-as.matrix(tcgs[,junes])

y <- as.numeric((tcgs[site,])/1000)
x <- 1:length(y)
#samp_time_series <- tcgjune[site,]
#par(mar=c(1,1,1,1))
plot(x, y, type = "l",
     xlab = "Growing Season Month",
     ylab = "Condition Score",
     )


# sample_site <- data.frame(
#   month = x,
#   score = y
# )
# colnames(sample_site) <- c("x","score")
# 
# #making a nice time series plot: COULD NOT GET IT TO PLOT
# time_series <- sample_site %>%
#   ggplot(aes(x="x", y="score")) +
#   geom_line() +
#   ylab("Condition Score") +
#   xlab("Growing Season Month")
# time_series

#testing subset sites:
distsites <- c(10,18,21,30,31,32,47,53,60,
               66,72,97,98,101,102,105,110,
               111,120,130,135,141,142,145,
               146,151,152,153,157,160,170,
               171,172,188,195,196,203)
distsdoub <- c(14,119,147,161,164,166,167,173,174,175,191,197)
fastrecov <- c(41,42,43,94,95)
slowrecov <- c(48,54,55,58,75,87)

ds <- tcgs[distsites, 21:79]
dd <- tcgs[distsdoub, 21:79]
fr <- tcgs[fastrecov, 21:79]
sr <- tcgs[slowrecov, 21:79]

ds_mr <- hf_data[distsites,]
dd_mr <- hf_data[distsdoub,]
fr_mr <- hf_data[fastrecov,]
sr_mr <- hf_data[slowrecov,]


