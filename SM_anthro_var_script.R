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

