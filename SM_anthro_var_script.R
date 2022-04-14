#This script is for wrangling the anthropogenic variables

#load datasets
#Human population density (county scale):
popdens<-read.csv("Anthro_Var_Data/County_Pop_Data_MA_CT_RI.csv")
#sample point counties:
pointcty<-read.csv("Anthro_Var_Data/point_counties.csv")

#match points to counties
cond.scores.mo<-cbind(cond.scores.mo,pointcty$NAME,pointcty$STATE_NAME)

