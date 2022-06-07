#script for disturbance probability exploratory analyses

#load condition scores (finding where prob of dist = 0,1)
cond.scores.mo<-read.csv("2020_07_10_sample_score_mean_MONTHLY.csv")
condscores<-cond.scores.mo[,2:131]

#load tcg data for column numbers:
distmagsdata<-read.csv("SM_distmagrecov_data.csv")


#adding 0,1 for disturbance occurrence:
for (i in nrow(condscores)){
  
}