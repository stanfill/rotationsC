library(rotations)

############
#Try the drill data of Rancourt's Disertation
############
data(drilldata)

#Remove NAs
drilldata<-na.omit(drilldata)
SJP<-paste(drilldata$Subject,drilldata$Joint,drilldata$Position)
length(unique(SJP))
qplot(Subject,Replicate,data=drilldata)+facet_grid(Joint~Position,labeller=label_both)


############
#Try the nickel data of Bingham's Disertation
############
library(ggplot2)
data(nickel)
nickel$location<-as.factor(paste(nickel$xpos,nickel$ypos))
rep1<-nickel[nickel$rep==1,]
qplot(xpos,ypos,data=rep1)