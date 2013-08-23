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
library(rotations)
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
#g_legend will strip the legend and add it back so they can share one legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

load("/Users/stanfill/Dropbox/Rotation matrices/Melissa data/datasetnickel.RData")

dat.out <- adply(data, .margins= c(1,3), function(x) {
  as.vector(x)
})
names(dat.out)[c(1,2)] <- c("rep", "location")
dat.out$xpos <- xpos[dat.out$location]
dat.out$ypos <- ypos[dat.out$location]
## convert data into a data frame including location and position information
head(dat.out)

## are the matrices actually rotations?
checks <- adply(dat.out, .margins=1, function(x) {
  is.SO3(unlist(x[3:11]))
})
dat.out$check <- checks$V1

## compute the projected estimators at each location
dat.ests <- dlply(dat.out, .(location), function(x) {
  res <- na.omit(x)
  res <- subset(res, check==TRUE)
  
  n <- nrow(res) 
  SE2  <- SE1  <- NULL
  if (n == 1) {
    R <- as.SO3(matrix(unlist(res[1,3:11]), 3, 3))
    SE2  <- SE1  <- R
  } else if (n > 0) {
    rots <- as.SO3(as.matrix(res[,3:11]))
    SE2 <- mean(rots)
    SE1 <- median(rots)
  }
  location <- as.numeric(as.character(unique(x$location)))
  return(list(location=location, n=n, SE2=SE2, SE1=SE1))
})

## find distances between estimators and angles to identity for each
loc.stats <- ldply(dat.ests, function(x) {  
  location <- as.numeric(as.character(unique(x$location)))
  if (x$n > 0)
    data.frame(location=x$location, n=x$n, 
               dE1=angle(x$SE1), dE2=angle(x$SE2),
               dE=dist(x$SE1, x$SE2, method="intrinsic", p=1)
    )
})
loc.stats$xpos <- xpos[loc.stats$location]
loc.stats$ypos <- ypos[loc.stats$location]

#Plot grain map
qplot(xpos,ypos,colour=dE,data=loc.stats,size=I(5))

#####################
#### Find grain boundries, locations where 
#### mean and median disagree the most
#####################

#Find the location where the mean and median disagree most, this is likely a "contaminated" sample
ex<-loc.stats[which.max(loc.stats$dE),]$location 

#Randomly select a location with enough spread to be interesting
possibles<-which(loc.stats$dE>.1 & loc.stats$dE<.15)
ex<-loc.stats[sample(possibles,1),]$location #1539 has two clusters of size 9 and 5
  
exRots<-as.SO3(data.matrix(dat.out[dat.out$location==ex,3:11]))
#exRots<-as.SO3(exRots[-4,]) #potentially remove the observation that isn't stricly a rotation
p1<-plot(exRots,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,m=300)
p2<-plot(exRots,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,m=300,col=2)
p3<-plot(exRots,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,m=300,col=3)

grid.arrange(p1+theme(legend.position='none'),p2+theme(legend.position='none'),p3+theme(legend.position='none'),ncol=3)

#manually zoom-in on the interesting area
legend <- g_legend(p1) #pull the legend off one so it can be added to grid.arrange seperately (if desired)
grid.arrange(p1+theme(legend.position='none'),p2+theme(legend.position='none'),p3+theme(legend.position='none'),
             legend,ncol=4)

lims<-0.75
p1<-p1+xlim(c(-lims,lims))+ylim(c(-lims,lims))
p1
p2<-p2+xlim(c(-lims,lims))+ylim(c(-lims,lims))
p2
p3<-p3+xlim(c(-lims,lims))+ylim(c(-lims,lims))
p3
setwd("/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
#Can't use ggsave, need to use 'Export' function

#####################
#### Randomly select locations within grains to
#### estimate witin grain precision
#####################

#Select a grain by finding locations with similar dE1 values
grain1<-loc.stats[loc.stats$dE1<2.55 & loc.stats$dE1>2.45 & loc.stats$n>10,]
grain1Locs<-grain1$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain1Locs,])#where one grain map is grain 1?


#That didn't work, try to find grains based on xpos ypos
grain1<-loc.stats[loc.stats$ypos>5 & loc.stats$ypos<6,]
grain1<-grain1[grain1$xpos>8.5 & grain1$xpos<11,]
grain1<-grain1[grain1$dE1>2.4 &grain1$dE1<2.5,]
grain1Locs<-grain1$location

#where one grain map is grain 1?
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain1Locs,],xlim=c(0,12.5),ylim=c(0,10))

#take the first rep because it is the most reliable
samp<-dat.out[dat.out$location%in%grain1Locs & dat.out$rep==1,]
sampRots<-as.SO3(data.matrix(samp[,3:11]))
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1) 
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1,col=2) 
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1,col=3) 

region(sampRots,method='moment',type='theory',estimator='median',alp=.1)*180/pi
region(sampRots,method='moment',type='bootstrap',estimator='median',alp=.1)*180/pi

plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1) 
plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1,col=2) 
plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1,col=3) 

region(sampRots,method='moment',type='theory',estimator='mean',alp=.1)*180/pi
region(sampRots,method='moment',type='bootstrap',estimator='mean',alp=.1)*180/pi

