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
dat.out<-dat.out[dat.out$check==T,]
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
qplot(xpos,ypos,colour=dE1,data=loc.stats,size=I(5))

#Identify grains based on their central direction dE1 or dE2
hist(loc.stats$dE1,breaks=100)
medianDF<-data.frame(location=as.factor(loc.stats$location),dE1=loc.stats$dE1)
#fit<-hclust(dist(medianDF),method='complete')
#plot(fit) #There appear to be 8 clusters before they get too small, use kmeans to define them
fit2<-kmeans(medianDF$dE1,8)
loc.stats$grain<-fit2$cluster
qplot(xpos,ypos,colour=grain,data=loc.stats,size=I(5)) #make sure this matches the dE1 grain-map

#randomly sample 20 observations from each grain

#####################
#### Estimate central direction based on entire
#### data set, average over scans first
#####################

avg.scans <- ldply(dat.ests, function(x) {  
    data.frame(SE1=matrix(x$SE1,1,9),
               SE2=matrix(x$SE2,1,9))
})

nE<-nrow(avg.scans)

MedSamp<-as.SO3(data.matrix(avg.scans[,2:10]))
median(MedSamp)
region(MedSamp,method='moment',type='theory',estimator='median',alp=.05)*180/pi
region(MedSamp,method='moment',type='bootstrap',estimator='median',alp=.05)*180/pi


MeanSamp<-as.SO3(data.matrix(avg.scans[,11:19]))
mean(MeanSamp)
region(MeanSamp,method='moment',type='theory',estimator='mean',alp=.05)*180/pi
region(MeanSamp,method='moment',type='bootstrap',estimator='mean',alp=.05)*180/pi

#####################
#### Estimate grain-specific central direction based on entire
#### data set, average over scans first
#####################

grain.ests <- dlply(dat.out, .(grain), function(x) {
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
  
  data.frame(SE2=SE2, SE1=SE1)
})
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

region(exRots,method='moment',type='bootstrap',estimator='median',alp=.1)*180/pi
region(exRots,method='moment',type='bootstrap',estimator='mean',alp=.1)*180/pi

#####################
#### Randomly select locations within grains to
#### estimate witin grain precision
#####################

qplot(xpos,ypos,colour=dE,data=loc.stats,size=I(4))

#Select a grain by finding locations with similar dE1 values
grain1<-loc.stats[loc.stats$dE1<2.55 & loc.stats$dE1>2.45 & loc.stats$n>10,]
grain1Locs<-grain1$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain1Locs,])#where one grain map is grain 1?


#That didn't work, try to find grains based on xpos ypos & dE1
grain1<-loc.stats[loc.stats$ypos>5 & loc.stats$ypos<6,]
grain1<-grain1[grain1$xpos>8.5 & grain1$xpos<11,]
grain1<-grain1[grain1$dE<.008,]  #I think this gives us the .5 deg threshold Bingham mentions
grain1Locs<-grain1$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain1Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp<-dat.out[dat.out$location%in%grain1Locs & dat.out$rep==1,]
nSamp<-nrow(samp)
#sampRots<-as.SO3(data.matrix(samp[sample(1:nSamp,10,replace=F),3:11])) #randomly select 10 locations
sampRots<-as.SO3(data.matrix(samp[,3:11])) #take them all
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1) 
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1,col=2) 
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1,col=3) 

region(sampRots,method='moment',type='theory',estimator='median',alp=.05)*180/pi
region(sampRots,method='moment',type='bootstrap',estimator='median',alp=.05,m=300)*180/pi

plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1) 
plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1,col=2) 
plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1,col=3) 

region(sampRots,method='moment',type='theory',estimator='mean',alp=.05)*180/pi
region(sampRots,method='moment',type='bootstrap',estimator='mean',alp=.05)*180/pi


####
#Try a different grain

#That didn't work, try to find grains based on xpos ypos & dE1
grain2<-loc.stats[loc.stats$ypos>5 & loc.stats$ypos<6,]
grain2<-grain2[grain2$xpos>6 & grain2$xpos<9,]
grain2<-grain2[grain2$dE1>2.5 & grain2$dE1<2.9,]
grain2Locs<-grain2$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain2Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp<-dat.out[dat.out$location%in%grain2Locs & dat.out$rep==2,]
nSamp<-nrow(samp)
#sampRots<-as.SO3(data.matrix(samp[sample(1:nSamp,10,replace=F),3:11])) #randomly select 10 locations
sampRots<-as.SO3(data.matrix(samp[,3:11])) #take them all
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1) 
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1,col=2) 
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1,col=3) 

region(sampRots,method='moment',type='theory',estimator='median',alp=.01)*180/pi
region(sampRots,method='moment',type='bootstrap',estimator='median',alp=.01,m=300)*180/pi

plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1) 
plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1,col=2) 
plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1,col=3) 

region(sampRots,method='moment',type='theory',estimator='mean',alp=.01)*180/pi
region(sampRots,method='moment',type='bootstrap',estimator='mean',alp=.01)*180/pi
