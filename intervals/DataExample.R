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

load("C:/Users/Brittney Ritchey/Dropbox/Rotation matrices/Melissa data/datasetnickel.RData")

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
  R1 <- SE2  <- SE1  <- NULL
  if (n == 1) {
    R <- as.SO3(matrix(unlist(res[1,3:11]), 3, 3))
    R1 <- SE2  <- SE1  <- R
  } else if (n > 0) {
    rots <- as.SO3(as.matrix(res[,3:11]))
    SE2 <- mean(rots)
    SE1 <- median(rots)
    R1 <- as.SO3(as.matrix(res[1,3:11]))
  }

  location <- as.numeric(as.character(unique(x$location)))
  return(list(location=location, n=n, SE2=SE2, SE1=SE1, R1=R1))
})

## find distances between estimators and angles to identity for each
loc.stats <- ldply(dat.ests, function(x) {  
  location <- as.numeric(as.character(unique(x$location)))
  if (x$n > 0)
    data.frame(location=x$location, n=x$n, 
               dE1=angle(x$SE1), dE2=angle(x$SE2),
               dE=dist(x$SE1, x$SE2, method="intrinsic", p=1),dR1=angle(x$R1)
    )
})
loc.stats$xpos <- xpos[loc.stats$location]
loc.stats$ypos <- ypos[loc.stats$location]

#Plot grain map
qplot(xpos,ypos,colour=dE1,data=loc.stats,size=I(5))

#Grain map based on median off all 14 scans
d <- ggplot(loc.stats, aes(xpos, ypos, color=dE1))
d2 <- d + geom_point(size=2.5) + scale_colour_gradient(expression(d[R](tilde(S)[E], I["3x3"])), low="grey99", high="grey10", limits=c(0, pi), breaks=c( pi/4, pi/2, 3*pi/4), labels=expression( pi/4, pi/2, 3*pi/4)) + 
  theme_bw() + xlab("") + ylab("") + coord_equal() + scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + 
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m)) + 
  #geom_point(shape="o", colour="yellow", size=5, data=loc.stats[idx,])  + 
  theme(plot.margin=unit(rep(0,4), "lines"))
d2
#ggsave(file="/Users/stanfill/Dropbox/Thesis/Intervals/Figures/grain-map.png", width=5.75, height=4.25)

#Grain map based on observation at scan1
d3 <- ggplot(loc.stats, aes(xpos, ypos, color=dR1))
d4 <- d3 + geom_point(size=2.5) + scale_colour_gradient(expression(r[gi]), low="grey99", high="grey10", limits=c(0, pi), breaks=c( pi/4, pi/2, 3*pi/4), labels=expression( pi/4, pi/2, 3*pi/4)) + 
  theme_bw() + xlab("") + ylab("") + coord_equal() + scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + 
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m)) + 
  #geom_point(shape="o", colour="yellow", size=5, data=loc.stats[idx,])  + 
  theme(plot.margin=unit(rep(0,4), "lines"))
d4
#ggsave(file="/Users/stanfill/Dropbox/Thesis/Intervals/Figures/grain-map-first-scan.png", width=5.75, height=4.25)

#Add labels for the identified grains
d4+annotate("text", x=c(3.5,7.5,.75,11,10,.75,7.5,.75),y=c(4,4.5,5,5,9,1,1,9),
            label=c(as.character(1:8)),colour='white',size=I(10))
#ggsave(file="C:/Users/Brittney Ritchey/Dropbox/Thesis/Intervals - Mean/Figures/grain-map-first-scan-labels.png", width=5.75, height=4.25)

#Grain map based on observation at scan1, degrees
d5 <- ggplot(loc.stats, aes(xpos, ypos, color=dR1*180/pi))
d6 <- d5 + geom_point(size=2.5) + scale_colour_gradient(expression(r[gi]), low="grey99", high="grey10", limits=c(0, 180), breaks=c(45,90,135), labels=expression(45^o,90^o,135^o)) + 
  theme_bw() + xlab("") + ylab("") + coord_equal() + scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + 
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m)) + 
  #geom_point(shape="o", colour="yellow", size=5, data=loc.stats[idx,])  + 
  theme(plot.margin=unit(rep(0,4), "lines"))+annotate("text", x=c(3.5,7.5,.75,11,10,.75,7.5,.75),y=c(4,4.5,5,5,9,1,1,9),
            label=c(as.character(1:8)),colour='white',size=I(10))
d6
#ggsave(file="C:/Users/Brittney Ritchey/Dropbox/Thesis/Intervals - Mean/Figures/grain-map-first-scan-labels.png", width=5.75, height=4.25)



##
#Identify grains based on their central direction dE1 or dE2
hist(loc.stats$dE1,breaks=100)
medianDF<-data.frame(location=as.factor(loc.stats$location),dE1=loc.stats$dE1)
#fit<-hclust(dist(medianDF),method='complete')
#plot(fit) #There appear to be 8 clusters before they get too small, use kmeans to define them
fit2<-kmeans(medianDF$dE1,8)
loc.stats$grain<-fit2$cluster
qplot(xpos,ypos,colour=grain,data=loc.stats,size=I(5)) #make sure this matches the dE1 grain-map

for(i in 1:length(dat.ests)){
  dat.ests[[i]]$grain<-loc.stats[loc.stats$location==dat.ests[[i]]$location,]$grain
}

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
region(MedSamp,method='moment',type='theory',estimator='median',alp=.01)*180/pi
region(MedSamp,method='moment',type='bootstrap',estimator='median',alp=.01)*180/pi


MeanSamp<-as.SO3(data.matrix(avg.scans[,11:19]))
mean(MeanSamp)
region(MeanSamp,method='moment',type='theory',estimator='mean',alp=.01)*180/pi
region(MeanSamp,method='moment',type='bootstrap',estimator='mean',alp=.01)*180/pi

#####################
#### Estimate grain-specific central direction based on entire
#### data set, average over scans first
#####################

#add grain identifier to dat.ests
avg.scans$grain<-1
avg.scans[avg.scans$location%in%loc.stats[loc.stats$grain==2,]$location,]$grain<-2
avg.scans[avg.scans$location%in%loc.stats[loc.stats$grain==3,]$location,]$grain<-3
avg.scans[avg.scans$location%in%loc.stats[loc.stats$grain==4,]$location,]$grain<-4
avg.scans[avg.scans$location%in%loc.stats[loc.stats$grain==5,]$location,]$grain<-5
avg.scans[avg.scans$location%in%loc.stats[loc.stats$grain==6,]$location,]$grain<-6
avg.scans[avg.scans$location%in%loc.stats[loc.stats$grain==7,]$location,]$grain<-7
avg.scans[avg.scans$location%in%loc.stats[loc.stats$grain==8,]$location,]$grain<-8


grain.ests <- ddply(avg.scans, .(grain), function(x) {
  
  SE1<-matrix(median(as.SO3(data.matrix(x[,2:10]))),1,9)
  SE2<-matrix(median(as.SO3(data.matrix(x[,11:19]))),1,9)
  
  data.frame(SE2=SE2, SE1=SE1)
})

gmeans<-as.SO3(data.matrix(grain.ests[,2:10]))
gmedians<-as.SO3(data.matrix(grain.ests[,11:19]))

region(gmeans,method='moment',type='theory',estimator='mean',alp=.01)*180/pi/sqrt(3383)
plot(gmeans,center=mean(gmeans),mean_regions='moment theory',alp=.1)

region(gmedians,method='moment',type='theory',estimator='median',alp=.01)*180/pi/sqrt(3383)
plot(gmedians,center=median(gmedians),median_regions='moment theory',alp=.1)

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

qplot(xpos,ypos,colour=dE2,data=loc.stats,size=I(4))

#Try to find grains based on xpos ypos & dE1
grain1<-loc.stats[loc.stats$xpos>0 & loc.stats$xpos<6,]
grain1<-grain1[grain1$ypos>2.5 & grain1$ypos<7.5,]
grain1<-grain1[grain1$dE1>2.5 & grain1$dE1<2.6,]  #I think this gives us the .5 deg threshold Bingham mentions
grain1Locs<-grain1$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain1Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp<-dat.out[dat.out$location%in%grain1Locs & dat.out$rep==1,]
samp<-samp[samp$V1<0 & samp$V2<0,]
nSamp<-nrow(samp)
nSamp
#sampRots<-as.SO3(data.matrix(samp[sample(1:nSamp,50,replace=F),3:11])) #randomly select 10 locations
sampRots<-as.SO3(data.matrix(samp[,3:11])) #take them all
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1) 
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1,col=2) 
plot(sampRots,center=median(sampRots),median_regions="all",alp=.1,col=3) 

region(sampRots,method='moment',type='theory',estimator='median',alp=.01)*180/pi*sqrt(nSamp)
region(sampRots,method='moment',type='bootstrap',estimator='median',alp=.01,m=500)*180/pi*sqrt(nSamp)

plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1) 
plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1,col=2) 
plot(sampRots,center=mean(sampRots),mean_regions="all",alp=.1,col=3) 

region(sampRots,method='moment',type='theory',estimator='mean',alp=.01)*180/pi*sqrt(nSamp)
region(sampRots,method='moment',type='bootstrap',estimator='mean',alp=.01)*180/pi*sqrt(nSamp)


####
#Try a different grain, grain 2

qplot(xpos,ypos,colour=dE2,data=loc.stats,size=I(4))

#That didn't work, try to find grains based on xpos ypos & dE1
grain2<-loc.stats[loc.stats$xpos>7.5 & loc.stats$xpos<12.5,]
grain2<-grain2[grain2$ypos>2.5 & grain2$ypos<8,]
grain2<-grain2[grain2$dE1>2.9 & grain2$dE1<3.0,]  #I think this gives us the .5 deg threshold Bingham mentions
grain2Locs<-grain2$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain2Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp2<-dat.out[dat.out$location%in%grain2Locs & dat.out$rep==1,]
samp2<-samp2[ samp2$V2<.5,]
nSamp2<-nrow(samp2)
nSamp2
#sampRots<-as.SO3(data.matrix(samp[sample(1:nSamp,10,replace=F),3:11])) #randomly select 10 locations
sampRots2<-as.SO3(data.matrix(samp2[,3:11])) #take them all
plot(sampRots2,center=median(sampRots2),median_regions="all",alp=.1) 
plot(sampRots2,center=median(sampRots2),median_regions="all",alp=.1,col=2) 
plot(sampRots2,center=median(sampRots2),median_regions="all",alp=.1,col=3) 

region(sampRots2,method='moment',type='theory',estimator='median',alp=.01)*180/pi*sqrt(nSamp2)
region(sampRots2,method='moment',type='bootstrap',estimator='median',alp=.01,m=300)*180/pi*sqrt(nSamp2)

plot(sampRots2,center=mean(sampRots2),mean_regions="all",alp=.1) 
plot(sampRots2,center=mean(sampRots2),mean_regions="all",alp=.1,col=2) 
plot(sampRots2,center=mean(sampRots2),mean_regions="all",alp=.1,col=3) 

region(sampRots2,method='moment',type='theory',estimator='mean',alp=.01)*180/pi*sqrt(nSamp2)
region(sampRots2,method='moment',type='bootstrap',estimator='mean',alp=.01)*180/pi*sqrt(nSamp2)

####Try another grain, grain 3

qplot(xpos,ypos,colour=dE2,data=loc.stats,size=I(4))

#Try to find grains based on xpos ypos & dE1
grain3<-loc.stats[loc.stats$xpos>0 & loc.stats$xpos<6,]
grain3<-grain3[grain3$ypos>2.5 & grain3$ypos<7.5,]
grain3<-grain3[grain3$dE1>1.5 & grain3$dE1<1.7,]  #I think this gives us the .5 deg threshold Bingham mentions
grain3Locs<-grain3$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain3Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp3<-dat.out[dat.out$location%in%grain3Locs & dat.out$rep==1,]
samp3<-samp3[samp3$V1>.4 & samp3$V1<.5,]
nSamp3<-nrow(samp3)
nSamp3
#sampRots3<-as.SO3(data.matrix(samp3[sample(1:nSamp3,50,replace=F),3:11])) #randomly select 50 locations
sampRots3<-as.SO3(data.matrix(samp3[,3:11])) #take them all


####Median
plot(sampRots3,center=median(sampRots3),median_regions="all",alp=.1) 
plot(sampRots3,center=median(sampRots3),median_regions="all",alp=.1,col=2) 
plot(sampRots3,center=median(sampRots3),median_regions="all",alp=.1,col=3) 

#confidence region radius for central direction
region(sampRots3,method='moment',type='theory',estimator='median',alp=.01)*180/pi
region(sampRots3,method='moment',type='bootstrap',estimator='median',alp=.01,m=500)*180/pi

#confidence region radius for single location
region(sampRots3,method='moment',type='theory',estimator='median',alp=.01)*180/pi*sqrt(nrow(sampRots3))
region(sampRots3,method='moment',type='bootstrap',estimator='median',alp=.01,m=1000)*180/pi*sqrt(nrow(sampRots3))


#standard error estimate of h vector
hTilNTH<-sqrt(region(sampRots3,method='moment',type='theory',estimator='median',alp=.01)^2*nSamp3/qchisq((1-.01),3))*180/pi
hTilNTH
hTilBOOT<-sqrt(region(sampRots3,method='moment',type='bootstrap',estimator='median',alp=.01)^2*nSamp3/qchisq((1-.01),3))*180/pi
hTilBOOT

#standard error estimate for r (very close to region radius / 2)
rMedSEnth<-sqrt(hTilNTH^2*3/nSamp3)
rMedSEnth
rMedSEboot<-sqrt(hTilBOOT^2*3/nSamp3)
rMedSEboot

####Mean
plot(sampRots3,center=mean(sampRots3),mean_regions="all",alp=.1) 
plot(sampRots3,center=mean(sampRots3),mean_regions="all",alp=.1,col=2) 
plot(sampRots3,center=mean(sampRots3),mean_regions="all",alp=.1,col=3) 

#confidence region radius (in SO(3)) for central direction
region(sampRots3,method='moment',type='theory',estimator='mean',alp=.01)*180/pi
region(sampRots3,method='moment',type='bootstrap',estimator='mean',alp=.01)*180/pi

#confidence region radius (in SO(3)) for single location!
region(sampRots3,method='moment',type='theory',estimator='mean',alp=.01)*180/pi*sqrt(nrow(sampRots3))
region(sampRots3,method='moment',type='bootstrap',estimator='mean',alp=.01,m=500)*180/pi*sqrt(nrow(sampRots3))


#standard error estimate of h vector (in R^3)
sqrt(region(sampRots3,method='moment',type='theory',estimator='mean',alp=.05)^2*nSamp3/qchisq((1-.05),3))*180/pi
sqrt(region(sampRots3,method='moment',type='bootstrap',estimator='mean',alp=.01)^2*nSamp3/qchisq((1-.01),3))*180/pi


####
#Try a different grain, grain 4

qplot(xpos,ypos,colour=dE2,data=loc.stats,size=I(4))

#That didn't work, try to find grains based on xpos ypos & dE1
grain4<-loc.stats[loc.stats$xpos>7.5 & loc.stats$xpos<12.5,]
grain4<-grain4[grain4$ypos>2.5 & grain4$ypos<8,]
grain4<-grain4[grain4$dE1>2.4 & grain4$dE1<2.5,]  #I think this gives us the .5 deg threshold Bingham mentions
grain4Locs<-grain4$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain4Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp4<-dat.out[dat.out$location%in%grain4Locs & dat.out$rep==1,]
nSamp4<-nrow(samp4)
nSamp4
#sampRots<-as.SO3(data.matrix(samp[sample(1:nSamp,10,replace=F),3:11])) #randomly select 10 locations
sampRots4<-as.SO3(data.matrix(samp4[,3:11])) #take them all
plot(sampRots4,center=median(sampRots4),median_regions="all",alp=.1) 
plot(sampRots4,center=median(sampRots4),median_regions="all",alp=.1,col=2) 
plot(sampRots4,center=median(sampRots4),median_regions="all",alp=.1,col=3) 

region(sampRots4,method='moment',type='theory',estimator='median',alp=.01)*180/pi*sqrt(nSamp4)
region(sampRots4,method='moment',type='bootstrap',estimator='median',alp=.01,m=500)*180/pi*sqrt(nSamp4)

plot(sampRots4,center=mean(sampRots4),mean_regions="all",alp=.1) 
plot(sampRots4,center=mean(sampRots4),mean_regions="all",alp=.1,col=2) 
plot(sampRots4,center=mean(sampRots4),mean_regions="all",alp=.1,col=3) 

region(sampRots4,method='moment',type='theory',estimator='mean',alp=.01)*180/pi*sqrt(nSamp4)
region(sampRots4,method='moment',type='bootstrap',estimator='mean',alp=.01)*180/pi*sqrt(nSamp4)

####
#Try a different grain, grain 5

qplot(xpos,ypos,colour=dE2,data=loc.stats,size=I(4))

#That didn't work, try to find grains based on xpos ypos & dE1
grain5<-loc.stats[loc.stats$xpos>7.5 & loc.stats$xpos<12.5,]
grain5<-grain5[grain5$ypos>7.5 & grain5$ypos<10,]
grain5<-grain5[ grain5$dE1<1.5,]  #I think this gives us the .5 deg threshold Bingham mentions
grain5Locs<-grain5$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain5Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp5<-dat.out[dat.out$location%in%grain5Locs & dat.out$rep==1,]
#samp5<-samp5[ samp4$V2<.5,]
nSamp5<-nrow(samp5)
nSamp5
#sampRots<-as.SO3(data.matrix(samp[sample(1:nSamp,10,replace=F),3:11])) #randomly select 10 locations
sampRots5<-as.SO3(data.matrix(samp5[,3:11])) #take them all
plot(sampRots5,center=median(sampRots5),median_regions="all",alp=.1) 
plot(sampRots5,center=median(sampRots5),median_regions="all",alp=.1,col=2) 
plot(sampRots5,center=median(sampRots5),median_regions="all",alp=.1,col=3) 

region(sampRots5,method='moment',type='theory',estimator='median',alp=.01)*180/pi*sqrt(nSamp5)
region(sampRots5,method='moment',type='bootstrap',estimator='median',alp=.01,m=500)*180/pi*sqrt(nSamp5)

plot(sampRots5,center=mean(sampRots5),mean_regions="all",alp=.1) 
plot(sampRots5,center=mean(sampRots5),mean_regions="all",alp=.1,col=2) 
plot(sampRots5,center=mean(sampRots5),mean_regions="all",alp=.1,col=3) 

region(sampRots5,method='moment',type='theory',estimator='mean',alp=.01)*180/pi*sqrt(nSamp5)
region(sampRots5,method='moment',type='bootstrap',estimator='mean',alp=.01,m=500)*180/pi*sqrt(nSamp5)

####
#Try a different grain, grain 6

qplot(xpos,ypos,colour=dE2,data=loc.stats,size=I(4))

#That didn't work, try to find grains based on xpos ypos & dE1
grain6<-loc.stats[ loc.stats$xpos<5,]
grain6<-grain6[ grain6$ypos<4,]
grain6<-grain6[ grain6$dE1>2.5 & grain6$dE1<2.7,]  #I think this gives us the .5 deg threshold Bingham mentions
grain6Locs<-grain6$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain6Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp6<-dat.out[dat.out$location%in%grain6Locs & dat.out$rep==1,]
samp6<-samp6[ samp6$V2>.5,]
nSamp6<-nrow(samp6)
nSamp6
#sampRots<-as.SO3(data.matrix(samp[sample(1:nSamp,10,replace=F),3:11])) #randomly select 10 locations
sampRots6<-as.SO3(data.matrix(samp6[,3:11])) #take them all
plot(sampRots6,center=median(sampRots6),median_regions="all",alp=.1) 
plot(sampRots6,center=median(sampRots6),median_regions="all",alp=.1,col=2) 
plot(sampRots6,center=median(sampRots6),median_regions="all",alp=.1,col=3) 

region(sampRots6,method='moment',type='theory',estimator='median',alp=.01)*180/pi*sqrt(nSamp6)
region(sampRots6,method='moment',type='bootstrap',estimator='median',alp=.01,m=500)*180/pi*sqrt(nSamp6)

plot(sampRots6,center=mean(sampRots6),mean_regions="all",alp=.1) 
plot(sampRots6,center=mean(sampRots6),mean_regions="all",alp=.1,col=2) 
plot(sampRots6,center=mean(sampRots6),mean_regions="all",alp=.1,col=3) 

region(sampRots6,method='moment',type='theory',estimator='mean',alp=.01)*180/pi*sqrt(nSamp6)
region(sampRots6,method='moment',type='bootstrap',estimator='mean',alp=.01,m=500)*180/pi*sqrt(nSamp6)

####
#Try a different grain, grain 7

qplot(xpos,ypos,colour=dE2,data=loc.stats,size=I(4))

#That didn't work, try to find grains based on xpos ypos & dE1
grain7<-loc.stats[ loc.stats$xpos>5 & loc.stats$xpos<10,]
grain7<-grain7[ grain7$ypos<4,]
grain7<-grain7[ grain7$dE1>2.75 & grain7$dE1<2.95,]  #I think this gives us the .5 deg threshold Bingham mentions
grain7Locs<-grain7$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain7Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp7<-dat.out[dat.out$location%in%grain7Locs & dat.out$rep==1,]
samp7<-samp7[ samp7$V1>.5,]
nSamp7<-nrow(samp7)
nSamp7
#sampRots<-as.SO3(data.matrix(samp[sample(1:nSamp,10,replace=F),3:11])) #randomly select 10 locations
sampRots7<-as.SO3(data.matrix(samp7[,3:11])) #take them all
plot(sampRots7,center=median(sampRots7),median_regions="all",alp=.1) 
plot(sampRots7,center=median(sampRots7),median_regions="all",alp=.1,col=2) 
plot(sampRots7,center=median(sampRots7),median_regions="all",alp=.1,col=3) 

region(sampRots7,method='moment',type='theory',estimator='median',alp=.01)*180/pi*sqrt(nSamp7)
region(sampRots7,method='moment',type='bootstrap',estimator='median',alp=.01,m=500)*180/pi*sqrt(nSamp7)

plot(sampRots7,center=mean(sampRots7),mean_regions="all",alp=.1) 
plot(sampRots7,center=mean(sampRots7),mean_regions="all",alp=.1,col=2) 
plot(sampRots7,center=mean(sampRots7),mean_regions="all",alp=.1,col=3) 

region(sampRots7,method='moment',type='theory',estimator='mean',alp=.01)*180/pi*sqrt(nSamp7)
region(sampRots7,method='moment',type='bootstrap',estimator='mean',alp=.01,m=500)*180/pi*sqrt(nSamp7)


####
#Try a different grain, grain 8

qplot(xpos,ypos,colour=dE2,data=loc.stats,size=I(4))

#That didn't work, try to find grains based on xpos ypos & dE1
grain8<-loc.stats[ loc.stats$xpos<2.5,]
grain8<-grain8[ grain8$ypos>7.5,]
grain8<-grain8[ grain8$dE1>2 & grain8$dE1<2.5,]  #I think this gives us the .5 deg threshold Bingham mentions
grain8Locs<-grain8$location
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain8Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))

#take the first rep because it is the most reliable
samp8<-dat.out[dat.out$location%in%grain8Locs & dat.out$rep==1,]
samp8<-samp8[ samp8$V1>.79,]
nSamp8<-nrow(samp8)
nSamp8
#sampRots<-as.SO3(data.matrix(samp[sample(1:nSamp,10,replace=F),3:11])) #randomly select 10 locations
sampRots8<-as.SO3(data.matrix(samp8[,3:11])) #take them all
plot(sampRots8,center=median(sampRots8),median_regions="all",alp=.1) 
plot(sampRots7,center=median(sampRots7),median_regions="all",alp=.1,col=2) 
plot(sampRots7,center=median(sampRots7),median_regions="all",alp=.1,col=3) 

region(sampRots8,method='moment',type='theory',estimator='median',alp=.01)*180/pi*sqrt(nSamp8)
region(sampRots8,method='moment',type='bootstrap',estimator='median',alp=.01,m=500)*180/pi*sqrt(nSamp8)

plot(sampRots7,center=mean(sampRots7),mean_regions="all",alp=.1) 
plot(sampRots7,center=mean(sampRots7),mean_regions="all",alp=.1,col=2) 
plot(sampRots7,center=mean(sampRots7),mean_regions="all",alp=.1,col=3) 

region(sampRots8,method='moment',type='theory',estimator='mean',alp=.01)*180/pi*sqrt(nSamp8)
region(sampRots8,method='moment',type='bootstrap',estimator='mean',alp=.01,m=500)*180/pi*sqrt(nSamp8)


#Slideshow to identify all the grains so far:
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain1Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain2Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain3Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain4Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain5Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain6Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain7Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))
qplot(xpos,ypos,colour=dE1,data=loc.stats[loc.stats$location%in%grain8Locs,],xlim=c(0,12),ylim=c(0,10.04591),size=I(4))


####################
#Formally identify grains and compute 4 CRs for each (doesn't work)

scan1<-dat.out[dat.out$check==T & dat.out$rep==1,]
hist(scan1$V2,breaks=100)
fit<-kmeans(scan1$V1,50)
scan1$cluster<-as.factor(fit$cluster)
qplot(xpos,ypos,colour=cluster,data=scan1,size=I(4))

rots<-as.SO3(data.matrix(scan1[scan1$cluster==1,3:11]))
n<-nrow(rots)
MeanMom<-region(rots,method='moment',type='theory',estimator='mean',alp=.01)
MeanMom<-MeanMom*180/pi*sqrt(n)
MeanMom

crs<-ddply(scan1,.(cluster),function(res) {
  
  rots <- as.SO3(as.matrix(res[,3:11]))
  n<-nrow(rots)
  MeanMom<-region(rots,method='moment',type='theory',estimator='mean',alp=.01)
  MeanMom<-MeanMom*180/pi*sqrt(n)

  data.frame(cluster=res$cluster[1], n=n, MeanMom=MeanMom)
})
crs