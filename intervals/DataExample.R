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

#Find the location where the mean and median disagree most, this is likely a "contaminated" sample
ex<-loc.stats[which.max(loc.stats$dE),]$location 
#Randomly select a location with enough spread to be interesting
possibles<-which(loc.stats$dE>.1)
ex<-loc.stats[sample(possibles,1),]$location
  
exRots<-as.SO3(data.matrix(dat.out[dat.out$location==ex,3:11]))
#exRots<-as.SO3(exRots[-4,]) #potentially remove the observation that isn't stricly a rotation
p1<-plot(exRots,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1)
p2<-plot(exRots,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,col=2)
p3<-plot(exRots,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,col=3)

legend <- g_legend(p1)

setwd("/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
grid.arrange(p1+theme(legend.position='none'),p2+theme(legend.position='none'),p3+theme(legend.position='none'),ncol=3)
grid.arrange(p1+theme(legend.position='none'),p2+theme(legend.position='none'),p3+theme(legend.position='none'),
             legend,ncol=4)

