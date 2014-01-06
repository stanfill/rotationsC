#This gets the dataset in a unable condition

library(rotations)
library(plyr)
library(reshape2)
library(gridExtra)
library(xtable)

#g_legend will strip the legend and add it back so they can share one legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#load("/Users/stanfill/Dropbox/Rotation matrices/Melissa data/datasetnickel.RData")
load("C:/Users/Brittney Ritchey/Dropbox/Rotation matrices/Melissa data/datasetnickel.RData")

dat.out <- adply(data, .margins= c(1,3), function(x) {
  as.vector(x)
})
names(dat.out)[c(1,2)] <- c("rep", "location")
dat.out$xpos <- xpos[dat.out$location]
dat.out$ypos <- ypos[dat.out$location]
## convert data into a data frame including location and position information
#head(dat.out)

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
               dE1=mis.angle(x$SE1), dE2=mis.angle(x$SE2),
               dE=rot.dist(x$SE1, x$SE2, method="intrinsic", p=1),dR1=mis.angle(x$R1)
    )
})
loc.stats$xpos <- xpos[loc.stats$location]
loc.stats$ypos <- ypos[loc.stats$location]