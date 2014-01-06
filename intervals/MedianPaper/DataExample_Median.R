############
#Try the nickel data of Bingham's Disertation
############
source("intervals/MedianPaper/PrepareData.R")

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

#####################
#### Draw eyeballs for worst sample to illustrate mean/median difference
#####################

#Find the location where the mean and median disagree most, this is likely a "contaminated" sample
ex<-loc.stats[which.max(loc.stats$dE),]$location #Should be 698
  
exRots<-as.SO3(data.matrix(dat.out[dat.out$location==ex,3:11]))
#exRots<-as.SO3(exRots[-4,]) #potentially remove the observation that isn't stricly a rotation

#Use new plot function to show all columns at once
plot(exRots,col=c(1,2,3),center=median(exRots),show_estimates=c('proj.mean','proj.median'))
plot(exRots,col=c(1,2,3),center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='theory',mean_regions='moment theory',alp=.1)
plot(exRots,col=c(1,2,3),center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,m=100)


#Manually put the three axes on the same plot
p1<-plot(exRots,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,m=300)
p2<-plot(exRots,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,m=300,col=2)
p3<-plot(exRots,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,m=300,col=3)
p4<-g_legend(p1)
p1<-p1+theme(legend.position='none')
p2<-p2+theme(legend.position='none')
p3<-p3+theme(legend.position='none')

grid.arrange(p1,p2,p3,p4,nrow=1,widths=c(2,2,2,2))
grid.arrange(p1,p2,p3,p4,nrow=2,widths=c(2,2,2,1))
#ggsave doesn't work, use "Export." I used "width=900" and "height=300"

#Manually put the three axes on the same plot, no regions
p1<-plot(exRots,col=2,center=median(exRots),show_estimates=c('proj.mean','proj.median'))+theme(legend.position='none')
p2<-plot(exRots,col=2,center=median(exRots),show_estimates=c('proj.mean','proj.median'),to_range=T)+theme(legend.position='none')
grid.arrange(p1,p2,nrow=1,widths=c(2,2))
#ggsave doesn't work, use "Export." I used "width=600" and "height=300"

#Just the data
p1<-plot(exRots,col=2,center=median(exRots))+theme(legend.position='none')
p2<-plot(exRots,col=2,center=median(exRots),to_range=T)+theme(legend.position='none')
grid.arrange(p1,p2,nrow=1,widths=c(2,2))

#With regions
p1<-plot(exRots,col=2,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,m=300)+theme(legend.position='none')
p2<-plot(exRots,col=2,center=median(exRots),show_estimates=c('proj.mean','proj.median'),median_regions='bootstrap',mean_regions='moment bootstrap',alp=.1,m=300,to_range=T)+theme(legend.position='none')
grid.arrange(p1,p2,nrow=1,widths=c(2,2))

#####################
#### Find all locations where mean/median differ substantially
#####################

#Randomly select a location with enough spread to be interesting
possibles<-which(loc.stats$dE>.1 & loc.stats$dE<.15)
possibles<-sort(c(possibles,which.max(loc.stats$dE)))
m<-length(possibles)
DataExDF<-data.frame(Location=rep(0,m),MeanNTH=rep(0,m),MedianNTH=rep(0,m),MeanBoot=rep(0,m),MedianBoot=rep(0,m))

for(i in 1:m){
  ex<-loc.stats[possibles[i],]$location 
  DataExDF$Location[i]<-ex
  exRots<-as.SO3(data.matrix(dat.out[dat.out$location==ex,3:11]))
  DataExDF$MeanNTH[i]<-region(exRots,method='moment',type='theory',estimator='mean',alp=.1)*180/pi
  DataExDF$MedianNTH[i]<-region(exRots,method='moment',type='theory',estimator='median',alp=.1)*180/pi
  DataExDF$MeanBoot[i]<-region(exRots,method='moment',type='bootstrap',estimator='mean',alp=.1,m=500)*180/pi
  DataExDF$MedianBoot[i]<-region(exRots,method='moment',type='bootstrap',estimator='median',alp=.1,m=500)*180/pi
}

DataExDF$Location<-as.factor(DataExDF$Location)
xtable(DataExDF,digits=3)
#####################
#### Find all locations where mean/median differ substantially, highlight them in grain map
#####################

#Grain map based on median off all 14 scans
d <- ggplot(loc.stats, aes(xpos, ypos, color=dE1))
d2 <- d + geom_point(size=4) + scale_colour_gradient(expression(d[R](tilde(S)[E], I["3x3"])), low="grey99", high="grey10", limits=c(0, pi), breaks=c( pi/4, pi/2, 3*pi/4), labels=expression( pi/4, pi/2, 3*pi/4)) + 
  theme_bw() + xlab("") + ylab("") + coord_equal() + scale_x_continuous(limits=c(0, 12.5), breaks=seq(0, 12.5, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m, 12.5*mu*m)) + 
  scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=2.5), labels=expression(0*mu*m, 2.5*mu*m, 5*mu*m, 7.5*mu*m, 10*mu*m)) + 
  geom_point(shape="o", colour="yellow", size=5, data=loc.stats[possibles,])  + 
  theme(plot.margin=unit(rep(0,4), "lines"))
d2
ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/Grain_map_with_circles.pdf",height=8,width=9)
