library(plyr)
library(reshape2)
library(ggplot2)

#Cayley
cayRes<-read.csv("intervals/BayesRegions/Results/CayleyResults.csv")[,-1]
caySum<-ddply(cayRes,.(Kappa,n),summarize,Perc=sum(cover)/length(cover),Width=mean(width),Sacc=mean(Saccept))
caySum
qplot(n,Perc,data=caySum,facets=.~Kappa,geom='line')+geom_hline(yintercept=0.9)

#Fisher
fishRes<-read.csv("intervals/BayesRegions/Results/FisherResults.csv")[,-1]
fishSum<-ddply(fishRes,.(Kappa,n),summarize,Perc=sum(cover)/length(cover),Width=mean(width),Sacc=mean(Saccept))
fishSum
qplot(n,Perc,data=fishSum,facets=.~Kappa,geom='line')+geom_hline(yintercept=0.9)

#Together
caySum$nu<-rep(c(0.75,.5,.25),each=4)
caySum$Dist="Cayley"
fishSum$nu<-rep(c(0.75,.5,.25),each=4)
fishSum$Dist<-"Fisher"
togSum<-rbind(caySum,fishSum)

qplot(n,Perc,data=togSum,facets=.~nu,geom='line',group=Dist,colour=Dist)+geom_hline(yintercept=0.9)

qplot(n,Width,data=togSum,facets=.~nu,geom='line',group=Dist,colour=Dist)
