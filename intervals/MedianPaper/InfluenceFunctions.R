IFMean<-function(r,kap,Fisher=T,SIF=F){
  
  #SIF - Standardize by Asymptotic variance?
  
  if(Fisher){
    d<-(((kap+1)/kap)*besselI(2*kap,1)-besselI(2*kap,0))/(3*(besselI(2*kap,0)-besselI(2*kap,1)))
    c<-(((kap+1)/kap)*besselI(2*kap,1)-besselI(2*kap,0))/(3*kap*(besselI(2*kap,0)-besselI(2*kap,1)))
  }else{
    d<-kap/(kap+2)
    c<-(4*kap+2)/(kap^2+5*kap+6)
  }
  
  if(SIF){
    return(2*d*sin(r)/c)
  }else{
    return(sin(r)/d)
  }
}

IFMedian<-function(r,kap,Fisher=T,SIF=F){
  fn<-sin(r)/sqrt(1-cos(r))
  
  if(Fisher){
    d<-exp(2*kap)*(6*sqrt(kap)-(3+8*kap)*dawson(2*sqrt(kap)))/(24*sqrt(2)*pi*kap^1.5*(besselI(2*kap,0)-besselI(2*kap,1)))
    c<-besselI(2*kap,1)/(12*kap*(besselI(2*kap,0)-besselI(2*kap,1)))
  }else{
    d<-sqrt(2/pi)*kap*gamma(kap+2)/(3*gamma(kap+2.5))
    c<-(2*kap+1)/(6*kap+12)
  }
  
  if(SIF){
    return(d*fn/c)
  }else{
    return(fn/(2*d))
  }
}

GESMean<-function(kap,Fisher=F,SIF=F,IF=F){
  
  if(Fisher){
    d<-(((kap+1)/kap)*besselI(2*kap,1)-besselI(2*kap,0))/(3*(besselI(2*kap,0)-besselI(2*kap,1)))
    c<-(((kap+1)/kap)*besselI(2*kap,1)-besselI(2*kap,0))/(3*kap*(besselI(2*kap,0)-besselI(2*kap,1)))
    #x<-2*d*d/c
    x<-(2/3)*((kap+1)*besselI(2*kap,1)-kap*besselI(2*kap,0))/(besselI(2*kap,0)-besselI(2*kap,1))
  }else{
    d<-kap/(kap+2)
    c<-(4*kap+2)/(kap^2+5*kap+6)
    x<-kap^2/(2*kap-1)
  }
  
  
  if(SIF){
    #Standardize
    if(IF){
      #Information standardized
      return(sqrt(x/d^2))
    }else{
      #Self Standardized
      return(sqrt(2/c))
    }
  }else{
    #Don't standardize
    return(1/d)
  }
}


GESMedian<-function(kap,Fisher=F,SIF=F,IF=F){
  
  if(Fisher){
    d<-exp(2*kap)*(6*sqrt(kap)-(3+8*kap)*dawson(2*sqrt(kap)))/(24*sqrt(2)*pi*kap^1.5*(besselI(2*kap,0)-besselI(2*kap,1)))
    c<-besselI(2*kap,1)/(12*kap*(besselI(2*kap,0)-besselI(2*kap,1)))
    #x<-2*d*d/c
    x<-(2/3)*((kap+1)*besselI(2*kap,1)-kap*besselI(2*kap,0))/(besselI(2*kap,0)-besselI(2*kap,1))
  }else{
    d<-sqrt(2/pi)*kap*gamma(kap+2)/(3*gamma(kap+2.5))
    c<-(2*kap+1)/(6*kap+12)
    x<-kap^2/(2*kap-1)
  }
  
  if(SIF){
    #Standardize
    if(IF){
      #Information standardized
      return(sqrt(x/(2*d^2)))
    }else{
      #Self Standardized
      return(1/sqrt(c))
    }
  }else{
    #Don't standardize
    return(1/(sqrt(2)*d))
  }
}

#############
#Compare IFs
library(gsl)
rs<-seq(-pi,pi,length=100)
kap<-20
plot(rs,IFMedian(rs,kap,T),main='Fisher Distribution',type='l',xlab='r',ylab="IF")
abline(v=0,h=0)
lines(rs,IFMean(rs,kap,T),col=2)

plot(rs,IFMedian(rs,kap,F),main='Cayley Distribution',type='l',xlab='r',ylab="IF")
abline(v=0,h=0)
lines(rs,IFMean(rs,kap,F),col=2)

#############
#Compare IFs as a function of kappa, Cayley
library(ggplot2)
library(reshape2)
library(plyr)
rs<-seq(-pi,pi,length=100)
kap<-c(1,5,20)

IFDF<-data.frame(r=rep(rs,3),kappa=as.factor(rep(kap,each=100)),Mean=c(IFMean(rs,kap[1],F),IFMean(rs,kap[2],F),IFMean(rs,kap[3],F)),
                     Median=c(IFMedian(rs,kap[1],F),IFMedian(rs,kap[2],F),IFMedian(rs,kap[3],F)))

#Plot mean and median seperately
qplot(r,Mean,data=IFDF,geom='line',colour=kappa,group=kappa,ylab="IF-Mean")+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')

qplot(r,Median,data=IFDF,geom='line',colour=kappa,group=kappa,ylab="IF-Median")+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')

#Put both estimators on one plot
IFDFm<-melt(IFDF,id=c("r","kappa"))
colnames(IFDFm)[3:4]<-c("Estimator","IF")
qplot(r,IF,data=IFDFm,geom='line',colour=kappa,linetype=Estimator)+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')

#Facet by estimator
qplot(r,IF,data=IFDFm,geom='line',colour=kappa,facets=.~Estimator)+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')

#Facet by concentration
qplot(r,IF,data=IFDFm,geom='line',colour=Estimator,ylab='IF(r,Cayley)',size=I(1.25))+theme_bw()+geom_vline(xintercept=0,colour='gray50')+
  geom_hline(yintercept=0,colour='gray50')+theme(aspect.ratio=1)+facet_grid(.~kappa, labeller = label_bquote(kappa==.(x)))+
  scale_x_continuous(breaks=c(-pi,-pi/2,0,pi/2,pi),labels=c(expression(-pi,-frac(pi,2),0,frac(pi,2),pi)))
#ggsave("C:/Users/Brittney Ritchey/Dropbox/Thesis/Intervals - Median/Figures/CayleyIF.pdf")

### Same with Fisher
kap<-c(1,2.5,5)
IFDFFish<-data.frame(r=rep(rs,3),kappa=as.factor(rep(kap,each=100)),Mean=c(IFMean(rs,kap[1],T),IFMean(rs,kap[2],T),IFMean(rs,kap[3],T)),
                 Median=c(IFMedian(rs,kap[1],T),IFMedian(rs,kap[2],T),IFMedian(rs,kap[3],T)))

IFDFFishm<-melt(IFDFFish,id=c("r","kappa"))
colnames(IFDFFishm)[3:4]<-c("Estimator","IF")

qplot(r,IF,data=IFDFFishm,geom='line',colour=Estimator,ylab='IF(r,Fisher)',size=I(1.25))+theme_bw()+geom_vline(xintercept=0,colour='gray50')+
  geom_hline(yintercept=0,colour='gray50')+theme(aspect.ratio=1)+facet_grid(.~kappa, labeller = label_bquote(kappa==.(x)))
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/FisherIF.pdf",width=4,height=6,units="in")


#############
#Compare SIFs as a function of kappa, Cayley
library(ggplot2)
library(reshape2)
library(plyr)
rs<-seq(-pi,pi,length=100)
kap<-c(1,5,20)

SIFDF<-data.frame(r=rep(rs,3),kappa=as.factor(rep(kap,each=100)),Mean=c(IFMean(rs,kap[1],F,T),IFMean(rs,kap[2],F,T),IFMean(rs,kap[3],F,T)),
                 Median=c(IFMedian(rs,kap[1],F,T),IFMedian(rs,kap[2],F,T),IFMedian(rs,kap[3],F,T)))
SIFDFm<-melt(IFDF,id=c("r","kappa"))
colnames(SIFDFm)[3:4]<-c("Estimator","SIF")
#Facet by concentration
qplot(r,SIF,data=SIFDFm,geom='line',colour=Estimator,ylab='SIF(r,Cayley)')+theme_bw()+geom_vline(xintercept=0,colour='gray50')+
  geom_hline(yintercept=0,colour='gray50')+theme(aspect.ratio=1)+facet_grid(.~kappa, labeller = label_bquote(kappa==.(x)))


#Same For Fisher
SIFDF<-data.frame(r=rep(rs,3),kappa=as.factor(rep(kap,each=100)),Mean=c(IFMean(rs,kap[1],T,T),IFMean(rs,kap[2],T,T),IFMean(rs,kap[3],T,T)),
                  Median=c(IFMedian(rs,kap[1],T,T),IFMedian(rs,kap[2],T,T),IFMedian(rs,kap[3],T,T)))
SIFDFm<-melt(IFDF,id=c("r","kappa"))
colnames(SIFDFm)[3:4]<-c("Estimator","SIF")
#Facet by concentration
qplot(r,SIF,data=SIFDFm,geom='line',colour=Estimator,ylab='SIF(r,Fisher)')+theme_bw()+geom_vline(xintercept=0,colour='gray50')+
  geom_hline(yintercept=0,colour='gray50')+theme(aspect.ratio=1)+facet_grid(.~kappa, labeller = label_bquote(kappa==.(x)))

#############
#Compare SIFs as a function of kappa standardized by KL divergence, Fisher
library(ggplot2)
library(reshape2)
library(plyr)
rs<-seq(-pi,pi,length=100)

kap<-c(1,2.5,5)
KF<-sqrt(kap*besselI(kap,0)/besselI(kap,1))
IFDFFish<-data.frame(r=rep(rs,3),kappa=as.factor(rep(kap,each=100)),Mean=c(IFMean(rs,kap[1],T),IFMean(rs,kap[2],T),IFMean(rs,kap[3],T)),
                     Median=c(IFMedian(rs,kap[1],T),IFMedian(rs,kap[2],T),IFMedian(rs,kap[3],T)))
IFDFFish

#############
#Compare estimator sensitivity as a function of kappa
library(ggplot2)
library(reshape2)
library(plyr)
library(gsl)
library(grid)
kap<-seq(.01,7,length=100)

#Use the unstandardized IF
GESDF<-data.frame(kappa=rep(kap,2),Dist=rep(c("Cayley","Fisher"),each=100),
                  Median=c(GESMedian(kap,F),GESMedian(kap,T)),Mean=c(GESMean(kap,F),GESMean(kap,T)))

GESdfm<-melt(GESDF,id=c("kappa","Dist"))
colnames(GESdfm)[3:4]<-c("Estimator","GES")
#Facet by Estimator
qplot(kappa,GES,data=GESdfm,geom='line',colour=Estimator,facets=.~Dist)+theme_bw()

#Use the self-standardized IF
SGESDF<-data.frame(kappa=rep(kap,2),Dist=rep(c("Cayley","Fisher"),each=100),
                  Median=c(GESMedian(kap,F,T),GESMedian(kap,T,T)),Mean=c(GESMean(kap,F,T),GESMean(kap,T,T)))

SGESdfm<-melt(SGESDF,id=c("kappa","Dist"))
colnames(SGESdfm)[3:4]<-c("Estimator","SGES")
#Facet by Estimator
qplot(kappa,SGES,data=SGESdfm,geom='line',colour=Estimator,facets=.~Dist)+theme_bw()

#Use the information-standardized IF
ISGESDF<-data.frame(kappa=rep(kap,2),Dist=rep(c("Cayley","matrix Fisher"),each=100),
                   Median=c(GESMedian(kap,F,T,T),GESMedian(kap,T,T,T)),Mean=c(GESMean(kap,F,T,T),GESMean(kap,T,T,T)))
ISGESDF$Ratio<-ISGESDF$Median/ISGESDF$Mean

ISGESdfm<-melt(ISGESDF,id=c("kappa","Dist"))
colnames(ISGESdfm)[3:4]<-c("Estimator","ISGES")
#Facet by Estimator
qplot(kappa,ISGES,data=ISGESdfm[ISGESdfm$Estimator!='Ratio',],size=I(1.25),geom='line',colour=Estimator,ylim=c(1,6),xlab=expression(kappa),ylab="SGES")+
  theme_bw()+facet_grid(.~Dist)+theme(aspect.ratio=1,panel.margin=unit(2, "lines"),legend.position=c(.9,.8))+geom_vline(xintercept=0,colour='gray')+
  geom_hline(yintercept=3*sqrt(2*pi)/4,colour='gray')
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/SGES.pdf",width=6.5,height=3.5)

qplot(kappa,ISGES,data=ISGESdfm[ISGESdfm$Estimator!='Ratio',],size=I(1.25),geom='line',linetype=Estimator,ylim=c(1,6),xlab=expression(kappa),ylab="SGES")+
  theme_bw()+facet_grid(.~Dist)+theme(aspect.ratio=1,panel.margin=unit(2, "lines"),legend.position=c(.9,.8))+geom_vline(xintercept=0,colour='gray')+
  geom_hline(yintercept=3*sqrt(2*pi)/4,colour='gray')
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/SGESbw.pdf",width=6.5,height=3.5)



qplot(kappa,ISGES,data=ISGESdfm[ISGESdfm$Estimator=='Ratio',],geom='line',facets=.~Dist)+theme_bw()

#################
#Unfair versus fair comparison

### Same with Fisher
kap<-c(20,5)
IFDFFish<-data.frame(r=rep(rs,2),kappa=as.factor(rep(kap,each=100)),Mean=c(IFMean(rs,kap[1],F),IFMean(rs,kap[2],F)),
                     Median=c(IFMedian(rs,kap[1],F),IFMedian(rs,kap[2],F)))

IFDFFishm<-melt(IFDFFish,id=c("r","kappa"))
colnames(IFDFFishm)[3:4]<-c("Estimator","IF")

qplot(r,IF,data=IFDFFishm,geom='line',colour=Estimator,ylab='IF(r,Cayley)',size=I(1.25))+theme_bw()+geom_vline(xintercept=0,colour='gray50')+
  geom_hline(yintercept=0,colour='gray50')+theme(aspect.ratio=1)+facet_grid(.~kappa, labeller = label_bquote(kappa==.(x)))+
  scale_x_continuous(breaks=c(-pi,-pi/2,0,pi/2,pi),labels=c(expression(-pi,-frac(pi,2),0,frac(pi,2),pi)))
#ggsave("C:/Users/Brittney Ritchey/Dropbox/Thesis/Defense/figure/CayleyIFComp.pdf",width=10,height=5,units="in")


