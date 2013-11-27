IFMean<-function(r,kap,Fisher=T){
  
  if(Fisher){
    d<-(((kap+1)/kap)*besselI(2*kap,1)-besselI(2*kap,0))/(3*(besselI(2*kap,0)-besselI(2*kap,1)))
  }else{
    d<-kap/(kap+2)
  }
  return(sin(r)/d)
  
}

IFMedian<-function(r,kap,Fisher=T){
  fn<-sin(r)/sqrt(1-cos(r))
  
  if(Fisher){
    d<-exp(2*kap)*(6*sqrt(kap)-(3+8*kap)*dawson(2*sqrt(kap)))/(24*sqrt(2)*pi*kap^1.5*(besselI(2*kap,0)-besselI(2*kap,1)))
  }else{
    d<-sqrt(2/pi)*kap*gamma(kap+2)/(3*gamma(kap+2.5))
  }
  
  return(fn/(2*d))
}


#############
#Compare IFs
library(gsl)
rs<-seq(-pi,pi,length=100)
kap<-100
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

IFDF<-data.frame(r=rep(rs,3),kappa=rep(kap,each=100),Mean=c(IFMean(rs,kap[1],F),IFMean(rs,kap[2],F),IFMean(rs,kap[3],F)),
                     Median=c(IFMedian(rs,kap[1],F),IFMedian(rs,kap[2],F),IFMedian(rs,kap[3],F)))

IFDF$kappa<-as.factor(IFDF$kappa)
qplot(r,Mean,data=IFDF,geom='line',colour=kappa,group=kappa,ylab="IF-Mean")+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')

qplot(r,Median,data=IFDF,geom='line',colour=kappa,group=kappa,ylab="IF-Median")+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')

#Put all on one plot
IFDFm<-melt(IFDF,id=c("r","kappa"))
colnames(IFDFm)[3:4]<-c("Estimator","IF")
qplot(r,IF,data=IFDFm,geom='line',colour=kappa,linetype=Estimator)+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')

qplot(r,IF,data=IFDFm,geom='line',colour=kappa,facets=.~Estimator)+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')

qplot(r,IF,data=IFDFm,geom='line',colour=Estimator)+facet_grid(.~kappa, labeller = label_bquote(kappa==.(x)))+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')+theme(aspect.ratio=1)

### Same with Fisher
kap<-c(1,2.5,5)
IFDFFish<-data.frame(r=rep(rs,3),kappa=rep(kap,each=100),Mean=c(IFMean(rs,kap[1],T),IFMean(rs,kap[2],T),IFMean(rs,kap[3],T)),
                 Median=c(IFMedian(rs,kap[1],T),IFMedian(rs,kap[2],T),IFMedian(rs,kap[3],T)))
IFDFFish$kappa<-as.factor(IFDFFish$kappa)
IFDFFishm<-melt(IFDFFish,id=c("r","kappa"))
colnames(IFDFFishm)[3:4]<-c("Estimator","IF")
qplot(r,IF,data=IFDFFishm,geom='line',colour=Estimator)+facet_grid(.~kappa, labeller = label_bquote(kappa==.(x)))+theme_bw()+
  geom_vline(xintercept=0,colour='gray50')+geom_hline(yintercept=0,colour='gray50')+theme(aspect.ratio=1)
