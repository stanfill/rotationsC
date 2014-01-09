######################################################
######################################################
#This code was originaly "CompareMeanMedianCoverage.R" and has been modified
######################################################
#compare regions for central orientation based on mean and median
#for contaminated distributions
library(gridExtra)
library(plyr)
library(reshape2)
library(ggplot2)
#library(rotations)
#source("intervals/MedianPaper/contaminationFunctions.R")
#sourceCpp("intervals/ZhangMethod.cpp")  

library(rotations,lib='../BayesCR')
source("contaminationFunctions.R")
sourceCpp("ZhangMethod.cpp") 

#ZhangMethod.cpp contains the functions that will compute c/d and perform the zhang bootstrap

alp<-.1
critVal<-qchisq(1-alp,3)
n<-c(10,50,100)

#These kappa values give a circular variance of 1/7 (arbitrary)
#kappa1<-20         #For Cayley kappa=20
kappa1<-20    #For Fisher kappa=5.393154      
kappa2<-kappa1      #Contamination distribution kappa


eps<-c(0,.1,.2) #Amount of contamination
B<-5000  			#Number of samples to use to estimate CDF
Distfn<-c(rcayley,rfisher)
Dist<-c("Fisher")

simSize<-length(eps)*length(n)*length(kappa1)*length(Dist)

CRcompare<-data.frame(expand.grid(eps=eps,kappa=kappa1,n=n,Dist=Dist),MeanDist=0,MeanNTH=0,MeanBoot=0,MedianDist=0,MedianNTH=0,MedianBoot=0)
CRcompare<-CRcompare[rep(seq_len(nrow(CRcompare)),each=B),]
dimS<-nrow(CRcompare)

Scont<-genR(pi/2)

for(j in 1:dimS){
  
  nj<-CRcompare$n[j]
  #rs<-rfisher(nj,kappa=CRcompare$kappa[j])
  #Rs<-genRCont(rs,p=CRcompare$eps[j],S=id.SO3,Scont=Scont)
  Rs<-ruarsCont(nj,rfisher,kappa1=CRcompare$kappa[j],p=CRcompare$eps[j],S=id.SO3,Scont=Scont)
  Qs<-newQ4(Rs)
  
  #Mean
  Qhat<-mean(Qs)
  #Record the bias of the mean
  CRcompare$MeanDist[j]<-angle(Qhat)
  
  cdHat<-cdfunsC(Qs,Qhat) #compute c and d hat using consistent estimators
  chat<-cdHat[1];  dhat<-cdHat[2]
  
  #compute the volume of the mean-based normal theory CR and save it under 'MeanNTH'
  CRcompare$MeanNTH[j]<-sqrt(chat*critVal/(2*nj*dhat^2))
  
  #compute the volume of the mean-based bootstrap theory CR and save it under 'MeanNTH'
  critValBoot<-as.numeric(quantile(zhangQ4(Qs,300),1-alp,na.rm=T))
  CRcompare$MeanBoot[j]<-sqrt(chat*critValBoot/(2*nj*dhat^2))  
  
  #Median
  Stilde<-median(Rs)
  #Record the bias of the median
  CRcompare$MedianDist[j]<-angle(Stilde)
  
  cdTilde<-cdfunsCSO3(Rs,Stilde) #compute c and d tilde using consistent estimators
  ctilde<-cdTilde[1]
  dtilde<-cdTilde[2]
  
  
  #compute the volume of the mean-based normal theory CR and save it under 'MeanNTH'
  CRcompare$MedianNTH[j]<-sqrt(ctilde*critVal/(2*nj*dtilde^2))
  
  #compute the volume of the mean-based bootstrap theory CR and save it under 'MeanNTH'
  critValBootMed<-as.numeric(quantile(zhangMedianC(Rs,300),1-alp,na.rm=T))
  CRcompare$MedianBoot[j]<-sqrt(ctilde*critValBootMed/(2*nj*dtilde^2))  
 
  if(j%%1000==0){
    write.csv(CRcompare,"Results/MeanMedianContCompFisherPi_4.csv")
  }
  
}

write.csv(CRcompare,"Results/MeanMedianContCompFisherPi_4.csv")
CRcompare<-read.csv("intervals/MedianPaper/Results/MeanMedianContCompCayleyPi_4.csv")[,-1]

#CRFisher<-read.csv("intervals/MedianPaper/Results/MeanMedianContCompFisherPart1.csv")[,-1]
#CRFisher2<-read.csv("intervals/MedianPaper/Results/MeanMedianContCompFisherPart2.csv")[,-1]
#CRcompare<-rbind(CRFisher,CRFisher2)

#Remove rows that haven't finished running yet
#CRcompare<-CRcompare[rowSums(CRcompare[,5:10])>0,]

CRcompSum<-ddply(CRcompare,.(eps,kappa,n,Dist),summarize,Mean=mean(MeanDist),Median=mean(MedianDist),
                 MeanNTH=mean((MeanNTH)),MedianNTH=mean((MedianNTH)),
                 MeanBoot=mean((MeanBoot)),MedianBoot=mean((MedianBoot)))

#Compare Estimator Bias
CRcompM<-melt(CRcompSum,id=c("eps","kappa","n","Dist"))
labs<-c("Mean","Median", "Mean (NTH)","Mean (Boot)","Median (NTH)","Median (Boot)")
CRcompM$variable<-factor(CRcompM$variable,levels=c("Mean","Median","MeanNTH","MeanBoot","MedianNTH","MedianBoot"),labels=labs)
labsN<-c(expression(n==10),expression(n==50), expression(n==100))
CRcompM$n<-factor(CRcompM$n,levels=c("10","50","100"),labels=labsN)


BiasDF<-CRcompM[CRcompM$variable%in%c("Mean","Median"),]
colnames(BiasDF)[5]<-"Estimator"
qplot(eps,value,data=BiasDF,geom='line',colour=Estimator,group=Estimator,size=I(1.25),xlab=expression(epsilon),ylab="Bias")+
  facet_grid(n~.,scales="free_y",labeller = label_parsed)+theme_bw()
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/BiasComp.pdf",width=7,height=5)

#Compare Confidence Region Volume (?)
VolumeDF<-CRcompM[!CRcompM$variable%in%c("Mean","Median"),]
VolumeDF$Boot<-0; VolumeDF[grep("Boot",VolumeDF$variable),]$Boot<-1
qplot(eps,value,data=VolumeDF,geom='line',size=I(1.25),colour=variable,group=variable,xlab=expression(epsilon),ylab="CR Size")+
  facet_grid(n~.,scales="free_y",labeller = label_parsed)+theme_bw()


#Median Boot with n=10 is bad, remove it
VolumeDFEdited<-VolumeDF[VolumeDF$value<pi,]
qplot(eps,value,data=VolumeDFEdited,geom='line',size=I(1.25),colour=variable,group=variable,xlab=expression(epsilon),ylab="Region Size")+
  facet_grid(n~.,scales="free_y",labeller = label_parsed)+theme_bw()+theme(aspect.ratio=1/2)+
  scale_x_continuous(breaks=c(0,.1,.2))+theme(legend.position="none")
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/VolumeCompFisher.pdf",width=4,height=6,units="in")
#ggsave("C:/Users/Brittney Ritchey/Dropbox/Thesis/Intervals - Median/Figures/VolumeCompFisherPi4.pdf",width=4,height=6,units="in")

#Same thing, but BW
qplot(eps,value,data=VolumeDFEdited,geom='line',size=I(1.25),linetype=variable,col=variable,group=variable,xlab=expression(epsilon),ylab="Region Size")+
  facet_grid(n~.,scales="free_y",labeller = label_parsed)+theme_bw()+theme(legend.position="none")+theme(aspect.ratio=1/2)+
  scale_x_continuous(breaks=c(0,.1,.2))+
  scale_linetype_manual(values = rep(rev(c("dashed","solid")),2))+
  scale_colour_manual(values = rep(rev(c("gray50","black")),each=2))
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/VolumeCompBWFisher.pdf",width=4,height=6,units="in")
#ggsave("C:/Users/Brittney Ritchey/Dropbox/Thesis/Intervals - Median/Figures/VolumeCompBWFisher.pdf",width=4,height=6,units="in")


#Compare region coverage rates
CRcoverage<-ddply(CRcompare,.(eps,kappa,n,Dist),summarize,MeanNTH=sum(MeanDist<MeanNTH)/length(eps),MedianNTH=sum(MedianDist<MedianNTH)/length(eps),
                  MeanBoot=sum(MeanDist<MeanBoot)/length(eps),MedianBoot=sum(MedianDist<MedianBoot)/length(eps))
CRcoverM<-melt(CRcoverage,id=c("eps","kappa","n","Dist"))
labs<-c("Mean","Median", "Mean (NTH)","Mean (Boot)","Median (NTH)","Median (Boot)")
CRcoverM$variable<-factor(CRcoverM$variable,levels=c("Mean","Median","MeanNTH","MeanBoot","MedianNTH","MedianBoot"),labels=labs)
labsN<-c(expression(n==10),expression(n==50), expression(n==100))
CRcoverM$n<-factor(CRcoverM$n,levels=c("10","50","100"),labels=labsN)


p1<-qplot(eps,100*value,data=CRcoverM,geom='line',colour=variable,group=variable,size=I(1.25),xlab=expression(epsilon),ylab="Coverage (%)")+
  facet_grid(n~.,scales='free_y',labeller = label_parsed)+geom_hline(yintercept=90)+theme_bw()+labs(colour="")+theme(legend.position='top')
qplot(eps,100*value,data=CRcoverM,geom='line',colour=variable,group=variable,size=I(1.25),xlab=expression(epsilon),ylab="Coverage (%)",ylim=c(0,100))+
  facet_grid(n~.,labeller = label_parsed)+geom_hline(yintercept=90)+theme_bw()+coord_fixed(.2/200)+theme(legend.position='none')+
  scale_x_continuous(breaks=c(0,.1,.2))
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/CoverageCompFisher.pdf",width=4,height=6,units="in")
3ggsave("C:/Users/Brittney Ritchey/Dropbox/Thesis/Intervals - Median/Figures/CoverageCompCayleyPi4.pdf",width=4,height=6,units="in")

p2<-qplot(eps,100*value,data=CRcoverM,geom='line',linetype=variable,colour=variable,group=variable,size=I(1.25),xlab=expression(epsilon),ylab="Coverage (%)")+
  facet_grid(n~.,labeller = label_parsed)+geom_hline(yintercept=90)+theme_bw()+coord_fixed(.2/200)+
  labs(colour="",linetype="")+  
  theme(legend.key.width=unit(2,"line"),legend.position='top')+
  scale_linetype_manual(values = rep(rev(c("dashed","solid")),2))+
  scale_colour_manual(values = rep(rev(c("gray50","black")),each=2))

qplot(eps,100*value,data=CRcoverM,geom='line',linetype=variable,colour=variable,group=variable,size=I(1.25),xlab=expression(epsilon),ylab="Coverage (%)")+
  facet_grid(n~.,labeller = label_parsed)+geom_hline(yintercept=90)+theme_bw()+coord_fixed(.2/200)+theme(legend.position='none')+
  scale_x_continuous(breaks=c(0,.1,.2))+
  scale_linetype_manual(values = rep(rev(c("dashed","solid")),2))+
  scale_colour_manual(values = rep(rev(c("gray50","black")),each=2))
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/CoverageCompBWFisher.pdf",width=4,height=6,units="in")
#ggsave("C:/Users/Brittney Ritchey/Dropbox/Thesis/Intervals - Median/Figures/CoverageCompBWFisher.pdf",width=4,height=6,units="in")


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend<-g_legend(p1)
#Need to use "Export" because it isn't a ggplot2 object
grid.newpage()
grid.draw(legend)

legend<-g_legend(p2)
#Need to use "Export" because it isn't a ggplot2 object
grid.newpage()
grid.draw(legend)

#Small sample Relative efficiency
#CRcompare$RE<-(CRcompare$MeanNTH^2)/(CRcompare$MedianNTH^2)
#RelativeE<-ddply(CRcompare,.(eps,kappa,n,Dist),summarize,SSRE=mean(RE))
#RelativeE$n<-as.factor(RelativeE$n)
#qplot(eps,SSRE,data=RelativeE,colour=n,group=n,geom='line')


#####
#Create Tables
#Coverage rate comparison
# library(xtable)
# CRcompare<-read.csv("intervals/MedianPaper/Results/MeanMedianContComp.csv")[,-1]
# CRFisher<-read.csv("intervals/MedianPaper/Results/MeanMedianContCompFisherPart1.csv")[,-1]
# CRFisher2<-read.csv("intervals/MedianPaper/Results/MeanMedianContCompFisherPart2.csv")[,-1]
# CRcompare<-rbind(CRcompare,CRFisher,CRFisher2)
# 
# ###
# #region coverage rate
# 
# CRcoverage<-ddply(CRcompare,.(eps,kappa,n,Dist),summarize,MeanNTH=sum(MeanDist<MeanNTH)/length(eps),MedianNTH=sum(MedianDist<MedianNTH)/length(eps),
#                   MeanBoot=sum(MeanDist<MeanBoot)/length(eps),MedianBoot=sum(MedianDist<MedianBoot)/length(eps))
# CRcoverM<-melt(CRcoverage,id=c("eps","kappa","n","Dist"))
# labs<-c("Mean","Median", "Mean (NTH)","Median (NTH)","Mean (Boot)","Median (Boot)")
# CRcoverM$variable<-factor(CRcoverM$variable,levels=c("Mean","Median","MeanNTH","MedianNTH","MeanBoot","MedianBoot"),labels=labs)
# CRcoverM$eps<-as.factor(CRcoverM$eps)
# 
# xtable(dcast(CRcoverM,Dist+n+eps~variable),caption="Coverage rate comparison.",digits=3)
# 
# ###
# #region size
# 
# CRcompSum<-ddply(CRcompare,.(eps,kappa,n,Dist),summarize,Mean=mean(MeanDist),Median=mean(MedianDist),
#                  MeanNTH=mean((MeanNTH)),MedianNTH=(mean(MedianNTH)),
#                  MeanBoot=mean((MeanBoot)),MedianBoot=(mean(MedianBoot)))
# 
# CRcompM<-melt(CRcompSum,id=c("eps","kappa","n","Dist"))
# CRcompM$variable<-factor(CRcompM$variable,levels=c("Mean","Median","MeanNTH","MedianNTH","MeanBoot","MedianBoot"),labels=labs)
# VolumeDF<-CRcompM[!CRcompM$variable%in%c("Mean","Median"),]
# VolumeDF$eps<-as.factor(VolumeDF$eps)
# 
# xtable(dcast(VolumeDF,Dist+n+eps~variable),caption="Radius comparison.",digits=3)


