######################################################
######################################################
#This code was originaly "CompareMeanMedianCoverage.R" and has been modified
######################################################
#compare regions for central orientation based on mean and median
#for contaminated distributions

library(plyr)
library(reshape2)
#library(rotations)
#source("intervals/MedianPaper/contaminationFunctions.R")
#sourceCpp("intervals/ZhangMethod.cpp")  

library(rotations,lib='../BayesCR')
source("contaminationFunctions.R")
sourceCpp("ZhangMethod.cpp") 

#ZhangMethod.cpp contains the functions that will compute c/d and perform the zhang bootstrap

alp<-.1
critVal<-qchisq(1-alp,3)
n<-c(25,50,100)
kappa1<-20      
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
  Rs<-ruarsCont(nj,rfisher,kappa1=CRcompare$kappa[j],p=CRcompare$eps[j],S=id.SO3,Scont=Scont)
  Qs<-Q4(Rs)
  
  #Mean
  Qhat<-mean(Qs)
  #Record the bias of the mean
  CRcompare$MeanDist[j]<-angle(Qhat)
  
  cdHat<-cdfunsC(Qs,Qhat) #compute c and d hat using consistent estimators
  chat<-cdHat[1];  dhat<-cdHat[2]
  hsqMean<-CRcompare$MeanDist[j]^2
  
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
    write.csv(CRcompare,"Results/MeanMedianContCompFisher.csv")
  }
  
}

write.csv(CRcompare,"Results/MeanMedianContCompFisher.csv")
#CRcompare<-read.csv("Results/MeanMedianContComp.csv")[,-1]

#Remove rows that haven't finished running yet
CRcompare<-CRcompare[rowSums(CRcompare[,5:10])>0,]

CRcompSum<-ddply(CRcompare,.(eps,kappa,n,Dist),summarize,Mean=mean(MeanDist),Median=mean(MedianDist),
                 MeanNTH=mean(MeanNTH),MedianNTH=min(mean(MedianNTH),pi),
                 MeanBoot=mean(MeanBoot),MedianBoot=min(mean(MedianBoot),pi))

#Compare Estimator Bias
CRcompM<-melt(CRcompSum,id=c("eps","kappa","n","Dist"))
BiasDF<-CRcompM[CRcompM$variable%in%c("Mean","Median"),]
colnames(BiasDF)[5]<-"Estimator"
qplot(eps,value,data=BiasDF,geom='line',colour=Estimator,group=Estimator,size=I(1.25),xlab=expression(epsilon),ylab="Bias")+
  facet_grid(n~.,scales="free_y")+theme_bw()
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/BiasComp.pdf",width=8,height=5)

#Compare Confidence Region Volume (?)
VolumeDF<-CRcompM[!CRcompM$variable%in%c("Mean","Median"),]
VolumeDF$Boot<-0; VolumeDF[grep("Boot",VolumeDF$variable),]$Boot<-1
qplot(eps,value,data=VolumeDF,geom='line',size=I(1.25),colour=variable,group=variable,xlab=expression(epsilon),ylab="CR Size")+
  facet_grid(n~.,scales="free_y")+theme_bw()


#Median Boot with n=10 is bad, remove it
VolumeDFEdited<-VolumeDF[VolumeDF$value<pi,]
qplot(eps,value,data=VolumeDFEdited,geom='line',size=I(1.25),colour=variable,group=variable,xlab=expression(epsilon),ylab="CR Size")+
  facet_grid(n~.,scales="free_y")+theme_bw()+labs(colour="")
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/VolumeComp.pdf",width=8,height=5)


#Compare region coverage rates
CRcoverage<-ddply(CRcompare,.(eps,kappa,n,Dist),summarize,MeanNTH=sum(MeanDist<MeanNTH)/length(eps),MedianNTH=sum(MedianDist<MedianNTH)/length(eps),
                  MeanBoot=sum(MeanDist<MeanBoot)/length(eps),MedianBoot=sum(MedianDist<MedianBoot)/length(eps))
CRcoverM<-melt(CRcoverage,id=c("eps","kappa","n","Dist"))
CRcoverM<-CRcoverM[CRcoverM$value>0,]#Remove rows that haven't finished running yet

qplot(eps,100*value,data=CRcoverM,geom='line',colour=variable,group=variable,size=I(1.25),xlab=expression(epsilon),ylab="Coverage (%)")+
  facet_grid(n~.,scales='free_y')+geom_hline(yintercept=90)+theme_bw()+labs(colour="")
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Median/Figures/CoverageComp.pdf",width=8,height=5)

#coverCompare<-read.csv("Results/MeanMedianComparison.csv")[,-1]

# ccrit<-qchisq(1-alp,3)
# coverCompare$MeanCoverC<-as.numeric(coverCompare$MeanStat<ccrit)
# coverCompare$MeanCoverZ<-as.numeric(coverCompare$MeanStat<coverCompare$MeanCrit)
# 
# coverCompare$MedianCoverC<-as.numeric(coverCompare$MedianStat<ccrit)
# coverCompare$MedianCoverZ<-as.numeric(coverCompare$MedianStat<coverCompare$MedianCrit)
# 
# compareRate<-ddply(coverCompare,.(Dist,nus,n),summarize,MeanCCover=100*sum(MeanCoverC)/B,MeanZCover=100*sum(MeanCoverZ)/B,
# 									 MedianCCover=100*sum(MedianCoverC)/B,MedianZCover=100*sum(MedianCoverZ)/B)
# 
# #####
# #Analyze the results
# 
# cRateM<-melt(compareRate,id=c('Dist','nus','n'))
# colnames(cRateM)[4]<-'Method'
#  
# levels(cRateM$Method)<-c("Mean(NTH)","Mean(Boot.)","Median(NTH)","Median(Boot.)")
#  
# levels(cRateM$Dist)<-c("Cayley","matrix~~Fisher")
# cRateM$nu<-factor(cRateM$nu,labels=c("nu == 0.25","nu == 0.5","nu == 0.75"))
# 
# require(grid) #necessary to use the unit() function
# 
# # Compare coverage rates by estimator and method
# qplot(n,value,data=cRateM,colour=Method,group=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
#  	facet_grid(Dist~nu,labeller=label_parsed)+
#  	geom_hline(yintercept=(1-alp)*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
#  	scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+theme(panel.margin=unit(0.5,'lines'))
# #ggsave("/Users/stanfill/Dropbox/Thesis/Intervals/Figures/MeanMedianCoverCompare.pdf",width=8,height=6)
# 
# qplot(n,value,data=cRateM,col=Method,linetype=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
#   scale_linetype_manual(values = rep(c("dashed","solid"),each=2))+
#   scale_colour_manual(values = rep(c("gray50","black"),2))+
#   facet_grid(Dist~nu,labeller=label_parsed)+
#   geom_hline(yintercept=(1-alp)*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
#   scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+
#   theme(panel.margin=unit(0.5,'lines'),legend.key.width=unit(3,"line"))  
# #ggsave("/Users/stanfill/Dropbox/Thesis/Intervals/Figures/MeanMedianCoverCompareBW.pdf",width=8,height=6)
# 
# 
# #####################
# #Plot test statistics against one another?
# ccrit<-qchisq(1-alp,3)
# 
# p<-ggplot(coverCompare[coverCompare$nus==0.25,],aes(MeanStat,MedianStat))
# p+xlim(0,30)+ylim(0,30)+facet_grid(n~Dist)+geom_hline(yintercept=ccrit,colour='red')+
#   geom_vline(xintercept=ccrit,colour='red')+geom_point(alpha=I(.1))+xlab(expression(hat(S)[E]))+
#   ylab(expression(tilde(S)[E]))+coord_fixed()
# 
# #ggsave("/Users/stanfill/Dropbox/Thesis/Intervals/Figures/MeanMedianStatCompare.pdf",width=5,height=10)
# 
# 
# qplot(MeanStat,MedianStat,data=coverCompare[coverCompare$nus==0.5,],facets=n~Dist,
# 	xlim=c(0,30),ylim=c(0,30))+geom_hline(yintercept=ccrit,colour='red')+
# 	geom_vline(xintercept=ccrit,colour='red')
# 
# p<-ggplot(coverCompare[coverCompare$nus==0.75,],aes(MeanStat,MedianStat))
# p+xlim(0,30)+ylim(0,30)+facet_grid(n~Dist)+geom_hline(yintercept=ccrit,colour='red')+
#   geom_vline(xintercept=ccrit,colour='red')+geom_point(alpha=I(.1))+xlab(expression(hat(S)[E]))+
#   ylab(expression(tilde(S)[E]))
# 
# 
# 
# #####################
# #plot Mean-test statistics against bootstrap critical value
# qplot(MeanStat,MeanCrit,data=coverCompare[coverCompare$nus==0.25,],facets=n~Dist)
# qplot(MeanStat,MeanCrit,data=coverCompare[coverCompare$nus==0.5,],facets=n~Dist)
# qplot(MeanStat,MeanCrit,data=coverCompare[coverCompare$nus==0.75,],facets=n~Dist)
# 
# #####################
# #plot Median-test statistics against bootstrap critical value
# qplot(MedianStat,MedianCrit,data=coverCompare[coverCompare$nus==0.25,],facets=n~Dist)
# qplot(MedianStat,MedianCrit,data=coverCompare[coverCompare$nus==0.5,],facets=n~Dist)
# qplot(MedianStat,MedianCrit,data=coverCompare[coverCompare$nus==0.75,],facets=n~Dist)
# 
# ####################
# #plot mean versus median critical value
# qplot(MeanCrit,MedianCrit,data=coverCompare[coverCompare$nus==0.25,],facets=n~Dist,ylim=c(0,100))
# 
# ####################
# #empirical asymptotic relative efficiency. 
# #Nope, the test statistic includes d(S,S-check)^2 which is going to zero making the ratio 1
# 
# cCover<-ddply(coverCompare,.(nus,n,Dist),summarize,ratio=mean(MedianCrit/MeanCrit))
# 
# qplot(n,ratio,data=cCover,facets=Dist~nus,geom='line',ylim=c(0,4))
