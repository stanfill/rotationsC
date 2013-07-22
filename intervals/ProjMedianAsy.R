#Note about these simulations:
#
#The C++ version of genR() seems to crash for the fisher or von Mises distributions but works
#with relative reliability for the Cayley distribution.  For now I've removed any call to the
#C++ version of genR(), namely "SO3defaultC()" from within "SO3.default()", which I thought I had
#already done.  Not it runs with little issue!

################################
#In this section I compare the empirical CDF of the test statistic to 
#the theoretical chi^2_3 limiting distribution

library(rotations2)
setwd("C:/Users/stanfill/Desktop/GitHub/rotationsC/intervals")
source("IntervalFuns.R")
B<-1000
n<-c(10,50,300)
numN<-length(n)
kap<-.1
tstats<-matrix(0,B,numN)

for(j in 1:numN){
	
	for(i in 1:B){
		rs<-rcayley(n[j],kappa=kap)
		Rs<-genR(rs)
		

		cosrs<-cos(rs)
		crs<-(cosrs+1)
		drs<-(1+3*cosrs)/(sqrt(1-cosrs))
		#cosrs2<-cos(rs/2)^2
		#cotrs<-cot(rs)

		c<-mean(crs)/6  #I think this is C for proj median according to notes from 7/16
		d<-mean(drs)/12

		
		Shat<-median(Rs)
	
		ShatMedian<-dist(Shat,method='intrinsic',p=2)
	
		tstats[i,j]<-2*n[j]*(d^2)*ShatMedian/c
	
	}
	
}
#hist(tstats,breaks=100,prob=T)

#ses<-seq(0,max(tstats),length=B)
#lines(ses,dchisq(ses,3))

xmax<-15

for(j in 1:numN){
	tstats[,j]<-sort(tstats[,j])
}
plot(tstats[,numN],ecdf(tstats[,numN]),type='l',xlim=c(0,xmax))

for(j in 1:(numN-1)){
	lines(tstats[,j],ecdf(tstats[,j]),col=(j+1))
}

seqChi<-seq(0,xmax,length=B)

lines(seqChi,pchisq(seqChi,3),lty=2)


###########################
#Make plots for proj.median similar to proj.mean

library(plyr)
library(reshape2)
library(rotations2)
setwd("C:/Users/stanfill/Desktop/GitHub/rotationsC/intervals")
source("IntervalFuns.R")

n<-c(10,50,100,300)
nus<-c(.25,.75)
B<-1000				#Number of samples to use to estimate CDF
Dist<-c('cayley','fisher','mises')

simSize<-length(n)*length(nus)*length(Dist)

tMat<-matrix(0,simSize,B)

cdfDF<-data.frame(expand.grid(nus=nus,n=n,Dist=Dist),tMat)

for(j in 1:simSize){
	
	if(cdfDF$Dist[j]=='cayley'){
		
		rangle<-rcayley
		kappaj <- cayley_kappa(cdfDF$nu[j])
		
	}else if(cdfDF$Dist[j]=='fisher'){
		
		rangle<-rfisher
		kappaj <- fisher_kappa(cdfDF$nu[j])
		
	}else{
		rangle<-rvmises
		kappaj <- vmises_kappa(cdfDF$nu[j])
	}
	
	for(i in 1:B){
		
		rs<-rangle(cdfDF$n[j],kappa=kappaj)
		
		Rs<-genR(rs)
		
		cosrs<-cos(rs)
		crs<-(cosrs+1)
		drs<-(1+3*cosrs)/(sqrt(1-cosrs))
		#cosrs2<-cos(rs/2)^2
		#cotrs<-cot(rs)
		
		c<-mean(crs)/6  #I think this is C for proj median according to notes from 7/16
		d<-mean(drs)/12
		
		Shat<-median(Rs)
		
		hsqMean<-dist(Shat,method='intrinsic',p=2)
		
		cdfDF[j,(3+i)]<-2*cdfDF$n[j]*d^2*hsqMean/c
		
	}
}

resM<-melt(cdfDF,id=c("Dist","nus","n"))
ss<-seq(0,15,length=B)
Probs<-pchisq(ss,3)

resM$n<-as.factor(resM$n)
resM$Prob<-0
resM$ID<-paste(resM$Dist,resM$nu,resM$n)
kns<-unique(resM$ID)

for(i in 1:length(unique(resM$ID))){
	resM[resM$ID==kns[i],]$value<-sort(resM[resM$ID==kns[i],]$value)
	resM[resM$ID==kns[i],]$Prob<-ecdf(resM[resM$ID==kns[i],]$value)
}

chiDF<-data.frame(Dist="All",nus=rep(nus,B),n='Chisq',variable='Tr',value=rep(ss,each=2),Prob=rep(pchisq(ss,3),each=2))
chiDF$ID<-paste(chiDF$Dist,chiDF$nus,chiDF$n)

fullDF<-rbind(resM,chiDF)

Newlabs<-c("Chisq","10","50","100","300")
fullDF$n<-factor(fullDF$n,levels=Newlabs)
fullDF$Stat<-1
fullDF[fullDF$n=='Chisq',]$Stat<-2
fullDF$Stat<-as.factor(fullDF$Stat)
fullDF$nus<-factor(fullDF$nus,labels=c("nu == 0.25","nu == 0.75"))
	

qplot(value,Prob,data=fullDF[fullDF$Dist%in%c("cayley",'All'),],colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,15))+
	scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
	facet_grid(.~nus,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
	scale_size_discrete(range=c(0.75,1.5),guide='none')+
	guides(colour=guide_legend(label.hjust=0))
	
#setwd("C:/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
#ggsave("CayleyECDFMedian.pdf",height=5,width=8)
#write.csv(fullDF,"CayleyECDF.csv")
	
	
qplot(value,Prob,data=fullDF[fullDF$Dist%in%c("fisher",'All'),],colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,15))+
	scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
	facet_grid(.~nus,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
	scale_size_discrete(range=c(0.75,1.5),guide='none')+
	guides(colour=guide_legend(label.hjust=0))
	
#setwd("C:/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
#ggsave("FisherECDFMedian.pdf",height=5,width=8)
#write.csv(fullDF,"FisherECDF.csv")
	

qplot(value,Prob,data=fullDF[fullDF$Dist%in%c("mises",'All'),],colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,15))+
	scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
	facet_grid(.~nus,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
	guides(colour=guide_legend(label.hjust=0))+
	scale_size_discrete("",range=c(0.75,1.5),guide='none')
	
#setwd("C:/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
#ggsave("vonMisesECDF.pdf",height=5,width=8)
#write.csv(fullDF,"vonMisesECDF.csv")
	

###########################
###########################
#How do c and d compare between proj.mean and proj.median

library(rotations2)
n<-100
kap<-1
B<-100
AvarHat<-rep(0,B)
AvarTilde<-rep(0,B)

for(i in 1:B){

	rs<-rvmises(n,kap)
	cosrs<-cos(rs)
	cos2rs<-cos(rs)^2

	crs<-(cosrs+1)
	drs<-(1+3*cosrs)/(sqrt(1-cosrs))

	chat<-2*mean(1-cos2rs)/3
	dhat<-mean(1+2*cosrs)/3
	AvarHat[i]<-chat/(2*dhat^2)

	ctilde<-mean(crs)/6  
	dtilde<-mean(drs)/12
	AvarTilde[i]<-ctilde/(2*dtilde^2)

}

#Empirically AvarTilde > AvarHat
plot(AvarHat,AvarTilde,pch=19)
abline(0,1)


###########################
###########################
#Compare theoretical to empirical c and d values for distributions

library(rotations2)
kap<-250

#Cayley distribution
rs<-rcayley(1000,kappa=kap)
crs<-cos(rs)

mean(1/sqrt(1-crs))
sqrt(2)*gamma(kap+2)/((2*kap+1)*gamma(kap+.5)*gamma(1.5)) #Pretty good

mean(crs/sqrt(1-crs))
(2*kap-1)*gamma(kap+2)/(sqrt(2*pi)*gamma(kap+2.5))  #Pretty good


#Fisher distribution
library(gsl)
kap<-100

rs<-rfisher(1000,kappa=kap)
crs<-cos(rs)
mean(1/sqrt(1-crs))

exp(2*kap)*sqrt(2)*dawson(2*sqrt(kap))/(pi*sqrt(kap)*(besselI(2*kap,0)-besselI(2*kap,1)))

mean(crs/sqrt(1-crs))
exp(2*kap)*(2*sqrt(kap)-(1+4*kap)*dawson(2*sqrt(kap)))/((2*kap)^(1.5)*pi*(besselI(2*kap,0)-besselI(2*kap,1)))  #Pretty good


######################################################
######################################################
#Find coverage rates for the parametric and bootstrap
#confidence regions


library(plyr)
library(reshape2)
library(rotations2)
sourceCpp("intervals/ZhangMethod.cpp")  

#ZhangMethod.cpp contains the functions that will compute c/d and perform the zhang bootstrap

alp<-.1
critVal<-qchisq(1-alp,3)
n<-c(10,20,50,100)
nus<-c(.25,.5,.75)
B<-10000  			#Number of samples to use to estimate CDF
Dist<-c('Cayley','matrix-Fisher')

simSize<-length(n)*length(nus)*length(Dist)

tMat<-matrix(0,simSize,B)

cdfDF<-data.frame(expand.grid(nus=nus,n=n,Dist=Dist),tMat)
coverRate<-data.frame(expand.grid(nus=nus,n=n,Dist=Dist),Chang=0,Zhang=0)
#coverRate<-read.csv("Results/MedianResultsM300.csv")[,-1]

for(j in 1:simSize){
  
  if(cdfDF$Dist[j]=='Cayley'){
    
    rangle<-rcayley
    kappaj <- cayley_kappa(cdfDF$nu[j])
    
  }else{
    
    rangle<-rfisher
    kappaj <- fisher_kappa(cdfDF$nu[j])
    
  }
  
  for(i in 1:B){
    
    rs<-rangle(cdfDF$n[j],kappa=kappaj)
    
    Rs<-genR(rs)
    
    Shat<-median(Rs)
    
    cdTilde<-cdfunsCSO3(Rs,Shat) #compute c and d tilde using consistent estimators
    c<-cdTilde[1]
    d<-cdTilde[2]
    
	  if(d>1000){
	  	break
	  }
    
    hsqMean<-dist(Shat,method='intrinsic',p=2)
    
    statIJ<-2*cdfDF$n[j]*d^2*hsqMean/c
    
    #dfDF[j,(3+i)]<-statIJ
    coverRate[j,]$Chang<-coverRate[j,]$Chang+as.numeric(statIJ<critVal)
    
    zhangIJ<-as.numeric(quantile(zhangMedianC(Rs,300),1-alp,na.rm=T))
    
    #cdfDF[j,(3+i)]<-statIJ
    coverRate[j,]$Zhang<-coverRate[j,]$Zhang+as.numeric(statIJ<zhangIJ)
    
  }
  coverRate$Chang[j]<-100*coverRate$Chang[j]/B
  coverRate$Zhang[j]<-100*coverRate$Zhang[j]/B
  write.csv(coverRate,"Results/MedianResultsM300.csv")
}

cRateM<-melt(coverRate,id=c('Dist','nus','n'))
colnames(cRateM)[4]<-'Method'

levels(cRateM$Method)<-c("NTH(C&R)","B(Z&N)")

levels(cRateM$Dist)<-c("Cayley","matrix~~Fisher")
cRateM$nu<-factor(cRateM$nu,labels=c("nu == 0.25","nu == 0.5","nu == 0.75"))


qplot(n,value,data=cRateM,colour=Method,group=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
	facet_grid(Dist~nu,labeller=label_parsed)+
	geom_hline(yintercept=(1-alp)*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
	scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+theme(panel.margin=unit(0.5,'lines'))

ggsave("C:/Users/stanfill/Dropbox/Thesis/Intervals/Figures/CoverRatesB1000Median.pdf",width=7,height=4.5)
#ggsave("CoverRatesB1000Median.pdf",width=5,height=4)

#library(xtable)
#xtable(coverRate)