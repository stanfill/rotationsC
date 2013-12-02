#Note about these simulations:
#
#The C++ version of genR() seems to crash for the fisher or von Mises distributions but works
#with relative reliability for the Cayley distribution.  For now I've removed any call to the
#C++ version of genR(), namely "SO3defaultC()" from within "SO3.default()", which I thought I had
#already done.  Not it runs with little issue!

################################
#In this section I compare the empirical CDF of the test statistic to 
#the theoretical chi^2_3 limiting distribution

library(rotations)
setwd("/Users/stanfill/Documents/GitHub/rotationsC/intervals")
source("IntervalFuns.R")
B<-1000
n<-c(10,50,300)
numN<-length(n)
kap<-2
tstats<-matrix(0,B,numN)

S<-as.SO3(matrix(genR(pi/4),3,3))
#S<-id.SO3

for(j in 1:numN){
	
	for(i in 1:B){
		rs<-rcayley(n[j],kappa=kap)
		Rs<-genR(rs,S=S)
		
		Shat<-median(Rs)
    
    rs2<-dist(Rs,Shat,method='intrinsic')
		cosrs<-cos(rs2)
		crs<-(cosrs+1)
		drs<-(1+3*cosrs)/(sqrt(1-cosrs))
		#cosrs2<-cos(rs/2)^2
		#cotrs<-cot(rs)

		c<-mean(crs)/6  #I think this is C for proj median according to notes from 7/16
		d<-mean(drs)/12

		ShatMedian<-dist(Shat,S,method='intrinsic',p=2)
	
		tstats[i,j]<-2*n[j]*(d^2)*ShatMedian/c
	
	}
	
}
#hist(tstats,breaks=100,prob=T)

#ses<-seq(0,max(tstats),length=B)
#lines(ses,dchisq(ses,3))

xmax<-10

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
library(rotations)
source("intervals/IntervalFuns.R")	#This is needed for the ecdf function

n<-c(10,50,100,300)
kappa<-c(2,8)
B<-1000				#Number of samples to use to estimate CDF
#Dist<-c('cayley','fisher','mises')
Dist<-c('cayley','fisher')
simSize<-length(n)*length(kappa)*length(Dist)

tMat<-matrix(0,simSize,B)

cdfDF<-data.frame(expand.grid(kappa=kappa,n=n,Dist=Dist),tMat)

for(j in 1:simSize){
	
	if(cdfDF$Dist[j]=='cayley'){
		
		rangle<-rcayley
		#kappaj <- cayley_kappa(cdfDF$nu[j])
		
	}else if(cdfDF$Dist[j]=='fisher'){
		
		rangle<-rfisher
		#kappaj <- fisher_kappa(cdfDF$nu[j])
		
	}else{
		rangle<-rvmises
		#kappaj <- vmises_kappa(cdfDF$nu[j])
	}
	
	for(i in 1:B){
		
		rs<-rangle(cdfDF$n[j],kappa=cdfDF$kappa[j])
		
		Rs<-genR(rs)
		
		Shat<-median(Rs)
    rs2<-dist(Rs,Shat,method='intrinsic')
    
		cosrs<-cos(rs2)
		crs<-(cosrs+1)
		drs<-(1+3*cosrs)/(sqrt(1-cosrs))
		#cosrs2<-cos(rs/2)^2
		#cotrs<-cot(rs)
		
		c<-mean(crs)/6  #I think this is C for proj median according to notes from 7/16
		d<-mean(drs)/12
		
		hsqMean<-dist(Shat,method='intrinsic',p=2)
		
		cdfDF[j,(3+i)]<-2*cdfDF$n[j]*d^2*hsqMean/c
		
	}
}

resM<-melt(cdfDF,id=c("Dist","kappa","n"))
ss<-seq(0,10,length=B)
Probs<-pchisq(ss,3)

resM$n<-as.factor(resM$n)
resM$Prob<-0
resM$ID<-paste(resM$Dist,resM$kappa,resM$n)
kns<-unique(resM$ID)

for(i in 1:length(unique(resM$ID))){
	resM[resM$ID==kns[i],]$value<-sort(resM[resM$ID==kns[i],]$value)
	resM[resM$ID==kns[i],]$Prob<-ecdf(resM[resM$ID==kns[i],]$value)
}

chiDF<-data.frame(Dist="All",kappa=rep(kappa,B),n='Chisq',variable='Tr',value=rep(ss,each=2),Prob=rep(pchisq(ss,3),each=2))
chiDF$ID<-paste(chiDF$Dist,chiDF$nus,chiDF$n)

fullDF<-rbind(resM,chiDF)

Newlabs<-c("Chisq","10","50","100","300")
fullDF$n<-factor(fullDF$n,levels=Newlabs)
fullDF$Stat<-1
fullDF[fullDF$n=='Chisq',]$Stat<-2
fullDF$Stat<-as.factor(fullDF$Stat)
fullDF$kappa<-factor(fullDF$kappa,labels=c("kappa == 2","kappa == 8"))
	

qplot(value,Prob,data=fullDF[fullDF$Dist%in%c("cayley",'All'),],colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,10))+
	scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
	facet_grid(.~kappa,labeller=label_parsed)+theme_bw()+
	scale_size_discrete(range=c(0.75,1.5),guide='none')+
	guides(colour=guide_legend(label.hjust=0))+coord_equal(10)
	
#setwd("/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
#ggsave("CayleyECDFMedian.pdf",height=5,width=8)

	
	
qplot(value,Prob,data=fullDF[fullDF$Dist%in%c("fisher",'All'),],colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,10))+
	scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
	facet_grid(.~kappa,labeller=label_parsed)+theme_bw()+
	scale_size_discrete(range=c(0.75,1.5),guide='none')+
	guides(colour=guide_legend(label.hjust=0))+coord_equal(10)
	
#setwd("/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
#ggsave("FisherECDFMedian.pdf",height=5,width=7)

	

qplot(value,Prob,data=fullDF[fullDF$Dist%in%c("mises",'All'),],colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,15))+
	scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
	facet_grid(.~nus,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
	guides(colour=guide_legend(label.hjust=0))+
	scale_size_discrete("",range=c(0.75,1.5),guide='none')
	
#setwd("C:/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
#ggsave("vonMisesECDF.pdf",height=5,width=8)

#write.csv(fullDF,"medianECDF.csv")

###########################
###########################
#How do c and d compare between proj.mean and proj.median

library(rotations2)
n<-100
kap<-c(.01,1)
B<-100
AvarHat<-rep(0,B)
AvarHat<-matrix(c(AvarHat,AvarHat),ncol=2)
AvarTilde<-rep(0,B)
AvarTilde<-matrix(c(AvarTilde,AvarTilde),ncol=2)

for(j in 1:2){
  for(i in 1:B){

	  rs<-rcayley(n,kap)
	  cosrs<-cos(rs)
	  cos2rs<-cos(rs)^2

	  crs<-(cosrs+1)
	  drs<-(1+3*cosrs)/(sqrt(1-cosrs))

	  chat<-2*mean(1-cos2rs)/3
	  dhat<-mean(1+2*cosrs)/3
	  AvarHat[i,j]<-chat/(2*dhat^2)

	  ctilde<-mean(crs)/6  
	  dtilde<-mean(drs)/12
	  AvarTilde[i,j]<-ctilde/(2*dtilde^2)

  }
}
colMeans(AvarHat)/colMeans(AvarTilde)
#Empirically AvarTilde > AvarHat
plot(AvarHat[,1],AvarTilde[,1],pch=19)
abline(0,1)
plot(AvarHat[,2],AvarTilde[,2],pch=19)
abline(0,1)

##
#Based on c,d forms compare AV(Stilde)/AV(Shat)
#Cayley distribution-requires stirling's formula to simplify Gamma function

kap<-.1
hatOVERtilde<-(8*exp(1)*(kap+2)^2*(kap+1)^(2*kap+3))/(3*pi*(kap+3)*(kap+1.5)^(2*kap+4))
hatOVERtilde

cayRatio<-function(kap){
	return((8*exp(1)*(kap+2)^2*(kap+1)^(2*kap+3))/(3*pi*(kap+3)*(kap+1.5)^(2*kap+4)))
}

ks<-seq(0,25,length=100)
plot(ks,cayRatio(ks),type='l')

#Fisher
library(gsl)

kap<-.1

fishRatio<-function(kap){
	return((exp(4*kap)*(6*sqrt(kap)-(3+8*kap)*dawson(2*sqrt(kap)))^2)/(32*pi^2*kap^2*besselI(2*kap,1)*((kap+1)*besselI(2*kap,1)-kap*besselI(2*kap,0))))
}


ks<-seq(0,25,length=100)
plot(ks,fishRatio(ks),type='l')

##Make a pretty plot comparing ARE of the two estimators
cayRatio2<-function(kap){
  #compue ARE of median for Cayley distribution without Stirings approximation of the 
  #gamma function
  return(8*(kap+2)^2*gamma(kap+2)^2/(3*pi*(kap+3)*gamma(kap+2.5)^2))
}


ARE<-data.frame(kap=c(ks,ks),Are=c(cayRatio2(ks),fishRatio(ks)),Distribution=c(rep("Cayley",100),rep("matrix Fisher",100)))
qplot(kap,Are,data=ARE,linetype=Distribution,group=Distribution,lwd=I(1.5),geom='line',
	xlab=expression(kappa),ylab="Asymptotic Relative Efficiency")+
  theme_bw()+geom_hline(yintercept=8/(3*pi),colour='gray50')

#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals/Figures/AREDist.pdf",width=6,height=4)

###########################
#Plot ARE for cayley with and without using Stirling's to simplify gamma function
ks<-seq(0,25,length=50)
plot(ks,cayRatio2(ks),type='l') #without stirling approximation
lines(ks,cayRatio(ks),lty=2) #Using stirling approximation


###########################
#Use Fisher information matrix (FIM) from Leon et al. (2006) and
#c and d to plot ARE for cayley with and without using Stirling's to simplify gamma function

cayAre<-function(kap){
  c<-(4*kap+2)/(kap^2+5*kap+6)
  d<-kap/(kap+2)
  Mean<-c/(2*d^2)
  
  FIM<-kap^2/(2*kap-1)
  
  c2<-(2*kap+1)/(6*kap+12)
  d2<-sqrt(2/pi)*kap*gamma(kap+2)/(3*gamma(kap+2.5))
  Median<-c2/(2*d2^2)
  return(list(Hat=(1/FIM)/Mean,Tilde=(1/FIM)/Median))
}

ks<-seq(.51,25,length=50)
plot(ks,cayAre(ks)$Hat,type='l')
lines(ks,cayAre(ks)$Tilde,lty=2)
lines(ks,cayAre(ks)$Tilde/cayAre(ks)$Hat,lty=3)
legend("bottomright",c("MLE/MOM","MLE/Median","MOM/Median"),lty=c(1,2,3))
###########################
#See if this is true for Beta distribution.  It is.

#median has AV equal to 1/[2f(mu)]^2 where mu=alpha/(alpha+beta)

betaARE<-function(kappa){
  alpha<-kappa+0.5
  beta<-1.5
  mu<-alpha/(alpha+beta)
  
  Median<-(2*dbeta(mu,alpha,beta))^-2
  Mean<-(alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
  return(Mean/Median)
}

ks<-seq(.51,25,length=50)
plot(ks,betaARE(ks),type='l')

x<-seq(0,1,length=100)
plot(x,dbeta(x,.51,1.5),type='l')
lines(x,dbeta(x,.6,1.5),col=2)
lines(x,dbeta(x,.7,1.5),col=3)
lines(x,dbeta(x,1,1.5),col=4)

normARE<-function(sigma){
  Median<-(2*dnorm(0,0,sigma))^-2
  Mean<-sigma^2
  return(Mean/Median)
}

ks<-seq(.51,25,length=50)
plot(ks,normARE(ks),type='l')

###########################
#Idea, compare f(0) for cayley and fisher as a function of kappa
kap<-seq(.1,2,length=20)
fcay<-dcayley(0.0001,kappa=kap,Haar=T)
ffisher<-dfisher(0.0001,kappa=kap,Haar=T)

plot(kap,ffisher,type='l',ylim=c(min(fcay,ffisher),max(fcay,ffisher)))
lines(kap,fcay,col=2)

#Nope, try d-cayley versus d-fisher for median
dcay<-sqrt(2/pi)*kap*gamma(kap+2)/(3*gamma(kap+2.5))
dfish<-exp(2*kap)*(6*sqrt(kap)-(3+8*kap)*dawson(2*sqrt(kap)))/(24*sqrt(2)*pi*kap^1.5*(besselI(2*kap,0)-besselI(2*kap,1)))
plot(kap,dfish,type='l',ylim=c(min(dcay,dfish),max(dcay,dfish)))
lines(kap,dcay,col=2)

#Not convinced, try d-hat verus d-tilde for cayley
dhat<-kap/(kap+2)
plot(kap,dhat,type='l',ylim=c(min(dcay,dhat),max(dcay,dhat)))
lines(kap,dcay,col=2)
###########################
###########################
#Compare theoretical to empirical c and d values for distributions

library(rotations2)

#####
#Cayley distribution
kap<-25
rs<-rcayley(1000,kappa=kap)
crs<-cos(rs)

#chat
mean(1-crs^2)*2/3
(4*kap+2)/(kap^2+5*kap+6)

#dhat
mean(1+2*crs)/3
kap/(kap+2)

#ctilde
mean(1+crs)/6
(2*kap+1)/(6*kap+12)

#dtilde
mean((1+3*crs)/(12*sqrt(1-crs)))
sqrt(2/pi)*kap*gamma(kap+2)/(3*gamma(kap+2.5))

#####
#Fisher distribution
library(gsl)
kap<-25

rs<-rfisher(1000,kappa=kap)
crs<-cos(rs)

#chat
mean(1-crs^2)*2/3
(((kap+1)/kap)*besselI(2*kap,1)-besselI(2*kap,0))/(3*kap*(besselI(2*kap,0)-besselI(2*kap,1)))

#dhat
mean(1+2*crs)/3
(((kap+1)/kap)*besselI(2*kap,1)-besselI(2*kap,0))/(3*(besselI(2*kap,0)-besselI(2*kap,1)))


#ctilde
mean(1+crs)/6
besselI(2*kap,1)/(12*kap*(besselI(2*kap,0)-besselI(2*kap,1)))

#dtilde
mean((1+3*crs)/(12*sqrt(1-crs)))
exp(2*kap)*(6*sqrt(kap)-(3+8*kap)*dawson(2*sqrt(kap)))/(24*sqrt(2)*pi*kap^1.5*(besselI(2*kap,0)-besselI(2*kap,1)))


#####
#von Mises
kap<-25
rs<-rvmises(1000,kappa=kap)
crs<-cos(rs)

#chat
mean(1-crs^2)*2/3
(1-besselI(kap,2)/besselI(kap,0))/3
(besselI(kap,0)-besselI(kap,2))/(3*besselI(kap,0))

#dhat
mean(1+2*crs)/3
(besselI(kap,0)+2*besselI(kap,1))/(3*besselI(kap,0))


#ctilde
mean(1+crs)/6
(besselI(kap,0)+besselI(kap,1))/(6*besselI(kap,0))

####################
#compare moments to theory moments

#Cayley dist
mean(1/sqrt(1-crs))
sqrt(2)*gamma(kap+2)/((2*kap+1)*gamma(kap+.5)*gamma(1.5)) #Pretty good

mean(crs/sqrt(1-crs))
(2*kap-1)*gamma(kap+2)/(sqrt(2*pi)*gamma(kap+2.5))  #Pretty good


#Fisher dist
exp(2*kap)*sqrt(2)*dawson(2*sqrt(kap))/(pi*sqrt(kap)*(besselI(2*kap,0)-besselI(2*kap,1)))

mean(crs/sqrt(1-crs))
exp(2*kap)*(2*sqrt(kap)-(1+4*kap)*dawson(2*sqrt(kap)))/((2*kap)^(1.5)*pi*(besselI(2*kap,0)-besselI(2*kap,1)))  #Pretty good

kap<-50
rs<-rfisher(1000,kappa=kap)

mean(cos(rs))
(besselI(2*kap,1)-.5*besselI(2*kap,0)-.5*besselI(2*kap,2))/(besselI(2*kap,0)-besselI(2*kap,1))
(((2*kap+1)/(2*kap))*besselI(2*kap,1)-besselI(2*kap,0))/(besselI(2*kap,0)-besselI(2*kap,1))


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

levels(cRateM$Method)<-c("Theory","Bootstrap")

levels(cRateM$Dist)<-c("Cayley","matrix~~Fisher")
cRateM$nu<-factor(cRateM$nu,labels=c("nu == 0.25","nu == 0.5","nu == 0.75"))


#if can't find "unit" run library(grid)
qplot(n,value,data=cRateM,colour=Method,group=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
	facet_grid(Dist~nu,labeller=label_parsed)+
	geom_hline(yintercept=(1-alp)*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
	scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+theme(panel.margin=unit(0.5,'lines'))

#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals/Figures/CoverRatesB10000Median.pdf",width=7,height=4.5)
#ggsave("CoverRatesB1000Median.pdf",width=5,height=4)

#library(xtable)
#xtable(coverRate)