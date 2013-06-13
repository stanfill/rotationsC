###############################################################
#Plot the test statistics cumulative dist function versus the 
#theoretical limiting dist
###############################################################
setwd("C:/Users/stanfill/Desktop/GitHub/rotationsC")
library(Rcpp)
Rcpp::sourceCpp('ZhangMethod.cpp')
#Rcpp::sourceCpp("FisherMethod.cpp")
Rcpp::sourceCpp('rotations2/src/estimators.cpp')
library(rotations)
library(reshape2)
library(plyr)
source("IntervalFuns.R")

n<-c(10,50,100,300)
B<-1000				#Number of samples to use to estimate CDF
Dist<-'cayley'

if(Dist=='cayley'){
	rangle<-rcayley
	ks<-c(0.5,8)
}else if(Dist=='fisher'){
	rangle<-rfisher
	ks<-c(2,8)
}else{
	rangle<-rvmises
	ks<-c(2,8)
}

simSize<-length(n)*length(ks)

tMat<-matrix(0,simSize,B)


cdfDFMean<-data.frame(expand.grid(kappa=ks,n=n),tMat)
#cdfDFMedian<-data.frame(expand.grid(kappa=ks,n=n),tMat)

for(j in 1:simSize){
	
	
	for(i in 1:B){
		
		rs<-rangle(cdfDFMean$n[j],kappa=cdfDFMean$kappa[j])
		ars<-abs(rs)
		#Qs<-genR(rs,space='Q4')
		#Rs<-SO3(Qs)
		Rs<-genR(rs)
		
		cosrs<-cos(rs)
		
		#ers<-mean(rs)
		#ers2<-mean(rs^2)
		#rfn<-mean((ars+2*(cos(ars)))/(sin(ars)))
		ecos<-mean(cosrs)
		ecos2<-mean(cosrs^2)
		
		#c<-(2/3)*(ers2)
		c<-(2/3)*(1-ecos2)
		#d<-(1/3)*(rfn)
		d<-(1/3)*(1+2*ecos)
		
		#Shat<-as.SO3(gmeanSO3C(Rs,200,1e-5))
		Shat<-as.SO3(meanSO3C(Rs))
		
		ShatMean<-rdistSO3C(Shat,id.SO3)^2
		
# 		if(is.na(ShatMean)){
# 			Shat<-median(Rs,type='geometric')
# 			ShatMean<-rdistSO3C(Shat,id.SO3)^2
# 			print(ShatMean)
# 		}
		
		#med<-as.SO3(medianSO3C(Rs,100,1e-5))
		#ShatMed<-dist(med,method='intrinsic')^2
		
		#ShatMean<-HartmedianSO3C(Rs,100,1e-5)
		#ShatMean<-Q4(as.SO3(matrix(ShatMean,nrow=1)))
		
		#hsqMean<-RdistC(ShatMean,id.Q4)^2
		
		cdfDFMean[j,(2+i)]<-2*cdfDFMean$n[j]*d^2*ShatMean/c
		#cdfDFMedian[j,(2+i)]<-2*cdfDF$n[j]*d^2*ShatMed/c

		
	}
}

resM<-melt(cdfDFMean,id=c("kappa","n"))
ss<-seq(0,max(resM$value),length=B)
Probs<-pchisq(ss,3)

resM$n<-as.factor(resM$n)
resM$Prob<-0
resM$ID<-paste(resM$kappa,resM$n)
kns<-unique(resM$ID)

for(i in 1:length(unique(resM$ID))){
	resM[resM$ID==kns[i],]$value<-sort(resM[resM$ID==kns[i],]$value)
	resM[resM$ID==kns[i],]$Prob<-ecdf(resM[resM$ID==kns[i],]$value)
}

chiDF<-data.frame(kappa=rep(ks,1000),n='Chisq',variable='Tr',value=rep(ss,each=2),Prob=rep(pchisq(ss,3),each=2))
chiDF$ID<-paste(chiDF$kappa,chiDF$n)

fullDF<-rbind(resM,chiDF)

Newlabs<-c("Chisq","10","50","100","300")
fullDF$n<-factor(fullDF$n,levels=Newlabs)
fullDF$Stat<-1
fullDF[fullDF$n=='Chisq',]$Stat<-2
fullDF$Stat<-as.factor(fullDF$Stat)

if(Dist=='cayley'){
	
	fullDF$kappa<-factor(fullDF$kappa,labels=c("kappa == 2.0","kappa == 8.0"))
	
	qplot(value,Prob,data=fullDF,colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,15))+
		scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
		facet_grid(.~kappa,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
		scale_size_discrete(range=c(0.75,1.5),guide='none')+
		guides(colour=guide_legend(label.hjust=0))
	
	#setwd("C:/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
	#ggsave("CayleyECDF.pdf",height=5,width=8)
	#write.csv(fullDF,"CayleyECDF.csv")
	
}else if(Dist=='fisher'){
	
	fullDF$kappa<-factor(fullDF$kappa,labels=c("kappa == 2.0","kappa == 8.0"))
	
	qplot(value,Prob,data=fullDF,colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,15))+
		scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
		facet_grid(.~kappa,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
		scale_size_discrete(range=c(0.75,1.5),guide='none')+
		guides(colour=guide_legend(label.hjust=0))
	
	#setwd("C:/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
	#ggsave("GeometricFisherECDF.pdf",height=5,width=8)
	#write.csv(fullDF,"FisherECDF.csv")
	
}else{
	
	fullDF$kappa<-factor(fullDF$kappa,labels=c("kappa == 2.0","kappa == 8.0"))
	
	qplot(value,Prob,data=fullDF,colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,15))+
		scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
		facet_grid(.~kappa,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
		guides(colour=guide_legend(label.hjust=0))+
		scale_size_discrete("",range=c(0.75,1.5),guide='none')
	
	#setwd("C:/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
	#ggsave("vonMisesECDF.pdf",height=5,width=8)
	#write.csv(fullDF,"vonMisesECDF.csv")
	
}