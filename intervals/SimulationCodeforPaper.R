###############################################################
#Use the Rcpp functions to perform the simulation study
#Compute coverage rates for the Rivest, Fisher and our methods
###############################################################

library(RcppArmadillo)
Rcpp::sourceCpp('ZhangMethod.cpp')
Rcpp::sourceCpp("FisherMethod.cpp")
library(rotations)
library(reshape2)
library(plyr)
library(grid)
source("IntervalFuns.R")


n<-c(10,20,50,100)

nu<-c(0.25,0.50,0.75)
cayKap<-c(10,4,2)
fishKap<-c(3.17,1.71,1.15)
misesKap<-c(2.40,1.16,0.52)

Dist<-c("Cayley","Fisher","Mises")


B<-5000				#Number of samples to use to estimate coverage probability (Zhang used 10,000)
m<-300				#Number of bootstrap samples used to estimate bootstrap test statistic (Zhang used 300)
alp<-0.9

#resultsDf<-data.frame(expand.grid(Dist=Dist,nu=nu,n=n),Rivest=0,Fisher=0,NormalMean=0,PivotMean=0)
resultsDf<-read.csv("Results/ResultsB5000M300Part2.csv")[,-1]
resultsDf
date()

for(p in 36:nrow(resultsDf)){

	distp<-resultsDf$Dist[p]
	np<-resultsDf$n[p]
	nup<-resultsDf$nu[p]
	
	if(distp=='Mises'){
		
		kapp<-misesKap[which(nu==nup)]
		rfn<-rvmises
		
	}else if(distp=='Cayley'){
		
		kapp<-cayKap[which(nu==nup)]
		rfn<-rcayley
		
	}else{
		
		kapp<-fishKap[which(nu==nup)]
		rfn<-rfisher
		
	}	
	Rivest<-Fisher<-PivotMean<-NormalMean<-0
		for(k in 1:B){
			
			Qs<-ruars(np,rfn,kappa=kapp,space='Q4')
			
			#Execute the Method in Rancourt 2000
						
			ti<-RivestCI(Qs)
			Rivest<-Rivest+as.numeric(ti<qchisq(alp,3))
			
			#Execute the Fisher, Hall, Jing and Wood Method
			fishC<-fisherBootC(Qs,m)
			Criticaltf<-quantile(fishC,alp,na.rm=T)
			testStat<-fisherAxisC(Qs,id.Q4)
			
			#Criticaltf
			#testStat
			#hist(fishC,breaks=100,prob=T)
			
			Fisher<-Fisher+as.numeric(testStat<Criticaltf)
			
			
			#Now for the methods in Zhang's MS for the Projected Mean
			zhangMean<-zhangQ4(Qs,m)		        #Bootstrap Pivotal cut point
			cutPt<-quantile(zhangMean,alp,na.rm=T)
			ShatE<-as.Q4(meanQ4C(Qs))													#Projected Mean of Shat
			RscdMean<-cdfunsC(Qs,ShatE)												#Compute c-hat and d-hat for sample quantity 
			
			#Use the exact test stat
			tStatNonMean<-RdistC(ShatE,id.Q4)^2
			tStatMean<-(2*np*RscdMean[2]^2*tStatNonMean/RscdMean[1])	#Pivotal sample test-stat
			
			NormalMean<-NormalMean+as.numeric(tStatMean<qchisq(alp,3)) #Cover true mean under normal assumption
			PivotMean<-PivotMean+as.numeric(tStatMean<cutPt)	 #Cover true mean with pivot Bootstrap
				 
			
		}
	
	resultsDf[p,]$Rivest<-100*Rivest/B
	resultsDf[p,]$Fisher<-100*Fisher/B
	
	resultsDf[p,]$NormalMean<-100*NormalMean/B
	resultsDf[p,]$PivotMean<-100*PivotMean/B
	print(resultsDf[p,])
	#write.csv(resultsDf,"Results/ResultsB5000M300Part2.csv")
}


date()
resultsDf

#write.csv(resultsDf,"Results/ResultsB5000M300Part2.csv")
#resultsDf<-read.csv("Results/ResultsB10000M300.csv")[,-1]
resM<-melt(resultsDf,id=c('Dist','nu','n'))
colnames(resM)[4]<-'Method'
#levels(resM$Method)[3:4]<-c("Nordman Normal","Nordman Bootstrap")
levels(resM$Method)<-c("Directional(LSA)","Directional(Boot)","Direct(LSA)","Direct(Boot)")

levels(resM$Dist)<-c("Cayley","matrix~~Fisher","circular-von~~Mises")
#Order alphabetically
resM$Dist<-factor(resM$Dist,levels=c("Cayley","circular-von~~Mises","matrix~~Fisher"))
#Order by tail weight
resM$Dist<-factor(resM$Dist,levels=c("Cayley","matrix~~Fisher","circular-von~~Mises"))


resM2<-resM
resM$nu<-factor(resM$nu,labels=c("nu == 0.25","nu == 0.50","nu == 0.75"))
#resM$n<-as.factor(resM$n)

qplot(n,value,data=resM,colour=Method,group=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
	facet_grid(Dist~nu,labeller=label_parsed)+
	geom_hline(yintercept=alp*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
	scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals/Figures/CoverRatesB10000.pdf",width=8,height=6)

qplot(n,value,data=resM,col=Method,linetype=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
  scale_linetype_manual(values = rep(c("dashed","solid"),2))+
  scale_colour_manual(values = rep(c("gray50","black"),each=2))+
  facet_grid(Dist~nu,labeller=label_parsed)+
  geom_hline(yintercept=(alp+c(-.01,.01,0))*100,colour=c('gray75','gray75','gray50'))+geom_line(lwd=I(1.25),alpha=I(.8))+
  scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+
  theme(legend.key.width=unit(3,"line"),legend.position='top',panel.margin = unit(.75, "lines"))
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Mean/Figures/CoverRatesB10000.pdf",width=7,height=6)

qplot(n,value,data=resM[resM$nu!="nu == 0.50",],col=Method,linetype=Method,ylab='Coverage Rate (%)',xlab='Sample Size',ylim=c(80,100),geom='blank')+
  scale_linetype_manual(values = rep(c("dashed","solid"),2))+
  scale_colour_manual(values = rep(c("gray50","black"),each=2))+
  facet_grid(Dist~nu,labeller=label_parsed)+
  geom_hline(yintercept=(alp+c(-.01,.01,0))*100,colour=c('gray75','gray75','gray50'))+geom_line(lwd=I(1),alpha=I(.8))+geom_point(size=I(1))+
  scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+coord_fixed(2.5)+
  theme(legend.key.width=unit(1.5,"line"),legend.position='top',panel.margin = unit(.75, "lines"))
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals - Mean/Figures/CoverRatesB10000No5.pdf",width=7,height=6)

qplot(n,value,data=resM[resM$nu!="nu == 0.50",],col=Method,linetype=Method,ylab='Coverage Rate (%)',xlab='Sample Size',ylim=c(80,100),geom='blank')+
  scale_linetype_manual(values = rep(c("dashed","solid"),2))+
  scale_colour_manual(values = rep(c("red","black"),each=2))+
  facet_grid(Dist~nu,labeller=label_parsed)+
  geom_hline(yintercept=(alp+c(-.01,.01,0))*100,colour=c('gray75','gray75','gray50'))+geom_line(lwd=I(1),alpha=I(.8))+geom_point(size=I(1))+
  scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+coord_fixed(2.5)+
  theme(legend.key.width=unit(1.5,"line"),legend.position='top',panel.margin = unit(.75, "lines"))
#ggsave("C:/Users/Brittney Ritchey/Dropbox/Thesis/Defense/figure/CoverRatesB10000No5Color.pdf",width=7,height=6)

#############
#Massage the plots to accentuate the differences as a function of n and nu

levels(resM2$Dist)<-c("Cayley","matrix Fisher","circular-von Mises")
resM2$LargeN<-"Small n"
resM2[resM2$n>30,]$LargeN<-"Large n"
resM2$LargeN<-factor(resM2$LargeN,levels=c("Small n","Large n"))
resM2$factorNu<-factor(resM2$nu,labels=c("nu == 0.25","nu == 0.50","nu == 0.75"))


qplot(n,value,data=resM2[resM2$nu==0.25,],col=Method,linetype=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
  scale_linetype_manual(values = rep(c("dashed","solid"),each=2))+
  scale_colour_manual(values = rep(c("gray50","black"),2))+
  facet_wrap(Dist~LargeN,scales='free',ncol=2)+
  geom_hline(yintercept=alp*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
  scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+
  theme(legend.key.width=unit(3,"line"),legend.position='top',panel.margin = unit(.75, "lines"))

qplot(n,value,data=resM2[resM2$nu==0.75,],col=Method,linetype=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
  scale_linetype_manual(values = rep(c("dashed","solid"),each=2))+
  scale_colour_manual(values = rep(c("gray50","black"),2))+
  facet_wrap(Dist~LargeN,scales='free',ncol=2)+
  geom_hline(yintercept=alp*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
  scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+
  theme(legend.key.width=unit(3,"line"),legend.position='top',panel.margin = unit(.75, "lines"))

resM2$Variability<-"High Var."
resM2[resM2$nu==0.5,]$Variability<-"Medium Var."
resM2[resM2$nu==0.25,]$Variability<-"Low Var."
resM2$Variability<-factor(resM2$Variability,levels=c("Low Var.","Medium Var.","High Var."))
levels(resM2$Dist)<-c("Cayley","matrix~~Fisher","circular-von~~Mises")
resM2$Dist<-factor(resM2$Dist,levels=c("Cayley","circular-von~~Mises","matrix~~Fisher"))

qplot(n,value,data=resM2[resM2$LargeN=="Small n" & resM2$nu!=0.5,],col=Method,linetype=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
  scale_linetype_manual(values = rep(c("dashed","solid"),each=2))+
  scale_colour_manual(values = rep(c("gray50","black"),2))+
  facet_grid(Dist~factorNu,labeller=label_parsed)+
  geom_hline(yintercept=alp*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
  scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+
  theme(legend.key.width=unit(3,"line"),legend.position='top',panel.margin = unit(.5, "lines"))

qplot(n,value,data=resM2[resM2$LargeN!="Small n" & resM2$nu!=0.5,],col=Method,linetype=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
  scale_linetype_manual(values = rep(c("dashed","solid"),each=2))+
  scale_colour_manual(values = rep(c("gray50","black"),2))+
  facet_grid(Dist~factorNu,labeller=label_parsed)+
  geom_hline(yintercept=alp*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
  scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+
  theme(legend.key.width=unit(3,"line"),legend.position='top',panel.margin = unit(.5, "lines"))

###############
#No circular von Mises for mean median comparison
resMnoCVM<-resM[resM$Dist!="circular-von~~Mises",]

qplot(n,value,data=resMnoCVM,col=Method,linetype=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
  scale_linetype_manual(values = rep(c("dashed","solid"),each=2))+
  scale_colour_manual(values = rep(c("gray50","black"),2))+
  facet_grid(Dist~nu,labeller=label_parsed)+
  geom_hline(yintercept=alp*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
  scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+
  theme(legend.key.width=unit(3,"line"))  
#ggsave("/Users/stanfill/Dropbox/Thesis/Intervals/Figures/CoverRatesB10000NoCVM.pdf",width=8,height=6)


###############################################################
#Plot the test statistics cumulative dist function versus the 
#theoretical limiting dist
###############################################################
setwd("/Users/stanfill/Documents/GitHub/rotationsC/intervals")
library(Rcpp)
Rcpp::sourceCpp('ZhangMethod.cpp')
Rcpp::sourceCpp("FisherMethod.cpp")
#Rcpp::sourceCpp('rotations2/src/estimators.cpp')
library(rotations)
library(reshape2)
library(plyr)
source("IntervalFuns.R")

n<-c(10,50,100,300)

ks<-c(2,8)
B<-1000				#Number of samples to use to estimate CDF
Dist<-'none'

if(Dist=='cayley'){
	rangle<-rcayley
}else if(Dist=='fisher'){
	rangle<-rfisher
}else{
	rangle<-rvmises
}

simSize<-length(n)*length(ks)

tMat<-matrix(0,simSize,B)

cdfDF<-data.frame(expand.grid(kappa=ks,n=n),tMat)

for(j in 1:simSize){

	for(i in 1:B){
		
		rs<-rangle(cdfDF$n[j],kappa=cdfDF$kappa[j])

		Rs<-genR(rs,space='Q4')
		#Rs<-genR(rs)
    
		ShatMean<-meanQ4C(Rs)
		#ShatMean<-HartmedianSO3C(Rs,100,1e-5)
		#ShatMean<-Q4(as.SO3(matrix(ShatMean,nrow=1)))
    
		rs2<-dist(Rs,ShatMean,method='intrinsic')
		cosrs<-cos(rs2)
    
		#cosrs<-cos(rs)
		
		ecos<-mean(cosrs)
		ecos2<-mean(cosrs^2)
		
		c<-(2/3)*(1-ecos2)
		d<-(1/3)*(1+2*ecos)
		
		hsqMean<-RdistC(ShatMean,id.Q4)^2

		cdfDF[j,(2+i)]<-2*cdfDF$n[j]*d^2*hsqMean/c
		
	}
}

resM<-melt(cdfDF,id=c("kappa","n"))
ss<-seq(0,max(resM$value,na.rm=T),length=B)
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
	
	fullDF$kappa<-factor(fullDF$kappa,labels=c("kappa == 2","kappa == 8"))
	
	qplot(value,Prob,data=fullDF,colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,10))+
		scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100",'n=300'))+
		facet_grid(.~kappa,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
    scale_size_discrete(range=c(1,2),guide='none')+theme(legend.position=c(.9,.3),panel.margin = unit(1, "lines"))+
	  guides(colour=guide_legend(label.hjust=0,override.aes=list(size=2)))+coord_equal(10)
	
	setwd("/Users/stanfill/Dropbox/Thesis/Intervals - Mean/Figures")
	ggsave("CayleyECDF.pdf",height=4,width=8)
	#write.csv(fullDF,"CayleyECDF.csv")
	
}else if(Dist=='fisher'){
	
	fullDF$kappa<-factor(fullDF$kappa,labels=c("kappa == 2","kappa == 8"))
	
	qplot(value,Prob,data=fullDF,colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,10))+
		scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
		facet_grid(.~kappa,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
	  scale_size_discrete(range=c(0.75,1.5),guide='none')+
	  guides(colour=guide_legend(label.hjust=0))+coord_equal(10)
	
	#setwd("/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
	#ggsave("FisherECDF.pdf",height=5,width=7)
	#write.csv(fullDF,"FisherECDF.csv")
	
}else{
	
	fullDF$kappa<-factor(fullDF$kappa,labels=c("kappa == 2","kappa == 8"))
	
	qplot(value,Prob,data=fullDF,colour=n,lwd=Stat,geom="line",xlab='x',ylab="F(x)",xlim=c(0,10))+
		scale_colour_grey("",labels=c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"))+
		facet_grid(.~kappa,labeller=label_parsed)+theme_bw()+coord_fixed(ratio=15/1)+
	  guides(colour=guide_legend(label.hjust=0))+
	  scale_size_discrete("",range=c(0.75,1.5),guide='none')+coord_equal(10)

	#setwd("/Users/stanfill/Dropbox/Thesis/Intervals/Figures")
	#ggsave("vonMisesECDF.pdf",height=5,width=8)
	#write.csv(fullDF,"vonMisesECDF.csv")
	
}
