###############################################################
#Use the Rcpp functions to perform the simulation study
#Compute coverage rates for the Rivest, Fisher and our methods
###############################################################

library(Rcpp)
Rcpp::sourceCpp('ZhangMethod.cpp')
Rcpp::sourceCpp("FisherMethod.cpp")
library(rotations)
library(reshape2)
library(plyr)
source("IntervalFuns.R")


n<-c(10,50,100)

nu<-c(0.25,0.50,0.75)
cayKap<-c(10,4,2)
fishKap<-c(3.17,1.71,1.15)
misesKap<-c(2.40,1.16,0.52)

Dist<-c("Cayley","Fisher","Mises")


B<-1000				#Number of samples to use to estimate coverage probability (Zhang used 10,000)
m<-300				#Number of bootstrap samples used to estimate bootstrap test statistic (Zhang used 300)
alp<-0.9

resultsDf<-data.frame(expand.grid(Dist=Dist,nu=nu,n=n),Rivest=0,Fisher=0,NormalMean=0,PivotMean=0)

date()

for(p in 18:nrow(resultsDf)){

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
			zhangMean<-bootQhat(Qs,m)		        #Bootstrap Pivotal cut point
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
}


date()
resultsDf

write.csv(resultsDf,"ResultsB1000M300.csv")

resM<-melt(resultsDf,id=c('Dist','nu','n'))
qplot(n,value,data=resM,facets=Dist~nu,colour=variable,
      geom='line',lwd=I(1.25))+geom_hline(yintercept=alp*100)


###############################################################
#Plot the test statistics cumulative dist function versus the 
#theoretical limiting dist
###############################################################

library(Rcpp)
Rcpp::sourceCpp('ZhangMethod.cpp')
Rcpp::sourceCpp("FisherMethod.cpp")
library(rotations)
library(reshape2)
library(plyr)
source("IntervalFuns.R")

n<-c(10,50,100,300)
ks<-c(2,8)
B<-1000				#Number of samples to use to estimate coverage probability (Zhang used 10,000)

c<-2*(1-(ks^2-ks+3)/(ks^2+5*ks+6))/3
d<-(1+2*(ks-1)/(ks+2))/3

simSize<-length(n)*length(ks)

tMat<-matrix(0,simSize,B)

resultsDfMean<-data.frame(expand.grid(kappa=ks,n=n),tMat)


for(j in 1:simSize){
	kapInd<-which(ks==resultsDfMean$kappa[j])
	cj<-c[kapInd]
	dj<-d[kapInd]
	
	for(i in 1:B){
		
		Rs<-ruars(resultsDfMean$n[j],rcayley,kappa=resultsDfMean$kappa[j],space='Q4')
		
		ShatMean<-meanQ4C(Rs)
		
		hsqMean<-RdistC(ShatMean,id.Q4)^2

		resultsDfMean[j,(2+i)]<-2*resultsDfMean$n[j]*dj^2*hsqMean/cj
		
	}
}

colNums<-ncol(resultsDfMean)

whichK<-1 # 1 for kappa=2 and 2 for kappa=8

tInt10Mean<-sort(as.numeric(as.character(resultsDfMean[whichK,3:colNums])))
tInt30Mean<-sort(as.numeric(as.character(resultsDfMean[(whichK+2),3:colNums])))
tInt100Mean<-sort(as.numeric(as.character(resultsDfMean[(whichK+4),3:colNums])))
tInt300Mean<-sort(as.numeric(as.character(resultsDfMean[(whichK+6),3:colNums])))

plot(tInt300Mean,ecdf(tInt300Mean),main=expression(kappa==2),type='l',col=2,xlab='x',ylab='F(x)',xlim=c(0,10))
lines(tInt100Mean,ecdf(tInt100Mean),col=3)
lines(tInt30Mean,ecdf(tInt30Mean),col=4)
lines(tInt10Mean,ecdf(tInt10Mean),col=5)
lines(tInt300Mean,pchisq(tInt300Mean,3))
legend(7,.6,c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"),col=c(1,5,4,3,2),lty=1,lwd=2,bty='n')

whichK<-2 # 1 for kappa=2 and 2 for kappa=8
tInt10Mean<-sort(as.numeric(as.character(resultsDfMean[whichK,3:colNums])))
tInt30Mean<-sort(as.numeric(as.character(resultsDfMean[(whichK+2),3:colNums])))
tInt100Mean<-sort(as.numeric(as.character(resultsDfMean[(whichK+4),3:colNums])))
tInt300Mean<-sort(as.numeric(as.character(resultsDfMean[(whichK+6),3:colNums])))

plot(tInt300Mean,ecdf(tInt300Mean),main=expression(kappa==8),type='l',col=2,xlab='x',ylab='F(x)',xlim=c(0,10))
lines(tInt100Mean,ecdf(tInt100Mean),col=3)
lines(tInt30Mean,ecdf(tInt30Mean),col=4)
lines(tInt10Mean,ecdf(tInt10Mean),col=5)
lines(tInt300Mean,pchisq(tInt300Mean,3))
legend(7,.6,c(expression(chi[3]^2),"n=10","n=50","n=100","n=300"),col=c(1,5,4,3,2),lty=1,lwd=2,bty='n')