#Use the Rcpp functions to perform the simulation study

library(Rcpp)
Rcpp::sourceCpp('ZhangMethod.cpp')
library(rotations)
#setwd("U:/Thesis/Intervals/Code")
#source("IntervalFuns.R")


n<-c(20,50,100)

nu<-c(0.25,0.50,0.75)
cayKap<-c(10,4,2)
fishKap<-c(3.17,1.71,1.15)
misesKap<-c(2.40,1.16,0.52)

Dist<-c("Cayley","Fisher","Mises")


B<-1000				#Number of samples to use to estimate coverage probability (Zhang used 10,000)
m<-300				#Number of bootstrap samples used to estimate bootstrap test statistic (Zhang used 300)
alp<-0.95

resultsDf<-data.frame(expand.grid(Dist=Dist,nu=nu,n=n),Rivest=0,Fisher=0,NormalMean=0,PivotMean=0)

date()

for(p in 1:nrow(resultsDf)){

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
			#qs<-Q4(Rs)			
			#ti<-RivestCI(qs)
			#Rivest<-Rivest+as.numeric(ti<qchisq(alp,3))
			
			#Execute the Fisher, Hall, Jing and Wood Method
			#Criticaltf<-fisherAxisBoot(qs,m,alp)
			#testStat<-fisherAxisCompute(qs,id.Q4)
			#Fisher<-Fisher+as.numeric(testStat<Criticaltf)
			
			
			#Now for the methods in Zhang's MS for the Projected Mean
			zhangMean<-bootQhat(Qs,m)		        #Bootstrap Pivotal cut point
			cutPt<-quantile(zhangMean,alp)
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
	
}


date()
resultsDf