###############################################################
#Use the Rcpp functions to perform the simulation study
#Compute coverage rates for the Rivest, Fisher and our methods
###############################################################

library(rotations2)
#library(RcppArmadillo)
#sourceCpp("genQ4only.cpp")
#source("rotations2/R/distributions.R")
#source("rotations2/R/preliminary.R")
#source("rotations2/R/parameterizations.R")
#library(rotations)
#source("IntervalFuns.R")

n<-c(10,20,50,100)

nu<-c(0.25,0.50,0.75)
cayKap<-c(10,4,2)
fishKap<-c(3.17,1.71,1.15)
misesKap<-c(2.40,1.16,0.52)

Dist<-c("Cayley","Fisher","Mises")

B<-1000				#Number of samples to use to estimate coverage probability (Zhang used 10,000)
alp<-0.9

resultsDf<-data.frame(expand.grid(Dist=Dist,nu=nu,n=n),Rivest=0)

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
	Rivest<-0
		for(k in 1:B){

      rs<-rfn(np,kappa=kapp)
      #print(rs)
     
      theta <- acos(runif(np, -1, 1))      
      phi <- runif(np, -pi, pi)
      u <- matrix(c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)),np,3)
      #Qs<-Q4defaultC(u,rs)
      #Qs<-genrC2(rs,space='Q4')
      Qs<-genrC2(rs,space='SO3')
      #Execute the Method in Rancourt 2000
			#ti<-RivestCI2(Qs)
			#Rivest<-Rivest+as.numeric(ti<qf(alp,3,np-3))	
      #ti<-RivestCI(Qs)
      #Rivest<-Rivest+as.numeric(ti<qchisq(alp,3))  
		}
	#resultsDf[p,]$Rivest<-100*Rivest/B
	#print(resultsDf[p,])
	#write.csv(resultsDf,"Results/ResultsB5000M300Part2.csv")
}
date()
resultsDf

#write.csv(resultsDf,"Results/ResultsB5000M300Part2.csv")

resM<-melt(resultsDf,id=c('Dist','nu','n'))
colnames(resM)[4]<-'Method'
levels(resM$Method)[3:4]<-c("Nordman Normal","Nordman Bootstrap")

levels(resM$Dist)<-c("Cayley","matrix~~Fisher","circular-von~~Mises")
resM$nu<-factor(resM$nu,labels=c("nu == 0.25","nu == 0.50","nu == 0.75"))
resM$n<-as.factor(resM$n)

qplot(n,value,data=resM,colour=Method,group=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
	facet_grid(Dist~nu,labeller=label_parsed)+
	geom_hline(yintercept=alp*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
	theme_bw()

