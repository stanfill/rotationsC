#Note about these simulations:
#
#The C++ version of genR() seems to crash for the fisher or von Mises distributions but works
#with relative reliability for the Cayley distribution.  For now I've removed any call to the
#C++ version of genR(), namely "SO3defaultC()" from within "SO3.default()", which I thought I had
#already done.  Not it runs with little issue!

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
		rs<-rfisher(n[j],kappa=kap)
		Rs<-genR(rs)
		
		#ars<-abs(rs)
		cosrs<-cos(rs)
		crs<-(cosrs+1)
		drs<-(1+2*cosrs-3*cosrs^2)/(1-cosrs)^1.5
		#cosrs2<-cos(rs/2)^2
		#cotrs<-cot(rs)


		c<-mean(crs)/6  #I think this is C for proj median according to notes from 7/16
		#c<-1

		d<-mean(drs)/12
		#d<-1
		
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



###########
#Compare 