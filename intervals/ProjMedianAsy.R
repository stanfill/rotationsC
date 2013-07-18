library(rotations2)
setwd("C:/Users/stanfill/Desktop/GitHub/rotationsC/intervals")
source("IntervalFuns.R")
B<-500
n<-c(10,50,100)
numN<-length(n)
kap<-.1
tstats<-matrix(0,B,numN)

for(j in 1:numN){
	
	for(i in 1:B){
		rs<-rcayley(n[j],kappa=kap)
		Rs<-genR(rs)
		
		#ars<-abs(rs)
		
		cos2rs<-cos(rs/2)^2
		#drs<-(3*cos(rs)+1)/(sin(rs/2)*6*sqrt(2))
		#cosrs2<-cos(rs/2)^2
		#cotrs<-cot(rs)

		ecos2<-mean(cos2rs)
		#ecos2<-mean(cosrs2)
		#ecot<-mean(cotrs)

		c<-(4/3)*(ecos2)  #I think this is C for proj median according to notes from 7/16
		#c<-1

		#d<--mean(drs)
		d<-1
		Shat<-median(Rs)
	
		ShatMedian<-dist(Shat,method='intrinsic',p=2)
	
		tstats[i,j]<-2*n[j]*(d^2)*ShatMedian/c
	
	}
	
}
#hist(tstats,breaks=100,prob=T)

#ses<-seq(0,max(tstats),length=B)
#lines(ses,dchisq(ses,3))

for(j in 1:numN){
	tstats[,j]<-sort(tstats[,j])
}
plot(tstats[,1],ecdf(tstats[,1]),type='l')

for(j in 2:numN){
	lines(tstats[,j],ecdf(tstats[,j]),col=(j+1))
}

lines(tstats[,3],pchisq(tstats[,3],3),lty=2)


