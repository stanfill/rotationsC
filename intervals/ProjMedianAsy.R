library(rotations2)

B<-500
n<-20
kap<-100
tstats<-rep(0,B)

for(i in 1:B){
	rs<-rcayley(n,kappa=kap)
	Rs<-genR(rs)
	
	cosrs<-cos(rs)
	#cotrs<-cot(rs)

	ecos<-mean(cosrs)
	#ecot<-mean(cotrs)

	c<-(1/6)*(ecos+1)  #I think this is C for proj median according to notes from 7/16

	d<--(1/sqrt(72))*(mean(cos(rs/2)/tan(rs/2)))
	

	Shat<-median(Rs)
	
	ShatMedian<-dist(Shat,method='intrinsic',p=2)
	

	tstats[i]<-2*n*d^2*ShatMedian/c
	
}

hist(tstats,breaks=100,prob=T)

ses<-seq(0,max(tstats),length=B)
lines(ses,dchisq(ses,3))