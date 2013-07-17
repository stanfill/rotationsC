library(rotations2)
source("IntervalFuns.R")
B<-500
n<-c(10,100)
kap<-10
tstats<-matrix(0,B,2)

for(j in 1:2){
	
	for(i in 1:B){
		rs<-rcayley(n[j],kappa=kap)
		Rs<-genR(rs)
	
		cos2rs<-cos(rs/2)^2
		cosrs<-cos(rs)
		#cosrs2<-cos(rs/2)^2
		#cotrs<-cot(rs)

		ecos2<-mean(cos2rs)
		#ecos2<-mean(cosrs2)
		#ecot<-mean(cotrs)

		c<-(1/3)*(ecos2)  #I think this is C for proj median according to notes from 7/16

		#d<--(1/sqrt(72))*(mean(cos(rs/2)/tan(rs/2)))
		d<-(1+mean(cosrs))#/(sqrt(2))

		Shat<-median(Rs)
	
		ShatMedian<-dist(Shat,method='intrinsic',p=2)
	
		tstats[i,j]<-2*n[j]*(d^2)*ShatMedian/c
	
	}
	
}
#hist(tstats,breaks=100,prob=T)

#ses<-seq(0,max(tstats),length=B)
#lines(ses,dchisq(ses,3))

tstats[,1]<-sort(tstats[,1])
tstats[,2]<-sort(tstats[,2])
plot(tstats[,1],ecdf(tstats[,1]),type='l')
lines(tstats[,2],ecdf(tstats[,2]),col=2)
lines(tstats[,1],pchisq(tstats[,1],3),col=3)


####################
#Just check c:

n=50
kap=50
rs<-rcayley(n,kappa=kap)
Rs<-genR(rs)

Shat<-median(Rs)
chat<-2*mean(cos(rs)^2)/3

dists<-dist(Rs,Shat,method='projected')

#chat and mean(dists^2) should be close
mean(dists^2)
chat


