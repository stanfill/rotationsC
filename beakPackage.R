Rs<-ruars(20,rvmises,kappa=.1)
medianC(Rs,type='geometric')
HartmedianSO3C(Rs,100,1e-5)


Rs2<-ruars(20,rvmises)
medianC(Rs2,type='geometric')
HartmedianSO3C(Rs2,100,1e-5);



tm0symm<-rep(0,1000)
tm0<-rep(0,1000)

for(i in 1:1000){
	
	Qs<-ruars(100,rcayley,space='Q4')
	
	tm0symm[i]<-fisherAxisCSymmetric(Qs,id.Q4)
	
	tm0[i]<-fisherAxisC(Qs,id.Q4)

}

hist(tm0symm,prob=T)
hist(tm0,prob=T)

ss<-seq(0,max(tm0symm),length=length(tm0symm))
plot(tm0symm,ecdf(tm0symm))
lines(ss,pchisq(ss,3))

ss<-seq(0,max(tm0),length=length(tm0))
plot(tm0,ecdf(tm0))
lines(ss,pchisq(ss,3))

######
#

n<-100
rs<-rcayley(n)
Qs<-genR(rs,space='Q4')

fisherAxisCSymmetric(Qs,id.Q4)

d<-mean((1+2*cos(rs))/3)
c<-mean(2*(1-cos(rs)^2)/3)

n*c/(2*d^2)



############
#

for(i in 1:5000){
	Qs<-ruars(100,rcayley,space='Q4')
	s<-mean(Qs)
	s<-mean(Qs,type='geometric')
	s<-median(Qs)
	s<-median(Qs,type='geometric')
}
