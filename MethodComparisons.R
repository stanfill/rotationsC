library(Rcpp)
library(rotations)
library(reshape2)
library(microbenchmark)
Rcpp::sourceCpp("ZhangMethod.cpp")
Rcpp::sourceCpp("FisherMethod.cpp")
source("IntervalFuns.R")

Qs<-ruars(20,rcayley,space='Q4',kappa=10)
cf<-fisherBootC(Qs,300)
rf<-fisherAxisBoot(Qs,300)
pMax<-max(c(cf,rf))
ss<-seq(0,pMax,length=1000)

par(mfrow=c(1,2))
hist(cf,breaks=100,xlim=c(0,pMax),prob=T)
lines(ss,dchisq(ss,3))
hist(rf,breaks=100,xlim=c(0,pMax),prob=T)
lines(ss,dchisq(ss,3))


#Check each element to see if the functions match
n<-20
Qs<-ruars(n,rcayley,space='Q4',kappa=10)
Qhat<-meanQ4C(Qs)

fish<-rep(0,2)

for(i in 1:100){
	Qstar<-Qs[sample(n,replace=T),]
	fish[1]<-fisherAxisC(Qstar,Qhat)
	fish[2]<-fisherAxisCompute(Qstar,as.Q4(matrix(Qhat,1,4)))
	print(diff(fish))
}

tim<-microbenchmark(
	fisherBootC(Qs,300),
	fisherAxisBoot(Qs,300)
	)

print(tim)
plot(tim)