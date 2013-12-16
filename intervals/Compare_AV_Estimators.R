library(rotations)
library(reshape2)
Rcpp::sourceCpp('intervals/ZhangMethod.cpp')

#For the Cayley distribution, Leon et al use a moment estimator for
#kappa to estimate the asymptotic variance for the projected mean

Rs<-ruars(20,rcayley,kappa=1)
Shat<-mean(Rs)

Ps<-matrix(Rs,20,9)
Phat<-matrix(colMeans(Ps),3,3)
decom<-svd(Phat)

Phat2<-decom$u%*%t(decom$v)
kHat<-(2)*mean(decom$d)/(1-mean(decom$d))

#This Leon's estimate of asymtptotic variance, Section 6
(kHat+2)*(2*kHat+1)/(kHat^2*(kHat+3))

#This is our estimate of asymptotic variance c/2d^2
Qs<-Q4(Rs)
RscdMean<-cdfunsC(Qs,mean(Qs))  

RscdMean[1]/(2*RscdMean[2]^2)

#Simulate a lot of them and see how they compare
k<-2.5
B<-1000
mine<-rep(0,B)
leon<-rep(0,B)
truth<-(k+2)*(2*k+1)/(k^2*(k+3))
AVcomp<-data.frame(mine=rep(0,B),leon=rep(0,B))

for(i in 1:B){
  Rs<-ruars(20,rcayley,kappa=k)
  Ps<-matrix(Rs,20,9)
  Phat<-matrix(colMeans(Ps),3,3)
  decom<-svd(Phat)

  kHat<-(2)*mean(decom$d)/(1-mean(decom$d))
  AVcomp$leon[i]<-(kHat+2)*(2*kHat+1)/(kHat^2*(kHat+3))
  Qs<-Q4(Rs)
  RscdMean<-cdfunsC(Qs,mean(Qs))  
  AVcomp$mine[i]<-RscdMean[1]/(2*RscdMean[2]^2)
}
mAVcomp<-melt(AVcomp)

qplot(value,data=mAVcomp,fill=variable,group=variable,geom='density',alpha=I(.4))+geom_vline(xintercept=truth,col=2)