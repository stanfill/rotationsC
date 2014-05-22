###########
#How good an approximation are the first 2,3,4 terms of the expansion of asin?


dr2Approx<-function(R,S,order=3){
  S<-matrix(S,3,3)
  tri<-rep(0,nrow(R))
  for(i in 1:nrow(R)){
    Ri<-matrix(R[i,],3,3)
    tri[i]<-sum(diag(t(Ri)%*%S))
  }
  if(order==4){
    dr2<-2349/560-279*tri/140+47*tri^2/168-41*tri^3/1260+tri^4/560
  }else{
    dr2<-81/20-9*tri/5+11*tri^2/60-tri^3/90
  }
  return(dr2)
}

Rs<-ruars(20,rcayley,kappa=1)
dr<-rot.dist(Rs,method='intrinsic',p=2)
adr<-dr2Approx(Rs,id.SO3)
adr2<-dr2Approx(Rs,id.SO3,order=4)
de<-rot.dist(Rs,p=2)

plot(adr,dr,pch=19)
abline(0,1)
points(adr,adr2,col=2,pch=19)

################

a3<-function(r){
  sin2<-sin(r)^2
  cos2<-cos(r)^2
  a33<-mean(sin2*cos2)+16*mean(sin2)
  return(a33)
}

a4<-function(r){
  cosr<-cos(r)
  cos2<-cosr^2
  a44<-4-10*mean(cosr)-3*mean(cos2)
  return(a44/3)
}


rs<-rcayley(500)
a3(rs)
a4(rs)

B<-500
n<-c(10,100)
numN<-length(n)
kap<-2
tstats<-matrix(c(rep(0,B),rep(0,B)),ncol=2)

#S<-as.SO3(matrix(genR(pi/4),3,3))
S<-id.SO3

for(j in 1:numN){
  
  for(i in 1:B){
    rs<-rcayley(n[j],kappa=kap)
    Rs<-genR(rs,S=S)
    
    Shat<-mean(Rs,type='geo')
    
    rs2<-rot.dist(Rs,Shat,method='intrinsic')
    a3i<-a3(rs2)
    a4i<-a4(rs2)
    
    ShatMedian<-rot.dist(Shat,S,method='intrinsic',p=2)
    
    tstats[i,j]<-27*n[j]*(a4i^2)*ShatMedian/(4*a3i)
    
  } 

}

tstats2<-tstats

plot(ecdf(tstats2[,1]),col=2)
lines(ecdf(tstats2[,2]))
seqChi<-seq(0,max(tstats2),length=100)
lines(seqChi,pchisq(seqChi,3),lty=2,col=3,lwd=3)

