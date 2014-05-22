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