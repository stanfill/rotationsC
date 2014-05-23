###########
#How good an approximation are the first 2,3,4 terms of the expansion of asin?


dr2Approx<-function(R,S,order=3){
  S<-matrix(S,3,3)
  tri<-rep(0,nrow(R))
  for(i in 1:nrow(R)){
    Ri<-matrix(R[i,],3,3)
    tri[i]<-sum(diag(t(Ri)%*%S))
  }
  
  dr2<-(3-tri)+(3-tri)^2/12
  
  if(order==3){
    dr2<-dr2+(3-tri)^3/90
  }else if(order==4){
    dr2<-dr2+(3-tri)^3/90+(3-tri)^4/560
  }
  
  #if(order==4){
    #dr2<-2349/560-279*tri/140+47*tri^2/168-41*tri^3/1260+tri^4/560
  #}else if(order==3){
    #dr2<-81/20-9*tri/5+11*tri^2/60-tri^3/90
  #}else if(order==2){
  #}
  return(dr2)
}

Rs<-ruars(20,rcayley,kappa=1)
dr<-rot.dist(Rs,method='intrinsic',p=2)
adr2<-dr2Approx(Rs,id.SO3,order=2)
adr3<-dr2Approx(Rs,id.SO3,order=3)
adr4<-dr2Approx(Rs,id.SO3,order=4)
de<-rot.dist(Rs,p=2)

plot(dr,adr2,pch=19)
abline(0,1)
points(dr,adr3,col=2,pch=19)
points(dr,adr4,col=3,pch=19)

####################
#How many terms are needed to approximate [asin(x)]^2?

asin2Approx<-function(x,k){
  #x - value at which to compue asin(x)^2
  #k - number of terms to take the expansion
  
  asin2<-0
  for(i in 1:k){
    num <- (2*x)^(2*i)
    denom <- 2*i^2*choose(2*i,i)
    
    asin2 <- asin2 + num/denom
  }
  return(asin2)
}

rs<-seq(0,1,length=100)
ars<-asin(rs)^2
plot(rs,ars,type='l')

for(i in 1:5){
  lines(rs,asin2Approx(rs,k=i),col=(i+1))
}

legend("topleft",c(1:5,"Inf"),col=c(2:6,1),lty=1)

############################
##How close is tr(R'hat(S))=1+2cos(r)

S<-genR(0)
Rs<-ruars(200,rcayley,S=S,kappa=1)
Shat<-mean(Rs,type='geo')
trRtS<-rep(0,nrow(Rs))

for(i in 1:nrow(Rs)){
  Ri<-matrix(Rs[i,],3,3)
  trRtS[i] <- sum(diag(t(Ri)%*%Shat))
}

act<-1+2*cos(mis.angle(Rs-S))

plot(trRtS,act,pch=19)
abline(0,1)
mean(trRtS-act)
