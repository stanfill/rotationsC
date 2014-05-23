################
library(rotations)

a3<-function(r,terms=2){
  sin2<-sin(r)^2
  cos2<-cos(r)^2
  
  if(terms==2){
    a33<-mean(sin2*cos2)+16*mean(sin2)#-8*mean(sin2*cos(r))
    a33<-a33*4/27
  }
  if(terms==3){
    cos4<-cos(r)^4
    

    #p1<-676*mean(cos2*sin2)    #Old calculation
    
    p1<-1212*mean(cos2*sin2)
    p2<-1936*mean(sin2)
    p3<-16*mean(cos4*sin2)
    p4<-0-1584*mean(cos(r)*sin2)-144*mean(cos(r)^3*sin2)
    a33<-(p1+p2+p3+p4)/675
  }
  return(a33)
}

a4<-function(r,terms=2){
  
  cosr<-cos(r)
  cos2<-cosr^2
  
  if(terms==2){

    #a44<-(4-10*mean(cosr)-3*mean(cos2))/3    #Old calculation
    
    a44<-(8+14*mean(cosr)-4*mean(cos2))/9
    
  }
  if(terms==3){
    
    cos3<-cos(r)^3
    p1<-140*mean(cosr)
    p2<--64*mean(cos2)
    p3<-16*mean(cos3)
    a44<-88+p1+p2+p3
    a44<-a44/90
    
  }
  return(a44)
}


rs<-rcayley(500)
a3(rs)
a4(rs)

####################
#Compare a3 and a4 with the different expansions
rs<-matrix(rcayley(500,kappa=100),ncol=50)

plot(apply(rs,2,a3),apply(rs,2,a3,terms=3),pch=19,asp=1)
abline(0,1)

plot(apply(rs,2,a4),apply(rs,2,a4,terms=3),pch=19,asp=1)
abline(0,1)

####################
#Plot test stat ECDF versus chi^2_3 CDF

B<-500
n<-c(10,100)
ts<-2
numN<-length(n)
kap<-5
tstats<-matrix(c(rep(0,B),rep(0,B)),ncol=2)

#S<-as.SO3(matrix(genR(pi/4),3,3))
S<-id.SO3

for(j in 1:numN){
  
  for(i in 1:B){
    rs<-rcayley(n[j],kappa=kap)
    Rs<-genR(rs,S=S)
    
    Shat<-mean(Rs,type='geo')
    
    #rs2<-rs
    rs2<-rot.dist(Rs,Shat,method='intrinsic')
    a3i<-a3(rs2,terms=3)
    a4i<-a4(rs2,terms=3)
    
    ShatMedian<-rot.dist(Shat,S,method='intrinsic',p=2)
    
    tstats[i,j]<-n[j]*((a4i)^2)*ShatMedian/(a3i)
    
  } 

}

tstats2<-tstats*1

seqChi<-seq(0,max(tstats2),length=100)

plot(ecdf(tstats2[,1]),main=paste("kappa=",kap))
lines(ecdf(tstats2[,2]),col=2)
lines(seqChi,pchisq(seqChi,2),lty=2,col=3,lwd=3)
legend("bottomright",c(paste("n=",n[1]),paste("n=",n[2])),col=c(1,2),lwd=1)

