################
#Projected mean
a2<-function(rs){
  return( (1+2*mean(cos(rs)))/3 )
}

IFMean<-function(Rs,Ri){
  Shat<-mean(Rs)
  
  RR <- as.SO3(rbind(Rs,Ri))
  Shatj<-mean(RR)
  
  empIF<-rot.dist(Shatj,Shat,method='i')*nrow(Rs)
  #rs <- mis.angle(Rs-Shat)
  #ri <- mis.angle(Ri-Shat)
  rs <- mis.angle(Rs)
  ri<-mis.angle(Ri)
  theIF<-sin(ri)/a2(rs)
  
  return(list(emp=empIF,theo=theIF))
}


Rs<-ruars(100,rcayley,kappa=10)
Ri<-genR(pi/8)
IFs<-IFMean(Rs,Ri)

IFs$theo/IFs$emp

#Plot the IFs
ris<-seq(0.1,pi,length=100)
Rs<-ruars(50,rcayley,kappa=5)
ifDF<-data.frame(Emp=rep(0,length(ris)),Theory=0)

for(i in 1:length(ris)){
  Ri<-genR(ris[i])
  ifi<-IFMean(Rs,Ri)
  ifDF$Emp[i]<-ifi$emp
  ifDF$Theory[i]<-ifi$theo
}

plot(ris,ifDF$Emp,type='l')
lines(ris,ifDF$Theory,col=2)
legend("topright",c("Empirical","Theory"),col=c(1,2),lty=1)

#Quaternion derivation

IFQuat<-function(Qs,Q){
  Shat<-mean(Qs)
  ShatM<-matrix(Shat,ncol=1)
  QM<-matrix(Q,ncol=1)
  the<-sqrt(sum((ShatM-QM)^2))
  
  aQs<-as.Q4(rbind(Qs,Q))
  Shat2<-mean(aQs)
  emp<-rot.dist(Shat,Shat2,method='intrinsic')*nrow(Qs)
  
  return(list(emp=emp,the=the))
}

#Plot the IFs
ris<-seq(0.1,pi,length=100)
Rs<-ruars(50,rcayley,kappa=50,space='Q4')
ifDF<-data.frame(Emp=rep(0,length(ris)),The=0)

for(i in 1:length(ris)){
  Ri<-genR(ris[i],space='Q4')
  ifi<-IFQuat(Rs,Ri)

  ifDF$The[i]<-ifi$the
  ifDF$Emp[i]<-ifi$emp
}

plot(ris,ifDF$The,pch=19)
lines(ris,ifDF$Emp,col=2)

plot(ris,ifDF$Emp)
################
#Projected median
a2med<-function(rs){
  
  num <- 1+3*cos(rs)
  denom <- 12*sqrt(1-cos(rs)) 
  
  return( mean(num/denom) )
}

IFMedian<-function(Rs,Ri){
  Shat<-median(Rs)
  
  RR <- as.SO3(rbind(Rs,Ri))
  Shatj<-median(RR)
  
  empIF<-rot.dist(Shatj,Shat,method='i')*nrow(Rs)
  #rs <- mis.angle(Rs-Shat)
  #ri <- mis.angle(Ri-Shat)
  rs<-mis.angle(Rs)
  ri<-mis.angle(Ri)
  
  theIF<-sin(ri)/(2*a2med(rs)*sqrt(1-cos(ri)))
  
  return(list(emp=empIF,theo=theIF))
}


Rs<-ruars(100,rcayley,kappa=10)
Ri<-genR(pi/2)
IFsMed<-IFMedian(Rs,Ri)
IFsMed$theo/IFsMed$emp

#Plot the IFs
ris<-seq(0.1,pi,length=100)
Rs<-ruars(50,rcayley,kappa=10)
ifDFMed<-data.frame(Emp=rep(0,length(ris)),Theory=0)

for(i in 1:length(ris)){
  Ri<-genR(ris[i])
  ifi<-IFMedian(Rs,Ri)
  ifDFMed$Emp[i]<-ifi$emp
  ifDFMed$Theory[i]<-ifi$theo
}

plot(ris,ifDFMed$Emp,type='l')
lines(ris,ifDFMed$Theory,col=2)
legend("topright",c("Empirical","Theory"),col=c(1,2),lty=1)

################
#Geometric mean
a4<-function(r,terms=2){
  
  cosr<-cos(r)
  cos2<-cos(r)^2
  
  if(terms==2){
    
    #a44<-(4-10*mean(cosr)-3*mean(cos2))/3    #Old calculation
    
    a44<-(8+14*mean(cosr)-4*mean(cos2))/9
    
  }else if(terms==3){
    
    cos3<-cos(r)^3
    
    p1<-140*mean(cosr)
    p2<--64*mean(cos2)
    p3<-16*mean(cos3)
    a44<-(88+p1+p2+p3)/90

  }else if(terms==4){
    
    cos3<-cos(r)^3
    cos4<-cos(r)^4
    
    p1 <- (-279/140+47/84-123/1260+4/560)*(-2/3)*(1+2*mean(cosr))
    p2 <- (94/84-492/1260+24/560)*(-2/3)*(mean(cosr)+2*mean(cos2))
    p3 <- (-492/1260+48/560)*(-2/3)*(mean(cos2)+2*mean(cos3))
    p4 <- (32/560)*(-2/3)*(mean(cos3)+2*mean(cos4))
    a44 <- p1+p2+p3+p4
    
  }
  return(a44)
}

IFGeoMean<-function(Rs,Ri){
  Shat<-mean(Rs,type='g')
  
  RR <- as.SO3(rbind(Rs,Ri))
  Shatj<-mean(RR,type='g')
  
  empIF<-rot.dist(Shatj,Shat,method='i')*nrow(Rs)
  #rs <- mis.angle(Rs-Shat)
  #ri <- mis.angle(Ri-Shat)
  rs <- mis.angle(Rs)
  ri <- mis.angle(Ri)
  
  theIF<-2*ri/(a4(rs,terms=3)) #The IF was derived using a different f(h), the factor of 4 corrects this
  
  return(list(emp=empIF,theo=theIF))
}


Rs<-ruars(50,rcayley,kappa=100)
rs<-mis.angle(Rs-mean(Rs))
Ri<-genR(pi/2)
GeoIFs<-IFGeoMean(Rs,Ri)

GeoIFs$theo/GeoIFs$emp


#Plot the IFs
ris<-seq(0.1,pi,length=100)
Rs<-ruars(250,rvmises,kappa=500)
ifDF<-data.frame(Emp=rep(0,length(ris)),Theory=0)

for(i in 1:length(ris)){
  Ri<-genR(ris[i])
  ifi<-IFGeoMean(Rs,Ri)
  ifDF$Emp[i]<-ifi$emp
  ifDF$Theory[i]<-ifi$theo
}

plot(ris,ifDF$Emp,type='l')
lines(ris,ifDF$Theory,col=2)
legend("topleft",c("Empirical","Theory"),col=c(1,2),lty=1)


################
#Geometric median

a5 <- function(r){
  sinr<-sin(r)^2
  return(16*mean(sinr)/(3))
}

IFGeoMedian<-function(Rs,Ri){
  Shat<-median(Rs,type='g')
  
  RR <- as.SO3(rbind(Rs,Ri))
  Shatj<-median(RR,type='g')
  
  empIF<-rot.dist(Shatj,Shat,method='i')*nrow(Rs)
  
  rs <- mis.angle(Rs)

  theIF <- a5(rs)
  
  return(list(emp=empIF,theo=theIF))
}

#Plot the IFs
ris<-seq(0.1,pi,length=100)
ifDF<-data.frame(Emp=rep(0,length(ris)),Theory=0)

Rs<-ruars(250,rcayley,kappa=50)

for(i in 1:length(ris)){
  Ri<-genR(ris[i])
  ifi<-IFGeoMedian(Rs,Ri)
  ifDF$Emp[i]<-ifi$emp
  ifDF$Theory[i]<-ifi$theo
}

sc<-(2)

plot(ris,ifDF$Emp*sc,type='l',ylim=c(0,max(c(ifDF$Theory,sc*ifDF$Emp))))
lines(ris,ifDF$Theory)
