library(rotations)
Qs<-ruars(1,rcayley,space='Q4')
y<-ruars(1,rcayley,space='Q4')

y%*%t(Qs)%*%Qs%*%t(y)

sum(diag(t(Qs)%*%y))^2


rs<-rcayley(20)
XXp<-diag(c(mean(cos(rs)^2),mean(sin(rs)^2)*c(1,1,1)/3)) 
eigen(XXp)

rst<-seq(-pi,pi,length=100)
plot(rst,cos(rst)^2,type='l')
abline(h=0,v=0)
lines(rst,sin(rst)^2/3)

##############
#
library(rotations)

Qs<-ruars(20,rcayley,kappa=10,space='Q4')
y<-ruars(1,rcayley,space='Q4')

QuatIF<-function(Qs,y){

  Shat<-mean(Qs)
  Shat2<-mean(as.Q4(rbind(Qs,y)))

  emp <- rot.dist(Shat,Shat2,method='i')*nrow(Qs)

  cQs<-Qs-Shat

  denom<-eigen((t(cQs)%*%cQs)/nrow(cQs))$values[1]

  m<-matrix(Shat)
  
  #if(rot.dist(Shat,y,method='i')<2){
  #  m<-matrix(c(0,m[2],0,0))
  #}else{
  #  m<-matrix(c(m[1],0,0,0))
  #}
  
  num<-2*t(y)%*%y%*%m
  theory<-num/(denom)
  #theory<-sqrt(sum(theory^2))
  
  return(list(Emp=emp,The=theory))
}

QuatIF(Qs,y)

########
Qs<-ruars(20,rcayley,kappa=100,space='Q4')
emp<-rep(0,50)
the<-matrix(0,50,4)
rs<-seq(0,pi,length=50)

for(i in 1:length(rs)){
  
  y<-genR(rs[i],space='Q4')
  QIF<-QuatIF(Qs,y)
  emp[i]<-QIF$Emp
  the[i,]<-QIF$The
}

plot(rs,sqrt(rowSums(the^2)))
lines(rs,emp)
abline(v=2*asin(sqrt(3/4)))

plot(rs,cos(rs/2)^2,type='l',col=2)
lines(rs,sin(rs/2)^2/3)
abline(v=2*asin(sqrt(3/4)))


