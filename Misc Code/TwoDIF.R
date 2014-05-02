r<-seq(-pi,pi,length=100)
u1<-seq(-1,1,length=100)

IFmean<-function(r,u1){
  return(sin(r)*u1)
}

z<-outer(r,u1,IFmean)

persp(r,u1,z,theta=45,phi=60)
