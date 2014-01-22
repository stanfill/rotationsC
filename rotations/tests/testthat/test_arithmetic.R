rs<-rvmises(20)
Rs<-genR(rs,space='SO3')
Qs<-as.Q4(Rs)
RsRs<-Rs+Rs
RsRs2<-matrix(NA,20,9)

rs2<-2*abs(rs)
for(i in 1:20){
  Rsi<-matrix(Rs[i,],3,3)
  RsRs2[i,]<-Rsi%*%Rsi
  if(rs2[i]>pi){
    rs2[i]<-2*pi-rs2[i]
  }
}
class(RsRs2)<-"SO3"
context("Arithmetic")

expect_equal(rs2,mis.angle(Rs+Rs))
expect_equal(RsRs,RsRs2)
expect_equal(sum(mean(Rs-mean(Rs))),3)

expect_equal(Qs+Qs,as.Q4(Rs+Rs))
