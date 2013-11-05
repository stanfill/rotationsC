library(rotations,lib='.')
#source('BayesRegionsFunctions.R')
Rcpp::sourceCpp('CppBayesFunctions.cpp')

n<-c(10,20,50,100)
#cayKap<-c(10,4,2)
#fishKap<-c(3.17,1.71,1.15)
B<-1000
Cay_rho<-rep(c(39,75,210,350),each=B)
Cay_sigma<-rep(c(.78,.6,.38,.28),each=B)

fish_rho<-rep(c(28,45,175,300),each=B)
fish_sigma<-rep(c(.65,.45,.28,.2),each=B)

alp<-0.9


FishresDf<-data.frame(Kappa=rep(1.71,4*B),n=rep(n,each=B),cover=rep(0,4*B),width=rep(0,4*B),Saccept=rep(0,4*B),Kaccept=rep(0,4*B))
CayresDf<-data.frame(Kappa=rep(4,4*B),n=rep(n,each=B),cover=rep(0,4*B),width=rep(0,4*B),Saccept=rep(0,4*B),Kaccept=rep(0,4*B))
id<-matrix(id.SO3,1,9)

for(i in 1:(4*B)){

  #Fisher simulations 
  
  Rs<-ruars(FishresDf$n[i],rfisher,kappa=FishresDf$Kappa[i])
  mcRes<-both_MCMC_CPP(Rs,mean(Rs),kappa0=FishresDf$Kappa[i],rho=fish_rho[i],sigma=fish_sigma[i],burnin=5000,B=5000,Cayley=FALSE)
  #mcRes$Saccept;mcRes$Kaccept
  
  Sres<-as.SO3(mcRes$S)
  Shat<-mean(Sres)
  ds<-afun_CPP(Sres,Shat)
  rad<-quantile(ds,alp)
  cover<-as.numeric(afun_CPP(id,Shat)<rad)
  FishresDf$cover[i]=cover
  FishresDf$width[i]=rad
  FishresDf$Saccept[i]=mcRes$Saccept
  FishresDf$Kaccept[i]=mcRes$Kaccept

  
  #Cayley Simulations
  
  Rs<-ruars(CayresDf$n[i],rcayley,kappa=CayresDf$Kappa[i])
  mcRes<-both_MCMC_CPP(Rs,mean(Rs),kappa0=CayresDf$Kappa[i],rho=Cay_rho[i],sigma=Cay_sigma[i],burnin=5000,B=5000,Cayley=TRUE)
  #mcRes$Saccept;mcRes$Kaccept
  
  Sres<-as.SO3(mcRes$S)
  Shat<-mean(Sres)
  ds<-afun_CPP(Sres,Shat)
  rad<-quantile(ds,alp)
  cover<-as.numeric(afun_CPP(id,Shat)<rad)
  CayresDf$cover[i]=cover
  CayresDf$width[i]=rad
  CayresDf$Saccept[i]=mcRes$Saccept
  CayresDf$Kaccept[i]=mcRes$Kaccept
  
  if(B%%1000==0){
    write.csv(FishresDf,"Results/FisherResults_nu25.csv")
    write.csv(CayresDf,"Results/CayleyResults_nu25.csv")
  }
  
}

write.csv(FishresDf,"Results/FisherResults_nu25.csv")
write.csv(CayresDf,"Results/CayleyResults_nu25.csv")