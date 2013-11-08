library(rotations,lib='.')
#source('BayesRegionsFunctions.R')
Rcpp::sourceCpp('CppBayesFunctions.cpp')

n<-rep(c(10,20,50,100),3)
cayKap<-c(2,4,10)
fishKap<-c(1.15,1.71,3.17)
B<-1000
Cay_rho<-rep(c(22,35,105,200,39,75,210,350,95,175,450,850),each=B)
Cay_sigma<-rep(c(.9,.75,.42,.3,.78,.6,.38,.28,.8,.55,.35,.225),each=B)

fish_rho<-rep(c(20,45,100,200,28,45,175,300,60,140,425,725),each=B)
fish_sigma<-rep(c(.65,.45,.3,.22,.65,.45,.28,.2,.65,.45,.28,.2),each=B)

alp<-0.9

FishresDf<-data.frame(Kappa=rep(fishKap,each=4*B),n=rep(n,each=B),cover=0,width=0,Saccept=0,Kaccept=0)
CayresDf<-data.frame(Kappa=rep(cayKap,each=4*B),n=rep(n,each=B),cover=0,width=0,Saccept=0,Kaccept=0)
id<-matrix(id.SO3,1,9)

for(i in 1:nrow(FishresDf)){

  #Fisher simulations 
  
  Rs<-ruars(FishresDf$n[i],rfisher,kappa=FishresDf$Kappa[i])
  mcRes<-both_MCMC_CPP(Rs,mean(Rs),kappa0=FishresDf$Kappa[i],rho=fish_rho[i],sigma=fish_sigma[i],burnin=5000,B=5000,f=gfUARSC)
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
  mcRes<-both_MCMC_CPP(Rs,mean(Rs),kappa0=CayresDf$Kappa[i],rho=Cay_rho[i],sigma=Cay_sigma[i],burnin=5000,B=5000,f=gcayUARSC)
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
    write.csv(FishresDf,"Results/FisherResults.csv")
    write.csv(CayresDf,"Results/CayleyResults.csv")
  }
  
}

write.csv(FishresDf,"Results/FisherResults.csv")
write.csv(CayresDf,"Results/CayleyResults.csv")