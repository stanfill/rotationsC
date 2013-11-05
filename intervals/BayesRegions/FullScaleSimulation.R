library(rotations,lib='.')
source('BayesRegionsFunctions.R')
Rcpp::sourceCpp('CppBayesFunctions.cpp')

n<-c(10,20,50,100)
cayKap<-c(10,4,2)
fishKap<-c(3.17,1.71,1.15)

B<-1000
alp<-0.9
kap<-1

FishresDf<-data.frame(cover=rep(0,B),width=rep(0,B))
CayresDf<-data.frame(cover=rep(0,B),width=rep(0,B))

for(i in 1:B){

  Rs<-ruars(100,rfisher,kappa=kap)

  #Table 2 in Bingham 2010 suggests phi=1000, sigma=1 when for sample n=100 and kappa=1
  mcRes<-both_MCMC_CPP(Rs,mean(Rs),kappa0=kap,rho=400,sigma=.8,burnin=5000,B=5000,Cayley=FALSE)
  #mcRes$Saccept;mcRes$Kaccept
  Sres<-as.SO3(mcRes$S)
  Shat<-mean(Sres)
  ds<-afun_CPP(Sres,Shat)
  rad<-quantile(ds,alp)
  cover<-as.numeric(afun(id.SO3,Shat)<rad)
  FishresDf$cover[i]=cover
  FishresDf$width[i]=rad
  

  Rs<-ruars(100,rcayley,kappa=kap)
  #Table 2 in Bingham 2010 suggests phi=1000, sigma=1 when for sample n=100 and kappa=1
  mcRes<-both_MCMC_CPP(Rs,mean(Rs),kappa0=kap,rho=900,sigma=.5,burnin=5000,B=5000,Cayley=TRUE)
  #mcRes$Saccept;mcRes$Kaccept
  Sres<-as.SO3(mcRes$S)
  Shat<-mean(Sres)
  ds<-afun_CPP(Sres,Shat)
  rad<-quantile(ds,alp)
  cover<-as.numeric(afun(id.SO3,Shat)<rad)
  CayresDf$cover[i]=cover
  CayresDf$width[i]=rad
  
  if(B%%100==0){
    write.csv(FishresDf,"Results/FisherResults_n100_k1.csv")
    write.csv(CayresDf,"Results/CayleyResults_n100_k1.csv")
  }
  
}

write.csv(FishresDf,"Results/FisherResults_n100_k1.csv")
write.csv(CayresDf,"Results/CayleyResults_n100_k1.csv")