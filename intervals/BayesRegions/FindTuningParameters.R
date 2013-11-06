library(rotations)
Rcpp::sourceCpp('intervals/BayesRegions/CppBayesFunctions.cpp')

f<-rfisher

kappa<-3.17
n<-100

#To decrease % for S, decrease rho
#To decrease % for kappa, increase sigma

Rs<-ruars(n,f,kappa=kappa)
rho<-725 ; sigma<- .2
mcRes<-both_MCMC_CPP(Rs,mean(Rs),kappa0=kappa,rho=rho,sigma=sigma,burnin=5000,B=5000,Cayley=FALSE)
mcRes$Saccept;mcRes$Kaccept