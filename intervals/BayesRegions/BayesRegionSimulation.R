source("intervals/BayesRegions/BayesRegionsFunctions.R")


kap<-1
Rs<-ruars(100,rvmises,kappa=kap)

#Table 1 in Bingham 2009 suggests rho=1000, sigma=1 when for sample n=100 and kappa=1
mcRes<-both_MCMC(Rs,mean(Rs),kappa0=1,rho=1000,sigma=1,50,100)

Sres<-as.SO3(mcRes$S)
plot(Sres)
ares<-angle(Sres)
plot(ares,type='l')
plot(mcRes$kappa,type='l')