source("intervals/BayesRegions/BayesRegionsFunctions.R")


#############
#von Mises Distribution
kap<-1
Rs<-ruars(100,rvmises,kappa=kap)

#Table 1 in Bingham 2009 suggests rho=1000, sigma=1 when for sample n=100 and kappa=1
mcRes<-both_MCMC(Rs,mean(Rs),kappa0=1,rho=1000,sigma=1,burnin=100,B=500,gfun=gvmUARS)

Sres<-as.SO3(mcRes$S)
plot(Sres)
ares<-angle(Sres)
plot(ares,type='l',ylab='r')
plot(mcRes$kappa,type='l',ylab=expression(kappa))

#Make CR, see if id.SO3 is in there
Shat<-mean(Sres)
ds<-afun(Sres,Shat)
rad<-quantile(ds,.95)
afun(id.SO3,Shat)<rad

#############
#matrix Fisher Distribution
kap<-1
Rs<-ruars(100,rfisher,kappa=kap)

#Table 2 in Bingham 2010 suggests phi=1000, sigma=1 when for sample n=100 and kappa=1
mcRes<-both_MCMC(Rs,mean(Rs),kappa0=kap,rho=1000,sigma=1,burnin=50,B=100,gfun=gfUARS)

Sres<-as.SO3(mcRes$S)
plot(Sres)
ares<-angle(Sres)
plot(ares,type='l',ylab='r')
plot(mcRes$kappa,type='l',ylab=expression(kappa))

#Make CR, see if id.SO3 is in there
Shat<-mean(Sres)
ds<-afun(Sres,Shat)
rad<-quantile(ds,.95)
afun(id.SO3,Shat)<rad
