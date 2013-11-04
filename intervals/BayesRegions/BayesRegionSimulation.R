source("intervals/BayesRegions/BayesRegionsFunctions.R")

#######################
####R-implementation
#######################
#############
#von Mises Distribution
kap<-1
Rs<-ruars(100,rvmises,kappa=kap)

#Table 1 in Bingham 2009 suggests rho=1000, sigma=1 when for sample n=100 and kappa=1
mcRes<-both_MCMC(Rs,mean(Rs),kappa0=1,rho=1000,sigma=1,burnin=1000,B=1000,gfun=gvmUARSC,rfun=rvmises)
mcRes$Sacc; mcRes$Kacc #Check acceptance rates

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
hist(ds);abline(v=c(rad,afun(id.SO3,Shat)),col=c(1,2))

#############
#matrix Fisher Distribution
kap<-1
Rs<-ruars(100,rfisher,kappa=kap)

#Table 2 in Bingham 2010 suggests phi=1000, sigma=1 when for sample n=100 and kappa=1
mcRes<-both_MCMC(Rs,mean(Rs),kappa0=kap,rho=1000,sigma=.5,burnin=1000,B=1000,gfun=gfUARSC,rfun=rcayley)
mcRes$Sacc; mcRes$Kacc #Check acceptance rates


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
hist(ds);abline(v=c(rad,afun(id.SO3,Shat)),col=c(1,2))

#############
#Cayley Distribution
kap<-1
Rs<-ruars(100,rcayley,kappa=kap)

#Table 2 in Bingham 2010 suggests phi=1000, sigma=1 when for sample n=100 and kappa=1
mcRes<-both_MCMC(Rs,mean(Rs),kappa0=kap,rho=150,sigma=.4,burnin=1000,B=5000,gfun=gcayUARSC,rfun=rcayley)
mcRes$Sacc; mcRes$Kacc #Check acceptance rates

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
hist(ds);abline(v=c(rad,afun(id.SO3,Shat)),col=c(1,2))


#######################
####C++-implementation
#######################
#############
#matrix Fisher Distribution
kap<-1
Rs<-ruars(100,rfisher,kappa=kap)

#Table 2 in Bingham 2010 suggests phi=1000, sigma=1 when for sample n=100 and kappa=1
mcRes<-both_MCMC_CPP(Rs,mean(Rs),kappa0=kap,rho=1000,sigma=.5,burnin=1000,B=1000,Cayley=FALSE)
mcRes$Sacc; mcRes$Kacc #Check acceptance rates


Sres<-as.SO3(mcRes$S)
plot(Sres)
ares<-angle(Sres)
plot(ares,type='l',ylab='r')
plot(mcRes$kappa,type='l',ylab=expression(kappa))

#Make CR, see if id.SO3 is in there
Shat<-mean(Sres)
ds<-afun_CPP(Sres,Shat)
rad<-quantile(ds,.95)
afun(id.SO3,Shat)<rad
hist(ds);abline(v=c(rad,afun(id.SO3,Shat)),col=c(1,2))

#############
#Cayley Distribution
kap<-1
Rs<-ruars(100,rcayley,kappa=kap)

#Table 2 in Bingham 2010 suggests phi=1000, sigma=1 when for sample n=100 and kappa=1
mcRes<-both_MCMC_CPP(Rs,mean(Rs),kappa0=kap,rho=2000,sigma=1,burnin=1000,B=5000,Cayley=TRUE)
mcRes$Sacc; mcRes$Kacc #Check acceptance rates

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
hist(ds);abline(v=c(rad,afun(id.SO3,Shat)),col=c(1,2))