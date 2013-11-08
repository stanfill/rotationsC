library(rotations)
library(microbenchmark)

Rs<-ruars(20,rcayley)
cRs<-center(Rs,mean(Rs))
cRs2<-centerCpp(Rs,mean(Rs))
microbenchmark(center(Rs,mean(Rs)),centerCpp(Rs,mean(Rs)))


sum(abs(matrix(cRs,20,9)-matrix(cRs2,20,9)))
mean(as.SO3(cRs2))

microbenchmark(center(Rs,mean(Rs)),centerCpp(Rs,mean(Rs)))

rowSums(cRs[,c(1,5,9)])

log(gvmUARS(Rs,mean(Rs),1))
gvmUARSC(Rs,mean(Rs),1)
microbenchmark(gvmUARS(Rs,mean(Rs),1),gvmUARSC(Rs,mean(Rs),1))

log(gfUARS(Rs,mean(Rs),2))
gfUARSC(Rs,mean(Rs),2)
microbenchmark(gcayUARS(Rs,mean(Rs),1),gcayUARSC(Rs,mean(Rs),1))

log(gcayUARS(Rs,mean(Rs),10))
gcayUARSC(Rs,mean(Rs),10)
microbenchmark(gcayUARS(Rs,mean(Rs),1),gcayUARSC(Rs,mean(Rs),1))

S_MCMC_CPP(Rs,mean(Rs),1000,1,TRUE)
S_MCMC(Rs,mean(Rs),1000,1,gcayUARSC,rcayley)
microbenchmark(S_MCMC(Rs,mean(Rs),1000,1,gcayUARSC,rcayley),S_MCMC_CPP(Rs,mean(Rs),1000,1))


kap_MCMC(Rs,1,1,mean(Rs),gcayUARSC)
kap_MCMC_CPP(Rs,1,1,mean(Rs),TRUE)
microbenchmark(kap_MCMC(Rs,1,1,mean(Rs),gcayUARSC),kap_MCMC_CPP(Rs,1,1,mean(Rs)))


afun(Rs,mean(Rs))
as.vector(afun_CPP(Rs,mean(Rs)))
microbenchmark(afun(Rs,mean(Rs)),afun_CPP(Rs,mean(Rs)))


kap<-1
Rs<-ruars(100,rcayley,kappa=kap)

#Table 2 in Bingham 2010 suggests phi=1000, sigma=1 when for sample n=100 and kappa=1
mcRes<-both_MCMC(Rs,mean(Rs),kappa0=kap,rho=150,sigma=.4,burnin=1000,B=5000,gfun=gcayUARSC,rfun=rcayley)
mcRes$Sacc; mcRes$Kacc #Check acceptance rates


mcResC<-both_MCMC_CPP(Rs,mean(Rs),kappa0=kap,rho=1000,sigma=.25,burnin=1000,B=5000,gfUARSC)
mcResC$Saccept; mcResC$Kaccept

SresC<-as.SO3(mcResC$S)
plot(SresC)
ares<-afun(SresC,mean(SresC))
plot(ares,type='l',ylab='r')
plot(mcResC$kappa,type='l',ylab=expression(kappa))

Rs<-ruars(20,rcayley)

res1<-both_MCMC_CPP_CPP( Rs, mean(Rs), 1, 1000, .4, 100, 1000, TRUE)
res1$Saccept;res1$Kaccept
Ss<-as.SO3(res1$Rs)
plot(Ss)


res2<-both_MCMC_CPP( Rs, mean(Rs), 1, 1000, .4, 100, 1000, TRUE)
res2$Sacc;res2$Kacc  

microbenchmark(both_MCMC_CPP_CPP( Rs, mean(Rs), 1, 1000, .4, 100, 100, TRUE),
               both_MCMC_CPP( Rs, mean(Rs), 1, 1000, .4, 100, 100, TRUE))




#Package version of bayesCR

Rs<-ruars(20,rcayley,kappa=2)
cr<-bayesCR(Rs,type='Cayley',S0=mean(Rs),kappa0=2,tuneS=22,tuneK=0.9,burn_in=100,m=500,alp=.9)
str(cr)