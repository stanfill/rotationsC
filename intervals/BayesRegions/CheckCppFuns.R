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

gvmUARS(Rs,mean(Rs),1)
gvmUARSC(Rs,mean(Rs),1)
microbenchmark(gvmUARS(Rs,mean(Rs),1),gvmUARSC(Rs,mean(Rs),1))

gfUARS(Rs,mean(Rs),1)
gfUARSC(Rs,mean(Rs),1)
microbenchmark(gcayUARS(Rs,mean(Rs),1),gcayUARSC(Rs,mean(Rs),1))

gcayUARS(Rs,mean(Rs),10)
gcayUARSC(Rs,mean(Rs),10)
microbenchmark(gcayUARS(Rs,mean(Rs),1),gcayUARSC(Rs,mean(Rs),1))

S_MCMC_CPP(Rs,mean(Rs),1000,1)
S_MCMC(Rs,mean(Rs),1000,1,gcayUARSC,rcayley)
microbenchmark(S_MCMC(Rs,mean(Rs),1000,1,gcayUARSC,rcayley),S_MCMC_CPP(Rs,mean(Rs),1000,1))
