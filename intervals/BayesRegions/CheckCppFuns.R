library(rotations)
library(microbenchmark)

Rs<-ruars(20,rcayley)
cRs<-center(Rs,mean(Rs))
cRs2<-centerCpp(Rs,mean(Rs))

sum(abs(cRs-cRs2))
mean(as.SO3(cRs2))

microbenchmark(center(Rs,mean(Rs)),centerCpp(Rs,mean(Rs)))

rowSums(cRs[,c(1,5,9)])
gvmUARSC(Rs,mean(Rs),10.1)



gvmUARS(Rs,mean(Rs),1)
gvmUARSC(Rs,mean(Rs),1)
microbenchmark(gvmUARS(Rs,mean(Rs),1),gvmUARSC(Rs,mean(Rs),1))

gfUARS(Rs,mean(Rs),1)
gfUARSC(Rs,mean(Rs),1)
microbenchmark(gcayUARS(Rs,mean(Rs),1),gcayUARSC(Rs,mean(Rs),1))

gcayUARS(Rs,mean(Rs),10)
gcayUARSC(Rs,mean(Rs),10)
microbenchmark(gcayUARS(Rs,mean(Rs),1),gcayUARSC(Rs,mean(Rs),1))