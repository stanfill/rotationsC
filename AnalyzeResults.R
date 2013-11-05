Cay<-read.csv("intervals/BayesRegions/Results/CayleyResults_n100_k1.csv")[,-1]

summary(Cay$Saccept)
summary(Cay$Kaccept)
sum(Cay$cover)/nrow(Cay)
summary(Cay$width)
hist(Cay$width,breaks=100)



Fish<-read.csv("intervals/BayesRegions/Results/FisherResults_n100_k1.csv")[,-1]
summary(Fish$Saccept)
summary(Fish$Kaccept)
sum(Fish$cover)/nrow(Fish)
summary(Fish$width)
hist(Fish$width,breaks=100)