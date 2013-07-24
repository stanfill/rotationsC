######################################################
######################################################
#Find coverage rates for the parametric and bootstrap
#confidence regions


library(plyr)
library(reshape2)
library(rotations2)
sourceCpp("intervals/ZhangMethod.cpp")  

#ZhangMethod.cpp contains the functions that will compute c/d and perform the zhang bootstrap

alp<-.1
critVal<-qchisq(1-alp,3)
n<-c(10,20,50,100)
nus<-c(.25,.5,.75)
B<-5  			#Number of samples to use to estimate CDF
Dist<-c('Cayley','matrix-Fisher')

simSize<-length(n)*length(nus)*length(Dist)

coverCompare<-data.frame(expand.grid(nus=nus,n=n,Dist=Dist),MeanStat=0,MeanCrit=0,MedianStat=0,MedianCrit=0)
coverCompare<-coverCompare[rep(seq_len(nrow(coverCompare)),each=B),]
dimS<-nrow(coverCompare)

#coverRate<-read.csv("Results/MedianResultsM300.csv")[,-1]

for(j in 1:dimS){
	
	if(coverCompare$Dist[j]=='Cayley'){
		
		rangle<-rcayley
		kappaj <- cayley_kappa(coverCompare$nu[j])
		
	}else{
		
		rangle<-rfisher
		kappaj <- fisher_kappa(coverCompare$nu[j])
		
	}
	
	rs<-rangle(coverCompare$n[j],kappa=kappaj)
		
	Rs<-genR(rs)
	Qs<-Q4(Rs)
	
	#Mean
	Qhat<-mean(Qs)
	
	cdHat<-cdfunsC(Qs,Qhat) #compute c and d hat using consistent estimators
	chat<-cdHat[1]
	dhat<-cdHat[2]
	
	
	hsqMean<-dist(Qhat,method='intrinsic',p=2)
	
	#compute the test statistic and save it under 'MeanStat'
	coverCompare$MeanStat[j]<-2*coverCompare$n[j]*dhat^2*hsqMean/chat
	
	#compute bootstrap critical value
	coverCompare$MeanCrit[j]<-as.numeric(quantile(zhangQ4(Qs,300),1-alp,na.rm=T))
	
	
	#Median
	Stilde<-median(Rs)
		
	cdTilde<-cdfunsCSO3(Rs,Stilde) #compute c and d tilde using consistent estimators
	ctilde<-cdTilde[1]
	dtilde<-cdTilde[2]
		
		
	hsqMedian<-dist(Stilde,method='intrinsic',p=2)
		
	#compute the test statistic and save it under 'MeanStat'
	coverCompare$MedianStat[j]<-2*coverCompare$n[j]*dtilde^2*hsqMedian/ctilde
	
	#compute bootstrap critical value
	coverCompare$MedianCrit[j]<-as.numeric(quantile(zhangMedianC(Rs,300),1-alp,na.rm=T))
		

	if(j%%B==0)
		write.csv(coverCompare,"Results/MeanMedianComparison.csv")
}

ccrit<-qchisq(1-alp,3)
coverCompare$MeanCoverC<-as.numeric(coverCompare$MeanStat<ccrit)
coverCompare$MeanCoverZ<-as.numeric(coverCompare$MeanStat<coverCompare$MeanCrit)

coverCompare$MedianCoverC<-as.numeric(coverCompare$MedianStat<ccrit)
coverCompare$MedianCoverZ<-as.numeric(coverCompare$MedianStat<coverCompare$MedianCrit)

compareRate<-ddply(coverCompare,.(Dist,nus,n),summarize,MeanCCover=100*sum(MeanCoverC)/B,MeanZCover=100*sum(MeanCoverZ)/B,
									 MedianCCover=100*sum(MedianCoverC)/B,MedianZCover=100*sum(MedianCoverZ)/B)

# cRateM<-melt(coverRate,id=c('Dist','nus','n'))
# colnames(cRateM)[4]<-'Method'
# 
# levels(cRateM$Method)<-c("NTH(C&R)","B(Z&N)")
# 
# levels(cRateM$Dist)<-c("Cayley","matrix~~Fisher")
# cRateM$nu<-factor(cRateM$nu,labels=c("nu == 0.25","nu == 0.5","nu == 0.75"))
# 
# 
# qplot(n,value,data=cRateM,colour=Method,group=Method,ylab='Coverage Rate (%)',xlab='Sample Size')+
# 	facet_grid(Dist~nu,labeller=label_parsed)+
# 	geom_hline(yintercept=(1-alp)*100,colour='gray50')+geom_line(lwd=I(1.25),alpha=I(.8))+
# 	scale_x_continuous(breaks=c(10,20,50,100))+theme_bw()+theme(panel.margin=unit(0.5,'lines'))
# 
# ggsave("C:/Users/stanfill/Dropbox/Thesis/Intervals/Figures/CoverRatesB1000Median.pdf",width=7,height=4.5)
#ggsave("CoverRatesB1000Median.pdf",width=5,height=4)

#library(xtable)
#xtable(coverRate)