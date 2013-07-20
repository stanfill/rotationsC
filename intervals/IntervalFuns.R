formatQ4<-function(Qs){
	
	if(length(Qs)%%4!=0)
		stop("Data needs to have length divisible by 4.")
	
	Qs<-matrix(Qs,length(Qs)/4,4)
	
	if (!all(apply(Qs, 1, is.Q4))) 
		warning("At least one of the given observations is not a unit quaternion.  Use result with caution.")
	
	
	if(length(Qs)==4)
		return(as.Q4(Qs))
	else
		return(as.Q4(Qs))
}


formatSO3<-function(Rs){
	len<-length(Rs)
	if(len%%9!=0)
		stop("Data needs to have length divisible by 9.")
	Rs<-matrix(Rs,len/9,9)
	if (!all(apply(Rs, 1, is.SO3))) 
		warning("At least one of the given observations is not in SO(3).  Use result with caution.")
	return(as.SO3(Rs))
}

log.SO3 <- function(R) {
	
	#R<-formatSO3(R)
	
	theta <- angle(R)  
	R<-matrix(R,3,3)
	
	if (abs(cos(theta)) >= 1) {
		return(diag(0, 3, 3))
	} else {
		
		c <- theta/(2 * sin(theta))
		mlog <- c * (R - t(R))
		return(mlog)
	}
}



ecdf<-function(rs){
	
	n<-length(rs)
	cdf<-rep(0,n)
	for(i in 1:n){
		cdf[i]<-length(which(rs<rs[i]))/n
	}
	
	return(cdf)
}

#These two functions do the exact same thing but cdfuns is faster than cdfunsAng
cdfuns<-function(Rs,Shat){
	Rs<-formatSO3(Rs)
	n<-nrow(Rs)
	c<-0
	d<-0

	StO<-centeringSO3(Rs,Shat)
	
	for(i in 1:n){
		StOi<-matrix(StO[i,],3,3)
		c<-c+(3-sum(diag(StOi%*%StOi)))
		d<-d+sum(diag(StOi))
	}
	
	c<-c/(6*n)
	d<-d/(3*n)
	return(list(c=c,d=d))
}

#Slower than the above equivilant function 'cdfuns'
cdfunsAng<-function(Rs,Shat){
	Rs<-formatSO3(Rs)
	n<-nrow(Rs)
	
	StO<-centeringSO3(Rs,Shat)
	StOAng<-angle(StO)
	
	c<-2*(1-mean(cos(StOAng)^2))/3
	d<-(1+2*mean(cos(StOAng)))/3
	return(list(c=c,d=d))
}

ZhangCI<-function(Rs,m,alpha,estimator='mean'){
	
	#Rs is a n-by-9 matrix where each row is an 3-by-3 rotation matrix
	#m is the number of resamples to find q_1-alpha
	#alpha is the level of confidence desired, e.g. 0.95 or 0.90
	#pivot logical; should the pivotal (T) bootstrap be used or nonpivotal (F)
	
	Rs<-formatSO3(Rs)
	n<-nrow(Rs)
	
	if(estimator=='mean'){
		Shat<-mean(Rs)
	}else{
		Shat<-median(Rs,type='intrinsic')
	}
	
	tstar<-rep(0,m)
	tstarPivot<-rep(0,m)
	
	for(i in 1:m){
		Ostar<-as.SO3(Rs[sample(1:n,replace=T),])
		
		if(estimator=='mean'){
			ShatStar<-mean(Ostar)
		}else{
			ShatStar<-median(Ostar,type='intrinsic')
		}
		
		#tstar[i]<-3-sum(diag(t(Shat)%*%ShatStar))      			#Approximation
		#print(tstar[i])
		tstar[i]<-dist(Shat,ShatStar,method='intrinsic')^2		#Exact
		#print(tstar[i])
		cd<-cdfuns(Ostar,ShatStar)
		
		tstarPivot[i]<-tstar[i]*2*n*cd$d^2/cd$c
		
	}
	
	return(list(ts=as.numeric(quantile(tstar,alpha)),tsPivot=as.numeric(quantile(tstarPivot,alpha))))
}

formatQ4<-function(Qs){
	
	if(length(Qs)%%4!=0)
		stop("Data needs to have length divisible by 4.")
	
	Qs<-matrix(Qs,length(Qs)/4,4)
	
	if (!all(apply(Qs, 1, is.Q4))) 
		warning("At least one of the given observations is not a unit quaternion.  Use result with caution.")
	
	
	if(length(Qs)==4)
		return(as.Q4(Qs))
	else
		return(as.Q4(Qs))
}

qMult<-function(q1,q2){
	#Forms quaternion product q1 x q2, i.e., rotate q2 by q1
	q1<-formatQ4(q1)
	q2<-formatQ4(q2)
	t0<-q2%*%matrix(c(q1[1],-q1[2:4]),4,1)
	t1<-q2%*%matrix(c(q1[2:1],-q1[4],q1[3]),4,1)
	t2<-q2%*%matrix(c(q1[c(3,4,1)],-q1[2]),4,1)
	t3<-q2%*%matrix(c(q1[4],-q1[3],q1[2:1]),4,1)
	return(formatQ4(c(t0,t1,t2,t3)))
}

pMat<-function(p){
	Pmat<-matrix(0,4,4)
	Pmat[,1]<-p
	Pmat[,2]<-p[c(2,1,4,3)]*c(-1,1,1,-1)
	Pmat[,3]<-c(-p[3:4],p[1:2])
	Pmat[,4]<-p[4:1]*c(-1,1,-1,1)
	return(Pmat)
}

RivestCI2<-function(qs,S=id.Q4){
  #This takes as input the dataset and true central direction S
  n<-nrow(qs)
  Shat<-meanC(qs)
  Phat<-pMat(Shat)
  
  Rhat<-qs%*%Phat
  resids<-matrix(0,n,3)
  VarShat<-matrix(0,3,3)
  
  resids<-2*Rhat[,1]*matrix(Rhat[,2:4],n,3)
  
  VarShat<-t(resids)%*%resids/(n-1)
  
  RtR<-t(Rhat)%*%Rhat
  Ahat<-(diag(RtR[1,1],3,3)-RtR[-1,-1])/n
  
  #St<-as.Q4(matrix(c(S[1],-S[2:4]),1,4))
  #StShat<-qMult(St,Shat)
  StShat<-as.Q4(matrix(c(Shat[1],-Shat[2:4]),1,4))
  avec<-matrix(axis2(StShat)*angle(StShat),1,3)
  
  Tm<-((n-3)/(3*n-3))*avec%*%Ahat%*%solve(VarShat)%*%Ahat%*%t(avec)
  
  return(Tm)
}

RivestCI<-function(qs,S=id.Q4){
	#This takes as input the dataset and true central direction S
	n<-nrow(qs)
	Shat<-meanQ4C(qs)
	Phat<-pMat(Shat)
	
	Rhat<-qs%*%Phat
	resids<-matrix(0,n,3)
	VarShat<-matrix(0,3,3)
	
	resids<-2*Rhat[,1]*matrix(Rhat[,2:4],n,3)
	
	VarShat<-t(resids)%*%resids/(n-1)
	
	RtR<-t(Rhat)%*%Rhat
	Ahat<-(diag(RtR[1,1],3,3)-RtR[-1,-1])/n
	
	#St<-as.Q4(matrix(c(S[1],-S[2:4]),1,4))
	#StShat<-qMult(St,Shat)
	StShat<-as.Q4(matrix(c(Shat[1],-Shat[2:4]),1,4))
	avec<-matrix(axis2(StShat)*angle(StShat),1,3)
	
	Tm<-n*avec%*%Ahat%*%solve(VarShat)%*%Ahat%*%t(avec)
	
	return(Tm)
}

proj<-function(u,v){
	#Project the vector v orthogonally onto the line spanned by the vector u
	num<-t(u)%*%v
	denom<-t(u)%*%u
	return(num*u/denom)
}

mean2<-function(Qs){
	mhat<-colSums(Qs)
	mhat<-mhat/sqrt(t(mhat)%*%mhat)
	return(matrix(mhat,1,4))
}

fisherBootCI<-function(Qs,mhat){
	n<-nrow(Qs)
	
	mvec<-mean2(Qs)
	Mmat<-diag(4)
	Mmat[1,]<-mvec
	
	#Find the orthogoanl matrix Mhat(d)
	Mmat[2,]<-Mmat[2,]-proj(Mmat[1,],Mmat[2,])
	Mmat[3,]<-Mmat[3,]-proj(Mmat[1,],Mmat[3,])-proj(Mmat[2,],Mmat[3,])
	Mmat[4,]<-Mmat[4,]-proj(Mmat[1,],Mmat[4,])-proj(Mmat[2,],Mmat[4,])-proj(Mmat[3,],Mmat[4,])
	
	Mmat[2,]<-Mmat[2,]/sqrt(t(Mmat[2,])%*%Mmat[2,])
	Mmat[3,]<-Mmat[3,]/sqrt(t(Mmat[3,])%*%Mmat[3,])
	Mmat[4,]<-Mmat[4,]/sqrt(t(Mmat[4,])%*%Mmat[4,])
	Mmat<-Mmat[-1,]
	
	#Construct G matrix, a consistent estimator of the mean
	xbarlen<-colMeans(Qs)
	xbarlenSq<-t(xbarlen)%*%xbarlen
	
	G<-Mmat%*%t(Qs)%*%Qs%*%t(Mmat)/(xbarlenSq[1,1]*n)

	
	Tm<-n*mhat%*%t(Mmat)%*%solve(G)%*%Mmat%*%t(mhat)
	return(Tm)
}

fisherAxisBoot<-function(Qs,m){
	
	n<-nrow(Qs)
	Tm<-rep(0,m)
	QsHat<-mean(Qs)
	
	for(i in 1:m){
		
		ni<-sample(n,replace=T)

		#If less than 4 unique obs from the data set are 
		#sampled then Q'Q will be rank deficient and G wong
		#have an inverse

		while(length(unique(ni))<4){
			ni<-sample(n,replace=T)
		}

		Qsi<-Qs[ni,]
		Tm[i]<-fisherAxisCompute(Qsi,QsHat)
		
	}
	
	return(Tm)
}

fisherAxisCompute<-function(Qs,Shat){
	
	n<-nrow(Qs)
	
	svdQs<-svd(t(Qs)%*%Qs/n)
	mhat<-svdQs$v[,1]
	Mhat<-t(svdQs$v[,-1])
	etad<-svdQs$d[1]
	etas<-svdQs$d[-1]
	G<-matrix(0,3,3)
		
	for(j in 1:3){
		for(k in j:3){
			denom<-1/(n*(etad-etas[j])*(etad-etas[k]))
			#print(denom)
			for(i in 1:n){
				G[j,k]<-G[k,j]<-G[j,k]+(Mhat[j,]%*%Qs[i,])*(Mhat[k,]%*%Qs[i,])*(mhat%*%Qs[i,])^2*denom
			}
		}
	}
	#print("etas:");print(etas)
	#print("G:");print(G)
	#print("Qhat:");print(Shat)
	#print("Mhat:");print(t(Mhat))
	
	Ginv<-solve(G)
	
	#Ginv<-try(solve(G),silent=T)
		
	#if(class(Ginv)!='matrix'){
	#	Ginv<-diag(1/diag(G))
	#}
		
	Tm<-n*Shat%*%t(Mhat)%*%Ginv%*%Mhat%*%t(Shat)
	
	return(Tm)
}

centeringQ4<-function(Qs,S){
	#This takes a set of observations in Q4 and centers them around S
	Qs<-formatQ4(Qs)
	S<-formatQ4(S)
	S[2:4]<--S[2:4]
	
	for(i in 1:nrow(Qs)){
		Qs[i,]<-qMult(S,Qs[i,])
	}
	
	return(Qs)
}

centeringSO3<-function(Rs,S){
	#This takes a set of observations in SO3 and centers them around S
	
	Rs<-formatSO3(Rs)
	S<-matrix(formatSO3(S),3,3)
	
	for(i in 1:nrow(Rs)){
		Rs[i,]<-t(S)%*%matrix(Rs[i,],3,3)
	}
	return(as.SO3(Rs))
}

zhangMedian<-function(Rs,alpha,m=300){
  
  n<-nrow(Rs)
  Shat<-median(Rs)
  hstar<-rep(0,m)
  
  for(i in 1:m){
    nstar<-sample(n,replace=T)
    Rstar<-as.SO3(Rs[nstar,])
    Sstar<-median(Rstar)
    cd<-cdfunsCSO3(Rstar,Sstar)
    hsqMean<-dist(Shat,Sstar,method='intrinsic',p=2)
    
    hstar[i]<-2*n*cd[2]^2*hsqMean/cd[1]
  
  }
  
  return(as.numeric(quantile(hstar,1-alpha,na.rm=T)))
}