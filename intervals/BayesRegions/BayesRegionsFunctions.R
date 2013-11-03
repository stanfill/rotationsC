#My attempt to code Dr. Bingham's credible regions from her 2009 paper:
#Bayes One-Sample and One-Way Random Effects Analyses for 3-D Orientations with
#Application to Materials Science

gvmUARS<-function(RS,S,kappa){
  #Compute equation (6) of Bingham et al. 2009b
  
  #Rs - n-by-9 matrix with rows corresponding to sample points
  #S - central orientation
  #kappa - concentration parameter
  
  cRs <- center(Rs,S) #Each row of cRs is S'R
  trcRs <- rowSums(cRs[,c(1,5,9)]) #each row is tr(S'R)
  
  n1<-exp(kappa*sum(trcRs-1)/2)
  I0k<-besselI(kappa,0)
  I1k<-besselI(kappa,1)
  n2<-sqrt(I0k^2-I0k*I1k/kappa-I1k^2)
  
  d1<-I0k^(nrow(Rs)+1)
  d2<-prod(3-trcRs)
  
  return(n1*n2/(d1*d2))
  
}

gfUARS<-function(RS,S,kappa){
  #Compute equation (13) of Bingham et al. 2010
  
  #Rs - n-by-9 matrix with rows corresponding to sample points
  #S - central orientation
  #kappa - concentration parameter
  
  cRs <- center(Rs,S) #Each row of cRs is S'R
  trcRs <- rowSums(cRs[,c(1,5,9)]) #each row is tr(S'R)
  
  n1<-exp(kappa*sum(trcRs-1))
  I02k<-besselI(2*kappa,0)
  I12k<-besselI(2*kappa,1)
  
  n2<-sqrt(2*I02k^2/kappa-2*I02k*I12k/(kappa^2)+((1/(kappa^2))-(2/kappa))*I12k^2)
  
  d1<-(I02k-I12k)^(nrow(Rs)+1)
  
  return((n1*n2)/(d1))
  
}

S_MCMC<-function(Rs,oldS,rho,kappa,gfun){
  
  #Rs - the sample
  #oldS - the previous draw from distribution on S
  #rho - tuning parameter
  #kappa - concentration for likelihood
  #gfun - the g() function, either gvmUARS or gfUARS
  
  Sstar <- ruars(1,rvmises,S=oldS,kappa=rho)
  
  rj1 <- gfun(Rs,Sstar,kappa)/gfun(Rs,oldS,kappa)
  Wj1 <-rbinom(1,1,min(1,rj1,na.rm=T))
 
  if(Wj1==1){
    return(Sstar)
  }else
    return(oldS)
}

kap_MCMC<-function(Rs,oldKappa,sigma,S,gfun){
  
  #Rs - the sample
  #oldKappa - the previous draw from distribution for kappa
  #sigma - tuning parameter
  #kappa - central orientation for likelihood
  #gfun - the g() function, either gvmUARS or gfUARS
  
  kappaStar <- exp(rnorm(1,log(oldKappa),sigma))
  
  rj2 <- (kappaStar*gfun(Rs,S,kappaStar))/(oldKappa*gfun(Rs,S,oldKappa))
  
  Wj2 <- rbinom(1,1,min(1,rj2,na.rm=T))
  
  if(Wj2==1){
    return(kappaStar)
  }else{
    return(oldKappa)
  }
  
}


both_MCMC<-function(Rs,S0,kappa0,rho,sigma,burnin,B,gfun){
  
  Sdraws<-matrix(0,B,9)
  Kdraws<-rep(0,B)
  
  Snew<-S0
  Knew<-kappa0
  
  for(i in 1:(burnin+1)){
    Snew<-S_MCMC(Rs,Snew,rho,Knew,gfun)
    Knew<-kap_MCMC(Rs,Knew,sigma,Snew,gfun)
  }
  
  Sdraws[1,]<-Snew
  Kdraws[1]<-Knew
  
  for(j in 2:B){
    Sold<-as.SO3(matrix(Sdraws[(j-1),]))
    Sdraws[j,]<-S_MCMC(Rs,Sold,rho,Kdraws[j-1],gfun)
    
    Snew<-as.SO3(matrix(Sdraws[(j),]))
    Kdraws[j]<-kap_MCMC(Rs,Kdraws[j-1],sigma,Snew,gfun)
  }
  
  return(list(S=Sdraws,kappa=Kdraws))
  
}

#Compute the angles between each axis of R1 and R2
afun<-function(R1,R2){
  
  n<-length(R1)/9
  R1<-matrix(R1,n,9)
  R2<-matrix(R2,3,3)
  as<-rep(0,n)
  
  for(i in 1:n){
    Ri<-matrix(R1[i,],3,3)
    as[i]<-max(acos(diag(t(Ri)%*%R2)))
  }
  return(as)
  
}
