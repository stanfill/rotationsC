#My attempt to code Dr. Bingham's credible regions from her 2009 paper:
#Bayes One-Sample and One-Way Random Effects Analyses for 3-D Orientations with
#Application to Materials Science

gvmUARS<-function(RS,S,kappa){
  
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

S_MCMC<-function(Rs,oldS,rho,kappa){
  
  #Rs - the sample
  #oldS - the previous draw from distribution on S
  #rho - tuning parameter
  #kappa - concentration for likelihood
  
  Sstar <- ruars(1,rvmises,S=oldS,kappa=rho)
  
  rj1 <- gvmUARS(Rs,Sstar,kappa)/gvmUARS(Rs,oldS,kappa)
 
  Wj1 <-rbinom(1,1,min(1,rj1))
 
  if(Wj1==1){
    return(Sstar)
  }else
    return(oldS)
}

kap_MCMC<-function(Rs,oldKappa,sigma,S){
  
  #Rs - the sample
  #oldKappa - the previous draw from distribution for kappa
  #sigma - tuning parameter
  #kappa - central orientation for likelihood
  
  kappaStar <- exp(rnorm(1,log(oldKappa),sigma))
  
  rj2 <- gvmUARS(Rs,S,kappaStar)/gvmUARS(Rs,S,oldKappa)
  
  Wj2 <- rbinom(1,1,min(1,rj2))
  
  if(Wj2==1){
    return(kappaStar)
  }else{
    return(oldKappa)
  }
  
}


both_MCMC<-function(Rs,S0,kappa0,rho,sigma,burnin,B){
  
  Sdraws<-matrix(0,B,9)
  Kdraws<-rep(0,B)
  
  Snew<-S0
  Knew<-kappa0
  
  for(i in 1:(burnin+1)){
    Snew<-S_MCMC(Rs,Snew,rho,Knew)
    Knew<-kap_MCMC(Rs,Knew,sigma,Snew)
  }
  
  Sdraws[1,]<-Snew
  Kdraws[1]<-Knew
  
  for(j in 2:B){
    Sold<-as.SO3(matrix(Sdraws[(j-1),]))
    Sdraws[j,]<-S_MCMC(Rs,Sold,rho,Kdraws[j-1])
    
    Snew<-as.SO3(matrix(Sdraws[(j),]))
    Kdraws[j]<-kap_MCMC(Rs,Kdraws[j-1],sigma,Snew)
  }
  
  return(list(S=Sdraws,kappa=Kdraws))
  
}
