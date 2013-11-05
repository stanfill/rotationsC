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

gcayUARS<-function(Rs,S,kappa){
  
  n<-nrow(Rs)
  cRs <- center(Rs,S) #Each row of cRs is S'R
  trcRs <- rowSums(cRs[,c(1,5,9)]) #each row is tr(S'R)
  
  p1<-(sqrt(pi)*gamma(kappa+2)/gamma(kappa+.5))^n
  p2<-sqrt(trigamma(kappa+.5)-trigamma(kappa+2))
  p3<-prod((.5+.25*(trcRs-1))^kappa)
  
  return(p1*p2*p3)
  
}

S_MCMC<-function(Rs,oldS,rho,kappa,gfun,rfun){
  
  #Rs - the sample
  #oldS - the previous draw from distribution on S
  #rho - tuning parameter
  #kappa - concentration for likelihood
  #gfun - the g() function, either gvmUARS or gfUARS
  
  Sstar <- matrix(ruars(n=1,rangle=rfun,S=oldS,kappa=rho),3,3)
  oldS<-matrix(oldS,3,3)
  
  rj1 <- gfun(Rs,Sstar,kappa)/gfun(Rs,oldS,kappa)
  
  if(is.nan(rj1)) rj1<-0  #Temporary fix until I think of something better
  
  
  Wj1 <-rbinom(1,1,min(1,rj1))
 
  if(Wj1==1){
    return(list(S=Sstar,Acc=1))
  }else
    return(list(S=oldS,Acc=0))
}

kap_MCMC<-function(Rs,oldKappa,sigma,S,gfun){
  
  #Rs - the sample
  #oldKappa - the previous draw from distribution for kappa
  #sigma - tuning parameter
  #kappa - central orientation for likelihood
  #gfun - the g() function, either gvmUARS or gfUARS
  
  kappaStar <- exp(rnorm(1,log(oldKappa),sigma))
  S<-matrix(S,3,3)
  rj2 <- (kappaStar*gfun(Rs,S,kappaStar))/(oldKappa*gfun(Rs,S,oldKappa))
  
  if(is.nan(rj2)) rj2<-0  #Temporary fix until I think of something better
  
  Wj2 <- rbinom(1,1,min(1,rj2))
  
  if(Wj2==1){
    return(list(kappa=kappaStar,Acc=1))
  }else{
    return(list(kappa=oldKappa,Acc=0))
  }
  
}


both_MCMC<-function(Rs,S0,kappa0,rho,sigma,burnin,B,gfun,rfun){
  
  Sdraws<-matrix(0,B,9)
  Kdraws<-rep(0,B)
  
  Snew<-list(S=S0,Acc=0)
  Knew<-list(kappa=kappa0,Acc=0)
  
  for(i in 1:(burnin+1)){
    Snew<-S_MCMC(Rs,Snew$S,rho,Knew$kappa,gfun,rfun)
    Knew<-kap_MCMC(Rs,Knew$kappa,sigma,Snew$S,gfun)
  }
  
  Sdraws[1,]<-Snew$S
  Kdraws[1]<-Knew$kappa
  Saccept<-1
  Kaccept<-1
  
  for(j in 2:B){
    Sold<-as.SO3(matrix(Sdraws[(j-1),]))
    Sd<-S_MCMC(Rs,Sold,rho,Kdraws[j-1],gfun,rfun)
    Sdraws[j,]<-Sd$S
    Saccept<-Saccept+Sd$Acc
    
    Snew<-as.SO3(matrix(Sdraws[(j),]))
    Kd<-kap_MCMC(Rs,Kdraws[j-1],sigma,Snew,gfun)
    Kdraws[j]<-Kd$kappa
    Kaccept<-Kaccept+Kd$Acc
  }
  
  return(list(S=Sdraws,Sacc=Saccept/B,kappa=Kdraws,Kacc=Kaccept/B))
  
}

# both_MCMC_CPP<-function(Rs,S0,kappa0,rho,sigma,burnin,B,Cayley){
#   #Valid for Cayley and matrix Fisher only right now
#   
#   Sdraws<-matrix(0,B,9)
#   Kdraws<-rep(0,B)
#   
#   Snew=S0
#   Knew=kappa0
#   
#   for(i in 1:(burnin+1)){
#     Snew<-S_MCMC_CPP(Rs,Snew,rho,Knew,Cayley)
#     Knew<-kap_MCMC_CPP(Rs,Knew,sigma,Snew,Cayley)
#   }
#   
#   Sdraws[1,]<-Snew
#   Kdraws[1]<-Knew
#   Saccept<-0
#   Kaccept<-0
#   
#   for(j in 2:B){
#     Sold<-matrix(Sdraws[(j-1),],3,3)
#     Sdraws[j,]<-S_MCMC_CPP(Rs,Sold,rho,Kdraws[j-1],Cayley)
#     Saccept<-Saccept+as.numeric(all(Sdraws[j-1,]==Sdraws[j,]))
#     
#     Snew<-matrix(Sdraws[(j),],3,3)
#     Kd<-kap_MCMC_CPP(Rs,Kdraws[j-1],sigma,Snew,Cayley)
#     Kdraws[j]<-Kd
#     Kaccept<-Kaccept+as.numeric(Kdraws[j]==Kdraws[j-1])
#   }
#   
#   return(list(S=Sdraws,Sacc=Saccept/B,kappa=Kdraws,Kacc=Kaccept/B))
#   
# }

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
