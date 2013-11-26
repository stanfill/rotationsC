ruarsCont<-function(n,rangle,kappa1,kappa2=kappa1,p,S=id.SO3,Scont,space='SO3'){
  
  #n   		- sample size
  #rangle - angular distribution from which to simulate
  #kappa1	- concentration parameter for F data
  #kappa2 - concentration for contaminated data
  #p			- percent of n that will be contaminated
  #S			- central direction of normal data
  #Scont	-	central direction of contaminated data
  #space  - SO3 (default) or quaternions("Q4")
    
  nCont<-floor(p*n)
  nNorm<-n-nCont
  rsNorm<-rangle(nNorm,kappa=kappa1)
  RsNorm<-genR(rsNorm,S)  #Simulate from the normal distribution
  
  if(nCont!=0){
    rsCont<-rangle(nCont,kappa=kappa2)
    RsCont<-genR(rsCont,Scont) #Simulated from the contaminated distribution
    Rs<-as.SO3(rbind(RsNorm,RsCont))
  }else{
    Rs<-RsNorm
  }
  

  if(space=='Q4')
    Rs<-Q4(Rs)
  
  return(Rs)
  
}

#I'm getting seg-faults when generating rotations using ruarsCont, so try to 
#generate angles and rotations seperately using rfisher then genRCont

genRCont<-function(rs,p,S=id.SO3,Scont,space='SO3'){
  
  #rs     - angle to generate rotations based on
  #p			- percent of n that will be contaminated
  #S			- central direction of normal data
  #Scont	-	central direction of contaminated data
  #space  - SO3 (default) or quaternions("Q4")
  n<-length(rs)
  nCont<-floor(p*n)
  nNorm<-n-nCont
  rsNorm<-rs[1:nNorm]
  RsNorm<-genR(rsNorm,S)  #Simulate from the normal distribution
  
  if(nCont!=0){
    rsCont<-rs[-c(1:nNorm)]
    RsCont<-genR(rsCont,Scont) #Simulated from the contaminated distribution
    Rs<-as.SO3(rbind(RsNorm,RsCont))
  }else{
    Rs<-RsNorm
  }
  
  
  if(space=='Q4')
    Rs<-Q4(Rs)
  
  return(Rs)
  
}

newQ4<-function(Rs){
  rs<-angle(Rs)
  us<-axis(Rs)
  Qs<-cbind(cos(rs/2),sin(rs/2)*us)
  Qs<-as.Q4(Qs)
  return(Qs)
}
