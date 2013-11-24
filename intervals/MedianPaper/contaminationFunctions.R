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
