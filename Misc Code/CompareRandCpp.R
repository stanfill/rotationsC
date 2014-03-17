#Code various 

############################################
#Cayley isn't worth coding in C++
rcayley2 <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- cayley.kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
    n<-lenn
  
  theta<-rcayleyCpp(n,kappa)
  return(theta)
}

rs<-rcayley(20,1)
rs2<-rcayley2(20,1)

microbenchmark(rcayley(2000,1),rcayley2(2000,1))

############################################
#Try von Mises
rvmises2 <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- vmises.kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
    n<-lenn
  
  theta<-rvmisesCPP(n,kappa)
  return(theta)
}