#' MCMC for rotation data
#'
#' Use non-informative Bayes to infer about the central orientation and concentration parameter for a sample of rotations.
#'
#' The procedures detailed in \cite{bingham2009b} and \cite{bingham2010} are implemented to get
#' draws from the posterior distribution for the central orientation and concentration parameters for 
#' a sample of 3D rotations.  A uniform prior on SO(3) is used for the central orientation and the appropriate
#' Jeffreys prior is used for the concentration parameter.  
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type Angular distribution assumed on R.  Options are \code{\link{Cayley}}, \code{\link{Fisher}} or \code{\link{Mises}}
#' @param S0 initial estimate of central orientation
#' @param kappa0 initial estimate of concentration parameter
#' @param tuneS central orientation tuning parameter, concentration of proposal distribution
#' @param tuneK concentration tuning parameter, standard deviation of proposal distribution
#' @param burn_in number of draws to use as burn-in
#' @param m number of draws to keep from posterior distribution
#' @return  list of \item{S}{Draws from the posterior distribution for central orientation S}
#'          \item{kappa}{Draws from the posterior distribution for concentration parameter kappa}
#'          \item{Saccept}{Acceptance rate for cenral orientaion draws}
#'          \item{Kaccept}{Acceptance rate for concentration draws}
#' @cite bingham2009b bingham2010
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=4)
#' draws<-MCMCSO3(Rs,type='Cayley',S0=mean(Rs),kappa0=2,tuneS=39,tuneK=.8,burn_in=100,m=5000)

MCMCSO3<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  UseMethod("MCMCSO3")
}


#' @rdname MCMCSO3
#' @method MCMCSO3 SO3
#' @S3method MCMCSO3 SO3

MCMCSO3.SO3<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  
  if(type %in% c("Cayley","cayley")){
    
    lpangle <- lpcayley
    
  }else if(type %in% c("Fisher","fisher")){
    
    lpangle <- lpfisher
    
  }else if(type %in% c("Mises","mises")){
    
    lpangle <- lpvmises
    
  }else{
    stop("Invalid choise of type: please choose Cayley, Fisher or Mises.")
  }
  
  listRes<-both_MCMC_CPP(x,S0, kappa0,tuneS,tuneK,burn_in,m, lpangle)
  listRes$S<-as.SO3(listRes$S)
  
  return(listRes)
}

#' @rdname MCMCSO3
#' @method MCMCSO3 Q4
#' @S3method MCMCSO3 Q4

MCMCSO3.Q4<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  
  Rs<-SO3(x)
  S0<-as.SO3(matrix(SO3(S0),3,3))
  SO3Res<-MCMCSO3(Rs,type,S0,kappa0,tuneS,tuneK,burn_in,m)
  Q4Res<-list(Q=Q4(SO3Res$S),kappa=SO3Res$kappa,Qaccept=SO3Res$Saccept,Kaccept=SO3Res$Kaccept)
  return(Q4Res)
  
}

#' Bayes credible regions
#'
#' Find the radius of a \eqn{100(1-\alpha)}\% credible region for the central orientation and concentration parameter using 
#' non-informative Bayes.
#'
#' Compute the radius of a \eqn{100(1-\alpha)}\% credible region for the central orientation and concentration parameter
#' as described in \cite{bingham2009b} and \cite{bingham2010}.  The posterior mode is returned along with the radius
#' of the credible region centered at the posterior mode.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type Angular distribution assumed on R.  Options are \code{\link{Cayley}}, \code{\link{Fisher}} or \code{\link{Mises}}
#' @param S0 initial estimate of central orientation
#' @param kappa0 initial estimate of concentration parameter
#' @param tuneS central orientation tuning parameter, concentration of proposal distribution
#' @param tuneK concentration tuning parameter, standard deviation of proposal distribution
#' @param burn_in number of draws to use as burn-in
#' @param m number of draws to keep from posterior distribution
#' @param alp alpha level desired, e.g. 0.05 or 0.10.
#' @return  list of \item{Shat,Qhat}{the posterior mode}
#'          \item{Radius}{the radius of the credible region centered at the posterior mode}
#' @seealso \code{\link{fisheretal}}, \code{\link{prentice}}, \code{\link{chang}}, \code{\link{zhang}}
#' @cite bingham2009b bingham2010
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=4)
#' region(Rs,type='Cayley',method='Bayes',estimator='mean',
#' S0=mean(Rs),kappa0=2,tuneS=39,tuneK=.8,burn_in=100,alp=.01)

bayesCR<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000,alp=0.1){
  UseMethod("bayesCR")
}


#' @rdname bayesCR
#' @method bayesCR SO3
#' @S3method bayesCR SO3

bayesCR.SO3<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000,alp=0.1){
  
  if(type %in% c("Cayley","cayley")){
    
    lpangle <- lpcayley
    
  }else if(type %in% c("Fisher","fisher")){
    
    lpangle <- lpfisher
    
  }else if(type %in% c("Mises","mises")){
    
    lpangle <- lpvmises
    
  }else{
    stop("Invalid choise of type: please choose Cayley, Fisher or Mises.")
  }
  
  listRes<-both_MCMC_CPP(x,S0, kappa0,tuneS,tuneK,burn_in,m, lpangle)
  Sdraws<-as.SO3(listRes$S)
  Shat<-mean(Sdraws)
  rs<-dist(Sdraws,Shat)
  
  return(list(Shat=Shat,Radius=quantile(rs,1-alp)))
}

#' @rdname bayesCR
#' @method bayesCR Q4
#' @S3method bayesCR Q4

bayesCR.Q4<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000,alp=0.1){
  
  Rs<-SO3(x)
  S0<-as.SO3(matrix(SO3(S0),3,3))
  SO3Res<-bayesCR(Rs,type,S0,kappa0,tuneS,tuneK,burn_in,m,alp)
  Q4Res<-list(Qhat=Q4(SO3Res$Shat),Radius=SO3Res$Radius)
  return(Q4Res)
  
}

#' Parameter estimates based on non-informative Bayes
#'
#' Use non-informative Bayes to estimate the central orientation and concentration parameter of a sample of rotations.
#'
#' The procedures detailed in \cite{bingham2009b} and \cite{bingham2010} are implemented to get
#' draws from the posterior distribution for the central orientation and concentration parameters for 
#' a sample of 3D rotations.  A uniform prior on SO(3) is used for the central orientation and the appropriate
#' Jeffreys prior is used for the concentration parameter.  
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type Angular distribution assumed on R.  Options are \code{\link{Cayley}}, \code{\link{Fisher}} or \code{\link{Mises}}
#' @param S0 initial estimate of central orientation
#' @param kappa0 initial estimate of concentration parameter
#' @param tuneS central orientation tuning parameter, concentration of proposal distribution
#' @param tuneK concentration tuning parameter, standard deviation of proposal distribution
#' @param burn_in number of draws to use as burn-in
#' @param m number of draws to keep from posterior distribution
#' @return  list of \item{Shat}{Draws from the posterior distribution for central orientation S}
#'          \item{kappa}{Draws from the posterior distribution for concentration parameter kappa}
#' @seealso \code{\link{mean.SO3}}, \code{\link{median.SO3}}
#' @cite bingham2009b bingham2010
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=4)
#' ests<-bayes.mean(Rs,type='Cayley',S0=mean(Rs),kappa0=4,tuneS=39,tuneK=.8,burn_in=100,m=5000)

bayes.mean<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  UseMethod("bayes.mean")
}


#' @rdname bayes.mean
#' @method bayes.mean SO3
#' @S3method bayes.mean SO3

bayes.mean.SO3<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  
  if(type %in% c("Cayley","cayley")){
    
    lpangle <- lpcayley
    
  }else if(type %in% c("Fisher","fisher")){
    
    lpangle <- lpfisher
    
  }else if(type %in% c("Mises","mises")){
    
    lpangle <- lpvmises
    
  }else{
    stop("Invalid choise of type: please choose Cayley, Fisher or Mises.")
  }
  
  listRes<-both_MCMC_CPP(x,S0, kappa0,tuneS,tuneK,burn_in,m, lpangle)
  postRes<-list(Shat=mean.SO3(listRes$S),kappa=mean(listRes$kappa))
  
  return(postRes)
}

#' @rdname bayes.mean
#' @method bayes.mean Q4
#' @S3method bayes.mean Q4

bayes.mean.Q4<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  
  Rs<-SO3(x)
  S0<-as.SO3(matrix(SO3(S0),3,3))
  SO3Res<-MCMCSO3(Rs,type,S0,kappa0,tuneS,tuneK,burn_in,m)
  Q4Res<-list(Qhat=Q4(mean(SO3Res$S)),kappa=mean(SO3Res$kappa))
  return(Q4Res)
  
}
