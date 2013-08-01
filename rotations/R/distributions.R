#' Angular distributions
#' 
#' Density and random variate generation for symmetric probability distributions in the rotations package
#' 
#' The functions for the density function and random variate generation are named in the usual form dxxxx, pxxxx and rxxxx 
#' respectively.        
#' \itemize{
#' 	\item See \code{\link{Cayley}} for the Cayley distribution.
#' 	\item See \code{\link{Fisher}} for the matrix Fisher distribution.
#' 	\item See \code{\link{Haar}} for the uniform distribution on the circle.
#' 	\item See \code{\link{Mises}} for the von Mises-Fisher distribution.
#' }
#' 
#' @name Angular-distributions

NULL


rar <- function(n, f, M, ...) {
  res <- vector("numeric", length = n)
  for (i in 1:n) res[i] <- arsample.unif(f, M, ...)
  return(res)
}


#' The Symmetric Cayley Distribution
#'
#' Density and random generation for the Cayley distribution with concentration kappa (\eqn{\kappa})
#'
#' The symmetric Cayley distribution with concentration kappa (or circular variance nu) had density 
#' \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r).}{C(r |\kappa)= \Gamma(\kappa+2)(1+cos r)^\kappa(1-cos r)/[\Gamma(\kappa+1/2)2^(\kappa+1)\sqrt\pi].}
#'
#' @name Cayley
#' @aliases Cayley rcayley dcayley
#' @param r,q vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @param kappa Concentration paramter
#' @param nu The circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @return  \item{dcayley}{gives the density}
#'          \item{pcayley}{gives the distribution function}
#'          \item{rcayley}{generates random deviates}
#' @seealso \link{Angular-distributions} for other distributions in the rotations package
#' @cite Schaeben97 leon06

NULL


#' @rdname Cayley
#' @aliases Cayley dcayley pcayley rcayley
#' @export

dcayley <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- cayley_kappa(nu)
  
 	den <- 0.5 * gamma(kappa + 2)/(sqrt(pi) * 2^kappa * gamma(kappa + 0.5)) * (1 + cos(r))^kappa * (1 - cos(r))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
}

#' @rdname Cayley
#' @aliases Cayley dcayley pcayley rcayley
#' @export

pcayley<-function(q,kappa=1,nu=NULL,lower.tail=TRUE){
  
  n<-length(q)
  cdf<-rep(NA,n)
  
  for(i in 1:n)
    cdf[i]<-max(min(integrate(dcayley,-pi,q[i],kappa,nu,Haar=F)$value,1),0)
  
  if(lower.tail) 
    return(cdf) else return((1-cdf))
}

#' @rdname Cayley
#' @aliases Cayley dcayley pcayley rcayley
#' @export

rcayley <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- cayley_kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  bet <- rbeta(n, kappa + 0.5, 3/2)
  theta <- acos(2 * bet - 1) * (1 - 2 * rbinom(n, 1, 0.5))
  return(theta)
}

#' The Matrix Fisher Distribution
#'
#' Density and random generation for the matrix Fisher distribution with concentration kappa (\eqn{\kappa})
#'
#' The matrix Fisher distribution with concentration kappa (or circular variance nu) has density
#' \deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}{C(r|\kappa)=exp[2\kappa cos(r)][1-cos(r)]/(2\pi[I0(2\kappa)-I1(2\kappa)])}
#' where \eqn{\mathrm{I_p}(\cdot)}{Ip()} denotes the Bessel function of order \eqn{p} defined as  
#' \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos r}dr}{Ip(\kappa)} is the modified Bessel function with parameters \eqn{p} and \eqn{kappa}.
#'
#' @name Fisher
#' @aliases Fisher dfisher rfisher pfisher
#' @param r,q vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @param kappa concentration paramter
#' @param nu circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @return  \item{dfisher}{gives the density}
#'          \item{pfisher}{gives the distribution function}
#'          \item{rfisher}{generates random deviates}
#' @seealso \link{Angular-distributions} for other distributions in the rotations package

NULL

#' @rdname Fisher
#' @aliases Fisher dfisher pfisher rfisher
#' @export

dfisher <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
  n<-length(r)
  den<-rep(0,n)
  
 	den <- exp(2 * kappa * cos(r)) * (1 - cos(r))/(2 * pi * (besselI(2 * kappa, 0) - besselI(2 * kappa, 1)))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
  
}

#' @rdname Fisher
#' @aliases Fisher dfisher pfisher rfisher
#' @export

pfisher<-function(q,kappa=1, nu= NULL, lower.tail=TRUE){
  
  n<-length(q)
  cdf<-rep(NA,n)
  
  for(i in 1:n)
    cdf[i]<-max(min(integrate(dfisher,-pi,q[i],kappa,nu,Haar=F)$value,1),0)
  
  if(lower.tail)
    return(cdf) else return((1-cdf))
}

#' @rdname Fisher
#' @aliases Fisher dfisher pfisher rfisher
#' @export


rfisher <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  M <- max(dfisher(seq(-pi, pi, length = 1000), kappa,Haar=F))
  return(rar(n, dfisher, M, kappa = kappa, Haar=F))
}

#' Haar Measure
#'
#' Uniform density on the circle
#' 
#' The uniform density on the circle  (also referred to as Haar measure)
#' has the density \deqn{C_U(r)=\frac{[1-cos(r)]}{2\pi}.}{C(r)=[1-cos(r)]/2\pi.}
#'
#' @name Haar
#' @aliases Haar dhaar phaar rhaar
#' @param r,q vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @return  \item{dhaar}{gives the density}
#'          \item{phaar}{gives the distribution function}
#'          \item{rhaar}{generates random deviates}
#' @seealso \link{Angular-distributions} for other distributions in the rotations package

NULL

#' @rdname Haar
#' @aliases Haar dhaar phaar rhaar
#' @export

dhaar <- function(r){
	
	den <-(1 - cos(r))/(2 * pi)
	
	return(den)
} 

#' @rdname Haar
#' @aliases Haar dhaar phaar rhaar
#' @export

phaar<-function(q,lower.tail=TRUE){
  
  cdf<-(q-sin(q)+pi)/(2*pi)
  
  ind<-which(cdf>1)
  cdf[ind]<-1
  
  ind2<-which(cdf<0)
  cdf[ind2]<-0
  
  if(lower.tail)
    return(cdf) else return((1-cdf))
}

#' @rdname Haar
#' @aliases Haar dhaar phaar rhaar
#' @export

rhaar<-function(n){
	
	lenn<-length(n)
	if(lenn>1)
		n<-lenn
	
  return(rar(n, dhaar, 1/pi))
}

#' The Circular-von Mises Distribution
#'
#' Density and random generation for the the circular von Mises distribution with concentration kappa
#' 
#' The circular von Mises-based distribution has the density
#' \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa cos(r)}.}{C(r|\kappa)=exp[\kappa cos(r)]/[2\pi I(\kappa)]}
#' where \eqn{\mathrm{I_0}(\kappa)}{I(\kappa)} is the modified bessel function of order 0.
#'
#' @name Mises
#' @aliases Mises dvmises rvmises pvmises
#' @param r,q vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required
#' @param kappa concentration paramter
#' @param nu The circular variance, can be used in place of kappa
#' @param Haar logical; if TRUE density is evaluated with respect to Haar
#' @param lower.tail logica; if TRUE probabilites are \eqn{P(X\leq x)}{P(X\le x)} otherwise, \eqn{P(X>x)}
#' @return  \item{dvmises}{gives the density}
#'          \item{pvmises}{gives the distribution function}
#'          \item{rvmises}{generates random deviates}
#' @seealso \link{Angular-distributions} for other distributions in the rotations package

NULL

#' @rdname Mises
#' @aliases Mises dvmises pvmises rvmises
#' @export

dvmises <- function(r, kappa = 1, nu = NULL, Haar = T) {
  
  if(!is.null(nu))
    kappa <- vmises_kappa(nu)
  
  den <- 1/(2 * pi * besselI(kappa, 0)) * exp(kappa * cos(r))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) {
    return(den/(1 - cos(r)))
  } else {
    return(den)
  }
}

#' @rdname Mises
#' @aliases Mises dvmises pvmises rvmises
#' @export

pvmises<-function(q,kappa=1,nu=NULL,lower.tail=TRUE){
  
  n<-length(q)
  cdf<-rep(NA,n)
  
  for(i in 1:n)
    cdf[i]<-max(min(integrate(dvmises,-pi,q[i],kappa,nu,Haar=F)$value,1),0)
  
  if(lower.tail) 
    return(cdf) else return((1-cdf))
}

#' @rdname Mises
#' @aliases Mises dvmises pvmises rvmises
#' @export

rvmises <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- vmises_kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  u <- runif(3, 0, 1)
  a <- 1 + sqrt(1 + 4 * kappa^2)
  b <- (a - sqrt(2 * a))/(2 * kappa)
  r <- (1 + b^2)/(2 * b)
  theta <- rep(10, n)
  
  for (i in 1:n) {
    
    while (theta[i] == 10) {
      # Step 1
      u <- runif(3, 0, 1)
      z <- cos(pi * u[1])
      f <- (1 + r * z)/(r + z)
      c <- kappa * (r - f)
      
      # Step 2
      u <- runif(3, 0, 1)
      if ((c * (2 - c) - u[2]) > 0) {
        
        theta[i] = sign(u[3] - 0.5) * acos(f)
        
      } else {
        
        if (log(c/u[2]) + 1 - c < 0) {
          u <- runif(3, 0, 1)
        } else {
          u <- runif(3, 0, 1)
          theta[i] = sign(u[3] - 0.5) * acos(f)
        }
      }
    }
  }
  return(theta)
}



#' Generic UARS Distribution
#'
#' Density and random generation for the the generic uniform-axis random-spin class of distributions
#' 
#' For the rotation R with central orientation S and concentration \eqn{kappa} the UARS density is given by 
#' \deqn{f(R|S,\kappa)=\frac{4\pi}{3-tr(S'R)}C(acos[tr(S'R)-1]/2|\kappa)}{f(R|S,\kappa)=4\pi C(acos[tr(S'R)-1]/2|\kappa)/[3-tr(S'R)]}
#' where \eqn{C(r|\kappa)} is one of the \link{Angular-distributions}.
#'
#' @name UARS
#' @aliases UARS duars ruars
#' @param R Value at which to evaluate the UARS density
#' @param dangle The function to evaulate the angles from: e.g. dcayley, dvmises, dfisher, dhaar
#' @param pangle The form of the angular density: e.g. pcayley, pvmises, pfisher, phaar
#' @param rangle The function from which to simulate angles: e.g. rcayley, rvmises, rhaar, rfisher
#' @param S principal direction of the distribution
#' @param kappa concentration of the distribution
#' @param space Indicates the desired representation: matrix (SO3) or quaternion (Q4)
#' @param ... additional arguments passed to the angular distribution
#' @return  \item{duars}{gives the density}
#'          \item{puars}{gives the distribution function}
#'          \item{ruars}{generates random deviates}
#' @seealso For more on the angular distribution options see \link{Angular-distributions}
#' @cite bingham09

NULL


#' @rdname UARS
#' @aliases UARS duars puars ruars
#' @export

duars<-function(R,dangle,S=diag(3),kappa=1,...){
	
	R<-formatSO3(R)
	rs<-angle(R)
	cr<-dangle(rs,kappa,...)	
	trStO<-2*cos(rs)+1
	
	den<-4*pi*cr/(3-trStO)
	
	return(den)
}

#' @rdname UARS
#' @aliases UARS duars puars ruars
#' @export

puars<-function(R,pangle,S=diag(3),kappa=1,...){
	
	#This is not a true CDF, but it will work for now
	R<-formatSO3(R)
	rs<-angle(R)
	
	if(is.null(pangle)){
		
		n<-length(rs)
		cr<-rep(0,n)
		
		for(i in 1:length(rs))
			cr[i]<-length(which(rs<=rs[i]))/n
		
	}else{		
		cr<-pangle(rs,kappa,...)
	}
	
	#trStO<-2*cos(rs)+1
	
	#den<-4*pi*cr/(3-trStO)
	
	return(cr)
	
}

#' @rdname UARS
#' @aliases UARS duars puars ruars
#' @export

ruars<-function(n,rangle,S=NULL,kappa=1,space="SO3",...){
  
  r<-rangle(n,kappa,...)
  Rs<-genR(r,S,space)
  
  return(Rs)
}