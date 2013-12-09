#'  Circular variance and concentration parameter 
#'  
#'  Return the concentration parameter that corresponds to a given circular variance.
#'  
#'  The concentration parameter \eqn{\kappa} does not translate across circular distributions.  A commonly used
#'  measure of spread in circular distributions that does translate
#'  is the circular variance defined as \eqn{\nu=1-E[\cos(r)]}{\nu=1-E[cos(r)]} where \eqn{E[\cos(r)]}{E[cos(r)]} 
#'  is the mean resultant length.  See \cite{mardia2000} for more details.
#'  This function translates the circular variance \eqn{\nu} into the corresponding concentration parameter \eqn{\kappa}
#'  for the Cayley distribution.
#'  
#'  @param nu circular variance
#'  @return Concentration parameter corresponding to nu.
#'  @cite mardia2000
#'  @seealso \code{\link{Cayley}}
#'  @export

cayley_kappa<-function(nu){
  (3/nu)-2
}

fisher_nu_kappa<-function(kappa,nu){
  (1-(besselI(2*kappa,1)-.5*besselI(2*kappa,2)-.5*besselI(2*kappa,0))/(besselI(2*kappa,0)-besselI(2*kappa,1))-nu)^2
}

#'  Circular variance and concentration parameter 
#'  
#'  Return the concentration parameter that corresponds to a given circular variance.
#'  
#'  The concentration parameter \eqn{\kappa} does not translate across circular distributions.  A commonly used
#'  measure of spread in circular distributions that does translate
#'  is the circular variance defined as \eqn{\nu=1-E[\cos(r)]}{\nu=1-E[cos(r)]} where \eqn{E[\cos(r)]}{E[cos(r)]} 
#'  is the mean resultant length.  See \cite{mardia2000} for more details.
#'  This function translates the circular variance \eqn{\nu} into the corresponding concentration parameter \eqn{\kappa}
#'  for the matrix-Fisher distribution.
#'  
#'  @param nu circular variance
#'  @return Concentration parameter corresponding to nu.
#'  @seealso \code{\link{Fisher}}
#'  @cite mardia2000
#'  @export

fisher_kappa<-function(nu){
  
  kappa<-rep(0,length(nu))
  
  for(i in 1:length(nu))
    kappa[i]<-optimize(fisher_nu_kappa,interval=c(0,10),tol=.00001,nu=nu[i])$minimum
  
  return(kappa)
}


mises_nu_kappa<-function(kappa,nu){
  (1-besselI(kappa,1)/besselI(kappa,0)-nu)^2
}

#'  Circular variance and concentration parameter
#'  
#'  Return the concentration parameter that corresponds to a given circular variance.
#'  
#'  The concentration parameter \eqn{\kappa} does not translate across circular distributions.  A commonly used
#'  measure of spread in circular distributions that does translate
#'  is the circular variance defined as \eqn{\nu=1-E[\cos(r)]}{\nu=1-E[cos(r)]} where \eqn{E[\cos(r)]}{E[cos(r)]} 
#'  is the mean resultant length.  See \cite{mardia2000} for more details.
#'  This function translates the circular variance \eqn{\nu} into the corresponding concentration parameter \eqn{\kappa}
#'  for the circular-von Mises distribution.
#'  
#'  @param nu circular variance
#'  @return Concentration parameter corresponding to nu.
#'  @seealso \code{\link{Mises}}
#'  @cite mardia2000
#'  @export

vmises_kappa<-function(nu){
  
  kappa<-rep(0,length(nu))
  
  for(i in 1:length(nu))
    kappa[i]<-optimize(mises_nu_kappa,interval=c(0,10),tol=.00001,nu=nu[i])$minimum
  
  return(kappa)
}
