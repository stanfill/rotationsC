#' Mean rotation
#'
#' Compute the sample geometric or projected mean.
#'
#' This function takes a sample of 3D rotations (in matrix or quaternion form) and returns the projected arithmetic mean denoted \eqn{\widehat{\bm S}_P}{S_P} or
#' geometric mean \eqn{\widehat{\bm S}_G}{S_G} according to the \code{type} option.
#' For a sample of \eqn{n} rotations in matrix form \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd^2(\bm{R}_i,\bm{S})}{argmin\sum d^2(Ri,S)} 
#' where \eqn{d} is the Riemannian or Euclidean distance.
#' For more on the projected mean see \cite{moakher02} and for the geometric mean see \cite{manton04}.
#' For the projected mean from a quaternion point of view see \cite{tyler1981}.
#' 
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix form (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type string indicating "projected" or "geometric" type mean estimator.
#' @param epsilon stopping rule for the geometric-mean.
#' @param maxIter maximum number of iterations allowed for geometric-mean.
#' @param ... additional arguments.
#' @return Estimate of the projected or geometric mean of the sample.
#' @aliases mean.Q4
#' @seealso \code{\link{median.SO3}}
#' @cite tyler1981, moakher02, manton04
#' @S3method mean SO3
#' @method mean SO3
#' @examples
#' Rs<-ruars(20,rvmises,kappa=0.01)
#' mean(Rs)
#' Qs<-Q4(Rs)
#' mean(Qs)

mean.SO3 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000, ...) {
	
	Rs<-formatSO3(x)	
	
	if(nrow(Rs)==1)
		return(Rs)
	
  if (!(type %in% c("projected", "geometric")))
    stop("type needs to be one of 'projected' or 'geometric'.")

	
	if(type=='projected'){
		R<-meanSO3C(Rs)
	}else{
		R<-gmeanSO3C(Rs,maxIter,epsilon)
	}
	
  class(R)<-"SO3"
  return(R)
}

#' @rdname mean.SO3
#' @aliases mean.SO3
#' @S3method mean Q4
#' @method mean Q4

mean.Q4 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000,...) {
	
	Qs<-formatQ4(x)
	
	if(nrow(Qs)==1)
		return(Qs)
	
	if(type=='projected'){
		
		R<-meanQ4C(Qs)
		
	}else{
		
		Rs<-SO3(Qs)
  	R<-gmeanSO3C(Rs,maxIter,epsilon)
		R<-Q4.SO3(R)
	}
	
	class(R)<-'Q4'
  return(R)
  
}


#' Median rotation
#' 
#' Compute the sample projected or geometric median.
#'
#' The median-type estimators are defined as \deqn{\widetilde{\bm{S}}=argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd(\bm{R}_i,\bm{S}).}{argmin\sum d(Ri,S).}  
#' If the choice of distance metric \eqn{d} is Riemannian then the estimator is called the geometric median, 
#' and if the distance metric in Euclidean then it is called the projected median.
#' The algorithm used in the geometric case is discussed in \cite{hartley11} 
#' and the projected case was written by the authors.
#'
#' @name median.SO3
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix form (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type string indicating "projected" or "geometric" type mean estimator.
#' @param epsilon stopping rule.
#' @param maxIter maximum number of iterations allowed before returning most recent estimate.
#' @param ... additional arguments.
#' @return An estimate of the projected or geometric mean.
#' @aliases median.Q4 median.SO3
#' @seealso \code{\link{mean.SO3}}
#' @cite hartley11
#' @export

median<-function(x,...){
  UseMethod("median")
}

#' @rdname median.SO3
#' @aliases median.Q4 median
#' @method median SO3
#' @S3method median SO3

median.SO3 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000,...) {
  
	Rs<-formatSO3(x)
	n<-nrow(Rs)
	
	if(nrow(Rs)==1)
		return(Rs)
	
  stopifnot(type %in% c("projected", "geometric"))
  
  if (type == "projected") {
		
  	S<-medianSO3C(Rs,maxIter,epsilon)
      
  } else {
      
		S<-HartmedianSO3C(Rs,maxIter,epsilon)
      
  }
    
	class(S)<-"SO3"
  return(S)
}


#' @rdname median.SO3
#' @aliases median.SO3 median
#' @method median Q4
#' @S3method median Q4

median.Q4 <- function(x, type = "projected", epsilon = 1e-05, maxIter = 2000,...) {
	
	Qs<-formatQ4(x)
	
	if(length(Qs)==4)
		return(Qs)

  Rs<-SO3(Qs)
  
  R<-median(Rs,type,epsilon,maxIter)
  
  return(Q4.SO3(R))
}


#' Weighted Mean Rotation
#'
#' Compute the weighted geometric or projected mean of a sample of rotations.
#'
#' This function takes a sample of 3D rotations (in matrix or quaternion form) and returns the weighted projected arithmetic mean denoted \eqn{\widehat{\bm S}_P}{S_P} or
#' geometric mean \eqn{\widehat{\bm S}_G}{S_G} according to the \code{type} option.
#' For a sample of \eqn{n} rotations in matrix form \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd^2(\bm{R}_i,\bm{S})}{argmin\sum d(bar(R),S)} where \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i}{bar(R)=\sum R_i/n} and the distance metric \eqn{d}
#' is the Riemannian or Euclidean.  For more on the projected mean see \cite{moakher02} and for the geometric mean see \cite{manton04}.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix form (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param w vector of weights the same length as the number of rows in x giving the weights to use for elements of x.
#' @param type string indicating "projectced" or "geometric" type mean estimator.
#' @param epsilon stopping rule for the geometric method.
#' @param maxIter maximum number of iterations allowed before returning most recent estimate.
#' @param ... only used for consistency with mean.default.
#' @return Weighted mean of the sample.
#' @seealso \code{\link{median.SO3}}, \code{\link{mean.SO3}}
#' @aliases weighted.mean.Q4
#' @cite moakher02
#' @S3method weighted.mean SO3
#' @method weighted.mean SO3
#' @examples
#' Rs<-ruars(20,rvmises,kappa=0.01)
#' wt<-abs(1/angle(Rs))
#' weighted.mean(Rs,wt)
#' Qs<-Q4(Rs)
#' weighted.mean(Qs,wt)

weighted.mean.SO3 <- function(x, w, type = "projected", epsilon = 1e-05, maxIter = 2000, ...) {
	
	Rs<-formatSO3(x)
	
	if(nrow(Rs)==1)
		return(Rs)
	
	if(length(w)!=nrow(Rs))
		stop("'Rs' and 'w' must have same length")
	
	if (!(type %in% c("projected", "geometric")))
		stop("type needs to be one of 'projected' or 'geometric'.")
	
	if(any(w<0))
		warning("Negative weights were given.  Their absolute value is used.")
	
	w<-abs(w/sum(w))
	
	wRs<-w*Rs
	
	R <- as.SO3(project.SO3(matrix(colSums(wRs), 3, 3)))
	
	if (type == "geometric") {
		n <- nrow(Rs)
		d <- 1
		iter <- 0
		s <- matrix(0, 3, 3)
		
		while (d >= epsilon) {
			
			R <- R %*% exp_skew(s)
			
			s <- matrix(colSums(w*t(apply(Rs, 1, tLogMat, S = R))), 3, 3)
			
			d <- norm(s, type = "F")
			
			iter <- iter + 1
			
			if (iter >= maxIter) {
				warning(paste("No convergence in ", iter, " iterations."))
				return(as.SO3(R))
			}
		}
		R<-as.SO3(R)	
	}
	
	return(R)
}

#' @rdname weighted.mean.SO3
#' @aliases weighted.mean.SO3
#' @S3method weighted.mean Q4
#' @method weighted.mean Q4


weighted.mean.Q4 <- function(x, w, type = "projected", epsilon = 1e-05, maxIter = 2000,...) {
	
	Qs<-formatQ4(x)
	
	if(nrow(Qs)==1)
		return(Qs)
	
	Rs<-SO3(Qs)
	
	R<-weighted.mean.SO3(Rs,w,type,epsilon,maxIter)
	
	return(Q4.SO3(R))
	
}
