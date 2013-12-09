#' Confidence and creible regions for the central orientation
#'
#' Find the radius of a \eqn{100(1-\alpha)}\% confidence or credible region for the central orientation based on the projected mean or median.
#' For more on the currently available methods see \code{\link{prentice}}, \code{\link{fisheretal}}, \code{\link{chang}},
#' \code{\link{zhang}} and \code{\link{bayesCR}}.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param method character string specifying which type of interval to report, "bayes", "eigen" or "moment" based theory.
#' @param type characted string, "bootstrap" or "theory" are available.  For Bayes regions, give the type of likelihood: "Cayley","Mises" or "Fisher."
#' @param estimator character string either "mean" or "median."  Note that not all method/type combinations are available for both estimators.
#' @param alp the alpha level desired, e.g. 0.05 or 0.10.
#' @param ... additional arguments that are method specific.
#' @return For frequentist regions only the radius of the confidence region centered at the specified estimator is returned.
#'   For Bayes regions the posterior mode and radius of the credible region centered at that mode is returned.
#' @seealso \code{\link{bayesCR}}, \code{\link{prentice}}, \code{\link{fisheretal}}, \code{\link{chang}}, \code{\link{zhang}}
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=100)
#' region(Rs,method='eigen',type='theory',estimator='mean',alp=0.1)
#' region(Rs,method='eigen',type='bootstrap',estimator='mean',alp=0.1,symm=TRUE)
#' region(Rs,method='moment',type='bootstrap',estimator='mean',alp=0.1,m=100)
#' region(Rs,method='moment',type='theory',estimator='mean',alp=0.1)
#' region(Rs,method='Bayes',type='Cayley',estimator='mean',
#' S0=mean(Rs),kappa0=2,tuneS=39,tuneK=.8,burn_in=100,alp=.01)

region<-function(x,method, type, estimator,alp,...){
	UseMethod("region")
}


#' @rdname region
#' @method region Q4
#' @S3method region Q4

region.Q4<-function(x,method, type, estimator,alp=NULL,...){
	
	Qs<-formatQ4(x)
	
	if(is.null(alp)){
		#Take a default alpha=0.1 if no level is specified
		alp<-.1
		warning("No alpha-level specified, 0.1 used by default.")
	}
	
	if(method%in%c('Eigen','eigen') & type%in%c("Theory","theory")){
		
		if(estimator!='mean'){
			stop("The method due to Prentice is only available for the mean estimator.")
		}
		
		r<-prentice.Q4(x=Qs,alp=alp)
		
		return(r)
		
	}else	if(method%in%c('Moment','moment') & type%in%c("Bootstrap","bootstrap")){
		
		r<-zhang.Q4(x=Qs,estimator=estimator,alp=alp,...)
		
		return(r)
		
	}else	if(method%in%c('Eigen','eigen') & type%in%c("Bootstrap","bootstrap")){
		
		if(estimator!='mean'){
			stop("The method due to Fisher et al. is only available for the mean estimator.")
		}
		
		r<-fisheretal.Q4(x=Qs,alp=alp,...)
		
		return(r)
		
	}else	if(method%in%c('Moment','moment') & type%in%c("Theory","theory")){
		
		r<-chang.Q4(x=Qs,estimator=estimator,alp=alp)
		
		return(r)
		
	}else if(method%in%c('Bayes','bayes')){
	  
	  if(estimator!='mean'){
	    stop("Bayes confidence regions are only available for the mean estimator.")
	  }
	  
	  r<-bayesCR.Q4(x=Qs,type=type,alp=alp,...)
	  
	  return(r)
	  
	}else{
		
		stop("Please choose a correct combination of method, type and estimator.  See help file.")
		
	}
	
}


#' @rdname region
#' @method region SO3
#' @S3method region SO3

region.SO3<-function(x,method,type,estimator,alp=NULL,...){
	
	Rs<-formatSO3(x)
	
	if(is.null(alp)){
		#Take a default alpha=0.1 if no level is specified
		alp<-.1
		warning("No alpha-level specified, 0.1 used by default.")
	}
	
	if(method%in%c('Eigen','eigen') & type%in%c("Theory","theory")){
		
		if(estimator!='mean'){
			stop("The method due to Prentice is only available for the mean estimator.")
		}
		
		r<-prentice.SO3(x=Rs,alp=alp)
		
		return(r)
		
	}else if(method%in%c('Moment','moment') & type%in%c("Bootstrap","bootstrap")){
		
		r<-zhang.SO3(x=Rs,estimator=estimator,alp=alp,...)
		
		return(r)
		
	}else if(method%in%c('Eigen','eigen') & type%in%c("Bootstrap","bootstrap")){
		
		if(estimator!='mean'){
			stop("The method due to Fisher et al. is only available for the mean estimator.")
		}
		
		r<-fisheretal.SO3(x=Rs,alp=alp,...)
		
		return(r)
		
	}else if(method%in%c('Moment','moment') & type%in%c("Theory","theory")){
		
		r<-chang.SO3(x=Rs,estimator=estimator,alp=alp)
		
		return(r)
		
	}else if(method%in%c('Bayes','bayes')){
	  
	  if(estimator!='mean'){
	    stop("Bayes confidence regions are only available for the mean estimator.")
	  }
	  
	  r<-bayesCR.SO3(x=Rs,type=type,alp=alp,...)
	  
	  return(r)
	  
	}else{
		
	  stop("Please choose a correct combination of method, type and estimator.  See help file.")
		
	}
	
}

#' Eigenvector theory confidence region
#'
#' Find the radius of a \eqn{100(1-\alpha)}\% confidence region for the projected mean based on eigenvector based result.
#'
#' Compute the radius of a \eqn{100(1-\alpha)}\% confidence region for the central orientation based on the projected mean
#' estimator using the method due to \cite{prentice1986}.  For a rotation specific version see \cite{rancourt2000}. The variablity
#' in each axis is different so each axis will have its own radius. 
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param alp alpha level desired, e.g. 0.05 or 0.10.
#' @return Radius of the confidence region centered at the projected mean for each of the x-, y- and z-axes.
#' @seealso \code{\link{bayesCR}}, \code{\link{fisheretal}}, \code{\link{chang}}, \code{\link{zhang}}
#' @cite prentice1986, rancourt2000
#' @export
#' @examples
#' Qs<-ruars(20,rcayley,kappa=100,space='Q4')
#' region(Qs,method='eigen',type='theory',alp=0.1,estimator='mean')

prentice<-function(x,alp){
	UseMethod("prentice")
}


#' @rdname prentice
#' @method prentice Q4
#' @S3method prentice Q4

prentice.Q4<-function(x,alp=NULL){
	#This takes a sample qs and returns the radius of the confidence region
	#centered at the projected mean

	if(is.null(alp)){
		#Take a default alpha=0.1 if no level is specified
		alp<-.1
		warning("No alpha-level specified, 0.1 used by default.")
	}
	Qs<-x
	n<-nrow(Qs)
	Shat<-mean(Qs)
	Phat<-pMat(Shat)
	
	Rhat<-Qs%*%Phat
	resids<-matrix(0,n,3)
	VarShat<-matrix(0,3,3)
	
	resids<-2*Rhat[,1]*matrix(Rhat[,2:4],n,3)
	
	VarShat<-t(resids)%*%resids/(n-1)
	
	RtR<-t(Rhat)%*%Rhat
	Ahat<-(diag(RtR[1,1],3,3)-RtR[-1,-1])/n
	
	Tm<-diag(n*Ahat%*%solve(VarShat)%*%Ahat)
	
	r<-sqrt(qchisq((1-alp),3)/Tm)
	return(r)
}


#' @rdname prentice
#' @method prentice SO3
#' @S3method prentice SO3

prentice.SO3<-function(x,alp=NULL){
	Qs<-Q4(x)
	r<-prentice.Q4(Qs,alp)
	return(r)
}

#' M-estimator theory pivotal boostrap confidence region
#'
#' Compute the radius of a \eqn{100(1-\alpha)}\% confidence region for the central orientation based on M-estiamtor theory.
#' 
#' Compute the radius of a \eqn{100(1-\alpha)}\% confidence region for the central orientation based on the projected mean
#' estimator using the method due to Zhang & Nordman (2009) (unpublished MS thesis).  By construction each axis will have the same
#' radius so the radius reported is for all three axis.  A normal theory version of this procedure uses the theoretical
#' chi-square limiting distribution and is given by the \code{\link{chang}} option.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param estimator character string either "mean" or "median."
#' @param alp alpha level desired, e.g. 0.05 or 0.10.
#' @param m number of replicates to use to estiamte the critical value.
#' @return Radius of the confidence region centered at the specified estimator.
#' @seealso \code{\link{bayesCR}}, \code{\link{prentice}}, \code{\link{fisheretal}}, \code{\link{chang}}
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=100)
#' region(Rs,method='moment',type='bootstrap',alp=0.1,estimator='mean')

zhang<-function(x,estimator,alp,m){
	UseMethod("zhang")
}


#' @rdname zhang
#' @method zhang SO3
#' @S3method zhang SO3

zhang.SO3<-function(x,estimator,alp=NULL,m=300){
	
	#Rs is a n-by-9 matrix where each row is an 3-by-3 rotation matrix
	#m is the number of resamples to find q_1-a
	#alp is the level of confidence desired, e.g. 0.95 or 0.90
	#pivot logical; should the pivotal (T) bootstrap be used or nonpivotal (F)
	
  Rs<-formatSO3(x)
  
  if(estimator=='median'){
  	
  	if(is.null(alp)){
  		#Take a default alpha=0.1 if no level is specified
  		alp<-.1
  		warning("No alpha-level specified, 0.1 used by default.")
  	}
  	
  	n<-nrow(Rs)
  	stats<-zhangMedianC(Rs,m)
  	cdtilde<-cdfunsCSO3(Rs,median(Rs))
  	rad<-sqrt(as.numeric(quantile(stats,1-alp,na.rm=T))*cdtilde[1]/(2*n*cdtilde[2]^2))
  	
  }else if(estimator=='mean'){
  
  	Qs<-Q4(Rs)
  	rad<-zhang.Q4(Qs,estimator,alp,m)
  	
  }else{
  	stop("Please choose an estimator mean or median.")
  }
  
	return(min(rad,pi))
}

#' @rdname zhang
#' @method zhang Q4
#' @S3method zhang Q4

zhang.Q4<-function(x,estimator,alp=NULL,m=300){

	if(estimator=='mean'){
	
		if(is.null(alp)){
			#Take a default alpha=0.1 if no level is specified
			alp<-.1
			warning("No alpha-level specified, 0.1 used by default.")
		}
	
		Qs<-formatQ4(x)
		n<-nrow(Qs)
  	stats<-zhangQ4(Qs,m)
		#Shat<-mean(Qs)
  	cdhat<-cdfuns(Qs,estimator)
		rad<-sqrt(as.numeric(quantile(stats,1-alp,na.rm=T))*cdhat$c/(2*n*cdhat$d^2))
		
	}else if(estimator=='median'){
		
		Rs<-SO3(Qs)
		rad<-zhang.SO3(Rs,estimator,alp,m)
		
	}else{
		stop("Please choose an estimator mean or median.")
	}
	return(min(rad,pi))
}


cdfuns<-function(Qs,estimator){
  
  
  if(estimator=='mean'){
  	
  	Shat<-mean(Qs)
		cd<-cdfunsC(Qs,Shat)
		
  }else if(estimator=='median'){
  	
  	Shat<-median(Qs)
  	cd<-cdfunsCMedian(Qs,Shat)
  	
  }else{
  	stop("Please choose an estimator mean or median.")
  }
  
	return(list(c=cd[1],d=cd[2]))
}

#' M-estimator theory confidence region
#'
#' Compute the radius of a \eqn{100(1-\alpha)}\% confidence region for the central orientation based on M-estimator theory.
#' 
#' Compute the radius of a \eqn{100(1-\alpha)}\% confidence region for the central orientation centered at the projected mean
#' or median based on a result due to \cite{chang2001} amongst others.  By construction each axis will have the same
#' radius so the radius reported is for all three axes.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param estimator character string either "mean" or "median."
#' @param alp alpha level desired, e.g. 0.05 or 0.10.
#' @return Radius of the confidence region centered at the specified estimator.
#' @cite chang2001
#' @seealso \code{\link{bayesCR}}, \code{\link{prentice}}, \code{\link{fisheretal}}, \code{\link{zhang}}
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=100)
#' region(Rs,method='moment',type='theory',alp=0.1,estimator='mean')

chang<-function(x,estimator,alp){
	UseMethod("chang")
}


#' @rdname chang
#' @method chang SO3
#' @S3method chang SO3

chang.SO3<-function(x,estimator,alp=NULL){
	
	#Rs is a n-by-9 matrix where each row is an 3-by-3 rotation matrix
	#alp is the level of confidence desired, e.g. 0.95 or 0.90
	#pivot logical; should the pivotal (T) bootstrap be used or nonpivotal (F)
	
	Rs<-formatSO3(x)
	Qs<-Q4(Rs)
	rad<-chang.Q4(Qs,estimator,alp)
	return(rad)
}

#' @rdname chang
#' @method chang Q4
#' @S3method chang Q4

chang.Q4<-function(x,estimator,alp=NULL){
	
	if(is.null(alp)){
		#Take a default alpha=0.1 if no level is specified
		alp<-.1
		warning("No alpha-level specified, 0.1 used by default.")
	}
	
	Qs<-formatQ4(x)
	n<-nrow(Qs)
	
	cdhat<-cdfuns(Qs,estimator)
	
	rad<-sqrt(as.numeric(qchisq(1-alp,3))*cdhat$c/(2*n*cdhat$d^2))
	
	return(min(rad,pi))
}


#' Eigenvector theory pivotal bootstrap confidence region
#'
#' Find the radius of a \eqn{100(1-\alpha)}\% confidence region for the central orientation based on eigenvector theory.
#'
#' Compute the radius of a \eqn{100(1-\alpha)}\% confidence region for the central orientation based on the projected mean
#' estimator using the method for the mean polar axis as proposed in \cite{fisher1996}.  To be able to reduce their method
#' to a radius requires the additonal assumption of rotational symmetry, equation (10) in \cite{fisher1996}. 
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param alp alpha level desired, e.g. 0.05 or 0.10.
#' @param boot should the bootstrap or normal theory critical value be used.
#' @param m number of bootstrap replicates to use to estimate critical value.
#' @param symm logical; if TRUE (default), a symmetric region is constructed.
#' @return Radius of the confidence region centered at the projected mean.
#' @seealso \code{\link{bayesCR}}, \code{\link{prentice}}, \code{\link{chang}}, \code{\link{zhang}}
#' @cite fisher1996
#' @export
#' @examples
#' Qs<-ruars(20,rcayley,kappa=100,space='Q4')
#' region(Qs,method='eigen',type='bootstrap',alp=0.1,symm=TRUE,estimator='mean')

fisheretal<-function(x,alp,boot,m,symm){
	UseMethod("fisheretal")
}


#' @rdname fisheretal
#' @method fisheretal Q4
#' @S3method fisheretal Q4

fisheretal.Q4<-function(x,alp=NULL,boot=T,m=300,symm=TRUE){
	
	if(is.null(alp)){
		#Take a default alpha=0.1 if no level is specified
		alp<-.1
		warning("No alpha-level specified, 0.1 used by default.")
	}
	
	Qs<-formatQ4(x)
	
	if(boot){
    
	  Tstats <- fisherBootC(Qs,m,symm)
    
		qhat<-as.numeric(quantile(Tstats,1-alp,na.rm=T))
		
	}else{
		
		qhat<-qchisq(1-alp,3)
		
	}
	
	rsym<-optim(.05,optimAxis,Qs=Qs,cut=qhat,symm=T,method='Brent',lower=0,upper=pi)$par
	
	return(min(rsym,pi))
}


optimAxis<-function(r,Qs,cut,symm){
	
	Shat<-Q4(axis(mean(Qs)),r)
	if(symm){
		Tm<-fisherAxisC(Qs,Shat)
	}else{
		Tm<-fisherAxisCSymmetric(Qs,Shat)
	}
	return((Tm-cut)^2)
}


#' @rdname fisheretal
#' @method fisheretal SO3
#' @S3method fisheretal SO3

fisheretal.SO3<-function(x,alp=NULL,boot=T,m=300,symm=T){
	
	Qs<-Q4(x)
	r<-fisheretal.Q4(Qs,alp,boot,m,symm)
	
	return(r)
}
