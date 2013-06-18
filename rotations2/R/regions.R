#' Confidence Region for Mean Rotation
#'
#' Find the radius of a \eqn{100(1-\alpha)%} confidence region for the central orientation based on the projected mean estimator.
#' The current methods available are due to \code{\link{prentice}}, \code{\link{fisher}}, \code{\link{chang}},
#' and \code{\link{zhang}}.
#'
#' @param Rs,Qs A \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion form (p=4)
#' @param method Character string specifying which type of interval is required
#' @param alpha The alpha level desired, e.g. 0.95 or 0.90
#' @param ... Additional arguments that are method specific
#' @return Radius of the confidence region centered at the projected mean
#' @seealso \code{\link{prentice}} \code{\link{fisher}} \code{\link{chang}} \code{\link{zhang}}
#' @cite prentice1986, fisher1996, rancourt2000, chang2001
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=100)
#' region(Rs,method='prentice',alpha=0.1)
#' region(Rs,method='fisher',alpha=0.1,symm=T)
#' region(Rs,method='zhang',alpha=0.1,m=100)
#' region(Rs,method='chang',alpha=0.1)

region<-function(Qs,method,alpha,...){
	UseMethod("region")
}


#' @rdname region
#' @method region Q4
#' @S3method region Q4

region.Q4<-function(Qs,method,alpha,...){
	
	Qs<-formatQ4(Qs)
	
	if(method%in%c('Prentice','prentice')){
		
		r<-prentice.Q4(Qs=Qs,a=alpha)
		
		return(r)
		
	}else	if(method%in%c('Zhang','zhang')){
		
		r<-zhang.Q4(Qs=Qs,a=alpha,...)
		
		return(r)
		
	}else	if(method%in%c('Fisher','fisher')){
		
		r<-fisher.Q4(Qs=Qs,a=alpha,...)
		
		return(r)
		
	}else	if(method%in%c('Chang','chang')){
		
		r<-chang.Q4(Qs=Qs,a=alpha)
		
		return(r)
		
	}else{
		
		stop("Only the Prentice, Zhang, Chang and Fisher options are currently available")
		
	}
	
}


#' @rdname region
#' @method region SO3
#' @S3method region SO3

region.SO3<-function(Rs,method,alpha,...){
	
	Rs<-formatSO3(Rs)
	
	if(method%in%c('Prentice','prentice')){
		
		r<-prentice.SO3(Rs=Rs,a=alpha)
		return(r)
		
	}else	if(method%in%c('Zhang','zhang')){
		
		r<-zhang.SO3(Rs=Rs,a=alpha,...)
		
		return(r)
		
	}else	if(method%in%c('Fisher','fisher')){
		
		r<-fisher.SO3(Rs=Rs,a=alpha,...)
		
		return(r)
		
	}else	if(method%in%c('Chang','chang')){
		
		r<-chang.SO3(Rs=Rs,a=alpha)
		
		return(r)
		
	}else{
		
		stop("Only the Prentice, Zhang, Chang and Fisher options are currently available")
		
	}
	
}

#' Prentice confidence region method
#'
#' Find the radius of a \eqn{100(1-\alpha)%} confidence region for the projected mean
#'
#' Compute the radius of a \eqn{100(1-\alpha)%} confidence region for the central orientation based on the projected mean
#' estimator using the method due to \cite{prentice1986}.  For a rotation specific version see \cite{rancourt2000}. The variablity
#' in each axis is different so each axis will have its own radius.  In \cite{bingham09} they take the largest radius and use it to
#' form regions that are symmetric about each axis.
#'
#' @param Rs,Qs A \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion form (p=4)
#' @param alpha The alpha level desired, e.g. 0.05 or 0.10
#' @return Radius of the confidence region centered at the projected mean for each of the x-, y- and z-axis
#' @seealso \code{\link{fisher}} \code{\link{chang}} \code{\link{zhang}}
#' @cite prentice1986, rancourt2000, bingham09
#' @export
#' @examples
#' Qs<-ruars(20,rcayley,kappa=100,space='Q4')
#' region(Qs,method='prentice',alpha=0.1)

prentice<-function(Qs,alpha){
	UseMethod("prentice")
}


#' @rdname prentice
#' @method prentice Q4
#' @S3method prentice Q4

prentice.Q4<-function(Qs,alpha){
	#This takes a sample qs and returns the radius of the confidence region
	#centered at the projected mean
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
	
	r<-sqrt(qchisq((1-alpha),3)/Tm)
	return(r)
}


#' @rdname prentice
#' @method prentice SO3
#' @S3method prentice SO3

prentice.SO3<-function(Rs,alpha){
	Qs<-Q4(Rs)
	r<-prentice.Q4(Qs,alpha)
	return(r)
}

#' Zhang confidence region method
#'
#' Compute the radius of a \eqn{100(1-\alpha)%} confidence region for the central orientation
#' 
#' Compute the radius of a \eqn{100(1-\alpha)%} confidence region for the central orientation based on the projected mean
#' estimator using the method due to Zhang & Nordman (2009) (unpublished MS thesis).  By construction each axis will have the same
#' radius so the radius reported is for all three axis.  A normal theory version of this procedure uses the theoretical
#' chi-square limiting distribution and is given by the \code{\link{chang}} option.
#'
#' @param Rs,Qs A \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion form (p=4)
#' @param alpha The alpha level desired, e.g. 0.05 or 0.10
#' @param m Number of replicates to use to estiamte cut point
#' @return Radius of the confidence region centered at the projected mean
#' @seealso \code{\link{prentice}} \code{\link{fisher}} \code{\link{chang}}
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=100)
#' region(Rs,method='zhang',alpha=0.1)

zhang<-function(Qs,alpha,m){
	UseMethod("zhang")
}


#' @rdname zhang
#' @method zhang SO3
#' @S3method zhang SO3

zhang.SO3<-function(Rs,alpha,m=300){
	
	#Rs is a n-by-9 matrix where each row is an 3-by-3 rotation matrix
	#m is the number of resamples to find q_1-a
	#alpha is the level of confidence desired, e.g. 0.95 or 0.90
	#pivot logical; should the pivotal (T) bootstrap be used or nonpivotal (F)
	
  Rs<-formatSO3(Rs)
  Qs<-Q4(Rs)
  rad<-zhang.Q4(Qs,alpha,m)
	return(rad)
}

#' @rdname zhang
#' @method zhang Q4
#' @S3method zhang Q4

zhang.Q4<-function(Qs,alpha,m=300){
	
	Qs<-formatQ4(Qs)
	n<-nrow(Qs)
  stats<-zhangQ4(Qs,m)
	Shat<-mean(Qs)
  cdhat<-cdfuns(Qs,Shat)
  
	rad<-sqrt(as.numeric(quantile(stats,1-alpha))*cdhat$c/(2*n*cdhat$d^2))
	
	return(rad)
}


cdfuns<-function(Qs,Shat){
  
  Shat<-matrix(Shat,4,1)
	cd<-cdfunsC(Qs,Shat)
	
	return(list(c=cd[1],d=cd[2]))
}


#' Fisher confidence region method
#'
#' Find the radius of a \eqn{100(1-\alpha)%} confidence region for the central orientation
#'
#' Compute the radius of a \eqn{100(1-\alpha)%} confidence region for the central orientation based on the projected mean
#' estimator using the method for the mean polar axis as proposed in \cite{fisher1996}.  To be able to reduce their method
#' to a radius requires the additonal assumption of rotational symmetry, equation (10) in \cite{fisher1996}. 
#'
#' @param Rs,Qs A \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion form (p=4)
#' @param alpha The alpha level desired, e.g. 0.05 or 0.10
#' @param boot Should the bootstrap or normal theory critical value be used
#' @param m number of bootstrap replicates to use to estimate critical value
#' @param symm true/false on if rotationally symmetric regions should be computed or not
#' @return radius of the confidence region centered at the projected mean
#' @seealso \code{\link{prentice}} \code{\link{chang}} \code{\link{zhang}}
#' @cite fisher1996
#' @export
#' @examples
#' Qs<-ruars(20,rcayley,kappa=100,space='Q4')
#' region(Qs,method='fisher',alpha=0.1,symm=T)

fisher<-function(Qs,alpha,boot,m,symm){
	UseMethod("fisher")
}


#' @rdname fisher
#' @method fisher Q4
#' @S3method fisher Q4

fisher.Q4<-function(Qs,alpha,boot=T,m=300,symm=T){
	
	Qs<-formatQ4(Qs)
	
	if(boot){
    
	  Tstats <- fisherBootC(Qs,m,symm)
    
		qhat<-as.numeric(quantile(Tstats,1-alpha))
		
	}else{
		
		qhat<-qchisq(a,3)
		
	}
	
	rsym<-optim(.05,optimAxis,Qs=Qs,cut=qhat,symm=T,method='Brent',lower=0,upper=pi)$par
	
	return(rsym)
}


optimAxis<-function(r,Qs,cut,symm){
	
	Shat<-Q4(axis2(mean(Qs)),r)
	if(symm){
		Tm<-fisherAxisC(Qs,Shat)
	}else{
		Tm<-fisherAxisCSymmetric(Qs,Shat)
	}
	return((Tm-cut)^2)
}


#' @rdname fisher
#' @method fisher SO3
#' @S3method fisher SO3

fisher.SO3<-function(Rs,alpha,boot=T,m=300,symm=T){
	
	Qs<-Q4(Rs)
	r<-fisher.Q4(Qs,alpha,boot,m,symm)
	
	return(r)
}

#' Chang and Rivest confidence region method
#'
#' Compute the radius of a \eqn{100(1-\alpha)%} confidence region for the central orientation
#' 
#' Compute the radius of a \eqn{100(1-\alpha)%} confidence region for the central orientation based on the projected mean
#' estimator based on a result due to \cite{chang2001}.  By construction each axis will have the same
#' radius so the radius reported is for all three axis.
#'
#' @param Rs,Qs A \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion form (p=4)
#' @param alpha The alpha level desired, e.g. 0.05 or 0.10
#' @return Radius of the confidence region centered at the projected mean
#' @cite chang2001
#' @seealso \code{\link{prentice}} \code{\link{fisher}} \code{\link{zhang}}
#' @export
#' @examples
#' Rs<-ruars(20,rcayley,kappa=100)
#' region(Rs,method='chang',alpha=0.1)

chang<-function(Qs,alpha){
	UseMethod("chang")
}


#' @rdname chang
#' @method chang SO3
#' @S3method chang SO3

chang.SO3<-function(Rs,alpha){
	
	#Rs is a n-by-9 matrix where each row is an 3-by-3 rotation matrix
	#alpha is the level of confidence desired, e.g. 0.95 or 0.90
	#pivot logical; should the pivotal (T) bootstrap be used or nonpivotal (F)
	
	Rs<-formatSO3(Rs)
	Qs<-Q4(Rs)
	rad<-chang.Q4(Qs,alpha)
	return(rad)
}

#' @rdname chang
#' @method chang Q4
#' @S3method chang Q4

chang.Q4<-function(Qs,alpha){
	
	Qs<-formatQ4(Qs)
	n<-nrow(Qs)
	Shat<-mean(Qs)
	cdhat<-cdfuns(Qs,Shat)
	
	rad<-sqrt(as.numeric(qchisq(1-alpha,3))*cdhat$c/(2*n*cdhat$d^2))
	
	return(rad)
}
