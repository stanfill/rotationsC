

arsample <- function(f, g, M, kappa, Haar, ...) {
  #generate a random observation from target density f
  found = FALSE
  while (!found) {
    x <- g(1, ...)
    y <- runif(1, min = 0, max = M)
    if (y < f(x, kappa, Haar)) 
      found = TRUE
  }
  return(x)
  # arsample(f, g, M, kappa, ...)
}



arsample.unif <- function(f, M, ...) {
  #generate a random observation from target density f assuming g is uniform
  found = FALSE
  while (!found) {
    x <- runif(1, -pi, pi)
    y <- runif(1, min = 0, max = M)
    if (y < f(x, ...)) 
      found = TRUE
  }
  return(x)
  # arsample.unif(f, M, ...)
}


#' @S3method print SO3
#' @method print SO3

print.SO3<-function(x,...){
  Rs<-x
  len<-length(Rs)
  
  if(len%%9!=0)
    stop("Input is not of the correct length.")
  
  if(len==9){
    print.default(as.SO3(matrix(Rs,3,3)),...)
  }else{
    print.default(Rs,...)
  }
}

#' @S3method print Q4
#' @method print Q4

print.Q4<-function(x,...){
  Qs<-x
  len<-length(Qs)
  
  if(len%%4!=0)
    stop("Input is not of the correct length.")
  
  if(len==4){
    
    negs<-length(which(Qs[2:4]<0))
    
    if(negs==0){ 
      
      print.default(bquote(.(Qs[1])+.(Qs[2])*i+.(Qs[3])*j+.(Qs[4])*k),...)
      
    }else if(negs==1){
      
      if(Qs[2]<0){
        Qs[2]<--Qs[2]
        print.default(bquote(.(Qs[1])-.(Qs[2])*i+.(Qs[3])*j+.(Qs[4])*k),...)
      }else if(Qs[3]<0){
        Qs[3]<--Qs[3]
        print.default(bquote(.(Qs[1])+.(Qs[2])*i-.(Qs[3])*j+.(Qs[4])*k),...)
      }else{
        Qs[4]<--Qs[4]
        print.default(bquote(.(Qs[1])+.(Qs[2])*i+.(Qs[3])*j-.(Qs[4])*k),...)
      }
      
    }else if(negs==2){
      
      if(all(Qs[2:3]<0)){
        
        Qs[2:3]<-abs(Qs[2:3])
        print.default(bquote(.(Qs[1])-.(Qs[2])*i-.(Qs[3])*j+.(Qs[4])*k),...)
        
      }else if(all(Qs[3:4]<0)){
        
        Qs[3:4]<-abs(Qs[3:4])
        print.default(bquote(.(Qs[1])+.(Qs[2])*i-.(Qs[3])*j-.(Qs[4])*k),...)
        
      }else{
        Qs[2:4]<-abs(Qs[2:4])
        print.default(bquote(.(Qs[1])-.(Qs[2])*i+.(Qs[3])*j-.(Qs[4])*k),...)
      }
    }else{ 
      Qs[2:4]<-abs(Qs[2:4])
      print.default(bquote(.(Qs[1])-.(Qs[2])*i-.(Qs[3])*j-.(Qs[4])*k),...)
    }
    
  }else{
    colnames(Qs)<-c("Real","i","j","k")
    print.default(Qs,...)
  }
}

#print.Q4<-function(Qs,...){
  
#  len<-length(Qs)
  
#  if(len%%4!=0)
#    stop("Input is not of the correct length.")
  
#  if(len==4){
#    print.default(sprintf("%f + %f*i+ %f*j+%f*k ",Qs[1],Qs[2],Qs[3],Qs[4]),...)
#  }else{
#    colnames(Qs)<-c("Real","i","j","k")
#    print.default(Qs,...)
#  }
#}

#' Compute the rotational distance
#'
#' Calculate the Euclidean or Riemannian distance between two rotations.
#'
#' This function will calculate the intrinsic (Riemannian) or projected (Euclidean) distance between two rotations.
#' \code{R2} and \code{Q2} are set to the identity rotations by default.  For rotations \eqn{R_1}{R1} and \eqn{R_2}{R2}
#' both in \eqn{SO(3)}, the Euclidean distance between them is \deqn{||R_1-R_2||_F}{||R1-R2||} where \eqn{||\cdot||_F}{|| ||} is the Frobenius norm.
#' The Riemannian distance is defined as \deqn{||Log(R_1^\top R_2)||_F}{||Log(R1'R2)||} where \eqn{Log} is the matrix logarithm, and it corresponds
#' to the misorientation angle of \eqn{R_1^\top R_2}{R1'R2}.
#' To compute the distance matrix use \code{stats::dist()}.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion (p=4) form.
#' @param R2,Q2 the second rotation in the same parameterization as x.
#' @param method string indicating "projected" or "intrinsic" method of distance. 
#' @param p the order of the distance.
#' @param ... additional arguments.
#' @return The rotational distance between each rotation in x and R2 or Q2.
#' @export

dist<-function(x,...){
  UseMethod("dist")
}


#' @rdname dist
#' @method dist SO3
#' @S3method dist SO3

dist.SO3 <- function(x, R2=id.SO3, method='projected' , p=1,...) {
  
  R1<-formatSO3(x)
  
  if(method=='projected'){
    
    R2<-matrix(R2,nrow(R1),9,byrow=T)
    so3dist<-sqrt(rowSums((R1-R2)^2))^p
    
  }else if(method=='intrinsic'){
    
    R2<-matrix(R2,3,3)
    
    thetas<-c(rdistSO3C(R1,R2))
    
    so3dist<-thetas^p
    
  }else{
    stop("Incorrect usage of method argument.  Please choose intrinsic or projected.")
  }
  
  return(so3dist)
  
}


#' @rdname dist
#' @method dist Q4
#' @S3method dist Q4

dist.Q4 <- function(x, Q2=id.Q4 ,method='projected', p=1,...) {

  Q1<-formatQ4(x)
  Q2<-formatQ4(Q2)
  
  if(method=='intrinsic'){
    
    q4dist<-c(RdistC(Q1,Q2))^p
    
  }else if(method=='projected'){
    
    q4dist<-c(EdistC(Q1,Q2))^p
    
  }else{
    stop("Incorrect usage of method argument.  Please choose intrinsic or projected.")
  }
  
  return(q4dist)
}


#' Misorientation angle
#' 
#' Compute the misorientation angle of a rotation.
#' 
#' Every rotation can be thought of as some reference coordinate system rotated about an axis through an angle.  These quantites
#' are referred to as the misorientation axis and misorientation angle, respectively, in the material sciences literature.
#' This function returns the misorentation angle associated with a rotation assuming the reference coordinate system
#' is the identity.
#'  
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion (p=4) form.
#' @return Angle of rotation.
#' @seealso \code{\link{axis}}
#' @export

angle<-function(x){
  UseMethod("angle")
}


#' @rdname angle
#' @method angle SO3
#' @S3method angle SO3

angle.SO3 <- function(x){
	
	Rs<-formatSO3(x)
	n<-nrow(Rs)
	theta<-c(rdistSO3C(Rs,id.SO3))
  return(theta)
}


#' @rdname angle
#' @method angle Q4
#' @S3method angle Q4

angle.Q4 <- function(x){
	
  Qs<-formatQ4(x)
	n<-nrow(Qs)
	theta<-rep(0,n)
	
	for(i in 1:n){
	  theta[i]<-2*acos(Qs[i,1])
	}
	 return(theta)
}


#' Misorientation axis
#' 
#' Determine the misorientation axis of a rotation.
#' 
#' Every rotation can be interpreted as some reference coordinate system rotated about an axis through an angle.  These quantites
#' are referred to as the misorientation axis and misorientation angle, respectively, in the material sciences literature.
#' This function returns the misorentation axis associated with a rotation assuming the reference coordinate system
#' is the identity.
#' 
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion (p=4) form.
#' @param ... additional arguements.
#' @return Axis in form of three dimensional vector of length one.
#' @seealso \code{\link{angle}}
#' @export

axis<-function(x,...){
  UseMethod("axis")
}

#' @rdname axis
#' @method axis SO3
#' @S3method axis SO3

axis.SO3<-function(x,...){
  
	R<-formatSO3(x)
  n<-nrow(R)
	u<-matrix(NA,n,3)
	
	for(i in 1:n){
		Ri<-matrix(R[i,],3,3)
  	X <- Ri - t(Ri)
  	u[i,] <- rev(X[upper.tri(X)])*c(-1,1,-1)
		u[i,]<-u[i,]/sqrt(sum(u[i,]^2))
	}
  return(u) # will be trouble, if R is symmetric, i.e. id,  .... 

}

#' @rdname axis
#' @method axis Q4
#' @S3method axis Q4

axis.Q4 <- function(x,...){
  
  q<-formatQ4(x)
  theta<-angle(q)
  
  u <- q[,2:4]/sin(theta/2)

	if(any(is.infinite(u)|is.nan(u))){
		infs<-which(is.infinite(u)|is.nan(u))
		u[infs]<-0
	}  
  
  u<-matrix(u,ncol=3)
  
  return(u)
}



eskew <- function(U) {
  
  ulen<-sqrt(sum(U^2))
  
  if(ulen!=0){
    U<-U/ulen
  }
  
  u <- U[1]
  v <- U[2]
  w <- U[3]
  
  res <- matrix((-1) * c(0, -w, v, w, 0, -u, -v, u, 0), ncol = 3)
  return(res)
}



#' Generate rotations
#'
#' Generate rotations according to Rodrigues' formula.
#'
#' Given a vector \eqn{u\in R^3}{u in R^3} of length one and angle of rotation \eqn{r}, a rotation can be formed using Rodrigues' formula
#' \deqn{\cos(r)I_{3\times 3}+\sin(r)\Phi(u)+(1-\cos(r))uu^\top}{cos(r)I+sin(r)\Phi(u)+(1-cos(r))uu'} 
#' where \eqn{I_{3\times 3}}{I} is the \eqn{3\times 3}{3-by-3} identity matrix, \eqn{\Phi(u)} is a \eqn{3\times 3}{3-by-3} skew-symmetric matirix
#' with upper triangular elements \eqn{-u_3}{-u3}, \eqn{u_2}{u2} and \eqn{-u_1}{-u1} in that order.
#'
#' @param r vector of angles.
#' @param S central orientation.
#' @param space indicates the desired representation: rotation matrix "SO3" or quaternions "Q4." 
#' @return A matrix where each row is a sample point in the desired space.
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' Rs<-genR(r,space="SO3")
#' Qs<-genR(r,space="Q4")

genR <- function(r, S = NULL, space='SO3') {
  
  if(!(space %in% c("SO3","Q4")))
    stop("Incorrect space argument.  Options are: SO3 and Q4. ")
  
  n<-length(r)
  
  theta <- acos(runif(n, -1, 1))
  
  # Generate angles phi from a uniform distribution from -pi to pi
  
  phi <- runif(n, -pi, pi)
  u <- matrix(c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)),n,3)
  
  if(space=="SO3"){
  	
  	#For now the C++ code is broken, use R functions
  	#S<-matrix(S,3,3)
  	#o<-SO3defaultC(u,r)
  	#o<-genrC(r,S,1,u)
  	
  	o<-SO3(u,r)
  	
  	if(is.null(S)){
  		
  		class(o) <- "SO3"
  		return(o)
  		
  	}else{

  	  S<-formatSO3(S)
      St<-t(matrix(S,3,3))
  	  o<-center.SO3(o,St)
  	  class(o) <- "SO3"
  	  return(o)
  	
  	}
  	
  }else{
  	
  	#S<-matrix(S,1,4)
  	#q<-Q4defaultC(u,r)
  	#q<-genrC(r,S,2,u)
  	
  	q<-matrix(c(cos(r/2),sin(r/2)*u),n,4)
  	
  	if(is.null(S)){
  		
  		class(q)<-"Q4"
  		return(q)
  		
  	}else{
  	
  		S<-formatQ4(S)
  		S[2:4]<--S[2:4]
  		q<-center.Q4(q,S)
  	
  		class(q)<-"Q4"
  		return(q)
  	}
  	
  }

}


#' Matrix exponential
#'
#' Compute the matrix exponent for skew-symmetric matrices according to the usual Taylor expansion.
#' The expansion is significantly simplified for skew-symmetric matrices, see \cite{moakher02}.
#' Maps a matrix belonging to the lie algebra \eqn{so(3)} into the lie group \eqn{SO(3)}.
#'
#' @param H single \eqn{3\times 3}{3-by-3} skew-symmetric matrix or \eqn{n\times 9}{n-by-9} sample of skew-symmetric matrices.
#' @return Matrix in \eqn{SO(3)} \eqn{e^{\bm A}}{e^A}.
#' @cite moakher02
#' @export

exp_skew <- function(H) {

  if(length(H)==9){
    
    H<-matrix(H,3,3)

    return(as.SO3(expskewC(H)))
    
  }else{
    return(as.SO3(expskewCMulti(H)))
  }
}


#' Rotation logarithm
#'
#' Compute the logarithm of a rotation matrix.  The result is a \eqn{3\times 3}{3-by-3} skew-symmetric matrix.  This function maps
#' the lie group \eqn{SO(3)} into its tangent space, which is the space of all \eqn{3\times 3}{3-by-3} skew symmetric matrices,
#' which is the lie algerbra \eqn{so(3)}.  For details see e.g. \cite{moakher02}.
#'
#' @param x \eqn{n\times 9}{n-by-9} matrix where each row corresponds to a random rotation matrix.
#' @param ... additional arguements.
#' @return Skew symmetric matrix \eqn{\log(R)}{log(R)}.
#' @cite moakher02
#' @S3method log SO3
#' @method log SO3

log.SO3 <- function(x,...) {
  if(length(x)==9){
	  x<-matrix(x,3,3)
	  return(logSO3C(x))
  }else{
    return(logSO3CMulti(x))
  }
  
}

#' Projection into SO(3)
#'
#' Project an arbitrary \eqn{3\times 3}{3-by-3} matrix into \eqn{SO(3)}.
#'
#' This function uses the process given in \cite{moakher02} to project an arbitrary \eqn{3\times 3}{3-by-3} matrix into \eqn{SO(3)}.
#' More specifically it finds the closest orthogonal 3-by-3 matrix with determinant one to the provided matrix.
#' 
#' @param M \eqn{3\times 3}{3-by-3} matrix to project into \eqn{SO(3)}.
#' @return Projection of \eqn{\bm M}{M} into \eqn{SO(3)}.
#' @seealso \code{\link{mean.SO3}}, \code{\link{median.SO3}}
#' @export
#' @examples
#' M<-matrix(rnorm(9),3,3)
#' project.SO3(M)

project.SO3 <- function(M) {
  
	M<-matrix(M,3,3)
  R<-projectSO3C(M)
  return(R)
}


#' Sample distance
#'
#' Compute the sum of the \eqn{p^{th}}{pth} order distances between Rs and S.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion (p=4) form.
#' @param S the individual matrix of interest, usually an estimate of the mean.
#' @param method type of distance used method in "projected" or "intrinsic"
#' @param p the order of the distances to compute.
#' @return The sum of the pth order distance between each sample in Rs and S.
#' @seealso \code{\link{dist.SO3}}, \code{\link{dist.Q4}}
#' @export
#' @examples
#' Rs<-ruars(20,rvmises,kappa=10)
#' Sp<-mean(Rs)
#' sum_dist(Rs,S=Sp,p=2)

sum_dist<-function(x, S = genR(0, space=class(x)), method='projected', p=1){
  
  UseMethod( "sum_dist" )

}

#' @rdname sum_dist
#' @method sum_dist SO3
#' @S3method sum_dist SO3

sum_dist.SO3 <- function(x, S = id.SO3, method='projected', p=1) {

  return(sum(dist(x,S, method=method, p=p)))
  
}

#' @rdname sum_dist
#' @method sum_dist Q4
#' @S3method sum_dist Q4

sum_dist.Q4 <- function(x, S = id.Q4, method='projected', p=1) {
  
  return(sum(dist(x,S, method=method, p=p)))
  
}

#' Center rotation data
#' 
#' This function will take the sample Rs and return the sample Rs centered at
#' S, i.e. the returned sample is \eqn{S^\top R}{S'R}.  If S is the true center then
#' the projected mean should be close to the 3-by-3 identity matrix 
#' 
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (p=9) or quaternion (p=4) form.
#' @param S the rotation about which to center x.
#' @return The centered sample.
#' @export
#' @examples
#' Rs<-ruars(5,rcayley)
#' cRs<-center(Rs,mean(Rs))
#' mean(cRs) #Should be close to identity matrix

center<-function(x,S){
  
  UseMethod( "center" )
  
}

#' @rdname center
#' @method center SO3
#' @S3method center SO3

center.SO3<-function(x,S){
	#This takes a set of observations in SO3 and centers them around S
	
	Rs<-formatSO3(x)
	S<-matrix(formatSO3(S),3,3)
	
	for(i in 1:nrow(Rs)){
		Rs[i,]<-t(S)%*%matrix(Rs[i,],3,3)
	}
	return(as.SO3(Rs))
}


#' @rdname center
#' @method center Q4
#' @S3method center Q4

center.Q4<-function(x,S){
	#This takes a set of observations in Q4 and centers them around S
	Qs<-formatQ4(x)
	S<-formatQ4(S)
	S[2:4]<--S[2:4]
	
	for(i in 1:nrow(Qs)){
		Qs[i,]<-qMult(S,Qs[i,])
	}

	return(as.Q4(Qs))
}



formatSO3<-function(Rs){
	#This function will take input and format it to work with our functions
	#It also checks that the data is actually SO3 and of appropriate dimension
	
	len<-length(Rs)
	if(len%%9!=0)
		stop("Data needs to have length divisible by 9.")
	
	Rs<-matrix(Rs,len/9,9)
	
	if (!all(apply(Rs, 1, is.SO3))) 
		warning("At least one of the given observations is not in SO(3).  Use result with caution.")
	
	
	return(as.SO3(Rs))

}

formatQ4<-function(Qs){
	
  if(length(Qs)%%4!=0)
    stop("Data needs to have length divisible by 4.")
  
  Qs<-matrix(Qs,length(Qs)/4,4)
  
  if (!all(apply(Qs, 1, is.Q4))) 
  	warning("At least one of the given observations is not a unit quaternion.  Use result with caution.")
  
  
  if(length(Qs)==4)
    return(as.Q4(Qs))
  else
    return(as.Q4(Qs))
}

pMat<-function(p){
	#Make the matrix P from quaternion p according to 3.1 of Rancourt, Rivest and Asselin (2000)
	#This is one way to multiply quaternions
	Pmat<-matrix(0,4,4)
	Pmat[,1]<-p
	Pmat[,2]<-p[c(2,1,4,3)]*c(-1,1,1,-1)
	Pmat[,3]<-c(-p[3:4],p[1:2])
	Pmat[,4]<-p[4:1]*c(-1,1,-1,1)
	return(Pmat)
}

qMult<-function(q1,q2){
	#Forms quaternion product q1 x q2, i.e., rotate q2 by q1
	#This functions utilizes the 
	q1<-formatQ4(q1)
	q2<-formatQ4(q2)
	q1q2<-pMat(q1)%*%matrix(q2,4,1)
	return(formatQ4(q1q2))
}

proj<-function(u,v){
	#Project the vector v orthogonally onto the line spanned by the vector u
	num<-t(u)%*%v
	denom<-t(u)%*%u
	return(num*u/denom)
}

tLogMat <- function(x, S) {
  tra <- log.SO3(t(S) %*% matrix(x, 3, 3))
  return(as.vector(tra))
}

tLogMat2 <- function(x, S) {
  tra <- log.SO3(matrix(x, 3, 3)%*%t(S))
  return(as.vector(tra))
}

vecNorm <- function(x, S, ...) {
  n <- sqrt(length(x))
  cenX <- x - as.vector(S)
  return(norm(matrix(cenX, n, n), ...))
}
