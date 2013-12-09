#' SO3 class.
#'
#' Class for \eqn{3\times 3}{3-by-3} matrix representation of rotations.
#'
#' @name SO3-class
#' @seealso See the \code{\link{SO3}} functions.
#'
#' @exportClass SO3
setOldClass("SO3")


#' Q4 class.
#'
#' Class for quaterion representation of rotations.
#' 
#' @name Q4-class
#' @seealso See the \code{\link{Q4}} functions.
#'
#' @exportClass Q4
setOldClass("Q4")

#' Quaternions
#' 
#' Creates or tests for objects of class "Q4."
#' 
#' Construct a unit quaternion to represent a rotation.  Each quaternion can be interpreted as a rotation of some reference frame 
#' about the axis \eqn{U} (of unit length) through the angle \eqn{\theta}.  For each axis and angle the quaternion is formed through
#' \deqn{q=[cos(\theta/2),sin(\theta/2)U]^\top.}{q=[cos(theta/2),sin(theta/2)U]'.}  If no angle is supplied then the 
#' length of each axis is taken to be the angle of rotation theta.  If an \code{\link{SO3}} object is given then this function will
#' return the quaternion equivalent.
#'
#' @export
#' @param q object to be coerced or tested.
#' @param theta vector of rotation angles.
#' @param ... additional arguments.
#' @format \code{id.Q4} is the identity rotation given by the matrix \eqn{[1,0,0,0]^\top}{[1,0,0,0]'}.
#' @return 	\item{as.Q4}{coerces its object into an Q4 type.} 
#' 					\item{is.Q4}{returns \code{TRUE} or \code{False} depending on whether its argument satifies the conditions to be an
#' 					quaternion; namely it must be four-dimensional and of unit length.}
#' 					\item{Q4.default}{returns an \eqn{n}-by-4 matrix where each row is a quaternion constructed from axis \eqn{U} and angle theta.}
#' 					\item{Q4.SO3}{returns \eqn{n}-by-4 matrix where each row is a quaternion constructed from the corresponding rotation matrix.}
#' @aliases Q4 as.Q4 is.Q4 id.Q4 Q4.default Q4.SO3 Q4.Q4

Q4<-function(q,...){
  UseMethod("Q4")
}

#' @rdname Q4
#' @method Q4 default
#' @S3method Q4 default
#' @aliases Q4 as.Q4 is.Q4 id.Q4 Q4.default Q4.SO3 Q4.Q4
#' @export

Q4.default <- function(q,theta=NULL,...){  
  
  U<-q
  
  n<-length(U)/3
  
  if(n%%1!=0)
    stop("Each axis must be in three-dimensions")	 
  
  U<-matrix(U,n,3)
  ulen<-sqrt(rowSums(U^2))
  
  if(is.null(theta)){ 
    theta<-ulen%%pi
  }
  
  ntheta<-length(theta)
  
  if(n!=ntheta)
    stop("Number of angles must match number of axes")
  
  #if(any(ulen!=1))
  #  U<-U/ulen
  
  #CPP version is causing seg faults, try just doing it in R
  #x <- Q4defaultC(U,theta)
  
  x <- cbind(cos(theta/2), sin(theta/2) * U)
  class(x)<-"Q4"
  return(x)
}

#' @rdname Q4
#' @method Q4 SO3
#' @S3method Q4 SO3
#' @aliases Q4 as.Q4 is.Q4 id.Q4 Q4.default Q4.SO3 Q4.Q4
#' @export

Q4.SO3 <- function(q,...) {
  
  R<-q
  R<-formatSO3(R)
  theta <- angle(R)
  u <- axis(R)
  x <- Q4(u,theta)
  
  return(x)
}

#' @rdname Q4
#' @method Q4 Q4
#' @S3method Q4 Q4
#' @aliases Q4 as.Q4 is.Q4 id.Q4 Q4.default Q4.SO3 Q4.Q4
#' @export

Q4.Q4 <- function(q,...) {
  
  return(q)
}

#' @rdname Q4
#' @aliases Q4 as.Q4 is.Q4 id.Q4 Q4.default Q4.SO3 Q4.Q4
#' @export

as.Q4<-function(q){
	class(q)<-"Q4"
	return(q)
}

#' @rdname Q4
#' @aliases Q4 as.Q4 is.Q4 id.Q4 Q4.default Q4.SO3 Q4.Q4
#' @export

is.Q4 <- function(q) {
	
	return(sum(q^2)-1<10e-10 & length(q)==4)
	
}

#' @rdname Q4
#' @aliases Q4 as.Q4 is.Q4 id.Q4 Q4.default Q4.SO3 Q4.Q4
#' @export
#' 
id.Q4 <- as.Q4(matrix(c(1,0,0,0),1,4))


#' Rotation matrices
#' 
#' Creates or tests for objects of class "SO3."
#' 
#' Construct a 3-by-3 matrix to represent a rotation.  Each rotation matrix can be interpreted as a rotation of some reference frame 
#' about the axis \eqn{U} (of unit length) through the angle \eqn{\theta}.  For each axis and angle the matrix is formed through
#' \deqn{R=\exp[\Phi(U\theta)]}{R=exp[\Phi(U\theta)].}  If no angle of rotation are supplied then the 
#' length of each axis is taken to be the angle of rotation theta.  If a \code{\link{Q4}} object is given then this function will
#' return the rotation matrix equivalent.
#'
#' @export
#' @param R object to be coerced or tested.
#' @param theta vector of rotation angles.
#' @param ... additional arguments.
#' @format \code{id.SO3} is the identity rotation given by the the 3-by-3 identity matrix.
#' @return 	\item{as.SO3}{coerces its object into an SO3 type.} 
#' 					\item{is.SO3}{returns \code{TRUE} or \code{False} depending on whether its argument satifies the conditions to be an
#' 					rotation matrix.  Namely, has determinant one and its transpose is its inverse.}
#' 					\item{SO3.default}{returns an \eqn{n}-by-9 matrix where each row is a rotation matrix constructed from axis \eqn{U} and angle theta.}
#' 					\item{SO3.Q4}{returns \eqn{n}-by-9 matrix where each row is a rotation matrix constructed from the corresponding quaternion.}
#' @aliases SO3 as.SO3 is.SO3 id.SO3 SO3.default SO3.Q4 SO3.SO3

SO3 <- function(R,...){
  UseMethod("SO3")
}

#' @rdname SO3
#' @method SO3 default
#' @S3method SO3 default
#' @aliases SO3 as.SO3 is.SO3 id.SO3 SO3.default SO3.Q4 SO3.SO3
#' @export


SO3.default <- function(R, theta=NULL,...) {
  
  U<-R
  n<-length(U)/3
  
  if(n%%1!=0)
    stop("This functions only works in three dimensions.")	
  
  U<-matrix(U,n,3)
  
  ulen<-sqrt(rowSums(U^2)) 
  
  if(is.null(theta)){ 
    theta<-ulen%%(pi)
    
    #if(theta>pi)
    #	theta<-2*pi-theta
  }
  
  R<-matrix(NA,n,9)
  #print(U)
  #print(ulen)
  for(i in 1:n){
    
    if(ulen[i]!=0)
      U[i,]<-U[i,]/ulen[i]
    
    P <- U[i,] %*% t(U[i,])
    
    R[i,] <- P + (diag(3) - P) * cos(theta[i]) + eskew(U[i,]) * sin(theta[i])
  }
  
  class(R) <- "SO3"
  return(R)
}

# C++ version still isn't working, comment out for now
# SO3.default <- function(U, theta=NULL) {
#   
# 	n<-length(U)/3
# 	
# 	if(n%%1!=0)
# 		stop("Each axis must be in three-dimensions")
# 	
# 	U<-matrix(U,n,3)
# 	ulen<-sqrt(rowSums(U^2)) 
# 	
# 	if(is.null(theta)){ 
# 		theta<-ulen%%pi
# 	}
# 	
# 	ntheta<-length(theta)	
# 	
# 	if(n!=ntheta)
# 		stop("Number of angles must match number of axes")
# 	
# 	if(any(ulen!=1))
# 		U<-U/ulen
# 
# 	R<-SO3defaultC(U,theta)
#  		
#  	class(R) <- "SO3"
#   return(R)
# }


#' @rdname SO3
#' @method SO3 Q4
#' @S3method SO3 Q4
#' @aliases SO3 as.SO3 is.SO3 id.SO3 SO3.default SO3.Q4 SO3.SO3
#' @export

SO3.Q4<-function(R,...){
  
  q<-R
  q<-formatQ4(q)
  
  if(any((rowSums(q^2)-1)>10e-10)){
    warning("Unit quaternions required.  Input was normalized.")
    nonq<-which((rowSums(q^2)-1)>10e-10)
    q[nonq,]<-as.Q4(q[nonq,]/sqrt(rowSums(q[nonq,]^2)))
  }else{
    q<-as.Q4(q)
  }
  
  theta<-angle(q)
  
  u<-axis(q)
  
  return(SO3(u, theta)) 
}

#' @rdname SO3
#' @method SO3 SO3
#' @S3method SO3 SO3
#' @aliases SO3 as.SO3 is.SO3 id.SO3 SO3.default SO3.Q4 SO3.SO3
#' @export

SO3.SO3<-function(R,...){
  return(R)
}

#' @rdname SO3
#' @aliases SO3 as.SO3 is.SO3 id.SO3 SO3.default SO3.Q4 SO3.SO3
#' @export

as.SO3<-function(R){
	class(R)<-"SO3"
	return(R)
}

#' @rdname SO3
#' @aliases SO3 as.SO3 is.SO3 id.SO3 SO3.default SO3.Q4 SO3.SO3
#' @export

is.SO3 <- function(R) {
	
	R <- matrix(R, 3, 3)
	if (any(is.na(R))) return(FALSE)
	
	return(all(sum(t(R) %*% R - diag(1, 3))<10e-10)) 
	
}

#' @rdname SO3
#' @aliases SO3 as.SO3 is.SO3 id.SO3 SO3.default SO3.Q4 SO3.SO3
#' @export

id.SO3 <- as.SO3(diag(c(1,1,1)))
