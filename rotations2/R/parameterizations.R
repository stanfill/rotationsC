#' "SO3" class
#'
#' @name SO3-class
#' @aliases SO3
#' @family SO3
#'
#' @exportClass SO3
setOldClass("SO3")


#' Q4 class
#'
#' Class for quaterion representation of rotations
#' 
#' @name Q4-class
#' @aliases Q4
#' @family Q4
#'
#' @exportClass Q4
setOldClass("Q4")

#' Quaternions
#' 
#' Create a unit quaternion
#' 
#' Construct a unit quaternion to represent a rotation.  Each quaternion can be interpreted as a rotation of some reference frame 
#' about the axis U (of unit length) through the angle theta.  Provided with a vector in three-dimensions
#' is provided then
#'
#' @export
#' @param U \eqn{n\times 3}{n-by-3} matrix where rows represent axes of rotation
#' @param theta vector of rotation angles
#' @param R matrix in SO(3) to be translated into quaternions
#' @param ... additional arguments
#' @return unit quaternion of class Q4
#' @family Q4

Q4<-function(U,...){
  UseMethod("Q4")
}


#' @rdname Q4
#' @method Q4 default
#' @S3method Q4 default
#' @family Q4

Q4.default <- function(U,theta=NULL){	
	
	n<-length(U)/3
	
	if(n%%1!=0)
		stop("This functions only works in three dimensions.")	
	
	U<-matrix(U,n,3)
	
	ulen<-sqrt(rowSums(U^2)) 
	
	if(any(ulen!=1)){
		U<-U/ulen
	}
	
	if(is.null(theta)){ 
		theta<-ulen%%pi
	}
	
	
	x <- Q4defaultC(U,theta)

	class(x)<-"Q4"
  return(x)
}

#' @rdname Q4
#' @method Q4 SO3
#' @S3method Q4 SO3
#' @family Q4

Q4.SO3 <- function(R) {
  
	R<-formatSO3(R)
 	theta <- angle(R)
 	u <- axis2(R)
 	x <- Q4(u,theta)

  return(x)
}


#' Convert anything into Q4 class
#' 
#' @param x can be anything
#' @return x with class "Q4"
#' @family Q4
#' @export

as.Q4<-function(x){
  class(x)<-"Q4"
  return(x)
}

#' Identity in Q4 space
#' @family Q4
#' @export
id.Q4 <- as.Q4(matrix(c(1,0,0,0),1,4))

#' A function to determine if a given object is in unit quaternion or not.
#'
#' @param x numeric 1-by-4  vector of length 4
#' @return logical T if the vector is a unit quaternion and false otherwise
#' @family Q4
#' @export

is.Q4 <- function(x) {

 return(sum(x^2)-1<10e-10 & length(x)==4)
	
}

#' Matrix in SO(3)
#' 
#' Rotation matrix construction
#' 
#' Construct a 3-by-3 non-singular orthogonal matrix with determinant one.  The construction is accodring to Rodrigues
#' formula.  It is interpreted as a rotation of all three axis of the identity matrix about the axis U (of unit length)
#' through the angle theta.  Alternatively, if U is not of unit length, the length of U is taken to be angle of 
#' rotation, theta.  Given a quaternion representation of a rotation this function will translate the quaternion
#' into a rotation matrix
#'
#' @export
#' @param U three-dimensional vector describing the axis of rotation
#' @param theta vector of angles to create the matrices
#' @param q quaternions to be translated into rotations
#' @param ... additional arguments
#' @return matrix of rotations in SO3 format
#' @family SO3

SO3 <- function(U,...){
  UseMethod("SO3")
}


#' @rdname SO3
#' @method SO3 default
#' @S3method SO3 default
#' @family SO3

SO3.default <- function(U, theta=NULL) {
  
	n<-length(U)/3
	
	if(n%%1!=0)
		stop("This functions only works in three dimensions.")	
	
	U<-matrix(U,n,3)
	
	ulen<-sqrt(rowSums(U^2)) 
  
	if(any(ulen!=1)){
		U<-U/ulen
	}
	
  if(is.null(theta)){ 
  	theta<-ulen%%(pi)
  }

	R<-SO3defaultC(U,theta)
 		
 	class(R) <- "SO3"
  return(R)
}


#' @rdname SO3
#' @method SO3 Q4
#' @S3method SO3 Q4
#' @family SO3

SO3.Q4<-function(q){
  
	q<-formatQ4(q)
	
  if(any((rowSums(q^2)-1)>10e-10)){
    warning("Unit quaternions required.  Input was normalized.")
    nonq<-which((rowSums(q^2)-1)>10e-10)
    q[nonq,]<-as.Q4(q[nonq,]/sqrt(rowSums(q[nonq,]^2)))
  }else{
  	q<-as.Q4(q)
  }
  
  theta<-angle(q)
  
  u<-axis2(q)
  
  return(SO3(u, theta)) 
}


#' Add SO3 class to object
#' 
#' @param x an R object
#' @return x with class "SO3"
#' @family SO3
#' @export

as.SO3<-function(x){
  class(x)<-"SO3"
  return(x)
}

#' Identity in SO(3) space
#' @family SO3
#' @export
id.SO3 <- as.SO3(diag(c(1,1,1)))


#' A function to determine if a given object is in \eqn{SO(3)} or not.
#'
#' @param x numeric 3-by-3 matrix or vector of length 9
#' @return logical T if the matrix is in SO(3) and false otherwise
#' @family SO3
#' @export
#' @examples
#' is.SO3(diag(1,3,3))
#' is.SO3(1:9)
is.SO3 <- function(x) {
  
  x <- matrix(x, 3, 3)
  if (any(is.na(x))) return(FALSE)
  
  return(all(sum(t(x) %*% x - diag(1, 3))<10e-10)) 
  
}

