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
#' @rdname Q4
#' @param q object to be coerced or tested.
#' @param theta vector of rotation angles.
#' @param ... additional arguments.
#' @format \code{id.Q4} is the identity rotation given by the matrix \eqn{[1,0,0,0]^\top}{[1,0,0,0]'}.
#' @return 	\item{as.Q4}{coerces its object into an Q4 type.} 
#' 					\item{is.Q4}{returns \code{TRUE} or \code{False} depending on whether its argument satifies the conditions to be an
#' 					quaternion; namely it must be four-dimensional and of unit length.}
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @examples
#' data(drill) #load the included drill data set
#' Qs <- subset(drill, Subject == '1' & Joint == 'Wrist') #Pull off subject 1's wrist measurements
#' 
#' ## The measurements are in columns 5:8
#' is.Q4(Qs[,5:8]) #TRUE, eventhough Qs is a data.frame, the rows satisfy the conditions necessary to be quaternions
#'                 #BUT, S3 methods (e.g. 'mean' or 'plot') for objects of class 'Q4' will not work until 'as.Q4' is used
#' Qs <- as.Q4(Qs[,5:8]) #Coerce measurements into 'Q4' type using as.Q4.data.frame
#' all(is.Q4(Qs)) #TRUE  
#' mean(Qs) #Estimate central orientation for subject 1's wrist, see ?mean.Q4
#' plot(Qs, col = c(1, 2, 3)) #Visualize the measuremenets, see ?plot.Q4 for more
#' Rs <- as.SO3(Qs) #Coerse a 'Q4' object into rotation matrix format, see ?as.SO3

as.Q4<-function(q,...){
  UseMethod("as.Q4")
}

#' @rdname Q4
#' @method as.Q4 default
#' @S3method as.Q4 default
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

as.Q4.default <- function(q,theta=NULL,...){  
  
  p<-ncol(q)
  n<-nrow(q)
  
  if(p==3){
    
    #If input is length 3, q is assumed to be vectors in R^3
    U<-q
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
  
    q <- cbind(cos(theta/2), sin(theta/2) * U)
    
  }else if(p==4){
    
    #If input has length divisible by 4, data are normalized and made into class "Q4"
    n<-n/4
    rowLens<-(rowSums(q^2))^0.5
    q<-q/rowLens
    
  }else if(p==9){
    
    #If input has 9 columns, q is assumed to be rotations
    q<-as.Q4.SO3(q)
    
  }else{
    
    stop("Unknown data type.  Pease see ?as.Q4 for more details.")
    
  }
  class(q)<-"Q4"
  return(q)
}

#' @rdname Q4
#' @method as.Q4 SO3
#' @S3method as.Q4 SO3
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

as.Q4.SO3 <- function(q,...) {
  
  R<-q
  R<-formatSO3(R)
  theta <- mis.angle(R)
  u <- mis.axis(R)
  x <- as.Q4.default(u,theta)
  
  return(x)
}

#' @rdname Q4
#' @method as.Q4 Q4
#' @S3method as.Q4 Q4
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

as.Q4.Q4 <- function(q,...) {
  
  return(q)
}

#' @rdname Q4
#' @method as.Q4 data.frame
#' @S3method as.Q4 data.frame
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

as.Q4.data.frame <- function(q,...) {
  n<-nrow(q)
  p<-ncol(q)
  q<-as.matrix(q,n,p)
  return(as.Q4.default(q))
}

#' @rdname Q4
#' @aliases Q4 as.Q4 is.Q4 id.Q4 Q4.default Q4.SO3 Q4.Q4 as.Q4.data.frame
#' @export

is.Q4 <- function(q) {

	apply(q,1,function(q){sum(q^2)-1<10e-10 & length(q)==4})
	
}

#' @rdname Q4
#' @aliases Q4 as.Q4 is.Q4 id.Q4 Q4.default Q4.SO3 Q4.Q4 as.Q4.data.frame
#' @export

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
#' @rdname SO3
#' @param R object to be coerced or tested.
#' @param theta vector of rotation angles.
#' @param ... additional arguments.
#' @format \code{id.SO3} is the identity rotation given by the the 3-by-3 identity matrix.
#' @return 	\item{as.SO3}{coerces provided data into an SO3 type.} 
#' 					\item{is.SO3}{returns \code{TRUE} or \code{False} depending on whether its argument satifies the conditions to be an
#' 					rotation matrix.  Namely, has determinant one and its transpose is its inverse.}
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @examples
#' data(nickel) #load included nickel dataset
#' Rs <- subset(nickel, location == 698)  #Select one location to focus on
#' 
#' Rs <- as.SO3(Rs[,5:13]) #Translate the Rs data.frame into an object of class 'SO3'
#' Rs <- Rs[is.SO3(Rs),] #Some observations are not rotations, remove them
#' mean(Rs) #Estimate the central orientation with the average
#' plot(Rs, col = c(1, 2, 3)) #Visualize the location, there appears to be two groups
#' median(Rs) #Resetimate central orientation robustly

as.SO3 <- function(R,...){
  UseMethod("as.SO3")
}

#' @rdname SO3
#' @method as.SO3 default
#' @S3method as.SO3 default
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

as.SO3.default <- function(R, theta=NULL,...) {

  p<-ncol(R)
  n<-nrow(R)
  
  if(is.null(p)){
    p<-length(R)
    n<-1
  }  
  
  if(n==3 && p==3 && is.SO3(R)){
    
    #If there are 3 rows and columns and the object is already a rotation matrix, the same rotation is returned
    class(R) <- "SO3"
    return(R)
    
  }else if(p==9){
    #If there are 9 columns, it's assumed the data are already rotation matrices so the SO3 class is appeneded and object returned
    class(R) <- "SO3"
    return(R)
  }else if(p==3){
    
  #If there are 3 columns, it's assumed the input R is the matrix of unit axes of rotations and the theta vector are the angles,
    #or the length of the axes is the angle of rotation
    
    U<-matrix(R,n,3)
  
    ulen<-sqrt(rowSums(U^2)) 
  
    if(is.null(theta)){ 
      theta<-ulen%%(pi)
    
      #if(theta>pi)
      #	theta<-2*pi-theta
    }
  
    R<-matrix(NA,n,9)

    for(i in 1:n){
    
      if(ulen[i]!=0)
        U[i,]<-U[i,]/ulen[i]
    
      P <- U[i,] %*% t(U[i,])
    
      R[i,] <- P + (diag(3) - P) * cos(theta[i]) + eskew(U[i,]) * sin(theta[i])
    }
    class(R) <- "SO3"
    return(R)
  }else if(p==4){
    
    #If there are 4 columns, it's assumed the input is an n-by-4 matrix with rows corresponding to quaternions 
    R<-as.Q4(R)
    return(as.SO3(R))
    
  }
  
  stop("Unknown data type.  Please see ?SO3 for more details.")

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
#' @method as.SO3 Q4
#' @S3method as.SO3 Q4
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

as.SO3.Q4<-function(R,...){
  
  q<-R
  q<-formatQ4(q)
  
  if(any((rowSums(q^2)-1)>10e-10)){
    warning("Unit quaternions required.  Input was normalized.")
    nonq<-which((rowSums(q^2)-1)>10e-10)
    q[nonq,]<-as.Q4(q[nonq,]/sqrt(rowSums(q[nonq,]^2)))
  }else{
    class(q)<-"Q4"
  }
  
  theta<-mis.angle(q)
  
  u<-mis.axis(q)
  
  return(as.SO3.default(u, theta)) 
}

#' @rdname SO3
#' @method as.SO3 SO3
#' @S3method as.SO3 SO3
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

as.SO3.SO3<-function(R,...){
  return(R)
}

#' @rdname SO3
#' @method as.SO3 data.frame
#' @S3method as.SO3 data.frame
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

as.SO3.data.frame <- function(R,...) {
  n<-nrow(R)
  p<-ncol(R)
  R<-as.matrix(R,n,p)
  return(as.SO3.default(R))
}

#' @rdname SO3
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

is.SO3 <- function(R) {
	
  if(length(R)==9){
    R <- matrix(R, 3, 3)
    if(any(is.na(R))) return(FALSE)
    if(abs(det(R)-1)>10e-10) return(FALSE)
    return(all(abs(t(R) %*% R - diag(1, 3))<10e-5))
  }else{
  
    apply(R,1,
	  function(R){R <- matrix(R, 3, 3)
	  if(any(is.na(R))) return(FALSE)
	  if(abs(det(R)-1)>10e-10) return(FALSE)
	  return(all(abs(t(R) %*% R - diag(1, 3))<10e-5))}) 
  }
	
}

#' @rdname SO3
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

id.SO3 <- genR(0)
