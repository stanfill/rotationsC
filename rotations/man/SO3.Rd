\docType{data}
\name{as.SO3}
\alias{SO3}
\alias{as.SO3}
\alias{as.SO3.Q4}
\alias{as.SO3.SO3}
\alias{as.SO3.data.frame}
\alias{as.SO3.default}
\alias{id.SO3}
\alias{is.SO3}
\title{Rotation matrices}
\format{\code{id.SO3} is the identity rotation given by the the 3-by-3 identity matrix.}
\usage{
  as.SO3(R, ...)

  \method{as.SO3}{default} (R, theta = NULL, ...)

  \method{as.SO3}{Q4} (R, ...)

  \method{as.SO3}{SO3} (R, ...)

  \method{as.SO3}{data.frame} (R, ...)

  is.SO3(R)

  id.SO3
}
\arguments{
  \item{R}{object to be coerced or tested.}

  \item{theta}{vector or single rotation angle; if
  \code{length(theta)==1} the same theta is used for all
  axes}

  \item{...}{additional arguments.}
}
\value{
  \item{as.SO3}{coerces provided data into an SO3 type.}
  \item{is.SO3}{returns \code{TRUE} or \code{False}
  depending on whether its argument satifies the conditions
  to be an rotation matrix.  Namely, has determinant one
  and its transpose is its inverse.}
}
\description{
  Creates or tests for objects of class "SO3."
}
\details{
  Construct a 3-by-3 matrix to represent a rotation.  Each
  rotation matrix can be interpreted as a rotation of some
  reference frame about the axis \eqn{U} (of unit length)
  through the angle \eqn{\theta}.  For each axis and angle
  the matrix is formed through
  \deqn{R=\exp[\Phi(U\theta)]}{R=exp[\Phi(U\theta)].} If no
  angle of rotation are supplied then the length of each
  axis is taken to be the angle of rotation theta.  If a
  \code{\link{Q4}} object is given then this function will
  return the rotation matrix equivalent.  If a
  \code{data.frame} is supplied with n rows and 9 columns
  where each row is a quaternion, then an object of class
  \code{'SO3'} is returned where each row represents a
  3-by-3 rotation matrix.  As demonstrated below,
  \code{is.SO3} may return \code{TRUE} for a data frame,
  but the functions defined for objects of class
  \code{'SO3'} will not be called until \code{as.SO3} has
  been used. If an n-by-9 matrix is supplied then rows are
  treated as 3-by-3 matrices; rows that don't form matrices
  in SO(3) are projected into SO(3) and those that are
  already in SO(3) are returned untouched.
}
\examples{
data(nickel)                   #Select one location to focus on
Rs <- subset(nickel, location == 698)
is.SO3(Rs[,5:13])              #Some of the rows are not rotations due to rounding or entry errors
                               #as.SO3 will project matrices not in SO(3) to the closest rotation in it

Rs <- as.SO3(Rs[,5:13])        #Translate the Rs data.frame into an object of class 'SO3'
mean(Rs)                       #Estimate the central orientation with the average
median(Rs)                     #Re-estimate central orientation robustly
Qs <- as.Q4(Rs)                #Coerse into "SO3" format, see ?as.SO3 for more

#Visualize the location, there appears to be two groups
\dontrun{
plot(Rs, col = c(1, 2, 3))}
}
\keyword{datasets}

