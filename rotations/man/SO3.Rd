\docType{data}
\name{as.SO3}
\alias{as.SO3}
\alias{as.SO3.data.frame}
\alias{as.SO3.default}
\alias{as.SO3.Q4}
\alias{as.SO3.SO3}
\alias{id.SO3}
\alias{is.SO3}
\alias{SO3}
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

  \item{theta}{vector of rotation angles.}

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
  return the rotation matrix equivalent.
}
\examples{
data(nickel) #load included nickel dataset
Rs <- subset(nickel, location == 698)  #Select one location to focus on

Rs <- as.SO3(Rs[,5:13]) #Translate the Rs data.frame into an object of class 'SO3'
Rs <- Rs[is.SO3(Rs),] #Some observations are not rotations, remove them
mean(Rs) #Estimate the central orientation with the average
plot(Rs, col = c(1, 2, 3)) #Visualize the location, there appears to be two groups
median(Rs) #Resetimate central orientation robustly
}
\keyword{datasets}

