\docType{data}
\name{as.Q4}
\alias{Q4}
\alias{as.Q4}
\alias{as.Q4.Q4}
\alias{as.Q4.SO3}
\alias{as.Q4.data.frame}
\alias{as.Q4.default}
\alias{id.Q4}
\alias{is.Q4}
\title{Quaternions}
\format{\code{id.Q4} is the identity rotation given by the matrix \eqn{[1,0,0,0]^\top}{[1,0,0,0]'}.}
\usage{
  as.Q4(q, ...)

  \method{as.Q4}{default} (x, theta = NULL, ...)

  \method{as.Q4}{SO3} (x, ...)

  \method{as.Q4}{Q4} (x, ...)

  \method{as.Q4}{data.frame} (x, ...)

  is.Q4(x)

  id.Q4
}
\arguments{
  \item{x}{object to be coerced or tested.}

  \item{theta}{vector or single rotation angle; if
  \code{length(theta)==1}, the same theta is used for all
  axes}

  \item{...}{additional arguments.}
}
\value{
  \item{as.Q4}{coerces its object into an Q4 type.}
  \item{is.Q4}{returns \code{TRUE} or \code{False}
  depending on whether its argument satifies the conditions
  to be an quaternion; namely it must be four-dimensional
  and of unit length.}
}
\description{
  Creates or tests for objects of class "Q4."
}
\details{
  Construct a unit quaternion to represent a rotation.
  Each quaternion can be interpreted as a rotation of some
  reference frame about the axis \eqn{U} (of unit length)
  through the angle \eqn{\theta}.  For each axis and angle
  the quaternion is formed through
  \deqn{q=[cos(\theta/2),sin(\theta/2)U]^\top.}{q=[cos(theta/2),sin(theta/2)U]'.}
  If no angle is supplied then the length of each axis is
  taken to be the angle of rotation theta.  If an
  \code{\link{SO3}} object is given then this function will
  return the quaternion equivalent.  If a \code{data.frame}
  is supplied with n rows and 4 columns where each row is a
  quaternion, then an object of class \code{'Q4'} is
  returned where each row represents a quaternion.  As
  demonstrated below, \code{is.Q4} may return \code{TRUE}
  for a data frame, but the functions defined for objects
  of class \code{'Q4'} will not be called until
  \code{as.Q4} has been used.  If an n-by-4 matrix is
  supplied then rows are treated as vectors in \eqn{R^4};
  rows that aren't unit vectors (quaternions) are
  normalized to have length one and those that are already
  of unit length are returned untouched.
}
\examples{
data(drill)            #Pull off subject 1's wrist measurements
Qs <- subset(drill, Subject == '1' & Joint == 'Wrist')

                       ## The measurements are in columns 5:8
is.Q4(Qs[,5:8])        #TRUE, even though Qs is a data.frame, the rows satisfy the
                       #conditions necessary to be quaternions BUT, S3 methods (e.g. 'mean' or
                       #'plot') for objects of class 'Q4' will not work until 'as.Q4' is used

Qs <- as.Q4(Qs[,5:8])  #Coerce measurements into 'Q4' type using as.Q4.data.frame
all(is.Q4(Qs))         #TRUE
mean(Qs)               #Estimate central orientation for subject 1's wrist, see ?mean.Q4
Rs <- as.SO3(Qs)       #Coerse a 'Q4' object into rotation matrix format, see ?as.SO3

#Visualize the measuremenets, see ?plot.Q4 for more
\dontrun{
plot(Qs, col = c(1, 2, 3))}
}
\keyword{datasets}

