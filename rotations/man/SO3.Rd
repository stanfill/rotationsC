\docType{data}
\name{SO3}
\alias{as.SO3}
\alias{id.SO3}
\alias{is.SO3}
\alias{SO3}
\alias{SO3.default}
\alias{SO3.Q4}
\title{Rotation matrices}
\format{\code{id.SO3} is the identity rotation given by the the 3-by-3 identity matrix.}
\usage{
  SO3(R, ...)

  \method{SO3}{default} (R, theta = NULL, ...)

  \method{SO3}{Q4} (R, ...)

  as.SO3(R)

  is.SO3(R)

  id.SO3
}
\arguments{
  \item{R}{object to be coerced or tested.}

  \item{theta}{vector of rotation angles.}

  \item{...}{additional arguments.}
}
\value{
  \item{as.SO3}{coerces its object into an SO3 type.}
  \item{is.SO3}{returns \code{TRUE} or \code{False}
  depending on whether its argument satifies the conditions
  to be an rotation matrix.  Namely, has determinant one
  and its transpose is its inverse.}
  \item{SO3.default}{returns an n-by-9 matrix where each
  row is a rotation matrix constructed from axis U and
  angle theta.} \item{SO3.Q4}{returns n-by-9 matrix where
  each row is a rotation matrix constructed from the
  corresponding quaternion.}
}
\description{
  Creates or tests for objects of class "SO3."
}
\details{
  Construct a 3-by-3 matrix to represent a rotation.  Each
  rotation matrix can be interpreted as a rotation of some
  reference frame about the axis U (of unit length) through
  the angle theta.  For each axis, U, and angle theta the
  matrix is formed through
  \deqn{R=exp[\Phi(Ur)]}{R=exp[\Phi(Ur)].} If the theta
  element is left empty then the length of each axis is
  taken to be the angle of rotation theta.  If a
  \code{\link{Q4}} object is given then this function will
  return the rotation matrix equivalent.
}
\keyword{datasets}

