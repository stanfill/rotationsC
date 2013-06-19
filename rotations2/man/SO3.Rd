\name{SO3}
\alias{as.SO3}
\alias{is.SO3}
\alias{SO3}
\alias{SO3.default}
\alias{SO3.Q4}
\title{Rotation Matrices}
\usage{
  SO3(R, ...)

  as.SO3(R)

  is.SO3(R)

  \method{SO3}{default} (U, theta = NULL)

  \method{SO3}{Q4} (q)
}
\arguments{
  \item{R}{object to be coerced or tested}

  \item{U}{n-by-3 matrix where rows represent axes of
  rotation}

  \item{theta}{vector of rotation angles}

  \item{q}{n-by-4 matrix of quaternions to be translated
  into rotations}

  \item{...}{additional arguments}
}
\value{
  matrix of rotations in SO3 format
}
\description{
  Creates or tests for objects of type "SO3"
}
\details{
  Construct a unit quaternion to represent a rotation.
  Each quaternion can be interpreted as a rotation of some
  reference frame about the axis U (of unit length) through
  the angle theta.  For each axis, U, and angle theta the
  quaternion is formed through
  \deqn{exp[\Phi(Ur)]}{q=[cos(theta/2),sin(theta/2)u]'.} If
  the theta element is left empty then the Frobenius-norm
  of each axis is taken to the the angle theta.  If an
  \code{\link{SO3}} object is given then this function will
  return the quaternion equivalent.
}
\seealso{
  Other SO3: \code{\link{id.SO3}}, \code{\link{SO3-class}}

  Other SO3: \code{\link{id.SO3}}, \code{\link{SO3-class}}

  Other SO3: \code{\link{id.SO3}}, \code{\link{SO3-class}}

  Other SO3: \code{\link{id.SO3}}, \code{\link{SO3-class}}

  Other SO3: \code{\link{id.SO3}}, \code{\link{SO3-class}}
}

