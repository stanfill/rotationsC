\docType{data}
\name{Q4}
\alias{as.Q4}
\alias{id.Q4}
\alias{is.Q4}
\alias{Q4}
\alias{Q4.default}
\alias{Q4.Q4}
\alias{Q4.SO3}
\title{Quaternions}
\format{\code{id.Q4} is the identity rotation given by the matrix \eqn{[1,0,0,0]^\top}{[1,0,0,0]'}.}
\usage{
  Q4(q, ...)

  \method{Q4}{default} (q, theta = NULL, ...)

  \method{Q4}{SO3} (q, ...)

  \method{Q4}{Q4} (q, ...)

  as.Q4(q)

  is.Q4(q)

  id.Q4
}
\arguments{
  \item{q}{object to be coerced or tested.}

  \item{theta}{vector of rotation angles.}

  \item{...}{additional arguments.}
}
\value{
  \item{as.Q4}{coerces its object into an Q4 type.}
  \item{is.Q4}{returns \code{TRUE} or \code{False}
  depending on whether its argument satifies the conditions
  to be an quaternion; namely it must be four-dimensional
  and of unit length.} \item{Q4.default}{returns an
  \eqn{n}-by-4 matrix where each row is a quaternion
  constructed from axis \eqn{U} and angle theta.}
  \item{Q4.SO3}{returns \eqn{n}-by-4 matrix where each row
  is a quaternion constructed from the corresponding
  rotation matrix.}
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
  return the quaternion equivalent.
}
\keyword{datasets}

