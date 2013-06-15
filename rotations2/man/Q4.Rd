\name{Q4}
\alias{Q4}
\alias{Q4.default}
\alias{Q4.SO3}
\title{Quaternions}
\usage{
  Q4(U, ...)

  \method{Q4}{default} (U, theta)

  \method{Q4}{SO3} (R)
}
\arguments{
  \item{U}{three-dimensional vector describing the axis of
  rotation}

  \item{theta}{vector of angles to create the matrices}

  \item{R}{matrix in SO(3) to be translated into
  quaternions}

  \item{...}{additional arguments}
}
\value{
  unit quaternion of class "Q4"
}
\description{
  Create a unit quaternion
}
\details{
  Construct a unit quaternion to represent a rotation.  It
  is interpreted as a rotation of all three axis of the
  identity matrix about the axis U (of unit length) through
  the angle theta.  Alternatively, if U is not of unit
  length, the length of U is taken to be angle of rotation,
  theta.  Given a matrix representation of a rotation this
  function will translate the matrix into a quaternion.
}
\seealso{
  Other Q4: \code{\link{as.Q4}}, \code{\link{id.Q4}},
  \code{\link{Q4-class}}

  Other Q4: \code{\link{as.Q4}}, \code{\link{id.Q4}},
  \code{\link{Q4-class}}

  Other Q4: \code{\link{as.Q4}}, \code{\link{id.Q4}},
  \code{\link{Q4-class}}
}

