\name{SO3}
\alias{SO3}
\alias{SO3.default}
\alias{SO3.Q4}
\title{Matrix in SO(3)}
\usage{
  SO3(U, ...)

  \method{SO3}{default} (U, theta = NULL)

  \method{SO3}{Q4} (q)
}
\arguments{
  \item{U}{three-dimensional vector describing the axis of
  rotation}

  \item{theta}{vector of angles to create the matrices}

  \item{q}{quaternions to be translated into rotations}

  \item{...}{additional arguments}
}
\value{
  matrix of rotations in SO3 format
}
\description{
  Rotation matrix construction
}
\details{
  Construct a 3-by-3 non-singular orthogonal matrix with
  determinant one.  The construction is accodring to
  Rodrigues formula.  It is interpreted as a rotation of
  all three axis of the identity matrix about the axis U
  (of unit length) through the angle theta.  Alternatively,
  if U is not of unit length, the length of U is taken to
  be angle of rotation, theta.  Given a quaternion
  representation of a rotation this function will translate
  the quaternion into a rotation matrix
}
\seealso{
  Other SO3: \code{\link{as.SO3}}, \code{\link{id.SO3}},
  \code{\link{is.SO3}}, \code{\link{SO3-class}}

  Other SO3: \code{\link{as.SO3}}, \code{\link{id.SO3}},
  \code{\link{is.SO3}}, \code{\link{SO3-class}}

  Other SO3: \code{\link{as.SO3}}, \code{\link{id.SO3}},
  \code{\link{is.SO3}}, \code{\link{SO3-class}}
}

