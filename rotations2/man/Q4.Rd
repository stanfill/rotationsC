\name{Q4}
\alias{as.Q4}
\alias{is.Q4}
\alias{Q4}
\alias{Q4.default}
\alias{Q4.SO3}
\title{Quaternions}
\usage{
  Q4(q, ...)

  as.Q4(q)

  is.Q4(q)

  \method{Q4}{default} (U, theta = NULL)

  \method{Q4}{SO3} (R)
}
\arguments{
  \item{q}{object to be coerced or tested}

  \item{U}{n-by-3 matrix where rows represent axes of
  rotation}

  \item{theta}{vector of rotation angles}

  \item{R}{n-by-9 matrix of rotation matrices to be
  translated into quaternions}

  \item{...}{additional arguments}
}
\value{
  unit quaternion of class Q4
}
\description{
  Creates or tests for objects of type "Q4"
}
\details{
  Construct a unit quaternion to represent a rotation.
  Each quaternion can be interpreted as a rotation of some
  reference frame about the axis U (of unit length) through
  the angle theta.  For each axis, U, and angle theta the
  quaternion is formed through
  \deqn{q=[cos(\theta/2),sin(\theta/2)U]^\top.}{q=[cos(theta/2),sin(theta/2)U]'.}
  If the theta element is left empty then the
  Frobenius-norm of each axis is taken to the the angle
  theta.  If an \code{\link{SO3}} object is given then this
  function will return the quaternion equivalent.
}
\seealso{
  Other Q4: \code{\link{id.Q4}}, \code{\link{Q4-class}}

  Other Q4: \code{\link{id.Q4}}, \code{\link{Q4-class}}

  Other Q4: \code{\link{id.Q4}}, \code{\link{Q4-class}}

  Other Q4: \code{\link{id.Q4}}, \code{\link{Q4-class}}

  Other Q4: \code{\link{id.Q4}}, \code{\link{Q4-class}}
}

