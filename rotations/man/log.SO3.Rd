\name{log.SO3}
\alias{log.SO3}
\title{Rotation logarithm}
\usage{
  \method{log}{SO3} (R)
}
\arguments{
  \item{R}{a single or sample of matrices in \eqn{SO(3)}}
}
\value{
  numeric matrix \eqn{\log(R)}{log(R)}
}
\description{
  Compute the matrix logarithm of a matrix in SO(3).  The
  result is a 3-by-3 skew-symmetric matrix.  This function
  maps the lie group SO(3) into its tangent space, which is
  the lie algerbra so(3).  For details see \cite{Moakher
  (2002)}
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.
}

