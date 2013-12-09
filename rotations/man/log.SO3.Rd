\name{log.SO3}
\alias{log.SO3}
\title{Rotation logarithm}
\usage{
  \method{log}{SO3} (x, ...)
}
\arguments{
  \item{x}{\eqn{n\times 9}{n-by-9} matrix where each row
  corresponds to a random rotation matrix.}

  \item{...}{additional arguements.}
}
\value{
  Skew symmetric matrix \eqn{\log(R)}{log(R)}.
}
\description{
  Compute the logarithm of a rotation matrix, which results
  in a \eqn{3\times 3}{3-by-3} skew-symmetric matrix.  This
  function maps the lie group \eqn{SO(3)} into its tangent
  space, which is the space of all \eqn{3\times 3}{3-by-3}
  skew symmetric matrices, the lie algerbra \eqn{so(3)}.
  For details see e.g. \cite{Moakher (2002)}.
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." SIAM Journal on Matrix Analysis and
  Applications, 24(1), pp. 1-16.
}

