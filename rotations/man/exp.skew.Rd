\name{exp.skew}
\alias{exp.skew}
\title{Matrix exponential}
\usage{
  exp.skew(H)
}
\arguments{
  \item{H}{3-by-3 skew-symmetric matrix, i.e. A such that
  \eqn{\bm H=-\bm H^\top}}
}
\value{
  numeric matrix \eqn{e^{\bm A}}{e^A}
}
\description{
  Compute the matrix exponent for skew-symmetric matrices
  according to the usual Taylor expansion. The expansion is
  significantly simplified for skew-symmetric matrices, see
  \cite{Moakher (2002)}.
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.
}

