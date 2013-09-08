\name{exp_skew}
\alias{exp_skew}
\title{Matrix exponential}
\usage{
  exp_skew(H)
}
\arguments{
  \item{H}{singular or sample of \eqn{3\times 3}{3-by-3}
  skew-symmetric matrices}
}
\value{
  matrix in SO(3) \eqn{e^{\bm A}}{e^A}
}
\description{
  Compute the matrix exponent for skew-symmetric matrices
  according to the usual Taylor expansion. The expansion is
  significantly simplified for skew-symmetric matrices, see
  \cite{Moakher (2002)}. Maps a matrix belonging to the lie
  algebra so(3) into the lie group SO(3).
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." _SIAM Journal on Matrix Analysis and
  Applications_, *24*(1), pp. 1-16.
}

