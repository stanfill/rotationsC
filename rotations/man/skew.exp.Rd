\name{skew.exp}
\alias{skew.exp}
\title{Matrix exponential}
\usage{
  skew.exp(x)
}
\arguments{
  \item{x}{single \eqn{3\times 3}{3-by-3} skew-symmetric
  matrix or \eqn{n\times 9}{n-by-9} sample of
  skew-symmetric matrices.}
}
\value{
  Matrix \eqn{e^{\bm H}}{e^H} in \eqn{SO(3)} .
}
\description{
  Compute the matrix exponential for skew-symmetric
  matrices according to the usual Taylor expansion. The
  expansion is significantly simplified for skew-symmetric
  matrices, see \cite{Moakher (2002)}. Maps a matrix
  belonging to the lie algebra \eqn{so(3)} into the lie
  group \eqn{SO(3)}.
}
\examples{
Rs <- ruars(20, rcayley)
lRs <- log(Rs)           #Take the matrix logarithm for rotation matrices
Rs2 <- skew.exp(lRs)     #Go back to rotation matrices
all.equal(Rs, Rs2)
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." SIAM Journal on Matrix Analysis and
  Applications, 24(1), pp. 1-16.
}

