\name{project.SO3}
\alias{project.SO3}
\title{Projection into SO(3)}
\usage{
  project.SO3(M)
}
\arguments{
  \item{M}{\eqn{3\times 3}{3-by-3} matrix to project into
  \eqn{SO(3)}.}
}
\value{
  Projection of \eqn{\bm M}{M} into \eqn{SO(3)}.
}
\description{
  Project an arbitrary \eqn{3\times 3}{3-by-3} matrix into
  \eqn{SO(3)}.
}
\details{
  This function uses the process detailed in Section 3.1 of
  \cite{Moakher (2002)} to project an arbitrary
  \eqn{3\times 3}{3-by-3} matrix into \eqn{SO(3)}. More
  specifically it finds the closest orthogonal 3-by-3
  matrix with determinant one to the provided matrix.
}
\examples{
M<-matrix(rnorm(9),3,3)
project.SO3(M)
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." SIAM Journal on Matrix Analysis and
  Applications, 24(1), pp. 1-16.
}
\seealso{
  \code{\link{mean.SO3}}, \code{\link{median.SO3}}
}

