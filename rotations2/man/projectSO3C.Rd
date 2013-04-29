\name{projectSO3C}
\alias{projectSO3C}
\title{Projection Procedure}
\usage{
  projectSO3C(M)
}
\arguments{
  \item{M}{\eqn{3\times 3}{3-by-3} matrix to project}
}
\value{
  projection of \eqn{\bm M}{M} into \eqn{SO(3)}
}
\description{
  Project an arbitrary \eqn{3\times 3}{3-by-3} matrix into
  SO(3)
}
\details{
  This function uses the process given in \cite{moakher02}
  to project an arbitrary \eqn{3\times 3}{3-by-3} matrix
  into \eqn{SO(3)}.
}
\examples{
M<-matrix(rnorm(9),3,3)
project.SO3(M)
}
\references{
  moakher02
}
\seealso{
  \code{\link{mean.SO3}}, \code{\link{median.SO3}}
}

