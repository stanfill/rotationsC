\name{chang}
\alias{chang}
\alias{chang.Q4}
\alias{chang.SO3}
\title{Chang and Rivest confidence region method}
\usage{
  chang(Qs, alpha)

  \method{chang}{SO3} (Rs, alpha)

  \method{chang}{Q4} (Qs, alpha)
}
\arguments{
  \item{Rs,Qs}{A \eqn{n\times p}{n-by-p} matrix where each
  row corresponds to a random rotation in matrix (p=9) or
  quaternion form (p=4)}

  \item{alpha}{The alpha level desired, e.g. 0.05 or 0.10}
}
\value{
  Radius of the confidence region centered at the projected
  mean
}
\description{
  Compute the radius of a \eqn{100(1-\alpha)%} confidence
  region for the central orientation
}
\details{
  Compute the radius of a \eqn{100(1-\alpha)%} confidence
  region for the central orientation based on the projected
  mean estimator based on a result due to \cite{Chang et
  al. (2001)}.  By construction each axis will have the
  same radius so the radius reported is for all three axis.
}
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Rs,method='chang',alpha=0.1)
}
\references{
  Chang T and Rivest L (2001). "M-estimation for location
  and regression parameters in group models: A case study
  using Stiefel manifolds." _Annals of statistics_, pp.
  784-814.
}
\seealso{
  \code{\link{prentice}} \code{\link{fisheretal}}
  \code{\link{zhang}}
}

