\name{zhang}
\alias{zhang}
\alias{zhang.Q4}
\alias{zhang.SO3}
\title{M-estimator theory pivotal boostrap confidence region}
\usage{
  zhang(x, estimator, alp, m)

  \method{zhang}{SO3} (x, estimator, alp = NULL, m = 300)

  \method{zhang}{Q4} (x, estimator, alp = NULL, m = 300)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}

  \item{estimator}{character string either "mean" or
  "median."}

  \item{alp}{alpha level desired, e.g. 0.05 or 0.10.}

  \item{m}{number of replicates to use to estiamte the
  critical value.}
}
\value{
  Radius of the confidence region centered at the specified
  estimator.
}
\description{
  Compute the radius of a \eqn{100(1-\alpha)}\% confidence
  region for the central orientation based on M-estiamtor
  theory.
}
\details{
  Compute the radius of a \eqn{100(1-\alpha)}\% confidence
  region for the central orientation based on the projected
  mean estimator using the method due to Zhang & Nordman
  (2009) (unpublished MS thesis).  By construction each
  axis will have the same radius so the radius reported is
  for all three axis.  A normal theory version of this
  procedure uses the theoretical chi-square limiting
  distribution and is given by the \code{\link{chang}}
  option.
}
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Rs,method='moment',type='bootstrap',alp=0.1,estimator='mean')
}
\seealso{
  \code{\link{bayesCR}}, \code{\link{prentice}},
  \code{\link{fisheretal}}, \code{\link{chang}}
}

