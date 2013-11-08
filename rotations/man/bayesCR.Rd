\name{bayesCR}
\alias{bayesCR}
\alias{bayesCR.SO3}
\title{Bayes credible regions}
\usage{
  bayesCR(x, lpangle, S0, kappa0, tuneS, tuneK, burn_in,
    m = 5000, alp = NULL)

  \method{bayesCR}{SO3} (x, lpangle, S0, kappa0, tuneS,
    tuneK, burn_in, m = 5000, alp = NULL)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (p=9) or
  quaternion (p=4) form.}

  \item{alp}{alpha level desired, e.g. 0.05 or 0.10.}
}
\value{
  Radius of the confidence region centered at the posterior
  mode.
}
\description{
  Find the radius of a \eqn{100(1-\alpha)}\% credible
  region for the central orientation and concentration
  parameter.
}
\details{
  Compute the radius of a \eqn{100(1-\alpha)}\% credible
  region for the central orientation and concentration
  parameter as proposed in \cite{bingham2009b} and
  \cite{bingham2010}.  To be able to reduce their method to
  a radius requires the additonal assumption of rotational
  symmetry, equation (10) in \cite{Fisher et al. (1996)}.
}
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Rs,method='eigen',type='bootstrap',alp=0.1,symm=TRUE,estimator='mean')
}
\references{
  bingham2009b

  bingham2010

  Fisher N, Hall P, Jing B and Wood A (1996). "Improved
  pivotal methods for constructing confidence regions with
  directional data." _Journal of the American Statistical
  Association_, *91*(435), pp. 1062-1070.
}
\seealso{
  \code{\link{prentice}}, \code{\link{chang}},
  \code{\link{zhang}}
}

