\name{region}
\alias{region}
\alias{region.Q4}
\alias{region.SO3}
\title{Confidence region for central orientation}
\usage{
  region(x, method, type, estimator, alp, ...)

  \method{region}{Q4} (x, method, type, estimator,
    alp = NULL, ...)

  \method{region}{SO3} (x, method, type, estimator,
    alp = NULL, ...)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}

  \item{method}{character string specifying which type of
  interval to report, "eigen" or "moment" based theory.}

  \item{type}{characted string, "bootstrap" or "theory" are
  available.  For Bayes regions, give the distribution:
  "Cayley","Mises" or "Fisher."}

  \item{estimator}{character string either "mean" or
  "median."}

  \item{alp}{the alpha level desired, e.g. 0.05 or 0.10.}

  \item{...}{additional arguments that are method
  specific.}
}
\value{
  Radius of the confidence region centered at the projected
  mean.
}
\description{
  Find the radius of a \eqn{100(1-\alpha)}\% confidence
  region for the central orientation based on the projected
  mean or median. For more on the currently available
  methods see \code{\link{prentice}},
  \code{\link{fisheretal}}, \code{\link{chang}},
  \code{\link{zhang}} and \code{\link{bayesCR}}.
}
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Rs,method='eigen',type='theory',estimator='mean',alp=0.1)
region(Rs,method='eigen',type='bootstrap',estimator='mean',alp=0.1,symm=TRUE)
region(Rs,method='moment',type='bootstrap',estimator='mean',alp=0.1,m=100)
region(Rs,method='moment',type='theory',estimator='mean',alp=0.1)
}
\references{
  Prentice MJ (1986). "Orientation statistics without
  parametric assumptions." _Journal of the Royal
  Statistical Society. Series B (Methodological)_, pp.
  214-222.

  Fisher N, Hall P, Jing B and Wood A (1996). "Improved
  pivotal methods for constructing confidence regions with
  directional data." _Journal of the American Statistical
  Association_, *91*(435), pp. 1062-1070.

  Rancourt D, Rivest L and Asselin J (2000). "Using
  orientation statistics to investigate variations in human
  kinematics." _Journal of the Royal Statistical Society:
  Series C (Applied Statistics)_, *49*(1), pp. 81-94.

  Chang T and Rivest L (2001). "M-estimation for location
  and regression parameters in group models: A case study
  using Stiefel manifolds." _Annals of statistics_, pp.
  784-814.
}
\seealso{
  \code{\link{bayesCR}}, \code{\link{prentice}},
  \code{\link{fisheretal}}, \code{\link{chang}},
  \code{\link{zhang}}
}

