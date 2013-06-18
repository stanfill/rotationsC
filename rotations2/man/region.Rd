\name{region}
\alias{region}
\alias{region.Q4}
\alias{region.SO3}
\title{Confidence Region for Mean Rotation}
\usage{
  region(Qs, method, alpha, ...)

  \method{region}{Q4} (Qs, method, alpha, ...)

  \method{region}{SO3} (Rs, method, alpha, ...)
}
\arguments{
  \item{Rs,Qs}{A \eqn{n\times p}{n-by-p} matrix where each
  row corresponds to a random rotation in matrix (p=9) or
  quaternion form (p=4)}

  \item{method}{Character string specifying which type of
  interval is required}

  \item{alpha}{The alpha level desired, e.g. 0.95 or 0.90}

  \item{...}{Additional arguments that are method specific}
}
\value{
  Radius of the confidence region centered at the projected
  mean
}
\description{
  Find the radius of a \eqn{100(1-\alpha)%} confidence
  region for the central orientation based on the projected
  mean estimator. The current methods available are due to
  \code{\link{prentice}}, \code{\link{fisher}}, and
  \code{\link{zhang}}.
}
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Qs,method='prentice',alpha=0.1)
region(Qs,method='fisher',alpha=0.1,symm=T)
region(Rs,method='zhang',alpha=0.1,m=100)
}
\references{
  Fisher N, Hall P, Jing B and Wood A (1996). "Improved
  pivotal methods for constructing confidence regions with
  directional data." _Journal of the American Statistical
  Association_, *91*(435), pp. 1062-1070.

  Prentice MJ (1986). "Orientation statistics without
  parametric assumptions." _Journal of the Royal
  Statistical Society. Series B (Methodological)_, pp.
  214-222.

  Rancourt D, Rivest L and Asselin J (2000). "Using
  orientation statistics to investigate variations in human
  kinematics." _Journal of the Royal Statistical Society:
  Series C (Applied Statistics)_, *49*(1), pp. 81-94.
}

