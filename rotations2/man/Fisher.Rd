\name{fisher}
\alias{fisher}
\alias{fisher.Q4}
\alias{fisher.SO3}
\title{Fisher confidence region method}
\usage{
  fisher(Qs, alpha, boot, m, symm)

  \method{fisher}{Q4} (Qs, alpha, boot = T, m = 300,
    symm = T)

  \method{fisher}{SO3} (Rs, alpha, boot = T, m = 300,
    symm = T)
}
\arguments{
  \item{Rs,Qs}{A \eqn{n\times p}{n-by-p} matrix where each
  row corresponds to a random rotation in matrix (p=9) or
  quaternion form (p=4)}

  \item{alpha}{The alpha level desired, e.g. 0.05 or 0.10}

  \item{boot}{Should the bootstrap or normal theory
  critical value be used}

  \item{m}{number of bootstrap replicates to use to
  estimate critical value}

  \item{symm}{true/false on if rotationally symmetric
  regions should be computed or not}
}
\value{
  radius of the confidence region centered at the projected
  mean
}
\description{
  Find the radius of a \eqn{100(1-\alpha)%} confidence
  region for the central orientation
}
\details{
  Compute the radius of a \eqn{100(1-\alpha)%} confidence
  region for the central orientation based on the projected
  mean estimator using the method for the mean polar axis
  as proposed in \cite{Fisher et al. (1996)}.  To be able
  to reduce their method to a radius requires the additonal
  assumption of rotational symmetry, equation (10) in
  \cite{Fisher et al. (1996)}.
}
\examples{
Qs<-ruars(20,rcayley,kappa=100,space='Q4')
region(Qs,method='fisher',alpha=0.1,symm=T)
}
\references{
  Fisher N, Hall P, Jing B and Wood A (1996). "Improved
  pivotal methods for constructing confidence regions with
  directional data." _Journal of the American Statistical
  Association_, *91*(435), pp. 1062-1070.
}
\seealso{
  \code{\link{prentice}}, \code{\link{zhang}}
}

