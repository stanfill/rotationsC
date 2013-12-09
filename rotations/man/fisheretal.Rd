\name{fisheretal}
\alias{fisheretal}
\alias{fisheretal.Q4}
\alias{fisheretal.SO3}
\title{Eigenvector theory pivotal bootstrap confidence region}
\usage{
  fisheretal(x, alp, boot, m, symm)

  \method{fisheretal}{Q4} (x, alp = NULL, boot = T,
    m = 300, symm = TRUE)

  \method{fisheretal}{SO3} (x, alp = NULL, boot = T,
    m = 300, symm = T)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}

  \item{alp}{alpha level desired, e.g. 0.05 or 0.10.}

  \item{boot}{should the bootstrap or normal theory
  critical value be used.}

  \item{m}{number of bootstrap replicates to use to
  estimate critical value.}

  \item{symm}{logical; if TRUE (default), a symmetric
  region is constructed.}
}
\value{
  Radius of the confidence region centered at the projected
  mean.
}
\description{
  Find the radius of a \eqn{100(1-\alpha)}\% confidence
  region for the central orientation based on eigenvector
  theory.
}
\details{
  Compute the radius of a \eqn{100(1-\alpha)}\% confidence
  region for the central orientation based on the projected
  mean estimator using the method for the mean polar axis
  as proposed in \cite{Fisher et al. (1996)}.  To be able
  to reduce their method to a radius requires the additonal
  assumption of rotational symmetry, equation (10) in
  \cite{Fisher et al. (1996)}.
}
\examples{
Qs<-ruars(20,rcayley,kappa=100,space='Q4')
region(Qs,method='eigen',type='bootstrap',alp=0.1,symm=TRUE,estimator='mean')
}
\references{
  Fisher N, Hall P, Jing B and Wood A (1996). "Improved
  pivotal methods for constructing confidence regions with
  directional data." Journal of the American Statistical
  Association, 91(435), pp. 1062-1070.
}
\seealso{
  \code{\link{bayesCR}}, \code{\link{prentice}},
  \code{\link{chang}}, \code{\link{zhang}}
}

