\name{Cayley}
\alias{Cayley}
\alias{dcayley}
\alias{pcayley}
\alias{rcayley}
\title{The Symmetric Cayley Distribution}
\usage{
  dcayley(r, kappa = 1, nu = NULL, Haar = TRUE)

  pcayley(q, kappa = 1, nu = NULL, lower.tail = TRUE)

  rcayley(n, kappa = 1, nu = NULL)
}
\arguments{
  \item{r,q}{vector of quantiles}

  \item{n}{number of observations.  If \code{length(n)>1},
  the length is taken to be the number required}

  \item{kappa}{Concentration paramter}

  \item{nu}{The circular variance, can be used in place of
  kappa}

  \item{Haar}{logical; if TRUE density is evaluated with
  respect to Haar}

  \item{lower.tail}{logical; if TRUE (default)
  probabilities are \eqn{P(X\leq x)}{P(X\le x)} otherwise,
  \eqn{P(X>x)}}
}
\value{
  \item{dcayley}{gives the density} \item{pcayley}{gives
  the distribution function} \item{rcayley}{generates a
  vector of random deviates}
}
\description{
  Density and random generation for the Cayley distribution
  with concentration kappa (\eqn{\kappa})
}
\details{
  The symmetric Cayley distribution with concentration
  kappa (or circular variance nu) had density \deqn{C(r
  |\kappa)=\frac{1}{\sqrt{\pi}}
  \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos
  r)^\kappa(1-\cos r).}{C(r |\kappa)=
  \Gamma(\kappa+2)(1+cos r)^\kappa(1-cos
  r)/[\Gamma(\kappa+1/2)2^(\kappa+1)\sqrt\pi].}
}
\references{
  Schaeben H (1997). "A Simple Standard Orientation Density
  Function: The Hyperspherical de la Vallee Poussin
  Kernel." _Phys. Stat. Sol. (B)_, *200*, pp. 367-376.

  Leon C, e JM and Rivest L (2006). "A statistical model
  for random rotations." _Journal of Multivariate
  Analysis_, *97*(2), pp. 412-430.
}
\seealso{
  \link{Angular-distributions} for other distributions in
  the rotations package
}

