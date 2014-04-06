\name{Cayley}
\alias{Cayley}
\alias{dcayley}
\alias{pcayley}
\alias{rcayley}
\title{The symmetric Cayley distribution}
\usage{
  dcayley(r, kappa = 1, nu = NULL, Haar = TRUE)

  pcayley(q, kappa = 1, nu = NULL, lower.tail = TRUE)

  rcayley(n, kappa = 1, nu = NULL)
}
\arguments{
  \item{r,q}{vector of quantiles.}

  \item{n}{number of observations.  If \code{length(n)>1},
  the length is taken to be the number required.}

  \item{kappa}{concentration parameter.}

  \item{nu}{circular variance, can be used in place of
  \code{kappa}.}

  \item{Haar}{logical; if TRUE density is evaluated with
  respect to the Haar measure.}

  \item{lower.tail}{logical; if TRUE (default)
  probabilities are \eqn{P(X\leq x)}{P(X\le x)} otherwise,
  \eqn{P(X>x)}.}
}
\value{
  \item{dcayley}{gives the density} \item{pcayley}{gives
  the distribution function} \item{rcayley}{generates a
  vector of random deviates}
}
\description{
  Density, distribution function and random generation for
  the Cayley distribution with concentration \code{kappa}
  \eqn{\kappa}.
}
\details{
  The symmetric Cayley distribution with concentration
  \eqn{\kappa} has density \deqn{C_C(r
  |\kappa)=\frac{1}{\sqrt{\pi}}
  \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos
  r)^\kappa(1-\cos r).}{C(r |\kappa)=
  \Gamma(\kappa+2)(1+cos r)^\kappa(1-cos
  r)/[\Gamma(\kappa+1/2)2^(\kappa+1)\sqrt\pi].} The Cayley
  distribution is equivalent to the de la Vallee Poussin
  distribution of \cite{Schaeben (1997)}.
}
\examples{
r <- seq(-pi, pi, length = 500)

#Visualize the Cayley density fucntion with respect to the Haar measure
plot(r, dcayley(r, kappa = 10), type = "l", ylab = "f(r)")

#Visualize the Cayley density fucntion with respect to the Lebesgue measure
plot(r, dcayley(r, kappa = 10, Haar = FALSE), type = "l", ylab = "f(r)")

#Plot the Cayley CDF
plot(r,pcayley(r,kappa = 10), type = "l", ylab = "F(r)")

#Generate random observations from Cayley distribution
rs <- rcayley(20, kappa = 1)
hist(rs, breaks = 10)
}
\references{
  Schaeben H (1997). "A Simple Standard Orientation Density
  Function: The Hyperspherical de la Vallee Poussin
  Kernel." Physica Status Solidi (B), 200, pp. 367-376.

  Leon C, e JM and Rivest L (2006). "A statistical model
  for random rotations." Journal of Multivariate
  Analysis, 97(2), pp. 412-430.
}
\seealso{
  \link{Angular-distributions} for other distributions in
  the rotations package.
}

