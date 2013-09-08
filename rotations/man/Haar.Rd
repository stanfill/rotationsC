\name{Haar}
\alias{dhaar}
\alias{Haar}
\alias{phaar}
\alias{rhaar}
\title{Uniform distribution on \eqn{SO(3)}}
\usage{
  dhaar(r)

  phaar(q, lower.tail = TRUE)

  rhaar(n)
}
\arguments{
  \item{r,q}{vector of quantiles.}

  \item{n}{number of observations.  If \code{length(n)>1},
  the length is taken to be the number required.}

  \item{lower.tail}{logical; if TRUE (default),
  probabilities are \eqn{P(X \le x)} otherwise, \eqn{P(X >
  x)}.}
}
\value{
  \item{dhaar}{gives the density} \item{phaar}{gives the
  distribution function} \item{rhaar}{generates random
  deviates}
}
\description{
  Density, distribution function and random generation for
  the uniform distribution on \eqn{SO(3)}.
}
\details{
  The uniform distribution on the space \eqn{SO(3)} (also
  referred to as the Haar measure) has density
  \deqn{C_U(r)=\frac{[1-cos(r)]}{2\pi}.}{C(r)=[1-cos(r)]/2\pi.}
}
\seealso{
  \link{Angular-distributions} for other distributions in
  the rotations package.
}

