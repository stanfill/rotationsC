\name{vmises_kappa}
\alias{vmises_kappa}
\title{Circular variance and concentration parameter}
\usage{
  vmises_kappa(nu)
}
\arguments{
  \item{nu}{circular variance.}
}
\value{
  Concentration parameter corresponding to nu.
}
\description{
  The concentration parameter \eqn{\kappa} does not
  translate across circular distributions.  A commonly used
  measure of spread in circular distributions that does
  translate is the circular variance defined as
  \eqn{\nu=1-E[\cos(r)]}{\nu=1-E[cos(r)]} where
  \eqn{E[\cos(r)]}{E[cos(r)]} is the mean resultant length.
  See \cite{Mardia et al. (2000)} for more details.  This
  function translates the circular variance \eqn{\nu} into
  the corresponding concentration parameter \eqn{\kappa}
  for the circular-von Mises distribution.
}
\references{
Mardia, K. V., & Jupp, P. E. Directional statistics. 2000.
}
\seealso{
  \code{\link{Mises}}
}

