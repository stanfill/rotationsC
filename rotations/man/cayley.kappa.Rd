\name{cayley.kappa}
\alias{cayley.kappa}
\title{Circular variance and concentration parameter}
\usage{
  cayley.kappa(nu)
}
\arguments{
  \item{nu}{circular variance}
}
\value{
  Concentration parameter corresponding to nu.
}
\description{
  Return the concentration parameter that corresponds to a
  given circular variance.
}
\details{
  The concentration parameter \eqn{\kappa} does not
  translate across circular distributions.  A commonly used
  measure of spread in circular distributions that does
  translate is the circular variance defined as
  \eqn{\nu=1-E[\cos(r)]}{\nu=1-E[cos(r)]} where
  \eqn{E[\cos(r)]}{E[cos(r)]} is the mean resultant length.
  See \cite{Mardia et al. (2000)} for more details.  This
  function translates the circular variance \eqn{\nu} into
  the corresponding concentration parameter \eqn{\kappa}
  for the Cayley distribution.
}
\examples{
#Find the concentration parameter for circular variances 0.25, 0.5, 0.75
 cayley.kappa(0.25)
 cayley.kappa(0.5)
 cayley.kappa(0.75)
}
\references{
  Mardia K and Jupp P (2000). Directional statistics.
  John Wiley & Sons, Chichester, England.
}
\seealso{
  \code{\link{Cayley}}
}

