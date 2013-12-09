\name{prentice}
\alias{prentice}
\alias{prentice.Q4}
\alias{prentice.SO3}
\title{Eigenvector theory confidence region}
\usage{
  prentice(x, alp)

  \method{prentice}{Q4} (x, alp = NULL)

  \method{prentice}{SO3} (x, alp = NULL)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}

  \item{alp}{alpha level desired, e.g. 0.05 or 0.10.}
}
\value{
  Radius of the confidence region centered at the projected
  mean for each of the x-, y- and z-axes.
}
\description{
  Find the radius of a \eqn{100(1-\alpha)}\% confidence
  region for the projected mean based on eigenvector based
  result.
}
\details{
  Compute the radius of a \eqn{100(1-\alpha)}\% confidence
  region for the central orientation based on the projected
  mean estimator using the method due to \cite{Prentice
  (1986)}.  For a rotation specific version see
  \cite{Rancourt et al. (2000)}. The variablity in each
  axis is different so each axis will have its own radius.
}
\examples{
Qs<-ruars(20,rcayley,kappa=100,space='Q4')
region(Qs,method='eigen',type='theory',alp=0.1,estimator='mean')
}
\references{
  Prentice MJ (1986). "Orientation statistics without
  parametric assumptions." Journal of the Royal
  Statistical Society. Series B (Methodological), pp.
  214-222.

  Rancourt D, Rivest L and Asselin J (2000). "Using
  orientation statistics to investigate variations in human
  kinematics." Journal of the Royal Statistical Society:
  Series C (Applied Statistics), 49(1), pp. 81-94.
}
\seealso{
  \code{\link{bayesCR}}, \code{\link{fisheretal}},
  \code{\link{chang}}, \code{\link{zhang}}
}

