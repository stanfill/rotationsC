\name{misaxis}
\alias{misaxis}
\title{Misorientation axis}
\usage{
  misaxis(x, ...)

  \method{misaxis}{SO3} (x, ...)

  \method{misaxis}{Q4} (x, ...)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}

  \item{...}{additional arguements.}
}
\value{
  Axis in form of three dimensional vector of length one.
}
\description{
  Determine the misorientation axis of a rotation.
}
\details{
  Every rotation can be interpreted as some reference
  coordinate system rotated about an axis through an angle.
  These quantites are referred to as the misorientation
  axis and misorientation angle, respectively, in the
  material sciences literature. This function returns the
  misorentation axis associated with a rotation assuming
  the reference coordinate system is the identity.
}
\seealso{
  \code{\link{misangle}}
}

