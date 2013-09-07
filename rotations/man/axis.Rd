\name{axis}
\alias{axis}
\alias{axis.Q4}
\alias{axis.SO3}
\title{Misorientation Axis}
\usage{
  axis(x, ...)

  \method{axis}{SO3} (x, ...)

  \method{axis}{Q4} (x, ...)
}
\arguments{
  \item{x}{3-by-3 matrix in SO3 or unit quaterion}

  \item{...}{additional arguements}
}
\value{
  axis in form of three dimensional vector of length one.
}
\description{
  Determine the misorientation axis of a rotation
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
  \code{\link{angle}}
}

