\name{misangle}
\alias{misangle}
\title{Misorientation angle}
\usage{
  misangle(x)

  \method{misangle}{SO3} (x)

  \method{misangle}{Q4} (x)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}
}
\value{
  Angle of rotation.
}
\description{
  Compute the misorientation angle of a rotation.
}
\details{
  Every rotation can be thought of as some reference
  coordinate system rotated about an axis through an angle.
  These quantites are referred to as the misorientation
  axis and misorientation angle, respectively, in the
  material sciences literature. This function returns the
  misorentation angle associated with a rotation assuming
  the reference coordinate system is the identity.
}
\seealso{
  \code{\link{misaxis}}
}

