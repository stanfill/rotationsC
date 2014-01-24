\name{mis.axis}
\alias{mis.axis}
\alias{mis.axis.Q4}
\alias{mis.axis.SO3}
\title{Misorientation axis}
\usage{
  mis.axis(x, ...)

  \method{mis.axis}{SO3} (x, ...)

  \method{mis.axis}{Q4} (x, ...)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}

  \item{...}{additional arguments.}
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
  These quantities are referred to as the misorientation
  axis and misorientation angle, respectively, in the
  material sciences literature. This function returns the
  misorentation axis associated with a rotation assuming
  the reference coordinate system is the identity.
}
\examples{
rs <- rcayley(20, kappa = 20)
Rs <- genR(rs, S = id.SO3)
mis.axis(Rs)
all.equal(Rs, as.SO3(mis.axis(Rs), mis.angle(Rs)))

Qs <- genR(rs, S = id.Q4, space = "Q4")
mis.axis(Qs)
all.equal(Qs, as.Q4(mis.axis(Qs), mis.angle(Qs)))
}
\seealso{
  \code{\link{mis.angle}}
}

