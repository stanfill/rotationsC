\name{sum.rotdist}
\alias{sum.dist.Q4}
\alias{sum.dist.SO3}
\alias{sum.rotdist}
\alias{sum.rotdist.Q4}
\alias{sum.rotdist.SO3}
\title{Sample distance}
\usage{
  sum.rotdist(x, S = genR(0, space = class(x)),
    method = "extrinsic", p = 1)

  \method{sum.rotdist}{SO3} (x, S = id.SO3,
    method = "extrinsic", p = 1)

  \method{sum.rotdist}{Q4} (x, S = id.Q4,
    method = "extrinsic", p = 1)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}

  \item{S}{the individual matrix of interest, usually an
  estimate of the mean.}

  \item{method}{type of distance used method in "extrinsic"
  or "intrinsic"}

  \item{p}{the order of the distances to compute.}
}
\value{
  The sum of the pth order distance between each row of x
  and S.
}
\description{
  Compute the sum of the \eqn{p^{th}}{pth} order distances
  between each row of x and S.
}
\examples{
Rs<-ruars(20, rvmises, kappa = 10)
Sp<-mean(Rs)
sum.rotdist(Rs, S = Sp, p = 2)
}
\seealso{
  \code{\link{rot.dist}}
}

