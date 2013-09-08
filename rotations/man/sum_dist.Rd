\name{sum_dist}
\alias{sum_dist}
\alias{sum_dist.Q4}
\alias{sum_dist.SO3}
\title{Sample distance}
\usage{
  sum_dist(x, S = genR(0, space = class(x)),
    method = "projected", p = 1)

  \method{sum_dist}{SO3} (x, S = id.SO3,
    method = "projected", p = 1)

  \method{sum_dist}{Q4} (x, S = id.Q4,
    method = "projected", p = 1)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (p=9) or
  quaternion (p=4) form.}

  \item{S}{the individual matrix of interest, usually an
  estimate of the mean.}

  \item{method}{type of distance used method in "projected"
  or "intrinsic"}

  \item{p}{the order of the distances to compute.}
}
\value{
  The sum of the pth order distance between each sample in
  Rs and S.
}
\description{
  Compute the sum of the \eqn{p^{th}}{pth} order distances
  between Rs and S.
}
\examples{
Rs<-ruars(20,rvmises,kappa=10)
Sp<-mean(Rs)
sum_dist(Rs,S=Sp,p=2)
}
\seealso{
  \code{\link{dist.SO3}}, \code{\link{dist.Q4}}
}

