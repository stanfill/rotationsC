\name{median.SO3}
\alias{median}
\alias{median.Q4}
\alias{median.SO3}
\title{Median rotation}
\usage{
  median(x, ...)

  \method{median}{SO3} (x, type = "projected",
    epsilon = 1e-05, maxIter = 2000, ...)

  \method{median}{Q4} (x, type = "projected",
    epsilon = 1e-05, maxIter = 2000, ...)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix form
  (\eqn{p=9}) or quaternion (\eqn{p=4}) form.}

  \item{type}{string indicating "projected" or "geometric"
  type mean estimator.}

  \item{epsilon}{stopping rule.}

  \item{maxIter}{maximum number of iterations allowed
  before returning most recent estimate.}

  \item{...}{additional arguments.}
}
\value{
  Estimate of the projected or geometric median in the same
  parameterization.
}
\description{
  Compute the sample projected or geometric median.
}
\details{
  The median-type estimators are defined as
  \deqn{\widetilde{\bm{S}}=argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd(\bm{R}_i,\bm{S}).}{argmin\sum
  d(Ri,S).} If the choice of distance metric \eqn{d} is
  Riemannian then the estimator is called the geometric
  median, and if the distance metric in Euclidean then it
  is called the projected median. The algorithm used in the
  geometric case is discussed in \cite{Hartley et al.
  (2011)} and the projected case is in \cite{Stanfill et
  al. (2013)}.
}
\references{
  Hartley R, Aftab K and Trumpf J (2011). "L1 rotation
  averaging using the Weiszfeld algorithm." In 2011 IEEE
  Conference on Computer Vision and Pattern Recognition
  (CVPR), pp. 3041-3048. IEEE.

  Stanfill B, Genschel U and Hofmann H (2013). "Point
  Estimation of the Central Orientation of Random
  Rotations." Technometrics, 55(4), pp. 524-535.
}
\seealso{
  \code{\link{mean.SO3}}, \code{\link{bayes.mean}},
  \code{\link{weighted.mean.SO3}}
}

