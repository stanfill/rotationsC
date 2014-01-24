\name{weighted.mean.SO3}
\alias{weighted.mean.Q4}
\alias{weighted.mean.SO3}
\title{Weighted mean rotation}
\usage{
  \method{weighted.mean}{SO3} (x, w, type = "projected",
    epsilon = 1e-05, maxIter = 2000, ...)

  \method{weighted.mean}{Q4} (x, w, type = "projected",
    epsilon = 1e-05, maxIter = 2000, ...)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix form
  (\eqn{p=9}) or quaternion (\eqn{p=4}) form.}

  \item{w}{vector of weights the same length as the number
  of rows in x giving the weights to use for elements of
  x.}

  \item{type}{string indicating "projected" or "geometric"
  type mean estimator.}

  \item{epsilon}{stopping rule for the geometric method.}

  \item{maxIter}{maximum number of iterations allowed
  before returning most recent estimate.}

  \item{...}{only used for consistency with mean.default.}
}
\value{
  Weighted mean of the sample in the same parameterization.
}
\description{
  Compute the weighted geometric or projected mean of a
  sample of rotations.
}
\details{
  This function takes a sample of 3D rotations (in matrix
  or quaternion form) and returns the weighted projected
  arithmetic mean \eqn{\widehat{\bm S}_P}{S_P} or geometric
  mean \eqn{\widehat{\bm S}_G}{S_G} according to the
  \code{type} option. For a sample of \eqn{n} rotations in
  matrix form \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in
  SO(3), i=1,2,\dots,n}, the weighted mean is defined as
  \deqn{\widehat{\bm{S}}=argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nw_id^2(\bm{R}_i,\bm{S})}{argmin\sum wi
  d^2(Ri,S)} where \eqn{d} is the Riemannian or Euclidean
  distance.  For more on the projected mean see
  \cite{Moakher (2002)} and for the geometric mean see
  \cite{Manton (2004)}.
}
\examples{
Rs <- ruars(20, rvmises, kappa = 0.01)
mean(Rs)                   #Find the equal-weight projected mean
wt <- abs(1/mis.angle(Rs)) #Use the rotation misorientation angle as weight
weighted.mean(Rs, wt)      #as weight
rot.dist(mean(Rs))
rot.dist(weighted.mean(Rs, wt)) #usually much smaller than unweighted mean

#Can do the same thing with quaternions
Qs <- as.Q4(Rs)
mean(Qs)
wt <- abs(1/mis.angle(Qs))
weighted.mean(Qs, wt)
rot.dist(mean(Qs))
rot.dist(weighted.mean(Qs, wt))
}
\references{
  Moakher M (2002). "Means and averaging in the group of
  rotations." SIAM Journal on Matrix Analysis and
  Applications, 24(1), pp. 1-16.

  Manton JH (2004). "A globally convergent numerical
  algorithm for computing the centre of mass on compact Lie
  groups." In Control, Automation, Robotics and Vision
  Conference, 2004. ICARCV 2004 8th, volume 3, pp.
  2211-2216. IEEE.
}
\seealso{
  \code{\link{median.SO3}}, \code{\link{mean.SO3}},
  \code{\link{bayes.mean}}
}

