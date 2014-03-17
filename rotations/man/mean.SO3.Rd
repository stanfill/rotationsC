\name{mean.SO3}
\alias{mean.Q4}
\alias{mean.SO3}
\title{Mean rotation}
\usage{
  \method{mean}{SO3} (x, type = "projected",
    epsilon = 1e-05, maxIter = 2000, ...)

  \method{mean}{Q4} (x, type = "projected",
    epsilon = 1e-05, maxIter = 2000, ...)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix form
  (\eqn{p=9}) or quaternion (\eqn{p=4}) form.}

  \item{type}{string indicating "projected" or "geometric"
  type mean estimator.}

  \item{epsilon}{stopping rule for the geometric-mean.}

  \item{maxIter}{maximum number of iterations allowed for
  geometric-mean.}

  \item{...}{additional arguments.}
}
\value{
  Estimate of the projected or geometric mean of the sample
  in the same parameterization.
}
\description{
  Compute the sample geometric or projected mean.
}
\details{
  This function takes a sample of 3D rotations (in matrix
  or quaternion form) and returns the projected arithmetic
  mean denoted \eqn{\widehat{\bm S}_P}{S_P} or geometric
  mean \eqn{\widehat{\bm S}_G}{S_G} according to the
  \code{type} option. For a sample of \eqn{n} rotations in
  matrix form \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in
  SO(3), i=1,2,\dots,n}, the mean-type estimator is defined
  as \deqn{\widehat{\bm{S}}=argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd^2(\bm{R}_i,\bm{S})}{argmin\sum
  d^2(Ri,S)} where \eqn{d} is the Riemannian or Euclidean
  distance. For more on the projected mean see
  \cite{Moakher (2002)} and for the geometric mean see
  \cite{Manton (2004)}. For the projected mean from a
  quaternion point of view see \cite{Tyler (1981)}.
}
\examples{
Rs <- ruars(20, rvmises, kappa = 0.01)
mean(Rs)                               #Projected mean
project.SO3(colMeans(Rs))              #Same as mean(Rs)

mean(Rs, type = "geometric")           #Geometric mean
rot.dist(mean(Rs))                     #Bias of the projected mean
rot.dist(mean(Rs, type = "geometric")) #Bias of the geometric mean

Qs <- as.Q4(Rs)
mean(Qs)                               #Projected mean
mean(Qs, type = "geometric")           #Geometric mean
rot.dist(mean(Qs))                     #Bias of the projected mean
rot.dist(mean(Qs, type = "geometric")) #Bias of the geometric mean
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

  Tyler DE (1981). "Asymptotic inference for eigenvectors."
  The Annals of Statistics, pp. 725-736.
}
\seealso{
  \code{\link{median.SO3}}, \code{\link{bayes.mean}},
  \code{\link{weighted.mean.SO3}}
}

