% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimators.R
\name{mean}
\alias{mean}
\alias{mean.SO3}
\alias{mean.Q4}
\title{Mean rotation}
\usage{
\method{mean}{SO3}(x, type = "projected", epsilon = 1e-05, maxIter = 2000, ...)

\method{mean}{Q4}(x, type = "projected", epsilon = 1e-05, maxIter = 2000, ...)
}
\arguments{
\item{x}{\eqn{n\times p}{n-by-p} matrix where each row corresponds to a
random rotation in matrix form (\eqn{p=9}) or quaternion (\eqn{p=4}) form.}

\item{type}{string indicating "projected" or "geometric" type mean estimator.}

\item{epsilon}{stopping rule for the geometric-mean.}

\item{maxIter}{maximum number of iterations allowed for geometric-mean.}

\item{...}{additional arguments.}
}
\value{
Estimate of the projected or geometric mean of the sample in the same
  parametrization.
}
\description{
Compute the sample geometric or projected mean.
}
\details{
This function takes a sample of 3D rotations (in matrix or quaternion form)
and returns the projected arithmetic mean denoted \eqn{\widehat{\bm
S}_P}{S_P} or geometric mean \eqn{\widehat{\bm S}_G}{S_G} according to the
\code{type} option. For a sample of \eqn{n} rotations in matrix form
\eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the
mean-type estimator is defined as \deqn{\widehat{\bm{S}}=argmin_{\bm{S}\in
SO(3)}\sum_{i=1}^nd^2(\bm{R}_i,\bm{S})}{argmin\sum d^2(Ri,S)} where \eqn{d}
is the Riemannian or Euclidean distance. For more on the projected mean see
\cite{moakher02} and for the geometric mean see \cite{manton04}. For the
projected mean from a quaternion point of view see \cite{tyler1981}.

tyler1981, moakher02, manton04
}
\examples{
Rs <- ruars(20, rvmises, kappa = 0.01)

# Projected mean
mean(Rs)

# Same as mean(Rs)
project.SO3(colMeans(Rs))

# Geometric mean
mean(Rs, type = "geometric")

# Bias of the projected mean
rot.dist(mean(Rs))

# Bias of the geometric mean
rot.dist(mean(Rs, type = "geometric"))

# Same thing with quaternion form
Qs <- as.Q4(Rs)
mean(Qs)
mean(Qs, type = "geometric")
rot.dist(mean(Qs))
rot.dist(mean(Qs, type = "geometric"))
}
\seealso{
\code{\link{median.SO3}}, \code{\link{bayes.mean}}, \code{\link{weighted.mean.SO3}}
}
