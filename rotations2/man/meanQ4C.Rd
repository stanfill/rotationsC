\name{meanQ4C}
\alias{meanQ4C}
\title{Quaternion Mean Rotation}
\usage{
  meanQ4C(Q)
}
\arguments{
  \item{Rs}{A \eqn{n\times 4}{n-by-4} matrix where each row
  corresponds to a random rotation in matrix form}
}
\value{
  Estimate of the projected or intrinsic mean of the sample
}
\description{
  Compute the intrinsic or projected mean of a sample of
  rotations
}
\details{
  This function takes a sample of \eqn{4\times 1}{4-by-1}
  rotations (in the form of a \eqn{n\times 4}{n-by-4}
  matrix where \eqn{n>1} is the sample size) and returns
  the projected arithmetic mean denoted \eqn{\widehat{\bm
  S}_P}{S_P} or intrinsic mean \eqn{\widehat{\bm S}_G}{S_G}
  according to the \code{type} option. For a sample of
  \eqn{n} random rotations \eqn{\bm{R}_i\in SO(3),
  i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the mean-type
  estimator is defined as
  \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})}{argmin
  d^2(bar(R),S)} where
  \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i}{bar(R)=\sum
  Ri/n} and the distance metric \eqn{d_D}{d} is the
  Riemannian or Euclidean.  For more on the projected mean
  see \cite{moakher02} and for the intrinsic mean see
  \cite{manton04}.
}
\examples{
Rs<-ruars(20,rvmises,kappa=0.01)
mean(Rs)
}
\references{
  moakher02

  manton04
}
\seealso{
  \code{\link{medianSO3C}}
}

