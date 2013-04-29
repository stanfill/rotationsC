\name{medianSO3C}
\alias{medianSO3C}
\title{Median Rotation}
\usage{
  medianSO3C(Rs)
}
\arguments{
  \item{Rs}{A \eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix form
  (\eqn{p=9}) or quaternion form (\eqn{p=4})}
}
\value{
  an estimate of the projected or intrinsic mean
}
\description{
  Compute the projected or intrinsic median of a sample of
  rotations
}
\details{
  The median-type estimators are defined as
  \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in
  SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}{argmin\sum
  d(Ri,S)}.  If the choice of distance metrid,
  \eqn{d_D}{d}, is Riemannian then the estimator is called
  the intrinsic, and if the distance metric in Euclidean
  then it projected. The algorithm used in the intrinsic
  case is discussed in \cite{hartley11} and the projected
  case was written by the authors.
}
\references{
  hartley11
}
\seealso{
  \code{\link{meanSO3C}}
}

