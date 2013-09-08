\name{dist}
\alias{dist}
\alias{dist.Q4}
\alias{dist.SO3}
\title{Compute the rotational distance}
\usage{
  dist(x, ...)

  \method{dist}{SO3} (x, R2 = id.SO3, method = "projected",
    p = 1, ...)

  \method{dist}{Q4} (x, Q2 = id.Q4, method = "projected",
    p = 1, ...)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (p=9) or
  quaternion (p=4) form.}

  \item{R2,Q2}{the second rotation in the same
  parameterization as x.}

  \item{method}{string indicating "projected" or
  "intrinsic" method of distance.}

  \item{p}{the order of the distance.}

  \item{...}{additional arguments.}
}
\value{
  The rotational distance between each rotation in x and R2
  or Q2.
}
\description{
  Calculate the Euclidean or Riemannian distance between
  two rotations.
}
\details{
  This function will calculate the intrinsic (Riemannian)
  or projected (Euclidean) distance between two rotations.
  \code{R2} and \code{Q2} are set to the identity rotations
  by default.  For rotations \eqn{R_1}{R1} and
  \eqn{R_2}{R2} both in \eqn{SO(3)}, the Euclidean distance
  between them is \deqn{||R_1-R_2||_F}{||R1-R2||} where
  \eqn{||\cdot||_F}{|| ||} is the Frobenius norm. The
  Riemannian distance is defined as \deqn{||Log(R_1^\top
  R_2)||_F}{||Log(R1'R2)||} where \eqn{Log} is the matrix
  logarithm, and it corresponds to the misorientation angle
  of \eqn{R_1^\top R_2}{R1'R2}. To compute the distance
  matrix use \code{stats::dist()}.
}

