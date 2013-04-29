\name{rdistSO3C}
\alias{rdistSO3C}
\title{Rotational Distance}
\usage{
  rdistSO3C(Rs, R2)
}
\arguments{
  \item{Rs}{(or Qs) a matrix of rotations in either matrix
  or quaternion form}

  \item{R2}{(or Q2) the second rotation in the same
  parameterization as R1}
}
\value{
  a vector of rotational distances from each row of Rs (or
  Qs) to R2
}
\description{
  Calculate the Euclidean or Riemannian distance between
  two rotations
}
\details{
  This function will calculate the intrinsic (Riemannian)
  or projected (Euclidean) distance between two rotations.
  If only one rotation is specified the other will be set
  to the identity and the distance between the two is
  returned.  For rotations \eqn{R_1}{R1} and \eqn{R_2}{R2}
  both in \eqn{SO(3)}, the Euclidean distance between them
  is \deqn{||R_1-R_2||_F}{||R1-R2||} where
  \eqn{||\cdot||_F}{|| ||} is the Frobenius norm. The
  intrinsic distance is defined as \deqn{||Log(R_1^\top
  R_2)||_F}{||Log(R1'R2)||} where \eqn{Log} is the matrix
  logarithm, and it corresponds to the misorientation angle
  of \eqn{R_1^\top R_2}{R1'R2}.
}

