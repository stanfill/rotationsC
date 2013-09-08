\name{genR}
\alias{genR}
\title{Generate rotations}
\usage{
  genR(r, S = NULL, space = "SO3")
}
\arguments{
  \item{r}{vector of angles.}

  \item{S}{central orientation.}

  \item{space}{indicates the desired representation:
  rotation matrix "SO3" or quaternions "Q4."}
}
\value{
  A matrix where each row is a sample point in the desired
  space.
}
\description{
  Generate rotations according to Rodrigues' formula.
}
\details{
  Given a vector \eqn{u\in R^3}{u in R^3} of length one and
  angle of rotation \eqn{r}, a rotation can be formed using
  Rodrigues' formula \deqn{\cos(r)I_{3\times
  3}+\sin(r)\Phi(u)+(1-\cos(r))uu^\top}{cos(r)I+sin(r)\Phi(u)+(1-cos(r))uu'}
  where \eqn{I_{3\times 3}}{I} is the \eqn{3\times
  3}{3-by-3} identity matrix, \eqn{\Phi(u)} is a
  \eqn{3\times 3}{3-by-3} skew-symmetric matirix with upper
  triangular elements \eqn{-u_3}{-u3}, \eqn{u_2}{u2} and
  \eqn{-u_1}{-u1} in that order.
}
\examples{
r<-rvmises(20,0.01)
Rs<-genR(r,space="SO3")
Qs<-genR(r,space="Q4")
}

