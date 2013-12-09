\name{Arithmetic}
\alias{-.Q4}
\alias{-.SO3}
\alias{"-.Q4"}
\alias{"-.SO3"}
\alias{"+.Q4"}
\alias{"+.SO3"}
\alias{+.Q4}
\alias{+.SO3}
\alias{Arithmetic}
\title{Arithmetic operators on SO(3)}
\usage{
  \method{+}{SO3} (x, y)

  \method{-}{SO3} (x, y = NULL)

  \method{+}{Q4} (x, y)

  \method{-}{Q4} (x, y = NULL)
}
\arguments{
  \item{x}{first arguement}

  \item{y}{second arguement (optional for subtraction)}
}
\value{
  \item{+}{the result of rotating the identity frame
  through x then y} \item{-}{the difference of the
  rotations, or the inverse rotation of only one arguement
  is provided}
}
\description{
  These binary operators perform arithmetic on rotations in
  quaternion or rotation matrix form (or objects which can
  be coerced into them).
}
\details{
  The rotation group SO(3) is a multiplicative group so
  ``adding" rotations \eqn{R_1}{R1} and \eqn{R_2}{R2}
  results in \eqn{R_1+R_2=R_2R_1}{R1+R2=R2R1}.  Similarly,
  the difference between rotations \eqn{R_1}{R1} and
  \eqn{R_2}{R2} is \eqn{R_1-R_2=R_2^\top R_1}{R1-R2=R2'R1}.
  With this definiton it is clear that
  \eqn{R_1+R_2-R_2=R_2^\top
  R_2R_1=R_1}{R1+R2-R2=R2'R2R1=R1}. If only one rotation is
  provided to subtraction then the inverse (transpose) it
  returned, e.g. \eqn{-R_2=R_2^\top}{-R2=R2'}.
}

