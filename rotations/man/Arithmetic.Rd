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
\title{Arithmetic Operators on SO(3)}
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
  \item{+}{the result of both rotations} \item{-}{the
  difference of the rotations, or the inverse rotation of
  only one arguement is provided}
}
\description{
  These binary operators perform arithmetic on rotations in
  quaternion or rotation matrix form (or objects which can
  be coerced to them).
}
\details{
  The rotation group SO(3) is a multiplicative group so
  addition of rotations \eqn{R1} and \eqn{R2} is
  \eqn{R1+R2=R2R1}.  The difference between rotations
  \eqn{R1} and \eqn{R2} is \eqn{R1-R2=R2'R1}.  With this
  definiton it is clear that \eqn{R1+R2-R2=R2'R2R1=R1}
  Finally, if only one rotation is provided to subtraction
  then the inverse (transpose) it returned, i.e.
  \eqn{-R2=R2'}.
}

