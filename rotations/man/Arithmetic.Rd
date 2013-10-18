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
  difference of the rotations}
}
\description{
  These binary operators perform arithmetic on rotations in
  quaternion or rotation matrix form (or objects which can
  be coerced to them).
}
\details{
  The rotation group SO(3) is a multiplicative group so
  addition of rotations \eqn{R1} and \eqn{R2} is
  \eqn{R1+R2=R2R1}.
}

