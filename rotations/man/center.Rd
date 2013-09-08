\name{center}
\alias{center}
\alias{center.Q4}
\alias{center.SO3}
\title{Center rotation data}
\usage{
  center(x, S)

  \method{center}{SO3} (x, S)

  \method{center}{Q4} (x, S)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (p=9) or
  quaternion (p=4) form.}

  \item{S}{the rotation about which to center x.}
}
\value{
  The centered sample.
}
\description{
  This function will take the sample Rs and return the
  sample Rs centered at S, i.e. the returned sample is
  \eqn{S^\top R}{S'R}.  If S is the true center then the
  projected mean should be close to the 3-by-3 identity
  matrix
}
\examples{
Rs<-ruars(5,rcayley)
cRs<-center(Rs,mean(Rs))
mean(cRs) #Should be close to identity matrix
}

