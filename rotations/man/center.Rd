\name{center}
\alias{center}
\alias{center.Q4}
\alias{center.SO3}
\title{Center rotation data}
\usage{
  center(Rs, S)

  \method{center}{SO3} (Rs, S)

  \method{center}{Q4} (Qs, S)
}
\arguments{
  \item{Rs,Qs}{the sample to be centered}

  \item{S}{the rotation to sample around}
}
\value{
  The centered sample
}
\description{
  This function will take the sample Rs and return teh
  sample Rs centered at S, i.e., the returned sample is
  S'Rs.  If S is the true center then the projected mean
  should be the 3-by-3 identity matrix
}
\examples{
Rs<-ruars(5,rcayley)
cRs<-center(Rs,mean(Rs))
mean(cRs) #Should be close to identity matrix
}

