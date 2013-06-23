\name{UARS}
\alias{duars}
\alias{puars}
\alias{ruars}
\alias{UARS}
\title{Generic UARS Distribution}
\usage{
  duars(os, S = diag(3), kappa = 1, dangle, ...)

  puars(os, S = diag(3), kappa = 1, pangle, ...)

  ruars(n, S = NULL, kappa = 1, rangle, space = "SO3", ...)
}
\arguments{
  \item{os}{Value at which to evaluate the UARS density}

  \item{S}{principal direction of the distribution}

  \item{kappa}{concentration of the distribution}

  \item{dangle}{The function to evaulate the angles from:
  e.g. dcayley, dvmises, dfisher, dhaar}

  \item{pangle}{The form of the angular density: e.g.
  pcayley, pvmises, pfisher, phaar}

  \item{rangle}{The function from which to simulate angles:
  e.g. rcayley, rvmises, rhaar, rfisher}

  \item{space}{Indicates the desired representation: matrix
  (SO3) or quaternion (Q4)}

  \item{...}{additional arguments passed to the angular
  distribution}
}
\value{
  \item{duars}{gives the density} \item{puars}{gives the
  distribution function} \item{ruars}{generates random
  deviates}
}
\description{
  Density and random generation for the the generic
  uniform-axis random-spin class of distributions
}
\details{
  Give the form of the UARS distribution
}
\seealso{
  For more on the angular distribution options see
  \link{Angular-distributions}
}

