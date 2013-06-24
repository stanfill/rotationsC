\name{UARS}
\alias{duars}
\alias{puars}
\alias{ruars}
\alias{UARS}
\title{Generic UARS Distribution}
\usage{
  duars(os, dangle, S = diag(3), kappa = 1, ...)

  puars(os, pangle, S = diag(3), kappa = 1, ...)

  ruars(n, rangle, S = NULL, kappa = 1, space = "SO3", ...)
}
\arguments{
  \item{os}{Value at which to evaluate the UARS density}

  \item{dangle}{The function to evaulate the angles from:
  e.g. dcayley, dvmises, dfisher, dhaar}

  \item{pangle}{The form of the angular density: e.g.
  pcayley, pvmises, pfisher, phaar}

  \item{rangle}{The function from which to simulate angles:
  e.g. rcayley, rvmises, rhaar, rfisher}

  \item{S}{principal direction of the distribution}

  \item{kappa}{concentration of the distribution}

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
  For the rotation R with central orientation S and
  concentration \eqn{kappa} the UARS density is given by
  \deqn{\frac{4\pi}{3-tr(S'R)}C(acos[tr(S'R)-1]/2|\kappa)}{4\pi
  C(acos[tr(S'R)-1]/2|\kappa)/(3-tr(S'R))} where
  \eqn{C(r|\kappa)} is one of the
  \link{Angular-distributions}.
}
\seealso{
  For more on the angular distribution options see
  \link{Angular-distributions}
}

