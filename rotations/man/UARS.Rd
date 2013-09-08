\name{UARS}
\alias{duars}
\alias{puars}
\alias{ruars}
\alias{UARS}
\title{Generic UARS Distribution}
\usage{
  duars(R, dangle, S = id.SO3, kappa = 1, ...)

  puars(R, pangle, S = id.SO3, kappa = 1, ...)

  ruars(n, rangle, S = NULL, kappa = 1, space = "SO3", ...)
}
\arguments{
  \item{R}{Value at which to evaluate the UARS density.}

  \item{n}{number of observations. If \code{length(n)>1},
  the length is taken to be the number required.}

  \item{dangle}{The function to evaluate the angles from
  e.g. dcayley, dvmises, dfisher, dhaar.}

  \item{pangle}{The form of the angular density e.g.
  pcayley, pvmises, pfisher, phaar.}

  \item{rangle}{The function from which to simulate angles
  e.g. rcayley, rvmises, rhaar, rfisher.}

  \item{S}{central orientation of the distribution.}

  \item{kappa}{concentration parameter.}

  \item{space}{indicates the desired representation: matrix
  ("SO3") or quaternion ("Q4").}

  \item{...}{additional arguments.}
}
\value{
  \item{duars}{gives the density} \item{puars}{gives the
  distribution function} \item{ruars}{generates random
  deviates}
}
\description{
  Density, distribution function and random generation for
  the the generic uniform axis-random spin (UARS) class of
  distributions.
}
\details{
  For the rotation R with central orientation S and
  concentration \eqn{kappa} the UARS density is given by
  \deqn{f(R|S,\kappa)=\frac{4\pi}{3-tr(S^\top
  R)}C(acos[tr(S^\top R)-1]/2|\kappa)}{f(R|S,\kappa)=4\pi
  C(acos[tr(S'R)-1]/2|\kappa)/[3-tr(S'R)]} where
  \eqn{C(r|\kappa)} is one of the
  \link{Angular-distributions}.
}
\references{
  Bingham M, Nordman D and Vardeman S (2009). "Modeling and
  Inference for Measured Crystal Orientations and a
  Tractable Class of Symmetric Distributions for Rotations
  in three Dimensions." Journal of the American
  Statistical Association, 104(488), pp. 1385-1397.
}
\seealso{
  For more on the angular distribution options see
  \link{Angular-distributions}.
}

