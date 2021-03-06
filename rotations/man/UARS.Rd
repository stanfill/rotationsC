% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{UARS}
\alias{UARS}
\alias{puars}
\alias{duars}
\alias{ruars}
\title{Generic UARS Distribution}
\usage{
duars(R, dangle, S = id.SO3, kappa = 1, ...)

puars(R, pangle = NULL, S = id.SO3, kappa = 1, ...)

ruars(n, rangle, S = NULL, kappa = 1, space = "SO3", ...)
}
\arguments{
\item{R}{Value at which to evaluate the UARS density.}

\item{dangle}{The function to evaluate the angles from, e.g. dcayley, dvmises, dfisher, dhaar.}

\item{S}{central orientation of the distribution.}

\item{kappa}{concentration parameter.}

\item{...}{additional arguments.}

\item{pangle}{The form of the angular density, e.g. pcayley, pvmises, pfisher, phaar.}

\item{n}{number of observations. If \code{length(n)>1}, the length is taken to be the number required.}

\item{rangle}{The function from which to simulate angles, e.g. rcayley, rvmises, rhaar, rfisher.}

\item{space}{indicates the desired representation: matrix ("SO3") or quaternion ("Q4").}
}
\value{
\item{duars}{gives the density}
         \item{puars}{gives the distribution function.  If pangle is left empty, the empirical CDF is returned.}
         \item{ruars}{generates random deviates}
}
\description{
Density, distribution function and random generation for the the generic uniform axis-random spin (UARS) class of distributions.
}
\details{
For the rotation R with central orientation S and concentration \eqn{\kappa} the UARS density is given by
\deqn{f(R|S,\kappa)=\frac{4\pi}{3-tr(S^\top R)}C(\cos^{-1}[tr(S^\top R)-1]/2|\kappa)}{f(R|S,\kappa)=4\pi C(acos[tr(S'R)-1]/2|\kappa)/[3-tr(S'R)]}
where \eqn{C(r|\kappa)} is one of the \link{Angular-distributions}.

bingham09
}
\examples{
#Generate random rotations from the Cayley-UARS distribution with central orientation
#rotated about the y-axis through pi/2 radians
S <- as.SO3(c(0, 1, 0), pi/2)
Rs <- ruars(20, rangle = rcayley, kappa = 1, S = S)

rs <- mis.angle(Rs-S)                          #Find the associated misorientation angles
frs <- duars(Rs, dcayley, kappa = 10, S = S)   #Compute UARS density evaluated at each rotations
plot(rs, frs)

cdf <- puars(Rs, pcayley, S = S)               #By supplying 'pcayley', it is used to compute the
plot(rs, cdf)                                  #the CDF

ecdf <- puars(Rs, S = S)                       #No 'puars' arguement is supplied so the empirical
plot(rs, ecdf)                                 #cdf is returned
}
\seealso{
For more on the angular distribution options see \link{Angular-distributions}.
}
