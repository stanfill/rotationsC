\name{bayes.mean}
\alias{bayes.mean}
\alias{bayes.mean.Q4}
\alias{bayes.mean.SO3}
\title{Parameter estimates based on non-informative Bayes}
\usage{
  bayes.mean(x, type, S0, kappa0, tuneS, tuneK, burn_in,
    m = 5000)

  \method{bayes.mean}{SO3} (x, type, S0, kappa0, tuneS,
    tuneK, burn_in, m = 5000)

  \method{bayes.mean}{Q4} (x, type, S0, kappa0, tuneS,
    tuneK, burn_in, m = 5000)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}

  \item{type}{Angular distribution assumed on R.  Options
  are \code{\link{Cayley}}, \code{\link{Fisher}} or
  \code{\link{Mises}}}

  \item{S0}{initial estimate of central orientation}

  \item{kappa0}{initial estimate of concentration
  parameter}

  \item{tuneS}{central orientation tuning parameter,
  concentration of proposal distribution}

  \item{tuneK}{concentration tuning parameter, standard
  deviation of proposal distribution}

  \item{burn_in}{number of draws to use as burn-in}

  \item{m}{number of draws to keep from posterior
  distribution}
}
\value{
  list of \item{Shat}{Draws from the posterior distribution for
  central orientation S} \item{kappa}{Draws from the
  posterior distribution for concentration parameter kappa}
}
\description{
  Use non-informative Bayes to estimate the central
  orientation and concentration parameter of a sample of
  rotations.
}
\details{
  The procedures detailed in \cite{Bingham et al. (2009)}
  and \cite{Bingham et al. (2010)} are implemented to get
  draws from the posterior distribution for the central
  orientation and concentration parameters for a sample of
  3D rotations.  A uniform prior on SO(3) is used for the
  central orientation and the appropriate Jeffreys prior is
  used for the concentration parameter.
}
\examples{
Rs<-ruars(20,rcayley,kappa=4)
ests<-bayes.mean(Rs,type='Cayley',S0=mean(Rs),kappa0=4,tuneS=39,tuneK=.8,burn_in=100,m=5000)
}
\references{
  Bingham MA, Vardeman SB and Nordman DJ (2009). "Bayes
  one-sample and one-way random effects analyses for 3-D
  orientations with application to materials science."
  Bayesian Analysis, 4(3), pp. 607-629.

  Bingham MA, Nordman DJ and Vardeman SB (2010).
  "Finite-sample investigation of likelihood and Bayes
  inference for the symmetric von Mises-Fisher
  distribution." Computational Statistics & Data
  Analysis, 54(5), pp. 1317-1327.
}
\seealso{
  \code{\link{mean.SO3}}, \code{\link{median.SO3}}
}

