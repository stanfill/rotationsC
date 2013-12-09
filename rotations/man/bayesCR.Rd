\name{bayesCR}
\alias{bayesCR}
\alias{bayesCR.Q4}
\alias{bayesCR.SO3}
\title{Bayes credible regions}
\usage{
  bayesCR(x, type, S0, kappa0, tuneS, tuneK, burn_in,
    m = 5000, alp = 0.1)

  \method{bayesCR}{SO3} (x, type, S0, kappa0, tuneS, tuneK,
    burn_in, m = 5000, alp = 0.1)

  \method{bayesCR}{Q4} (x, type, S0, kappa0, tuneS, tuneK,
    burn_in, m = 5000, alp = 0.1)
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

  \item{alp}{alpha level desired, e.g. 0.05 or 0.10.}
}
\value{
  list of \item{Shat,Qhat}{the posterior mode} \item{Radius}{the
  radius of the credible region centered at the posterior
  mode}
}
\description{
  Find the radius of a \eqn{100(1-\alpha)}\% credible
  region for the central orientation and concentration
  parameter using non-informative Bayes.
}
\details{
  Compute the radius of a \eqn{100(1-\alpha)}\% credible
  region for the central orientation and concentration
  parameter as described in \cite{Bingham et al. (2009)}
  and \cite{Bingham et al. (2010)}.  The posterior mode is
  returned along with the radius of the credible region
  centered at the posterior mode.
}
\examples{
Rs<-ruars(20,rcayley,kappa=4)
region(Rs,type='Cayley',method='Bayes',estimator='mean',
S0=mean(Rs),kappa0=2,tuneS=39,tuneK=.8,burn_in=100,alp=.01)
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
  \code{\link{fisheretal}}, \code{\link{prentice}},
  \code{\link{chang}}, \code{\link{zhang}}
}

