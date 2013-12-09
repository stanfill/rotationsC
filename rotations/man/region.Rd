\name{region}
\alias{region}
\alias{region.Q4}
\alias{region.SO3}
\title{Confidence and creible regions for the central orientation}
\usage{
  region(x, method, type, estimator, alp, ...)

  \method{region}{Q4} (x, method, type, estimator,
    alp = NULL, ...)

  \method{region}{SO3} (x, method, type, estimator,
    alp = NULL, ...)
}
\arguments{
  \item{x}{\eqn{n\times p}{n-by-p} matrix where each row
  corresponds to a random rotation in matrix (\eqn{p=9}) or
  quaternion (\eqn{p=4}) form.}

  \item{method}{character string specifying which type of
  interval to report, "bayes", "eigen" or "moment" based
  theory.}

  \item{type}{characted string, "bootstrap" or "theory" are
  available.  For Bayes regions, give the type of
  likelihood: "Cayley","Mises" or "Fisher."}

  \item{estimator}{character string either "mean" or
  "median."  Note that not all method/type combinations are
  available for both estimators.}

  \item{alp}{the alpha level desired, e.g. 0.05 or 0.10.}

  \item{...}{additional arguments that are method
  specific.}
}
\value{
  For frequentist regions only the radius of the confidence
  region centered at the specified estimator is returned.
  For Bayes regions the posterior mode and radius of the
  credible region centered at that mode is returned.
}
\description{
  Find the radius of a \eqn{100(1-\alpha)}\% confidence or
  credible region for the central orientation based on the
  projected mean or median. For more on the currently
  available methods see \code{\link{prentice}},
  \code{\link{fisheretal}}, \code{\link{chang}},
  \code{\link{zhang}} and \code{\link{bayesCR}}.
}
\examples{
Rs<-ruars(20,rcayley,kappa=100)
region(Rs,method='eigen',type='theory',estimator='mean',alp=0.1)
region(Rs,method='eigen',type='bootstrap',estimator='mean',alp=0.1,symm=TRUE)
region(Rs,method='moment',type='bootstrap',estimator='mean',alp=0.1,m=100)
region(Rs,method='moment',type='theory',estimator='mean',alp=0.1)
region(Rs,method='Bayes',type='Cayley',estimator='mean',
S0=mean(Rs),kappa0=2,tuneS=39,tuneK=.8,burn_in=100,alp=.01)
}
\seealso{
  \code{\link{bayesCR}}, \code{\link{prentice}},
  \code{\link{fisheretal}}, \code{\link{chang}},
  \code{\link{zhang}}
}

