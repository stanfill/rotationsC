\name{nickel}
\alias{nickel}
\docType{data}
\title{
Data from Bingham et. al. (2009)
}
\description{
ESBD scans from a piece of nickel
}
\usage{data(nickel)}
\format{
  A data frame with 48286 rows and the following 12 columns
  \describe{
    \item{\code{xpos}}{location x position}
    \item{\code{ypos}}{location y position}
    \item{\code{rep}}{Which scan that observation was collected during}
    \item{\code{v1}\code{v2}\code{v3}}{x-axis of crystal orientation at that location/rep}
    \item{\code{v4}\code{v5}\code{v6}}{y-axis of crystal orientation at that location/rep}
    \item{\code{v7}\code{v8}\code{v9}}{z-axis of crystal orientation at that location/rep}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Bingham M, Nordman D and Vardeman S (2009). "Modeling and Inference for Measured Crystal Orientations and a Tractable Class of Symmetric Distributions for Rotations in three Dimensions." Journal of the American Statistical Association, 104 (488), pp. 1385-1397.
}
\examples{
data(nickel)
Rep1<-nickel[nickel$rep==1,]
qplot(xpos,ypos,data=Rep1,colour=v1,size=I(2))
}
\keyword{datasets}
