\name{drill}
\alias{drill}
\docType{data}
\title{
Drill data set
}
\description{
The \code{drill} data set was collected to assess variation in human movement while performing a task (Rancourt 1995).  Eight subjects drilled into a metal plate while being monitored by infared cameras. Quaternions are used to represent the orientation of each subjects' wrist, elbow and shoulder in one of six positions.  For some subjects several replicates are available.  See Rancourt et al. (2000) for one approach to analyzing these data.
}
\usage{data(drill)}
\format{
  A data frame with 720 observations on the following 8 variables.
  \describe{
    \item{\code{Subject}}{Subject number (1-8)}
    \item{\code{Joint}}{Joint name (Wrist, elbow, shoulder)}
    \item{\code{Position}}{Drilling position (1-6)}
    \item{\code{Replicate}}{Replicate number (1-5)}
    \item{\code{Q1}}{First element of orientation (quaternion)}
    \item{\code{Q2}}{Second element of orientation (quaternion)}
    \item{\code{Q3}}{Third element of orientation (quaternion)}
    \item{\code{Q4}}{Fourth element of orientation (quaternion)}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
http://www.mat.ulaval.ca/departement-et-professeurs/direction-et-personnel/professeurs/louis-paul-rivest/publications/
}
\references{
  Rancourt D (1995). "Arm posture and hand mechanical impedance in the control
  of a hand-held power drill." Ph.D. Thesis, MIT.

  Rancourt D, Rivest L and Asselin J (2000). "Using
  orientation statistics to investigate variations in human
  kinematics." Journal of the Royal Statistical Society:
  Series C (Applied Statistics), 49(1), pp. 81-94.
}
\examples{
#Load dataset into workspace
data(drill)

#Estimate central orientation of the first subect's wrist
Subject1Wrist <- subset(drill, Subject == 1 & Joint == "Wrist")
Qs <- as.Q4(Subject1Wrist[,5:8])
mean(Qs)

#Plot Subject 1's wrist measurements using the connection to rotation matrices
\dontrun{
plot(Qs, col = c(1, 2, 3))}

#Translate the quaternion measurements into rotations and estimate 
#the central orientation in terms of rotations
Rs <- as.SO3(Qs)
mean(Rs)
}
\keyword{datasets}
