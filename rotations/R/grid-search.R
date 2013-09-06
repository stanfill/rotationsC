phi <- seq(0, 2*pi, length=9)
theta <- seq(0, pi, length=9)

angles <- data.frame(expand.grid(list(phi=phi, theta=theta)))
angles <- subset(angles, (theta != 0) & (phi != 0))

require(plyr)
us <- ldply(1:nrow(angles), function(x) {
	angle <- unlist(angles[x,])
	theta <- angle[2]
	phi <- angle[1]
	U <- c(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
	return(cbind(phi=phi, theta=theta, U1=U[1], U2=U[2], U3=U[3]))
})

sphereA <- function(A, theta) {
	## gives a grid of rotation matrices with the same distance from rotation A
  R <- SO3(as.matrix(us[,3:5]), theta=rep(theta, length=nrow(us)))
# multiplication isn't right
  X <- matrix(A, nrow=3)
  for (i in 1:nrow(R)) {
    R[i,] <- as.SO3(X %*% matrix(R[i,], nrow=3))   
  }
  R
}



L2.error <- function(sample, Shat) {
  sum(dist(sample, Shat, method="intrinsic", p=2))
}

error.grid <- function(sample, Shat, theta=1, error) {
	rShat <- sphereA(Shat, theta)
	err <- vector(length=nrow(rShat))
  
	for (i in 1:nrow(rShat)) {
		R <- matrix(unlist(rShat[i, 1:9]), ncol=3)
		err[i] <- error(sample, R)
	}
	return(err)
}

#' Grid based optimization for user defined main direction of a rotation sample
#' 
#' @param sample sample of rotations
#' @param error user defined function to observed distance between sample and estimate, has to have parameters for the sample and the estimate
#' @param minerr minimal distance to consider for convergence
#' @param start starting value for the estimation
#' @param theta size of the grid considered  
#' @return list of 
#' \itemize{
#' \item Shat estimate of the main direction
#' \item iter number of iterations necessary for convergence
#' \item theta final size of the grid
#' \item minerr error used for convergence
#' \item error numeric value of total sample's distance from main direction
#' }
#' @export
#' @examples 
#' # minimize L1 norm:
#' L1.error <- function(sample, Shat) {
#'   sum(dist(sample, Shat, method="intrinsic", p=1))
#' }
#' 
#' cayley.sample <- ruars(n = 10, rangle = rcayley, nu = 1, space = 'SO3')
#' SL1 <- grid.search(cayley.sample, L1.error)
#' plot(cayley.sample, center=SL1$Shat, show_estimates="all")

grid.search <- function(sample, error, minerr =1e-5, start = mean(sample), theta=NULL) {
# 	if (length(start) == 1)
# 		Shat <- as.SO3(sample[start,])
# 	else if (all(dim(start) == 3))
# 		Shat <- as.SO3(start)
  Shat <- start
  
  err <- error(sample, Shat)
	if (is.null(theta)) theta <- 0.5*err/nrow(sample)
  
	iter <- 0
	while (theta > minerr) {
	  rA <- sphereA(Shat, theta)
	  rA.error <- error.grid(sample, Shat, theta, error)
	  if (min(rA.error, na.rm=TRUE) < err) {
	    Shat <- as.SO3(matrix(unlist(rA[which.min(rA.error),]), ncol=3))
	  } else {
	    theta <- theta/2
	  }
	  err <- error(sample, Shat)
	  iter <- iter+1
    if (is.na(theta) | is.na(minerr)) browser()
	}

#	print(paste("iterations:",iter))
	return(list(Shat=Shat, iter=iter, theta=theta, minerr=minerr, error=err))
}

