library(ggplot2)
require(grid)

#this had to be added because the .Call version of SO3.default can't be called from insider the function...I think
oldSO3 <- function(U, theta=NULL) {
	n<-length(U)/3
	if(n%%1!=0)
			stop("This functions only works in three dimensions.")	
	U<-matrix(U,n,3)
	ulen<-sqrt(rowSums(U^2)) 
	if(is.null(theta)){ 
			theta<-ulen%%(pi)
					
					#if(theta>pi)
					#	theta<-2*pi-theta
			}
	R<-matrix(NA,n,9)
 	for(i in 1:n){
					
				if(ulen[i]!=0)
						U[i,]<-U[i,]/ulen[i]
					
				P <- U[i,] %*% t(U[i,])
					
				R[i,] <- P + (diag(3) - P) * cos(theta[i]) + eskew(U[i,]) * sin(theta[i])
	}
	class(R) <- "SO3"
	return(R)
}

suppressMessages(library(ggplot2))
require(grid)

# set origin of concentric circles
origin <- matrix(oldSO3(c(1,-1,0), pi/16),3,3)

# construct helper grid lines for sphere

theta <- seq(0,pi, by=pi/8)
phi <- seq(0,2*pi, by=0.005)
df <- data.frame(expand.grid(theta=theta, phi=phi))

#qplot(theta,phi, geom="point", data=df) + coord_polar()

x <- with(df, sin(theta)*cos(phi))
y <- with(df, sin(theta)*sin(phi))
z <- with(df, cos(theta))
circles <- data.frame(cbind(x,y,z))
circles$ID <- as.numeric(factor(df$theta))

theta <- seq(0,pi, by=0.005)
phi <- seq(0,2*pi, by=pi/8) 
df <- data.frame(expand.grid(theta=theta, phi=phi))

x <- with(df, sin(theta)*cos(phi))
y <- with(df, sin(theta)*sin(phi))
z <- with(df, cos(theta))
circles.2 <- data.frame(cbind(x,y,z))
circles.2$ID <- as.numeric(factor(df$phi))+9

circles <- rbind(circles, circles.2)


setOrigin <- function(origin = matrix(oldSO3(c(1,-1,0), pi/8),3,3)) {
	origin <<- origin
	pcircles <- data.frame(as.matrix(circles[,1:3]) %*% origin)
	pcircles
}


# this is the coordinate system and should be fixed, no matter what column of the rotation matrices is shown

base <- ggplot(aes(x=X1, y=X2), data=setOrigin(matrix(oldSO3(c(1,-1,0), pi/16),3,3))) + 
	coord_equal() + 
	geom_point(aes(alpha=X3), size=0.6, colour="grey65") + 
	scale_alpha(range=c(0,0.8),  guide="none") + 
	theme(panel.background=element_blank(),
				panel.grid.minor=element_blank(),
				panel.grid.major=element_blank(),
				axis.title.x=element_blank(),
				axis.title.y=element_blank(),
				axis.text.x=element_blank(),
				axis.text.y=element_blank(),
				axis.ticks=element_blank(), 
				plot.margin = unit(rep(0, 4), "lines"))


roteye <- function(origin, center, column=1) {
	R <- list(matrix(oldSO3(c(0,1,0), pi/2),3,3), matrix(oldSO3(c(1,0,0), -pi/2),3,3), diag(c(1,1,1)))[[column]]
	rot <- center %*% R %*% origin 
}


#' Project rotation data onto sphere
#' 
#' Projection of rotation matrices onto sphere with given center.
#'
#' @param data data frame of rotation matrices in \eqn{3\times 3}{3-by-3} matrix representation.
#' @param center point about which to center the observations.
#' @param column integer 1 to 3 indicating which column to display.
#' @return  Data frame with columns X, Y, Z standing for the respective coordinates in 3D space.
#' @export
#' 
pointsXYZ <- function(data, center, column=1) {
	rot <- roteye(origin, center, column)
	idx <- list(1:3,4:6, 7:9)[[column]]
	data <- as.matrix(data[,idx])
	
	psample1 <- data.frame(data %*% rot)
	names(psample1) <- c("X","Y","Z")
	
	#  psample1 <- data.frame(psample1, data)
	#  psample1 <- psample1[order(psample1$Z, decreasing=FALSE),]
	psample1  
}


#' Visualizing random rotations.
#'
#' This function produces a three-dimensional globe onto which  one of the  columns of the provided sample of rotations is drawn.  The data are centered around a provided
#' matrix and the user can choose to display this center or not.  Based on \code{ggplot2} package by \cite{wickham09}.
#'
#' @param x n rotations in \code{SO3} format.
#' @param center rotation about which to center the observations.
#' @param col integer 1 to 3 indicating which column to display.
#' @param to_range show only part of the globe that is in range of the data?
#' @param show_estimates character vector to specify  which of the four estimates of the principal direction to show. Possibilities are "all", "proj.mean", "proj.median", "geom.mean", "geom.median."
#' @param label_points  vector of labels.
#' @param mean_regions character vector to specify which of the three confidence regions to show for the projected mean.  Possibilities are "all", "eigen theory","eigen bootstrap, "moment theory", "moment bootstrap."
#' @param median_regions character vector to specify which of the three confidence regions to show for the projected median.  Possibilities are "all", "theory", "bootstrap."
#' @param alp alpha level to be used for confidence regions.
#' @param m number of bootstrap replicates to use in Zhang confidence region.
#' @param ... parameters passed onto the points layer.
#' @return  A \code{ggplot2} object with the data displayed on spherical grid.
#' @S3method plot SO3
#' @method plot SO3
#' @cite wickham09
#' @export
#' @examples
#' r<-rvmises(200,1.0)
#' Rs<-genR(r)
#' plot(Rs,center=mean(Rs),show_estimates=NULL,shape=4)
#' # Z is computed internally and contains information on depth
#' plot(Rs,center=mean(Rs),show_estimates=c("proj.mean", "geom.mean"), 
#'  label_points=sample(LETTERS, 200, replace=TRUE)) + aes(size=Z, alpha=Z) + 
#'  scale_size(limits=c(-1,1), range=c(0.5,2.5))

plot.SO3 <- function(x, center=mean(x), col=1, to_range=FALSE, show_estimates=NULL, label_points=NULL, mean_regions=NULL, median_regions=NULL, alp=NULL, m=300,  ...) {

  Rs <- as.SO3(x)
	xlimits <- c(-1,1)
	ylimits <- c(-1,1)
	
	X <- Y <- Est <- NULL
  center<-matrix(center,3,3)
	proj2d <- pointsXYZ(Rs, center=center, column=col)
	if(to_range) {
		xlimits <- range(proj2d$X)
		ylimits <- range(proj2d$Y)
		xbar <- mean(xlimits)
		xlimits <- xbar + 1.1*(xlimits-xbar)
		ybar <- mean(ylimits)
		ylimits <- ybar + 1.1*(ylimits-ybar)
	}
	
	estimates <- NULL
  regs<-NULL
	regsMed<-NULL
	
	if (!is.null(show_estimates)) {
		ShatP <- StildeP <- ShatG <- StildeG <- NA
		if(any(show_estimates%in%c('all','All'))) show_estimates<-c("proj.mean","proj.median","geom.mean","geom.median")
		if (length(grep("proj.mean", show_estimates)) > 0) ShatP<-mean(Rs, type="projected")
		if (length(grep("proj.median", show_estimates)) >0)    StildeP<-median(Rs, type="projected")
		if (length(grep("geom.mean", show_estimates)) > 0)    ShatG<-mean(Rs, type="geometric")
		if (length(grep("geom.median", show_estimates)) > 0)    StildeG<-median(Rs, type="geometric")
		
		Shats<-data.frame(rbind(as.vector(ShatP),as.vector(StildeP),as.vector(ShatG),as.vector(StildeG)),Est=1:4)
		Shats$Est <- factor(Shats$Est)
		Estlabels <- c(expression(hat(S)[E]), expression(tilde(S)[E]), expression(hat(S)[R]), expression(tilde(S)[R]))
		
		
		levels(Shats$Est) <- Estlabels
		
		
		rmNA<-which(!is.na(Shats$X1))
		NAs<-c(1:4)[-rmNA]
		Shats<-na.omit(Shats)
		
		#Shats <- Shats[rmNA,]
		Estlabels<-Estlabels[c(rmNA,NAs)]
		
		if(!is.null(mean_regions) || !is.null(median_regions)){
			vals<-3:(2+nrow(Shats)) #Make the shapes noticable, 15:18
			estimates <- list(geom_point(aes(x=X, y=Y, shape=Est),size=3.5, data=data.frame(pointsXYZ(Shats, center=center, column=col), Shats)),
												scale_shape_manual(name="Estimates", labels=Estlabels,values=vals))
		}else{
			estimates <- list(geom_point(aes(x=X, y=Y, colour=Est),size=3.5, data=data.frame(pointsXYZ(Shats, center=center, column=col), Shats)),
												scale_colour_brewer(name="Estimates", palette="Paired", labels=Estlabels))
		}
	}
  
	if (!is.null(mean_regions)) {
	  prentr <- fishr <- changr <- zhangr  <- NA
	  if(any(mean_regions%in%c('all','All'))) mean_regions<-c("eigen theory","eigen bootstrap","moment theory","moment bootstrap")
	  if (length(grep("eigen theory", mean_regions)) > 0) prentr<-region(Rs,estimator='mean',method='eigen',type='theory',alp=alp)[col]
	  if (length(grep("eigen bootstrap", mean_regions)) > 0) fishr<-region(Rs,estimator='mean',method='eigen',type='bootstrap',alp=alp,m=m)
    if (length(grep("moment theory", mean_regions)) >0)    changr<-region(Rs,estimator='mean',method='moment',type='theory',alp=alp)
	  if (length(grep("moment bootstrap", mean_regions)) > 0)    zhangr<-region(Rs,estimator='mean',method='moment',type='bootstrap',alp=alp,m=m)

	  Regions<-data.frame(X1=c(prentr,fishr,changr,zhangr),Meth=c('Mean\nEigen Theory','Mean\nEigen Bootstrap','Mean\nMoment Theory','Mean\nMoment Bootstrap'))
	  Regions <- na.omit(Regions)
	  
    cisp.boot<-NULL
    
    for(i in 1:nrow(Regions)){
      if(col==1)
        cisp.boot <- rbind(cisp.boot,t(replicate(500, oldSO3(c(0,runif(2,-1,1)), Regions$X1[i]),simplify="matrix")))
      
      if(col==2)
        cisp.boot <- rbind(cisp.boot,t(replicate(500, oldSO3(c(runif(1,-1,1),0,runif(1,-1,1)), Regions$X1[i]),simplify="matrix")))
      
      if(col==3)
	      cisp.boot <- rbind(cisp.boot,t(replicate(500, oldSO3(c(runif(2,-1,1),0), Regions$X1[i]),simplify="matrix")))
    }
	  
	  regs <- geom_point(aes(x=X, y=Y,colour=Regions), data=data.frame(pointsXYZ(cisp.boot, center=t(mean(Rs))%*%center, column=col),Regions=rep(Regions$Meth,each=500)))

	}
  
	if (!is.null(median_regions)) {
		changr <- zhangr  <- NA
		if(any(median_regions%in%c('all','All'))) median_regions<-c("bootstrap","theory")
		if (length(grep("heory", median_regions)) >0)    changr<-region(Rs,method='moment',type='theory',estimator='median',alp=alp)
		if (length(grep("ootstrap", median_regions)) > 0)    zhangr<-region(Rs,method='moment',type='bootstrap',estimator='median',alp=alp,m=m)
		
		MedRegions<-data.frame(X1=c(changr,zhangr),Meth=c('Median\nMoment Theory','Median\nMoment Bootstrap'))
		MedRegions <- na.omit(MedRegions)
		
		cisp.boot<-NULL
		
		for(i in 1:nrow(MedRegions)){
			if(col==1)
				cisp.boot <- rbind(cisp.boot,t(replicate(500, oldSO3(c(0,runif(2,-1,1)), MedRegions$X1[i]),simplify="matrix")))
			
			if(col==2)
				cisp.boot <- rbind(cisp.boot,t(replicate(500, oldSO3(c(runif(1,-1,1),0,runif(1,-1,1)), MedRegions$X1[i]),simplify="matrix")))
			
			if(col==3)
				cisp.boot <- rbind(cisp.boot,t(replicate(500, oldSO3(c(runif(2,-1,1),0), MedRegions$X1[i]),simplify="matrix")))
		}
		
		regsMed <- geom_point(aes(x=X, y=Y,colour=Regions), data=data.frame(pointsXYZ(cisp.boot, center=t(median(Rs))%*%center, column=col),Regions=rep(MedRegions$Meth,each=500)))
		
	}
	
	labels <- NULL
	if (!is.null(label_points)) {
		proj2d$labels <- label_points
		labels <- geom_text(aes(x=X+0.05, y=Y, label=labels), size=3.25, data=proj2d, ...) 
	}
	base + geom_point(aes(x=X, y=Y), data=proj2d, ...) + 
		labels + 
		estimates +
    regs+
		regsMed+
		xlim(xlimits) + ylim(ylimits)
}

