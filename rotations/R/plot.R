# set origin of concentric circles
origin <- matrix(as.SO3(c(1,-1,0), pi/16),3,3)

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


setOrigin <- function(origin = matrix(as.SO3(c(1,-1,0), pi/8),3,3)) {
	origin <<- origin
	pcircles <- data.frame(as.matrix(circles[,1:3]) %*% origin)
	pcircles
}


# this is the coordinate system and should be fixed, no matter what column of the rotation matrices is shown

base <- ggplot(aes(x=X1, y=X2), data=setOrigin(matrix(as.SO3(c(1,-1,0), pi/16),3,3))) +
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
	R <- list(matrix(as.SO3(c(0,1,0), pi/2),3,3), matrix(as.SO3(c(1,0,0), -pi/2),3,3), diag(c(1,1,1)))[[column]]
	rot <- center %*% R %*% origin
}


#' Project rotation data onto sphere
#'
#' Projection of rotation matrices onto sphere with given center.
#'
#' @param data data frame of rotation matrices in \eqn{3\times 3}{3-by-3} matrix representation.
#' @param center rotation matrix about which to center the observations.
#' @param column integer 1 to 3 indicating which column to display.
#' @return  Data frame with columns X, Y, Z standing for the respective coordinates in 3D space.
#' @export
#' @examples
#' Rs<-ruars(20, rcayley)
#'
#' #Project the sample's 3 axes onto the 3-shere centered at the identity rotation
#'
#' pointsXYZ(Rs, center = id.SO3, column = 1)  #x-axis
#' pointsXYZ(Rs, center = id.SO3, column = 2)  #y-axis
#' pointsXYZ(Rs, center = id.SO3, column = 3)  #z-axis

pointsXYZ <- function(data, center=id.SO3, column=1) {

  data<-as.SO3(data)
  #data<-matrix(data,ncol=9)
  center<-as.SO3(center)
  data<-data-center
  data<-matrix(data,ncol=9)

  idx <- list(1:3,4:6,7:9)[[column]]
  data <- matrix(data[,idx],ncol=3)

  psample1 <- data.frame(data)
  names(psample1) <- c("X","Y","Z")

  #  psample1 <- data.frame(psample1, data)
  #  psample1 <- psample1[order(psample1$Z, decreasing=FALSE),]
  psample1
}

pointsXYZ_plot <- function(data, center=id.SO3, column=1) {

  data<-as.SO3(data)
  data<-matrix(data,length(data)/9,9)
  center<-matrix(as.SO3(center),3,3)

	rot <- roteye(origin, center, column)
	idx <- list(1:3,4:6,7:9)[[column]]
	data <- matrix(data[,idx],ncol=3)

	psample1 <- data.frame(data %*% rot)
	names(psample1) <- c("X","Y","Z")

	#  psample1 <- data.frame(psample1, data)
	#  psample1 <- psample1[order(psample1$Z, decreasing=FALSE),]
	psample1
}



#' Visualizing random rotations
#'
#' This function produces a static three-dimensional globe onto
#' which  one of the  columns of the provided sample of rotations is projected.
#' The data are centered around a user-specified rotation matrix.  The static
#' plot uses \code{ggplot2}.  Interactive plots are no longer supported.
#'
#' @name plot
#'
#' @param x n rotations in \code{SO3} or \code{Q4} format.
#' @param center rotation about which to center the observations.
#' @param col integer or vector comprised of 1, 2, 3 indicating which column(s)
#'   to display.  If \code{length(col)>1} then each eyeball is labelled with the
#'   corresponding axis.
#' @param to_range logical; if \code{TRUE} only part of the globe relevant to
#'   the data is displayed
#' @param show_estimates character vector to specify  which of the four
#'   estimates of the principal direction to show. Possibilities are "all",
#'   "proj.mean", "proj.median", "geom.mean", "geom.median".
#' @param label_points  vector of labels.
#' @param mean_regions character vector to specify which of the three confidence
#'   regions to show for the projected mean.  Possibilities are "all",
#'   "trans.theory","trans.bootstrap, "direct.theory", "direct.bootstrap".
#' @param median_regions character vector to specify which of the three
#'   confidence regions to show for the projected median.  Possibilities are
#'   "all", "theory", "bootstrap."
#' @param alp alpha level to be used for confidence regions.  See
#'   \code{\link{region}} for more details.
#' @param m number of bootstrap replicates to use in bootstrap confidence
#'   regions.
#' @param interactive deprecated; \code{sphereplot} was set to be removed from CRAN
#'  and was going to take this package down with it
#' @param ... parameters passed onto the points layer.
#'
#' @return  A visualization of rotation data.
#'
#' @examples
#' r <- rvmises(200, kappa = 1.0)
#' Rs <- genR(r)
#' plot(Rs, center = mean(Rs), show_estimates = "proj.mean", shape = 4)
#'
#' \donttest{
#'   # Z is computed internally and contains information on depth
#'   plot(
#'     Rs,
#'     center = mean(Rs),
#'     show_estimates = c("proj.mean", "geom.mean"),
#'     label_points = sample(LETTERS, 200, replace = TRUE)
#'  ) +
#'    ggplot2::aes(size = Z, alpha = Z) +
#'    ggplot2::scale_size(limits = c(-1, 1), range = c(0.5, 2.5))
#' }
NULL

#' @rdname plot
#' @export
plot.SO3 <- function(x,
                     center = mean(x),
                     col = 1,
                     to_range = FALSE,
                     show_estimates = NULL,
                     label_points = NULL,
                     mean_regions = NULL,
                     median_regions = NULL,
                     alp = NULL,
                     m = 300,
                     interactive = FALSE,
                     ...) {


  if (length(col) > 1)
    mplotSO3(
      x = x,
      center = center,
      col = col,
      to_range = to_range,
      show_estimates = show_estimates,
      label_points = label_points,
      mean_regions = mean_regions,
      median_regions = median_regions,
      alp = alp,
      m = m,
      ...
    )
  else {
  Rs <- as.SO3(x)
	xlimits <- c(-1,1)
	ylimits <- c(-1,1)

	X <- Y <- Est <- NULL
  center<-matrix(center,3,3)

	proj2d <- pointsXYZ_plot(Rs, center=center, column=col)

	if(to_range) {
		xlimits <- range(proj2d$X)
		ylimits <- range(proj2d$Y)
		xbar <- mean(xlimits)
		xlimits <- xbar + 1.1*(xlimits-xbar)
		ybar <- mean(ylimits)
		ylimits <- ybar + 1.1*(ylimits-ybar)
	}

  #Static plot objects
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
		Shats<-stats::na.omit(Shats)

		if(nrow(Shats)==0){
      warning("Incorrect input to show_estimates")
      show_estimates<-NULL
		}

		Estlabels<-Estlabels[c(rmNA,NAs)]

		if(!is.null(mean_regions) || !is.null(median_regions)){
			vals<-3:(2+nrow(Shats)) #Make the shapes noticable, 15:18

      
      estDF<-pointsXYZ_plot(Shats[,1:9], center=center, column=col)
			estimates <- list(geom_point(aes(x=X, y=Y, shape=Est),size=3.5, data=data.frame(estDF, Shats)),
												scale_shape_manual(name="Estimates", labels=Estlabels,values=vals))
      
		}else{
      
		  estDF<-pointsXYZ_plot(Shats[,1:9], center=center, column=col)
		  estimates <- list(geom_point(aes(x=X, y=Y, colour=Est),size=3.5, data=data.frame(estDF, Shats)),
												scale_colour_brewer(name="Estimates", palette="Paired", labels=Estlabels))
      
		}
	}

	if (!is.null(mean_regions)) {
	  prentr <- fishr <- changr <- zhangr  <- NA
	  if(any(mean_regions%in%c('all','All'))) mean_regions<-c("trans.theory","trans.bootstrap","direct.theory","direct.bootstrap")
	  if (length(grep("trans.theory", mean_regions)) > 0) prentr<-region(Rs,estimator='mean',method='trans',type='theory',alp=alp)[col]
	  if (length(grep("trans.bootstrap", mean_regions)) > 0) fishr<-region(Rs,estimator='mean',method='trans',type='bootstrap',alp=alp,m=m)
    if (length(grep("direct.theory", mean_regions)) >0)    changr<-region(Rs,estimator='mean',method='direct',type='theory',alp=alp)
	  if (length(grep("direct.bootstrap", mean_regions)) > 0)    zhangr<-region(Rs,estimator='mean',method='direct',type='bootstrap',alp=alp,m=m)

	  Regions<-data.frame(X1=c(prentr,fishr,changr,zhangr),Meth=c('Mean\nTrans. Theory','Mean\nTrans. Bootstrap','Mean\nDirect Theory','Mean\nDirect Bootstrap'),
                        Meth2=c('Mean Trans. Theory','Mean Trans. Boot.','Mean Direct Theory','Mean Direct Boot.'))
	  Regions <- stats::na.omit(Regions)

	  if(nrow(Regions)==0){
	    warning("Incorrect input to mean_regions")
	    mean_regions<-NULL
	  }

    cisp.boot<-NULL

    for(i in 1:nrow(Regions)){
      if(col==1)
        cisp.boot <- rbind(cisp.boot,t(replicate(500, as.SO3(c(0,stats::runif(2,-1,1)), Regions$X1[i]),simplify="matrix")))

      if(col==2)
        cisp.boot <- rbind(cisp.boot,t(replicate(500, as.SO3(c(stats::runif(1,-1,1),0,stats::runif(1,-1,1)), Regions$X1[i]),simplify="matrix")))

      if(col==3)
	      cisp.boot <- rbind(cisp.boot,t(replicate(500, as.SO3(c(stats::runif(2,-1,1),0), Regions$X1[i]),simplify="matrix")))
    }
    
	  meanregDF<-pointsXYZ_plot(cisp.boot, center=t(mean(Rs))%*%center, column=col)
	  regs <- geom_point(aes(x=X, y=Y,colour=Regions), data=data.frame(meanregDF,Regions=rep(Regions$Meth,each=500)))
    
	}

	if (!is.null(median_regions)) {
		changr <- zhangr  <- NA
		if(any(median_regions%in%c('all','All'))) median_regions<-c("bootstrap","theory")
		if (length(grep("heory", median_regions)) >0)    changr<-region(Rs,method='direct',type='theory',estimator='median',alp=alp)
		if (length(grep("ootstrap", median_regions)) > 0)    zhangr<-region(Rs,method='direct',type='bootstrap',estimator='median',alp=alp,m=m)

		MedRegions<-data.frame(X1=c(changr,zhangr),Meth=c('Median Theory','Median Bootstrap'))
		MedRegions <- stats::na.omit(MedRegions)

		if(nrow(MedRegions)==0){
		  warning("Incorrect input to median_regions")
		  median_regions<-NULL
		}

		cisp.boot<-NULL

		for(i in 1:nrow(MedRegions)){
			if(col==1)
				cisp.boot <- rbind(cisp.boot,t(replicate(500, as.SO3(c(0,stats::runif(2,-1,1)), MedRegions$X1[i]),simplify="matrix")))

			if(col==2)
				cisp.boot <- rbind(cisp.boot,t(replicate(500, as.SO3(c(stats::runif(1,-1,1),0,stats::runif(1,-1,1)), MedRegions$X1[i]),simplify="matrix")))

			if(col==3)
				cisp.boot <- rbind(cisp.boot,t(replicate(500, as.SO3(c(stats::runif(2,-1,1),0), MedRegions$X1[i]),simplify="matrix")))
		}
		medianregDF<-pointsXYZ_plot(cisp.boot, center=t(median(Rs))%*%center, column=col)
		regsMed <- geom_point(aes(x=X, y=Y,colour=Regions), data=data.frame(medianregDF,Regions=rep(MedRegions$Meth,each=500)))
    
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
}

#' @rdname plot
#' @export
plot.Q4 <- function(x,
                    center = mean(x),
                    col = 1,
                    to_range = FALSE,
                    show_estimates = NULL,
                    label_points = NULL,
                    mean_regions = NULL,
                    median_regions = NULL,
                    alp = NULL,
                    m = 300,
                    interactive = FALSE,
                    ...) {
  Rs <- as.SO3(x)
  center <- as.SO3(center)
  graphics::plot(
    x = Rs,
    center = center,
    col = col,
    to_range = to_range,
    show_estimates = show_estimates,
    label_points = label_points,
    mean_regions = mean_regions,
    median_regions = median_regions,
    alp = alp,
    m = m,
    interactive = interactive,
    ...
  )
}

# Function written by Luciano Selzer and published on Stackoverflow on Aug 9 2012
# and edited by user "sebastian-c". It removes the guide from a ggplot2 object
# that can then be drawn by calling "grid.draw()" on what is returned.
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) {
    legend <- tmp$grobs[[leg]]
    return(legend)
  } else
    return(NULL)
}

# If more then one column is called in plot.SO3 then this is called to
# independently create an eyeball for each column and print them in a single row
# with the legend at the end, if applicable.
mplotSO3 <- function(x,
                     center = mean(x),
                     col = 1,
                     to_range = FALSE,
                     show_estimates = NULL,
                     label_points = NULL,
                     mean_regions = NULL,
                     median_regions = NULL,
                     alp = NULL,
                     m = 300,
                     ...) {
  p4 <- NULL

  if (1 %in% col) {
    p1 <- graphics::plot(
      x,
      center = center,
      col = 1,
      to_range = to_range,
      show_estimates = show_estimates,
      label_points = label_points,
      mean_regions = mean_regions,
      median_regions = median_regions,
      alp = alp,
      m = m,
      ...
    )
    p1 <- p1 +
      theme(axis.title.x = element_text(size = rel(1.5))) +
      labs(x = "x-axis")
    p4 <- g_legend(p1)
    p1 <- p1 + theme(legend.position = 'none')
  } else
    p1 <- NULL

  if (2 %in% col) {
    p2 <- graphics::plot(
      x,
      center = center,
      col = 2,
      to_range = to_range,
      show_estimates = show_estimates,
      label_points = label_points,
      mean_regions = mean_regions,
      median_regions = median_regions,
      alp = alp,
      m = m,
      ...
    )
    p2 <- p2 +
      theme(axis.title.x = element_text(size = rel(1.5))) +
      labs(x = "y-axis")
    p4 <- g_legend(p2)
    p2 <- p2 + theme(legend.position = 'none')
  } else
    p2 <- NULL

  if (3 %in% col) {
    p3 <- graphics::plot(
      x,
      center = center,
      col = 3,
      to_range = to_range,
      show_estimates = show_estimates,
      label_points = label_points,
      mean_regions = mean_regions,
      median_regions = median_regions,
      alp = alp,
      m = m,
      ...
    )
    p3 <- p3 +
      theme(axis.title.x = element_text(size = rel(1.5))) +
      labs(x = "z-axis")
    p4 <- g_legend(p3)
    p3 <- p3 + theme(legend.position = 'none')
  } else
    p3 <- NULL

  ps <- list(p1, p2, p3, p4)
  ps <- !sapply(ps, is.null)

  if (all(ps == c(TRUE, TRUE, TRUE, TRUE)))
    gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2, widths = c(2, 2, 2, 1))
  else if (all(ps == c(TRUE, TRUE, FALSE, TRUE)))
    gridExtra::grid.arrange(p1, p2, p4, nrow = 1, widths = c(2, 2, 1))
  else if (all(ps == c(TRUE, FALSE, TRUE, TRUE)))
    gridExtra::grid.arrange(p1, p3, p4, nrow = 1, widths = c(2, 2, 1))
  else if (all(ps == c(FALSE, TRUE, TRUE, TRUE)))
    gridExtra::grid.arrange(p2, p3, p4, nrow = 1, widths = c(2, 2, 1))
  else if (all(ps == c(TRUE, TRUE, TRUE, FALSE)))
    gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
  else if (all(ps == c(TRUE, TRUE, FALSE, FALSE)))
    gridExtra::grid.arrange(p1, p2, nrow = 1)
  else if (all(ps == c(TRUE, FALSE, TRUE, FALSE)))
    gridExtra::grid.arrange(p1, p3, nrow = 1)
  else if (all(ps == c(FALSE, TRUE, TRUE, FALSE)))
    gridExtra::grid.arrange(p2, p3, nrow = 1)
  else
    stop("Specify the columns correctly.")
}
