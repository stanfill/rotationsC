theta <- acos(runif(1, -1, 1))
phi <- runif(1, -pi, pi)
u<- c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
r<-rvmises(1)

context("Conversions")
expect_equal(as.Q4(as.SO3(u,r)),as.Q4(u,r))
expect_equal(as.SO3(as.Q4(u,r)),as.SO3(u,r))