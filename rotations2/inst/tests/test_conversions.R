theta <- acos(runif(1, -1, 1))
phi <- runif(1, -pi, pi)
u<- c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
r<-rvmises(1)

context("Conversions")
expect_equal(Q4(SO3(u,r)),Q4(u,r))
expect_equal(SO3(Q4(u,r)),SO3(u,r))