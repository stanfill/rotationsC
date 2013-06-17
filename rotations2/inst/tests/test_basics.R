theta <- acos(runif(1, -1, 1))
phi <- runif(1, -pi, pi)
u<- matrix(c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)),1,3)
r<-rvmises(1)

R<-SO3(u,r)
Q<-Q4(u,r)

context("Basics")
#Does the angle function extract the correct angle from the rotation
expect_equal(angle(R),abs(r))
expect_equal(angle(Q),abs(r))

#Does the axis function extract the correct axis from the rotation
expect_equal(abs(axis2(R)),abs(u))
expect_equal(abs(axis2(Q)),abs(u))

#Can we recreate the rotation matrix from the axis and angle functions
expect_equal(SO3(axis2(R),angle(R)),as.SO3(matrix(R,1,9)))
expect_equal(Q4(axis2(Q),angle(Q)),Q)

#Make sure we can use angle-axis represnetation to generate a single rotation
expect_equal(SO3(u*r),SO3(u,r))
expect_equal(Q4(u*r),Q4(u,r))

#Check it works for a vector
rs<-rvmises(20)
theta <- acos(runif(20, -1, 1))
phi <- runif(20, -pi, pi)
us<- matrix(c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)),20,3)
Rs<-SO3(us,rs)
Qs<-Q4(us,rs)

#Does the angle function extract the correct angle from a sample of rotations
expect_equal(angle(Rs),abs(rs))
expect_equal(angle(Qs),abs(rs))

#Does the axis function extract the correct axis from a sample of rotations
expect_equal(SO3(axis2(Rs),angle(Rs)),Rs)
expect_equal(Q4(axis2(Qs),angle(Qs)),Qs)


#Make sure we can use angle-axis represnetation to generate a sample of rotations
expect_equal(SO3(us*rs),Rs)
expect_equal(Q4(us*rs),Qs)