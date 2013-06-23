rs<-rcayley(20)
Rs<-genR(rs)

context("distance metrics")
expect_equal(dist(Rs,method='intrinsic'),abs(rs))
expect_equal(dist(Rs,method='projected'),sqrt(8)*sin(abs(rs)/2))