Rs<-ruars(20,rcayley)

context('plot')

expect_error(plot(Rs),"'center' is missing")
expect_error(plot(Rs,center=id.SO3,show_estimates=T),"undefined columns selected")