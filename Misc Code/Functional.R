library(rotations)
Qs<-ruars(1,rcayley,space='Q4')
y<-ruars(1,rcayley,space='Q4')

y%*%t(Qs)%*%Qs%*%t(y)

sum(diag(t(Qs)%*%y))^2
