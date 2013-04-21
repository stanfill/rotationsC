library(Rcpp)

cppFunction('
	int add(int x, int y, int z){
		int sum = x+y+z;
		return sum;
	}'
)

add(1,2,3)

Rcpp::sourceCpp("testFns.cpp")

convolveCpp(c(1,4,3),c(2,3,4),2)

Rcpp::sourceCpp('svd.c')

