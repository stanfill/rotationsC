/* To source run- Rcpp::sourceCpp('cayley.cpp') */

#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcayleyC(int n, double kappa){
	
	double temp = 0;	
	NumericVector mat(n);
	for(int i = 0; i<n ;i++){
		
		temp = rbeta(1,kappa+0.5,1.5)[0];
		mat[i] = acos(2*temp-1)*(1-2*rbinom(1,1,0.5)[0]);
		
	}
	return mat;
}


/*** R

#Apparently C is not required to speed rcayley up any

library(rotations)
library(microbenchmark)

res<-microbenchmark(
rcayleyC(100,10),
rcayley(100,10)
)
print(res)
boxplot(res)
*/
