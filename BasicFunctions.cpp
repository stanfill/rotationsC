#include <RcppArmadillo.h>   
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]   
arma::vec meanQ4C(arma::mat Q) { 
	/*Compute the projected mean of the sample Q*/
	arma::mat Qsq=Q.t()*Q;
	arma::mat eigvec;
	arma::vec eigval;
  arma::eig_sym(eigval,eigvec,Qsq);   
  arma::vec qhat=eigvec.col(3);
  
  if(qhat[0]<0){
  	qhat = -qhat;
  }
  
  return qhat;
} 


// [[Rcpp::export]]
NumericVector bootQhat(NumericMatrix Q){
	
	int n=Q.nrow();
	NumericVector samp = runif(n);
	samp = floor(samp * n);
	
	arma::mat Qstar(n,4);
	arma::mat QSamp = as<arma::mat>(Q);
	
	for(int i=0;i<n;i++){
		Qstar.row(i) = QSamp.row(samp[i]);		//Copying a matrix row by row produces a bunch of junk messages
	}																				//so I do it with arma instead of standard Rcpp
	
	arma::vec QhatStar = meanQ4C(Qstar);
	
	return as<Rcpp::NumericVector>(wrap(QhatStar));
}

/*** R

#library(Rcpp)
#library(microbenchmark)
library(rotations)
Qs<-ruars(10,rcayley,space='Q4')
bootQhat(Qs)

#tim<-microbenchmark(          
#mean(Qs),
#meanQ4C(Qs)
#)

#print(tim)
#plot(tim)

#Q1<-matrix(rnorm(16),4,4)
#Q1%*%Q1
#qhatC(Q1)

*/

