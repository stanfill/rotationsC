#include <RcppArmadillo.h>   
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]  
NumericVector fisherAxisC(NumericMatrix Qs){
	
	arma::mat armaQs=as<arma::mat>(Qs);
	
	arma::mat Qsq=armaQs.t()*armaQs;
	arma::mat eigvec;
	arma::vec eigval;
  arma::eig_sym(eigval,eigvec,Qsq);   
  arma::vec qhat=eigvec.col(3);
	
  arma::vec m=eigvec.col(3);
  
  if(m[0]<0){
  	m = -m;
  }
	
}