#include <RcppArmadillo.h>   
#include <Rcpp.h>
//#include "../inst/include/rotations2.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::mat eskewC(arma::rowvec U) {
  
  double ulen=norm(U,2);
  
  U = U/ulen;
  
  double u = U(0);
  double v = U(1);
  double w = U(2);
  
  arma::mat33 res;
  res.zeros();
  res(0,1) = -w;
  res(1,0) = w;
  res(0,2) = v;
  res(2,0) = -v;
  res(1,2) = -u;
  res(2,1) = u;
  
  /*res <- matrix((-1) *  c(0, -w, v,
                          w, 0, -u,
                          -v, u, 0), ncol = 3)*/
  return res;
}

// [[Rcpp::export]]
arma::mat SO3defaultC(arma::mat U, arma::vec theta) {
  
  //This function expects U to be n-by-3 and theta to be a vector of length n 
  //each row of U needs to be length 1
  
  int n=U.n_rows, i=0;	
  //int n=1, i=0;
	arma::mat Ri(3,3), Rs(n,9), I(3,3), SS(3,3);
  Ri.zeros(); Rs.zeros(); I.eye();
  arma::rowvec Rir;

  for(i=0;i<n;i++){
  	
 		Ri = U.row(i).t() * U.row(i);
    SS = eskewC(U.row(i));
  	Ri = Ri + (I - Ri) * cos(theta[i]) +  SS * sin(theta[i]);
    Rs.row(i) = as<arma::rowvec>(wrap(Ri));

  }
 		
  return Rs;
}

// [[Rcpp::export]]
arma::mat genrC(arma::vec r, arma::mat S , int SO3) {

  int n=r.size(), i=0;
  
  NumericVector theta = runif(n,-1,1);
  theta = acos(theta);
    
  NumericVector phi = runif(n, -M_PI, M_PI);
  arma::mat u(n,3);
  
  for(i=0;i<n;i++){
    u(i,0)=sin(theta[i]) * cos(phi[i]);
    u(i,1)=sin(theta[i]) * sin(phi[i]);
    u(i,2)=cos(theta[i]);
  }  
  
  if(SO3==1){
    
    int j;
    arma::mat Rs(n,9);
    arma::mat33 Rsi;

    Rs = SO3defaultC(u, r);
    
    
    for(i=0;i<n;i++){
      
      for(j=0;j<9;j++){
        Rsi(j) = Rs(i,j);
      }
      
      Rsi = S * Rsi;
      Rs.row(i) = as<arma::rowvec>(wrap(Rsi));
      
    }
      
    return Rs;
  }else{
    return u;
  }

}