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

/*NumericMatrix genrC <- function(NumericVector r, NumericMatrix S , int SO3) {

  
  int n=r.size(), i=0;
  
  NumericVector theta = runif(n,-1,1);
  theta = acos(theta);
    
  NumericVector phi = runif(n, -pi, pi);
  NumericMatrix u(n,3);
  
  for(i=0;i<n;i++){
    u[i,0]=sin(theta) * cos(phi);
    u[i,1]=sin(theta) * sin(phi);
    u[i,2]=cos(theta);
  }  
  if(SO3==1){
    
  	if(is.null(S))
  		S<-id.SO3
  	S<-formatSO3(S)
  	o<-SO3(u,r)
  	o<-centeringSO3(o,t(S))
  	
  	class(o) <- "SO3"
  	return(o)
  	
  }else{
  	
  	if(is.null(S))
  		S<-id.Q4
  	
  	S<-formatQ4(S)
  	
  	S[2:4]<--S[2:4]
  	q<-matrix(c(cos(r/2),sin(r/2)*u),n,4)
  	q<-centeringQ4(q,S)
  	
  	class(q)<-"Q4"
  	return(q)
  }

}*/