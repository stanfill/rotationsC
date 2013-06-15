#include <RcppArmadillo.h>   
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]
arma::mat Q4defaultC(arma::mat U, arma::vec theta){
	
	int n1 = U.n_rows, n = theta.size();
	
	arma::mat q;
	q.zeros(n,4);

	if(n1 != n){
		printf("Error, u and theta different length");
		return q;
	}
	
	arma::vec stheta = sin(theta/2);
	//q.row(0).print("Q0: ");
	q.col(0) = cos(theta/2);
	q.col(1) = U.col(0) % stheta;
	q.col(2) = U.col(1) % stheta;
	q.col(3) = U.col(2) % stheta;
 
	return q;
	
}

// [[Rcpp::export]]   
arma::rowvec meanQ4C(arma::mat Q) { 
	//Compute the projected mean of the sample Q
	
	NumericMatrix Qss = as<NumericMatrix>(wrap(Q));
	
	arma::mat Qsq=Q.t()*Q;
	arma::mat eigvec;
	arma::vec eigval;
  arma::eig_sym(eigval,eigvec,Qsq);   
  arma::vec qhat=eigvec.col(3);
  
  if(qhat[0]<0){
  	qhat = -qhat;
  }
  
  return qhat.t(); //Want to return it in a row vector so transpose it
}

/*// [[Rcpp::export]]
arma::mat SO3defaultC(arma::mat U, arma::vec theta) {
  
  //U is an n-by-3 matrix, each row is a misorentation axis
  //theta is a vector of length n, each item is a misorientation angle
  
  int n=U.n_rows, i=0;	
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
}*/
