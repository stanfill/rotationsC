#include <RcppArmadillo.h>   
#include <Rcpp.h>
//#include "../inst/include/rotations2.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]  
double fisherAxisC(arma::mat Qs, arma::rowvec Qhat){

	NumericMatrix Qss = as<NumericMatrix>(wrap(Qs));
	
  int n=Qs.n_rows;
  
	arma::mat Qsq=(Qs.t()*Qs)/n;
	arma::mat eigvec, Mhat(4,3);
	arma::vec eigval, eta(3);
  
  arma::eig_sym(eigval,eigvec,Qsq);   

	arma::mat G(3,3);
	G.zeros();

  int i, j, k;
  double Tm, denom;

  
  for(i=0;i<3;i++){
    Mhat.col(i)=eigvec.col(i);
    eta(i)=eigval(i);
  }
  
  for(j=0;j<3;j++){
		for(k=j;k<3;k++){
      
			denom = pow((n*(eigval(3)-eta(j))*(eigval(3)-eta(k))),-1);
			
			for(i=0;i<n;i++){
				G(j,k) = G(j,k) + arma::as_scalar(Qs.row(i)*Mhat.col(j)*Qs.row(i)*Mhat.col(k)*pow(Qs.row(i)*eigvec.col(3),2));
			}
      
      G(j,k) = G(j,k)*denom;
      G(k,j) = G(j,k);
		}
	}
  
  arma::mat Ginv = G.i();  
  
  Tm = arma::as_scalar(n*Qhat*Mhat*Ginv*Mhat.t()*Qhat.t());
  
  return Tm;
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


//[[Rcpp::export]]
arma::vec fisherBootC(arma::mat Qs, int m){

  int n = Qs.n_rows;
  int i , j , numUn=0, maxSamp=0;
  
  //arma::rowvec qhat = rotations2::meanQ4C(Qs);
	arma::rowvec qhat = meanQ4C(Qs);

  arma::vec Tm(m);
  NumericVector unSamp;
  IntegerVector samp(n);
  arma::mat Qstar(n,4);
  Qstar.zeros();
  
  for(i=0;i<m;i++){
    
    samp = floor(runif(n,0,n));			//Bootstrap sample of size n, with replacement
	  unSamp = unique(samp);
    numUn = unSamp.size();
    maxSamp = max(samp);
    
    while(numUn<4 || maxSamp>n-1){
      samp = floor(runif(n,0,n));	 //If bootstrap samp is less than 4 obs then							
	    unSamp = unique(samp);       //draw a new sample
      numUn = unSamp.size();
      maxSamp = max(samp);
    }
    
    for(j=0;j<n;j++){
      Qstar.row(j) = Qs.row(samp[j]);
    }
    
    Tm(i)=fisherAxisC(Qstar,qhat);
    
  }
  return Tm;
}

