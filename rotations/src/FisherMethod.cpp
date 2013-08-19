#include <RcppArmadillo.h>   
#include <Rcpp.h>
#include "../inst/include/rotations.h"
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
      
			denom = pow((n*(eigval[3]-eta[j])*(eigval[3]-eta[k])),-1);
			
			for(i=0;i<n;i++){
				G(j,k) = G(j,k) + arma::as_scalar(Qs.row(i)*Mhat.col(j)*Qs.row(i)*Mhat.col(k)*pow(Qs.row(i)*eigvec.col(3),2));
			}
      
      G(j,k) = G(j,k)*denom;
      G(k,j)=G(j,k);
		}
	}
  
  arma::mat Ginv = G.i();  
  
  Tm = arma::as_scalar(n*Qhat*Mhat*Ginv*Mhat.t()*Qhat.t());
  
  return Tm;
}


// [[Rcpp::export]]  
double fisherAxisCSymmetric(arma::mat Qs, arma::rowvec Qhat){

	//This is the same as fisherAxisC but the test statistic is much reduced
	//See equation 10 of Fisher et. al. (1996)

	NumericMatrix Qss = as<NumericMatrix>(wrap(Qs));

  int n=Qs.n_rows;
  
	arma::mat Qsq=(Qs.t()*Qs)/n;
	arma::mat eigvec, Mhat(4,3);
	arma::vec eigval, eta(3);
	
	//arma::mat33 unM;
  
  arma::eig_sym(eigval,eigvec,Qsq);   

  int i, j;
  double Tm, denom, trGhat=0.0, Gjj=0.0, Tnum=0.0;

  
  for(i=0;i<3;i++){
    Mhat.col(i)=eigvec.col(i);
    eta(i)=eigval(i);
  }
  
  for(j=0;j<3;j++){

			denom = pow((n*(eigval[3]-eta[j])*(eigval[3]-eta[j])),-1);
			
			for(i=0;i<n;i++){
				Gjj = Gjj + arma::as_scalar(Qs.row(i)*Mhat.col(j)*Qs.row(i)*Mhat.col(j)*pow(Qs.row(i)*eigvec.col(3),2));
			}
      
      Gjj = Gjj*denom;
      trGhat = trGhat + Gjj;
      Gjj = 0.0;
      
      Tnum += arma::as_scalar(pow(Qhat*Mhat.col(j),2));
		
	}
  
  trGhat = pow(trGhat,-1);
  
  //unM = 3*n*Mhat.t()*Mhat*trGhat;
  //unM.print("stat: ");
  
  Tm = arma::as_scalar(trGhat*(3*n)*Tnum);
  
  return Tm;
}



//[[Rcpp::export]]
arma::vec fisherBootC(arma::mat Qs, int m, bool symm){

  int n = Qs.n_rows;
  int i , j , numUn;
  
  arma::rowvec qhat = rotations::meanQ4C(Qs);
	//arma::rowvec qhat = meanQ4C(Qs);

  arma::vec Tm(m);
  NumericVector unSamp;
  IntegerVector samp(n);
  arma::mat Qstar(n,4);
  Qstar.zeros();
  
  for(i=0;i<m;i++){
    
    samp = floor(runif(n,0,n));			//Bootstrap sample of size n, with replacement
	  unSamp = unique(samp);
    numUn = unSamp.size();
    
    while(numUn<4){
      samp = floor(runif(n,0,n));	 //If bootstrap samp is less than 4 obs then							
	    unSamp = unique(samp);       //draw a new sample
      numUn = unSamp.size();
    }
    
    for(j=0;j<n;j++){
      Qstar.row(j) = Qs.row(samp[j]);
    }
    
    if(symm){
    	Tm[i]=fisherAxisCSymmetric(Qstar,qhat);
    }else{
    	Tm[i]=fisherAxisC(Qstar,qhat);
    }
  }
  return Tm;
}
