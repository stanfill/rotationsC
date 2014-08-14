
#include <RcppArmadillo.h>   

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]
arma::rowvec HnCpp(arma::mat Qs){
  //Compute the Hn tests statistics
  
  int n = Qs.n_rows, i=0;
  arma::mat T = Qs.t()*Qs;
  arma::mat eigvec, eigvecJ;
  arma::vec eigval, eigvalJ;
  arma::eig_sym(eigval,eigvec,T);
  arma::rowvec Hn(n);
  arma::rowvec Qj;
  arma::mat Tj;

  for(i = 0;i<n; i++){
    Qj = Qs.row(i);
    
    Tj = T-Qj.t()*Qj;
    arma::eig_sym(eigvalJ,eigvecJ,Tj);
    Hn(i)=(n-2)*(1+eigvalJ(3)-eigval(3))/(n-1-eigvalJ(3));
    
  }
  return Hn;
}

// [[Rcpp::export]]
arma::rowvec HnCppBloc(arma::mat Qs, arma::mat Cs){
  //Compute the Hn tests statistics
  
  int n = Qs.n_rows, i = 0, j = 0, rowNum = 0;
  int t = Cs.n_rows, nc = Cs.n_cols;
  //printf("nc: %i",nc);  
  //printf("t: %i",t);
  //Cs.print("Cs: ");
  
  arma::mat T = Qs.t()*Qs;
  arma::mat eigvec, eigvecJ;
  arma::vec eigval, eigvalJ;
  arma::eig_sym(eigval,eigvec,T);
  arma::rowvec Hn(nc);
  arma::mat Qj;
  Qj.zeros(t,4);
  arma::mat::fixed<4,4> Tj;

  for(i = 0;i<nc; i++){
    for(j = 0; j<t ; j++){
      
      rowNum = Cs(j,i) - 1;
      Qj.row(j) = Qs.row(rowNum);
    
      Tj = T-Qj.t()*Qj;
      arma::eig_sym(eigvalJ,eigvecJ,Tj);
      Hn(i)=(n-t-1)*(t+eigvalJ(3)-eigval(3))/(t*(n-t-eigvalJ(3)));
    }
    
  }
  return Hn;
}

