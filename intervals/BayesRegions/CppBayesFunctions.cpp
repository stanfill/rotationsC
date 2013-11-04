// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>  

using namespace Rcpp;


// [[Rcpp::export]]
arma::mat centerCpp(arma::mat Rs, arma::mat S){
  //Center the dataset Rs around S, so each row of Rs denoted R is replaces with S'R
  
  int n = Rs.n_rows;
  arma::mat Rsi(3,3), cRs(n,9);
  int i = 0, j=0;
  
  for(i=0;i<n;i++){
    
    for(j=0;j<9;j++){
     Rsi[j] = Rs(i,j); 
    } 
    
    Rsi = S.t()*Rsi;
    
    for(j=0;j<9;j++){
      cRs(i,j) = Rsi[j];
    }
    
  }
  return cRs;
}

// [[Rcpp::export]]
double gvmUARSC(arma::mat Rs, arma::mat S, double kappa){
  
  int n = Rs.n_rows;
  arma::mat cRs = centerCpp(Rs,S);
  arma::mat trcRs(n,3);
  
  trcRs.col(0)=cRs.col(0);
  trcRs.col(1)=cRs.col(4);
  trcRs.col(2)=cRs.col(8);
  
  arma::colvec traces = sum(trcRs,1);

  double n1 = exp(kappa*sum(traces-1)/2);
  
  double I0k = R::bessel_i(kappa,0,1);
  double I1k = R::bessel_i(kappa,1,1);
  
  double n2 = sqrt(pow(I0k,2)-(I0k*I1k/kappa)-pow(I1k,2));
  double d1 = pow(I0k,(n+1));
  arma::colvec d2 = 3-traces;
  double d3 = 1;
  for(int i=0; i<n ; i++){
    d3 *= d2[i];
  }
  
  return (n1*n2)/(d1*d3);
}


// [[Rcpp::export]]
double gfUARSC(arma::mat Rs, arma::mat S, double kappa){
  
  int n = Rs.n_rows;
  arma::mat cRs = centerCpp(Rs,S);
  arma::mat trcRs(n,3);
  
  trcRs.col(0)=cRs.col(0);
  trcRs.col(1)=cRs.col(4);
  trcRs.col(2)=cRs.col(8);
  
  arma::colvec traces = sum(trcRs,1);
  
  double n1 = exp(kappa*sum(traces)-n);
  double I02k = R::bessel_i(2*kappa,0,1);
  double I12k = R::bessel_i(2*kappa,1,1);
  double n2 = sqrt(2*pow(I02k,2)/kappa-2*I02k*I12k/pow(kappa,2)+((1/pow(kappa,2))-(2/kappa))*pow(I12k,2)); 
  double d1 = pow(I02k-I12k,n+1);
  return (n1*n2)/d1;
}



// [[Rcpp::export]]
double gcayUARSC(arma::mat Rs, arma::mat S, double kappa){
 
   int n = Rs.n_rows;
  arma::mat cRs = centerCpp(Rs,S);
  arma::mat trcRs(n,3);
  
  trcRs.col(0)=cRs.col(0);
  trcRs.col(1)=cRs.col(4);
  trcRs.col(2)=cRs.col(8);
  
  arma::colvec traces = sum(trcRs,1);

  double p1 = pow(sqrt(PI)*R::gammafn(kappa+2)/R::gammafn(kappa+0.5),n);
  double p2 = sqrt(R::trigamma(kappa+0.5)-R::trigamma(kappa+2));
  arma::colvec p3 = pow(0.5+0.25*(traces-1),kappa);
  double p4 = 1;
  
  for(int i=0;i<n;i++){
    p4 *= p3[i]; 
  }

  return p1*p2*p4;
 
}

