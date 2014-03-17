#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericVector rcayleyCpp(int n, double kappa){
  RNGScope scope;
  NumericVector bet(n), alp(n), theta(n);
    
  bet = rbeta(n,kappa+0.5,1.5);
  alp = rbinom(n,1,0.5);
  
  for(int i=0; i<n; i++){
    theta[i] = acos(2*bet[i]-1)*(1-2*alp[i]);
  }
  
  return theta;
}

// [[Rcpp::export]]
int sign(double x){
  if(x < 0){
    return(-1);
  }else{
    return(1);
  }
}

// [[Rcpp::export]]
NumericVector rvmisesCPP(int n, double kappa){
  RNGScope scope;
  NumericVector u(3), theta(n, 10.0);
    
  u = runif(3, 0, 1);
  double a = 1 + sqrt(1 + 4 * pow(kappa,2));
  double b = (a - sqrt(2 * a))/(2 * kappa);
  double r = (1 + pow(b,2))/(2 * b);
  double z = 0, f = 0, c = 0;
  
  for (int i = 0; i<n; i++) {
    
    while (theta[i] > 4) {
      // Step 1
      u = runif(3, 0, 1);
      z = cos(PI * u[0]);
      f = (1 + r * z)/(r + z);
      c = kappa * (r - f);
      
      // Step 2
      u = runif(3, 0, 1);
      if ((c * (2 - c) - u[1]) > 0) {
        
        theta[i] = (sign(u[2] - 0.5)) * acos(f);
        
      } else {
        
        if ((log(c/u[1]) + 1 - c) < 0) {
          u = runif(3, 0, 1);
        } else {
          u = runif(3, 0, 1);
          theta[i] = (sign(u[2] - 0.5)) * acos(f);
        }
      }
    }
  }
  
  return theta;
}