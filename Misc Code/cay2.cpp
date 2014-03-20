#include <Rcpp.h>
using namespace Rcpp;

/////////////////////////////////////////////////////////////
// Generate Cayley random deviates using C++
/////////////////////////////////////////////////////////////
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

/////////////////////////////////////////////////////////////
// Generate von Mises random deviates using C++
/////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////
// Generate matrix Fisher random deviates using C++
/////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double dfisherCpp(double r, double kappa) {
    
  double den;
  double I02k = R::bessel_i(2*kappa,0,1);
  double I12k = R::bessel_i(2*kappa,1,1);
  
  den = exp(2 * kappa * cos(r)); 
  den *= (1 - cos(r));
  
  den /= (2 * PI * (I02k - I12k));
  
  return den;
}

// [[Rcpp::export]]
double arsample_unifCpp(double M, double kappa) {
  RNGScope scope;
  //generate a random observation from target density f assuming g is uniform
  int found = 0; //FALSE
  NumericVector y(1);
  double x, evalF = 0.0;
  
  while (!found) {
    x = as<double>(runif(1, -PI, PI));
    y = runif(1, 0, M);
    
    evalF = dfisherCpp(x,kappa);
    
    if (y[0] < evalF) 
      found = 1;
  }
  return x;

}

// [[Rcpp::export]]
NumericVector rarCpp(int n, double kappa, double M) {
  
  NumericVector res(n);
  for (int i=0;i<n;i++){
    res[i] = arsample_unifCpp(M,kappa);
  } 
  return res;
}

// [[Rcpp::export]]

NumericVector rfisherCpp(int n, double kappa) {
  double step = 0.01, prog = -PI;
  double M = 0.0, Mi=0.0;
  NumericVector res(n);
  
  while(prog < 0){
    Mi = dfisherCpp(prog,kappa);
    if(M<Mi){
      M = Mi;
    }
    prog += step;
  }
  
  res = rarCpp(n, kappa ,M);
  return res;  
}