#ifndef ZHANGMETHOD_H
#define ZHANGMETHOD_H

#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
Rcpp::NumericVector RdistC(Rcpp::NumericMatrix Q1, Rcpp::NumericVector Q2);

// [[Rcpp::export]]
arma::rowvec rdistSO3C(arma::mat Rs, arma::mat R2);

// [[Rcpp::export]]
Rcpp::NumericVector EdistC(Rcpp::NumericMatrix Q1, Rcpp::NumericVector Q2);

// [[Rcpp::export]]
double oneRdistC(Rcpp::NumericMatrix Q1, Rcpp::NumericVector Q2);

// [[Rcpp::export]]
Rcpp::NumericVector cdfunsC(Rcpp::NumericMatrix Qs, Rcpp::NumericVector Qhat);

// [[Rcpp::export]]
Rcpp::NumericVector cdfunsCMedian(Rcpp::NumericMatrix Qs, Rcpp::NumericVector Qhat);

// [[Rcpp::export]]
Rcpp::NumericVector zhangQ4(Rcpp::NumericMatrix Q, int m);

// [[Rcpp::export]]
Rcpp::NumericVector cdfunsCSO3(arma::mat Rs, arma::mat Rhat);

// [[Rcpp::export]]
Rcpp::NumericVector zhangMedianC(arma::mat Rs, int m);

#endif /* ZHANGMETHOD_H */
