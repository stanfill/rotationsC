#ifndef ESTIMATORS_H
#define ESTIMATORS_H

#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
int checkQ4(Rcpp::NumericMatrix Q);

// [[Rcpp::export]]
int checkSO3(arma::mat Rs);

// [[Rcpp::export]]
arma::mat expskewC(arma::mat M);

// [[Rcpp::export]]
arma::mat expskewCMulti(arma::mat M);

// [[Rcpp::export]]
arma::mat logSO3C(arma::mat R);

// [[Rcpp::export]]
arma::mat logSO3CMulti(arma::mat R);

// [[Rcpp::export]]
arma::mat projectSO3C(arma::mat M);

// [[Rcpp::export]]
arma::mat meanSO3C(arma::mat Rs);

// [[Rcpp::export]]
arma::rowvec meanQ4C(arma::mat Q);

//[[Rcpp::export]]
arma::mat medianSO3C(arma::mat Rs, int maxIterations, double maxEps);

//[[Rcpp::export]]
arma::mat HartmedianSO3C(arma::mat Rs, int maxIterations, double maxEps);

// [[Rcpp::export]]
arma::mat gmeanSO3C(arma::mat Rs, int maxIterations, double maxEps);

#endif /* ESTIMATORS_H */
