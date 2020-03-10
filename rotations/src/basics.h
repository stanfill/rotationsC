#ifndef BASICS_H
#define BASICS_H

#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::mat eskewC(arma::rowvec U);

// [[Rcpp::export]]
arma::mat SO3defaultC(arma::mat U, arma::vec theta);

// [[Rcpp::export]]
arma::mat Q4defaultC(arma::mat U, arma::vec theta);

// [[Rcpp::export]]
arma::mat pMatC(arma::mat p);

// [[Rcpp::export]]
arma::mat genrC(arma::vec r, arma::mat S , int SO3, arma::mat u);

#endif /* BASICS_H */
