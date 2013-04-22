#include <RcppArmadillo.h>   
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]   
arma::vec meanQ4C(arma::mat Q) { 
	/*Compute the projected mean of the sample Q*/
	arma::mat Qsq=Q.t()*Q;
	arma::mat eigvec;
	arma::vec eigval;
  arma::eig_sym(eigval,eigvec,Qsq);   
  arma::vec qhat=eigvec.col(3);
  
  if(qhat[0]<0){
  	qhat = -qhat;
  }
  
  return qhat;
} 

// [[Rcpp::export]]
double RdistC(NumericVector Q1, NumericVector Q2){
	/*Compute the geodesic distance between quaternions Q1 and Q2*/
	double cp = sum(Q1*Q2);
	return acos(2*cp*cp-1);
	
}

// [[Rcpp::export]]
NumericVector cdfunsC(NumericMatrix Qs, NumericVector Qhat){
	int n = Qs.nrow();
	double crs;
	NumericVector cds(2);
	NumericVector rs(n);
	
	for(int i=0; i<n; i++){
		rs[i] = 2*acos(Qs(i,0));
		
		crs = cos(rs[i]);
		
		cds[0] += 1-pow(crs,2);
		cds[1] += 1+2*crs;
	}
	
	cds[0] = (2*cds[0])/(3*n);
	cds[1] = (cds[1]/(3*n));
	
	return cds;
}

// [[Rcpp::export]]
NumericVector bootQhat(NumericMatrix Q){
	
	int n=Q.nrow();
	NumericVector samp = runif(n);
	samp = floor(samp * n);
	
	arma::mat Qstar(n,4);
	arma::mat QSamp = as<arma::mat>(Q);
	
	for(int i=0;i<n;i++){
		Qstar.row(i) = QSamp.row(samp[i]);		//Copying a matrix row by row produces a bunch of junk messages
	}																				//so I do it with arma instead of standard Rcpp
	
	arma::vec QhatStar = meanQ4C(Qstar);
	
	return as<Rcpp::NumericVector>(wrap(QhatStar));
}

/*** R

library(microbenchmark)
library(rotations)
rs<-rcayley(20)
Qs<-genR(rs,space='Q4')
rs
angle(Qs)
cdfunsC(Qs,id.Q4)

#Qs<-ruars(10,rcayley,space='Q4')

#Qhat<-as.Q4(meanQ4C(Qs))
#Qhat
#mean(Qs)

#dist(Qhat,method='intrinsic')
#RdistC(Qhat,id.Q4)

#tim<-microbenchmark(          
#dist(Qhat,method='intrinsic'),
#RdistC(Qhat,id.Q4)
#)

#print(tim)
#plot(tim)

*/

