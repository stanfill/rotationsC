
#include <RcppArmadillo.h>   
#include "../inst/include/rotations.h"  

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

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


arma::rowvec RdistCArma(arma::mat Q1, arma::rowvec Q2){
  /*Compute the geodesic distance between quaternions Q1 and Q2*/
  /* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single (1x4) quaternion*/
	
	int n = Q1.n_rows, i=0; 
	double cp;
	arma::rowvec rs(n);
	
	for(i=0;i<n;i++){
		
		cp = sum(Q1.row(i)*Q2.t());
		rs(i) = acos(2*cp*cp-1);
		
	}
	
	return rs;
}


// [[Rcpp::export]]
arma::rowvec HnCppIntrinsic(arma::mat Qs){
  
  //Compute the intrinsic Hn tests statistics
  
  int n = Qs.n_rows, i=0, j=0;

  //Get T matrix of whole sample to make it easier later on
  arma::mat T = Qs.t()*Qs;
  arma::rowvec Qhat = rotations::meanQ4C(Qs);
  arma::rowvec dists(n);

  //Sum of squared geometric distances between proj. mean and each obs
  double SSE = 0.0, SSEJ=0.0;
  dists = square(RdistCArma(Qs,Qhat));
  SSE = sum(dists);

  arma::rowvec Hn(n);
  arma::rowvec Qhatj;
  arma::mat QsJ(n,4);
  
  //Variables for reduced sample mean
  arma::rowvec Qj, distsJ(n-1);
  arma::mat Tj(4,4), eigvecJ(4,4);

  for(i = 0;i<n; i++){
    
    QsJ.resize(n,4);
    
    for(j = 0; j<n; j++){
      if(j!=i){
        QsJ.row(j) = Qs.row(j);
      }
    }
    
    QsJ.shed_row(i);
    
    //Compute projected mean when jth row is cut out
    Qhatj = rotations::meanQ4C(QsJ);
    distsJ = square(RdistCArma(QsJ,Qhatj));
    SSEJ = sum(distsJ);
    
    //Rcpp::Rcout << "SSEJ: " << SSEJ << std::endl;
    
    Hn(i)=(n-2)*(SSE-SSEJ)/(SSEJ);
    
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

