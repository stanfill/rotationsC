#include <RcppArmadillo.h>   
#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
arma::mat eskewC(arma::rowvec U) {
  
  double ulen=norm(U,2);
  
  U = U/ulen;
  
  double u = U(0);
  double v = U(1);
  double w = U(2);
  
  arma::mat33 res;
  res.zeros();
  res(0,1) = -w;
  res(1,0) = w;
  res(0,2) = v;
  res(2,0) = -v;
  res(1,2) = -u;
  res(2,1) = u;
  
  return res;
}

// [[Rcpp::export]]
arma::mat SO3defaultC(arma::mat U, arma::vec theta) {
  
  //U is an n-by-3 matrix, each row is a misorentation axis
  //theta is a vector of length n, each item is a misorientation angle
  
  int n=U.n_rows, i=0;	
	arma::mat Ri(3,3), Rs(n,9), I(3,3), SS(3,3);
  Ri.zeros(); Rs.zeros(); I.eye();
  arma::rowvec Rir;

  for(i=0;i<n;i++){
 		Ri = U.row(i).t() * U.row(i);
    SS = eskewC(U.row(i));
  	Ri = Ri + (I - Ri) * cos(theta[i]) +  SS * sin(theta[i]);
    Rs.row(i) = as<arma::rowvec>(wrap(Ri));
  }
 		
  return Rs;
}


// [[Rcpp::export]]
arma::mat Q4defaultC(arma::mat U, arma::vec theta){
	
	int n1 = U.n_rows, n = theta.size();
	arma::mat q(n,4);
	q.zeros();
	
	if(n1 != n){
		//printf("Error, u and theta different length");
		return q;
	}
	
	arma::vec stheta = sin(theta/2);
	
	q.col(0) = cos(theta/2);
	q.col(1) = U.col(0) % stheta;
	q.col(2) = U.col(1) % stheta;
	q.col(3) = U.col(2) % stheta;
	
	return q;
	
}


// [[Rcpp::export]]
arma::mat pMatC(arma::mat p){

	arma::mat Pmat(4,4);
	Pmat.zeros();
	arma::mat revI(4,4);
	revI.zeros();
	
	//int n1=p.n_rows, n2=p.n_cols;
	
	//if(n1!=4 && n2!=4){}
	
	p.reshape(4,1);
	Pmat.col(0)=p;

	revI(0,1) = -1;revI(1,0) = 1;	revI(2,3) = 1; revI(3,2) = -1;
	Pmat.col(1) = revI*p;

	revI.zeros();
	revI(0,2) = -1;revI(1,3) = -1;	revI(2,0) = 1; revI(3,1) = 1;
	Pmat.col(2) = revI*p;

	revI.zeros();
	revI(0,3) = -1;revI(1,2) = 1;	revI(2,1) = -1; revI(3,0) = 1;
	Pmat.col(3)=revI*p;
	

	return Pmat;
}

// [[Rcpp::export]]
arma::mat genrC(arma::vec r, arma::mat S , int SO3, arma::mat u) {
	RNGScope scope;
	// r is a vector of angles
	// S is the central direction
	// SO3 is an integer, 1 means SO3, anything else gives
  int n=r.size(), i=0,n1 = u.n_rows, n2 = u.n_cols;
  
  //GetRNGstate();PutRNGstate();
  
  //NumericVector theta = runif(n,-1,1);
  //theta = acos(theta);
    
  //NumericVector phi = runif(n, -M_PI, M_PI);
  
  //int n1 = phi.size(), n2 = theta.size();
    
  if(n1 != n || n2!=3){
    //printf("u is wrong size");
    arma::mat q(n,4);
    q.zeros();
    return q;
  }
  
  /*arma::mat u(n,3);
  
  for(i=0;i<n;i++){
    u(i,0)=sin(theta[i]) * cos(phi[i]);
    u(i,1)=sin(theta[i]) * sin(phi[i]);
    u(i,2)=cos(theta[i]);
  }  */
  
  if(SO3==1){
    
    int j;
    arma::mat Rs(n,9);
    arma::mat33 Rsi;

    Rs = SO3defaultC(u, r);
    
    
    for(i=0;i<n;i++){
      
      for(j=0;j<9;j++){
        Rsi(j) = Rs(i,j);
      }
      
      Rsi = S * Rsi;
      Rs.row(i) = as<arma::rowvec>(wrap(Rsi));
      
    }
      
    return Rs;
    
  }else{
  	
  	arma::mat q;
  	q.zeros(n,4);
  	//arma::mat Smat;
    
    /*int ssize = S.n_rows;
    int ssize2 = S.n_cols;
    
    if(ssize!=4 && ssize2!=4){
      printf("S isn't big enough");
      q.zeros(n,4);
      return q;
    }*/
    
  	//Smat = pMatC(S);

  	q = Q4defaultC(u,r);
 		
 		//q = q*Smat.t();
  	
    return q;
    
  }

}
