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

//' generate a rotation matrix with axis and angles of rotation u and theta, respectively
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


//' A function to create a rotation in quaternion form with axis U and angle theta
// [[Rcpp::export]]
arma::mat Q4defaultC(arma::mat U, arma::vec theta){
	
	int n = U.n_rows, i=0;
  int n2 = theta.n_elem, n3 = U.n_cols;
    
	arma::mat q(n,4);
	q.zeros();
  
  if(n!=n2 || n3<3){
    throw Rcpp::exception("Error in Q4defaultC, u and theta not same length.");
    return q;
  }
  
	arma::rowvec stheta;
	
	for(i=0;i<n;i++){
		
		//U.row(i) = U.row(i)/norm(U.row(i),2);
		
		q(i,0) = cos(theta(i)/2);
		stheta = sin(theta(i)/2)*U.row(i);
      
		q(i,1) = stheta(0);
		q(i,2) = stheta(1);
		q(i,3) = stheta(2);
		
	}
	
  return q;
}

// [[Rcpp::export]]
arma::mat pMatC(arma::mat p){

	arma::mat Pmat(4,4);
	Pmat.zeros();
	arma::mat revI(4,4);
	
	p.reshape(4,1);
	
	Pmat.col(0)=p;

	revI.zeros();
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

//' a function to generate UARS rotations with angles of rotations r and central direction S
// [[Rcpp::export]]
arma::mat genrC(arma::vec r, arma::mat S , int SO3) {
	// r is a vector of angles
	// S is the central direction
	// SO3 is an integer, 1 means SO3, anything else gives
  int n=r.size(), i=0;
  
  NumericVector theta = runif(n,-1,1);
  theta = acos(theta);
    
  NumericVector phi = runif(n, -M_PI, M_PI);
  
  int n1 = phi.size(), n2 = theta.size();
  
  if(n1 != n && n2 != n ){
    printf("runif screwed me");
    arma::mat33 q;
    q.zeros();
    return q;
  }
  
  arma::mat u(n,3);
  
  for(i=0;i<n;i++){
    u(i,0)=sin(theta[i]) * cos(phi[i]);
    u(i,1)=sin(theta[i]) * sin(phi[i]);
    u(i,2)=cos(theta[i]);
  }  
  
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
  	arma::mat Smat;
    
    int ssize = S.n_rows;
    int ssize2 = S.n_cols;
    
    if(ssize!=4 && ssize2!=4){
      printf("S isn't big enough");
      q.zeros(3,3);
      return q;
    }
    
  	Smat = pMatC(S);

  	q = Q4defaultC(u,r);
 		
 		q = q*Smat.t();
  	
    return q;
    
  }

}