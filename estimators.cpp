#include <RcppArmadillo.h>   
#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]
int checkSO3(arma::mat Rs){
	
	/*This function will check that the rows in the matrix Rs are proper rotations*/
	int n = Rs.n_rows, p = Rs.n_cols, i, j;
	double deti, inv;
	arma::mat Ri(3,3), I(3,3);
	I.eye();
	
	if(n!=9 && p!=9){
		throw Rcpp::exception("The data are not each of length 9.");	
		return 1;
	}
	for(i=0;i<n;i++){
		
		for(j = 0; j<9; j++){
			Ri(j) = Rs(i,j);
		}
		deti = det(Ri); //Check that each determinant is one
		
		if(deti > 1.1 || deti < 0.9 ){
			throw Rcpp::exception("The data are not all determinant 1, so they are not rotations.");
			return 1;
		}
		
		inv = sum(sum(Ri*Ri.t()-I)); //Check that each inverse is the transpose
		
		if(inv > 0.001 || inv < -0.001 ){
			throw Rcpp::exception("Atleast on observation's transpose is not its invers, so they are not rotations.");
			return 1;
		}
	}
	return 0;
}


// [[Rcpp::export]]
arma::mat projectSO3C(arma::mat M){
	
	/*This function will project the an arbitrary 3-by-3 matrix M in M(3) into SO(3)
	It is expecting a 3-by-3 matrix*/
	
	arma::mat Msq = M.t()*M;
	arma::mat eigvec;
	arma::vec eigval;
  arma::eig_sym(eigval,eigvec,Msq); 
  arma::mat dMat(3,3);
  arma::mat u = fliplr(eigvec);
  dMat.zeros();
  
  int sign = 1;
  
  if(det(M)<0){
  	sign = -1;
  }
  
  dMat(0,0) = pow(eigval[2],-0.5);
  dMat(1,1) = pow(eigval[1],-0.5);
  dMat(2,2) = sign*pow(eigval[0],-0.5);

  return M * u * dMat * u.t();
	 
}

// [[Rcpp::export]]
arma::mat meanSO3C(arma::mat Rs){
	
	/*Compute the projected mean for a sample of n roations, Rs.  
	This function expects Rs to be a n-by-9 matrix where each row
	represents an observations in SO(3)*/
	
	int cSO3 = checkSO3(Rs);
	if(cSO3){
		throw Rcpp::exception("The data are not in SO(3).");
	}
	
	int i;
	arma::mat Rbarels = mean(Rs);
	arma::mat Rbar(3,3);
	
	for(i=0;i<9;i++){
			Rbar[i] = Rbarels[i];
	}
	
	return projectSO3C(Rbar);
}

/*// [[Rcpp::export]]
arma::mat medianSO3C(arma::mat Rs){
	
}*/

/*** R
library(rotations)
library(microbenchmark)

Rs<-ruars(200,rcayley)
tim<-microbenchmark(
meanSO3C(Rs),
mean(Rs))

plot(tim)
print(tim)

*/

