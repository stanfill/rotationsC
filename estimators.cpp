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
arma::mat expskewC(arma::mat M){
  /*This function takes a 3-by-3 skew symmetric matrix (in so(3)) and
  returs the exponential, a 3-by-3 roations (in SO(3))*/
  
  double MMt = sum(sum(M-M.t()));
  
  if(abs(MMt)>0.01){
    throw Rcpp::exception("The expskewC function is expecting a 3-by-3 skew symmetric matrix");
  }
  
  arma::mat expM(3,3);
  expM.eye();
  
  double a = pow(0.5*trace(M.t()*M),0.5);
  
  if(abs(a)<0.001){
    
    return expM;
  
  }
   
  expM = expM + (sin(a)/a) * M + (1-cos(a))*pow(a,-2)*M*M;
  return expM;
  
}


// [[Rcpp::export]]
arma::rowvec rdistSO3C(arma::mat Rs, arma::mat R2){
  /* This function takes the matrix of matrices Rs and returns
  each of the observatiosn distances from R2, which is assumed to be a
  3-by-3 rotation matrix.*/
  
  int cSO3 = checkSO3(Rs);
	if(cSO3){
		throw Rcpp::exception("The data are not in SO(3).");
	}
  
  int n = Rs.n_rows, i,j;
  
  arma::rowvec theta(n);
  theta.zeros();
  arma::mat33 Rsi;
  
  for(i=0; i<n ; i++){
    
    for(j = 0; j<9 ;j++){
      Rsi(j)=Rs(i,j);
    }
    
    Rsi = Rsi * R2.t();
    theta(i) = acos(0.5*trace(Rsi)-0.5);
    
  }
  return theta;
}

// [[Rcpp::export]]
arma::mat logSO3C(arma::mat R){
  
  /*int cSO3 = checkSO3(R);
  if(cSO3){
		throw Rcpp::exception("The data are not in SO(3).");
	}*/
  
  arma::mat I(3,3), logR(3,3);
  I.eye();
  
  double theta = acos(0.5*trace(R)-0.5);
  
  if(theta < 0.0001){
    logR.zeros();
    return logR;
  }
  
  logR = (R-R.t())*theta/(2*sin(theta));
  
  return logR;
  
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

//[[Rcpp::export]]
arma::mat medianSO3C(arma::mat Rs){
  
  /*Estimate the central direction with the projected median according
  to our algorithm in point estimation paper*/
  
  int cSO3 = checkSO3(Rs);
	if(cSO3){
		throw Rcpp::exception("The data are not in SO(3).");
	}
  int n = Rs.n_rows, i,j,iterations=0;
  arma::mat S = meanSO3C(Rs), RsCopy = Rs, Snew;
  arma::mat33 delta;
  arma::rowvec vnInv(n), deltaV(9);
  double denom,epsilon = 1.0;
  arma::rowvec Svec(9);
  
  while(epsilon > 0.0005 && iterations < 1500){
  
    for(j=0;j<9;j++){
      Svec(j) = S(j); 
    }
  
    denom = 0;
    for(i=0;i<n;i++){
      
      vnInv(i) = pow(norm(Rs.row(i)-Svec,2),-1);
      RsCopy.row(i) = Rs.row(i)*vnInv(i);
      denom += vnInv(i);
      
    }
  
    deltaV = sum(RsCopy)/denom;

    for(j=0;j<9;j++){
      delta(j) = deltaV(j);
    }

    Snew = projectSO3C(delta);
    
    iterations += 1;
    epsilon = norm(Snew-S,2);
    S = Snew;
  }

  //printf(" %i ",iterations);
  
  return S;
}

/*** R
library(rotations)
library(microbenchmark)
rs<-rcayley(100)
Rs<-genR(rs)

median(Rs)
medianSO3C(Rs)

tim<-microbenchmark(
medianSO3C(Rs),
median(Rs))

plot(tim)
print(tim)

*/

