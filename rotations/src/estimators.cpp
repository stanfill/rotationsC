#include <RcppArmadillo.h>   
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
int checkQ4(NumericMatrix Q){
	/*This function will check that the rows in the matrix Q are unit quaternions*/
	int n = Q.nrow(), p = Q.ncol(), i;
	double len;
	
	if(n!=4 && p!=4){
		throw Rcpp::exception("The data are not of length 4 each.");	
		return 1;
	}
	
	for(i=0;i<n;i++){
		
		len = sum(Q(i,_)*Q(i,_));
		if(len > 1.1 || len < 0.9){
			
			throw Rcpp::exception("The data are not all unit length so are not quaternions.");
			return 1;
			
		}		
	}
		
	return 0;
}

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
  returns the exponential, a 3-by-3 rotations (in SO(3))*/
  
  double MMt = sum(sum(M-M.t()));
  
  if(fabs(MMt)>0.01){
    throw Rcpp::exception("The exp.skew function is expecting a 3-by-3 skew symmetric matrix");
  }
  
  arma::mat expM(3,3);
  expM.eye();
  
  double a = pow(0.5*trace(M.t()*M),0.5);
  
  if(a < 0.000001 && a > -0.000001){
    return expM;
  }
   
  expM = expM + (sin(a)/a) * M + (1-cos(a))*pow(a,-2)*M*M;
  
  return expM;
  
}

// [[Rcpp::export]]
arma::mat expskewCMulti(arma::mat M){
  /*This function takes a sample of 3-by-3 skew symmetric matrices (in so(3)) and
  returs the exponential, a sample of 3-by-3 rotations (in SO(3))*/
  
  int n = M.n_rows, i, j;
  arma::mat eMi(3,3), Mi(3,3);
  arma::mat expM;
  expM = M; /*take dimensionality of input matrix*/
  expM.zeros(); /*make it all zeros*/
  
  for(i = 0 ; i<n ; i++ ){
    
    for(j = 0; j<9 ; j++ ){
      Mi(j)=M(i,j);
    }
    
    eMi = expskewC(Mi);
    
    for(j=0; j<9 ; j++ ){
      expM(i,j) = eMi(j);
    }
    
  }
  return expM;
  
}

// [[Rcpp::export]]
arma::mat logSO3C(arma::mat R){
  
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
arma::mat logSO3CMulti(arma::mat R){
  //This is a version of logSO3C that allows for multiple rows in Rs
  
  int n = R.n_rows, i, j;
  arma::mat I(3,3), logR(n,9), Ri(3,3), logRi(3,3);
  logR.zeros();
  I.eye();
  Ri.zeros();
  logRi.zeros();
  double theta;
  
  for(i = 0; i<n ; i++ ){
    
    for(j = 0; j<9; j++){
			Ri(j) = R(i,j);
		}
    
    theta = acos(0.5*trace(Ri)-0.5);
    
    //If theta<0.0001 leave that row as zeros
    if(theta > 0.00001){ 
  
      logRi = (Ri-Ri.t())*theta/(2*sin(theta));
    
      for(j = 0; j<9; j++){
  	  	logR(i,j) = logRi(j);
		  }
      
    }
    
  }
  
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
	
	int i;
	arma::mat Rbarels = mean(Rs);
	arma::mat Rbar(3,3);
	
	for(i=0;i<9;i++){
			Rbar[i] = Rbarels[i];
	}
	
	return projectSO3C(Rbar);
}



// [[Rcpp::export]]   
arma::rowvec meanQ4C(arma::mat Q) { 
	//Compute the projected mean of the sample Q
	
	NumericMatrix Qss = as<NumericMatrix>(wrap(Q));
	int cq4 = checkQ4(Qss);
	if(cq4){
		throw Rcpp::exception("The data are not in Q4.");
	}
	
	arma::mat Qsq=Q.t()*Q;
	arma::mat eigvec;
	arma::vec eigval;
  arma::eig_sym(eigval,eigvec,Qsq);   
  arma::vec qhat=eigvec.col(3);
  
  if(qhat[0]<0){
  	qhat = -qhat;
  }
  
  return qhat.t(); //Want to return it in a row vector so transpose it
}



//[[Rcpp::export]]
arma::mat medianSO3C(arma::mat Rs, int maxIterations, double maxEps){
  
  /*Estimate the central direction with the projected median according
  to our algorithm in point estimation paper*/
  
  int n = Rs.n_rows, i,j,iterations=0;
  arma::mat S = meanSO3C(Rs), RsCopy = Rs, Snew;
  arma::mat33 delta;
  arma::rowvec vnInv(n), deltaV(9), Svec(9);
  double denom,epsilon = 1.0;
  
  while(epsilon > maxEps && iterations < maxIterations){
  
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

  return S;
}

 
//[[Rcpp::export]]
arma::mat HartmedianSO3C(arma::mat Rs, int maxIterations, double maxEps){
  
  /*Estimate the central direction with the projected median according
  to our algorithm in point estimation paper*/
  
  int n = Rs.n_rows, i,j,iterations=0;
  arma::mat S = meanSO3C(Rs), RsCopy = Rs, Snew;
  arma::mat33 delta, Rsi, vi;
  arma::rowvec vnInv(n);
  double denom,epsilon = 1.0;
  double vin;
  
  while(epsilon > maxEps && iterations < maxIterations){
    
    delta.zeros();
    denom = 0;
    
    //printf("vi: ");
    
    for(i=0;i<n;i++){
      
      for(j=0;j<9;j++){
        Rsi(j)=Rs(i,j);
      }
      
      vi = logSO3C(Rsi*S.t());
      vin = std::max(norm(vi,2),1e-5);
      
      //printf(" %lf ",vin);
      
      vnInv(i) = pow(vin,-1);
      delta+=vi*vnInv(i);
      denom += vnInv(i);
      
    }
  	//printf("denom: %lf\n",denom);
    delta = delta/denom;

    Snew = expskewC(delta)*S;
    
    iterations += 1;
    epsilon = norm(Snew-S,2);
    S = Snew;
  }
  return S;
}


// [[Rcpp::export]]
arma::mat gmeanSO3C(arma::mat Rs, int maxIterations, double maxEps){
  int n = Rs.n_rows, i=0, j=0, iterations=0;
  arma::mat33 Rsi, r;
  arma::mat S = meanSO3C(Rs);
  double eps=1.0;
  Rsi.zeros();

  while(eps > maxEps && iterations < maxIterations){
  	
 	 	r.zeros();
  	
    for(i = 0; i<n ; i++){
    
      for(j=0;j<9;j++){
        Rsi(j) = Rs(i,j); 
      }
    
      r = r + logSO3C(S.t()*Rsi);
    
    }

    r = r/n;

 	  S = S*expskewC(r);
    eps = norm(r,"fro");
    iterations++;
  }
  return S;
}

