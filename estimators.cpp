#include <RcppArmadillo.h>   
#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 
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


//' Rotational Distance
//'
//' Calculate the Euclidean or Riemannian distance between two rotations
//'
//' This function will calculate the intrinsic (Riemannian) or projected (Euclidean) distance between two rotations.  If only one rotation is specified
//' the other will be set to the identity and the distance between the two is returned.  For rotations \eqn{R_1}{R1} and \eqn{R_2}{R2}
//' both in \eqn{SO(3)}, the Euclidean distance between them is \deqn{||R_1-R_2||_F}{||R1-R2||} where \eqn{||\cdot||_F}{|| ||} is the Frobenius norm.
//' The intrinsic distance is defined as \deqn{||Log(R_1^\top R_2)||_F}{||Log(R1'R2)||} where \eqn{Log} is the matrix logarithm, and it corresponds
//' to the misorientation angle of \eqn{R_1^\top R_2}{R1'R2}.
//'
//' @param Rs (or Qs) a matrix of rotations in either matrix or quaternion form
//' @param R2 (or Q2) the second rotation in the same parameterization as R1
//' @return a vector of rotational distances from each row of Rs (or Qs) to R2
//' @export
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


//' Natural Logarithm in SO(3)
//'
//' For details see \cite{moakher02}
//'
//' @param R numeric matrix in \eqn{SO(n)}
//' @return mlog numeric matrix \eqn{\log(R)}{log(R)}
//' @cite moakher02
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


//' Projection Procedure
//'
//' Project an arbitrary \eqn{3\times 3}{3-by-3} matrix into SO(3)
//'
//' This function uses the process given in \cite{moakher02} to project an arbitrary \eqn{3\times 3}{3-by-3} matrix into \eqn{SO(3)}.
//' 
//' @param M \eqn{3\times 3}{3-by-3} matrix to project
//' @return projection of \eqn{\bm M}{M} into \eqn{SO(3)}
//' @seealso \code{\link{mean.SO3}}, \code{\link{median.SO3}}
//' @export
//' @examples
//' M<-matrix(rnorm(9),3,3)
//' project.SO3(M)
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

//' Rotation Mean Rotation
//'
//' Compute the intrinsic or projected mean of a sample of rotations
//'
//' This function takes a sample of \eqn{3\times 3}{3-by-3} rotations (in the form of a \eqn{n\times 9}{n-by-9} matrix where \eqn{n>1} is the sample size) and returns the projected arithmetic mean denoted \eqn{\widehat{\bm S}_P}{S_P} or
//' intrinsic mean \eqn{\widehat{\bm S}_G}{S_G} according to the \code{type} option.
//' For a sample of \eqn{n} random rotations \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})}{argmin d^2(bar(R),S)} where \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i}{bar(R)=\sum Ri/n} and the distance metric \eqn{d_D}{d}
//' is the Riemannian or Euclidean.  For more on the projected mean see \cite{moakher02} and for the intrinsic mean see \cite{manton04}.
//'
//' @param Rs A \eqn{n\times 9}{n-by-9} matrix where each row corresponds to a random rotation in matrix form
//' @return Estimate of the projected or intrinsic mean of the sample
//' @seealso \code{\link{medianSO3C}}
//' @cite moakher02, manton04
//' @examples
//' Rs<-ruars(20,rvmises,kappa=0.01)
//' mean(Rs)
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


//' Quaternion Mean Rotation
//'
//' Compute the intrinsic or projected mean of a sample of rotations
//'
//' This function takes a sample of \eqn{4\times 1}{4-by-1} rotations (in the form of a \eqn{n\times 4}{n-by-4} matrix where \eqn{n>1} is the sample size) and returns the projected arithmetic mean denoted \eqn{\widehat{\bm S}_P}{S_P} or
//' intrinsic mean \eqn{\widehat{\bm S}_G}{S_G} according to the \code{type} option.
//' For a sample of \eqn{n} random rotations \eqn{\bm{R}_i\in SO(3), i=1,2,\dots,n}{Ri in SO(3), i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})}{argmin d^2(bar(R),S)} where \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i}{bar(R)=\sum Ri/n} and the distance metric \eqn{d_D}{d}
//' is the Riemannian or Euclidean.  For more on the projected mean see \cite{moakher02} and for the intrinsic mean see \cite{manton04}.
//'
//' @param Rs A \eqn{n\times 4}{n-by-4} matrix where each row corresponds to a random rotation in matrix form
//' @return Estimate of the projected or intrinsic mean of the sample
//' @seealso \code{\link{medianSO3C}}
//' @cite moakher02, manton04
//' @examples
//' Rs<-ruars(20,rvmises,kappa=0.01)
//' mean(Rs)
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


//' Median Rotation
//' 
//' Compute the projected or intrinsic median of a sample of rotations
//'
//' The median-type estimators are defined as \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}{argmin\sum d(Ri,S)}.  If the choice of distance metrid, \eqn{d_D}{d}, is Riemannian then the estimator is called the intrinsic, and if the distance metric in Euclidean then it projected.
//' The algorithm used in the intrinsic case is discussed in \cite{hartley11} and the projected case was written by the authors.
//'
//' @param Rs A \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix form (\eqn{p=9}) or quaternion form (\eqn{p=4})
//' @return an estimate of the projected or intrinsic mean
//' @seealso \code{\link{meanSO3C}}
//' @cite hartley11
//' @export
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
  arma::rowvec vnInv(n), deltaV(9), Svec(9);
  double denom,epsilon = 1.0;
  
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

  return S;
}

//[[Rcpp::export]]
arma::mat HartmedianSO3C(arma::mat Rs){
  
  /*Estimate the central direction with the projected median according
  to our algorithm in point estimation paper*/
  
  int cSO3 = checkSO3(Rs);
  if(cSO3){
		throw Rcpp::exception("The data are not in SO(3).");
	}
  int n = Rs.n_rows, i,j,iterations=0;
  arma::mat S = meanSO3C(Rs), RsCopy = Rs, Snew;
  arma::mat33 delta, Rsi, vi;
  arma::rowvec vnInv(n);
  double denom,epsilon = 1.0;
  
  while(epsilon > 0.0001 && iterations < 1500){
    
    delta.zeros();
    denom = 0;
    for(i=0;i<n;i++){
      
      for(j=0;j<9;j++){
        Rsi(j)=Rs(i,j);
      }
      
      vi = logSO3C(Rsi*S.t());
      
      vnInv(i) = pow(norm(vi,2),-1);
      delta+=vi*vnInv(i);
      denom += vnInv(i);
      
    }
  
    delta = delta/denom;

    Snew = expskewC(delta)*S;
    
    iterations += 1;
    epsilon = norm(Snew-S,2);
    S = Snew;
  }
  return S;
}

// [[Rcpp::export]]
arma::mat gmeanSO3C(arma::mat Rs){
  int n = Rs.n_rows, i=0, j=0, iterations=0;
  arma::mat33 Rsi, r;
  arma::mat S = meanSO3C(Rs);
  double eps=1.0;
  
  while(eps > 0.0001 && iterations < 2000){
    
    r.zeros();
    S = S*expskewC(r);
  
    for(i = 0; i<n ; i++){
    
      for(j=0;j<9;j++){
        Rsi(j) = Rs(i,j); 
      }
    
      r += logSO3C(Rsi.t()*S);
    
    }
  
    r = r/n;

    eps = norm(r,2);
    iterations++;
  }
  return S;
}

