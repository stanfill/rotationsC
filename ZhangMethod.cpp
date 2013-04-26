#include <RcppArmadillo.h>   
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

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
arma::rowvec meanQ4C(arma::mat Q) { 
	/*Compute the projected mean of the sample Q.*/
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
  
  return qhat.t(); /*Want to return it in a row vector so transpose it*/
} 

// [[Rcpp::export]]
NumericVector RdistC(NumericMatrix Q1, NumericVector Q2){
	/*Compute the geodesic distance between quaternions Q1 and Q2*/
	/* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single quaternion*/
	
	int cq4 = checkQ4(Q1);
	
	if(cq4){
		throw Rcpp::exception("The data are not in Q4.");
	}
	
	int n = Q1.nrow(), i; 
	double cp;
	NumericVector rs(n);
	
	for(i=0;i<n;i++){
		
		cp = sum(Q1(i,_)*Q2);
		rs[i] = acos(2*cp*cp-1);
		
	}
	
	return rs;
}

// [[Rcpp::export]]
NumericVector cdfunsC(NumericMatrix Qs, NumericVector Qhat){
	
	//Compute the values c and d to form the pivotal test statistic
	
	int n = Qs.nrow(), i;
	double crs;
	
	NumericVector cds(2);
	cds[0]=0;
	cds[1]=0;
	
	NumericVector rs(n);
	
	rs = RdistC(Qs,Qhat);
	
	for(i=0; i<n; i++){
		
		crs = cos(rs[i]);
		
		cds[0] += pow(crs,2);				//c=2E[1-cos(r)^2]/3
		cds[1] += crs;							//d=E[1+2cos(r)]/3
	}
	
	cds[0] = 2*(1-cds[0]/n)/3;
	cds[1] = (1+2*cds[1]/n)/3;
	
	return cds;
}

// [[Rcpp::export]]
NumericVector bootQhat(NumericMatrix Q, int m){
	
	int cq4 = checkQ4(Q);
	
	if(cq4){
		throw Rcpp::exception("The data are not in Q4.");
	}
	
	int n=Q.nrow(), i=0, j=0;
	NumericVector samp, cdstar;
	
	NumericVector testStat(m), sqrth;
	
	arma::mat Qstar(n,4);
	NumericVector QhatStar;
  NumericMatrix QhatStarMat(1,4);
	
	arma::mat QSamp = as<arma::mat>(Q); //Convert the sample into armadillo mode
	
	NumericMatrix QstarRcpp;
	
	NumericVector Qhat = as<NumericVector>(wrap(meanQ4C(QSamp)));
	
	for(j=0;j<m;j++){
		
		samp = floor(runif(n,0,n));			//Bootstrap sample of size n, with replacement
    
    
		for(i=0;i<n;i++){
			Qstar.row(i) = QSamp.row(samp[i]);		//Copying a matrix row by row produces a bunch of junk messages
		}																				//so I do it with arma instead of standard Rcpp
	
		QhatStar = as<NumericVector>(wrap(meanQ4C(Qstar))); //Both of these functinos return arma variables so
		QstarRcpp = as<NumericMatrix>(wrap(Qstar));					//They need to be converted to Rcpp type
		
		cdstar = cdfunsC(QstarRcpp,QhatStar);
		
		QhatStarMat = as<NumericMatrix>(QhatStar); /*QhatStar needs to be a matrix to be used in RdistC*/
		sqrth = RdistC(QhatStarMat,Qhat);
		
		testStat[j] = 2*n*pow(cdstar[1],2)*pow(sqrth[0],2)/cdstar[0];
		
	}
	
	return testStat;
	
}

/*** R

#library(microbenchmark)
library(rotations)
#source("U:/Thesis/Intervals/Code/IntervalFuns.R")
n<-10
rs<-rcayley(n)
#Rs<-genR(rs,space='SO3')
Qs<-genR(rs,space='Q4')
meanQ4C(Qs)
bootQhat(Qs,300)


*/

