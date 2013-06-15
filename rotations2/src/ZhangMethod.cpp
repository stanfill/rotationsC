#include <RcppArmadillo.h>   
#include <Rcpp.h>
#include "../inst/include/rotations2.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 
// [[Rcpp::interfaces(r, cpp)]]

// [[Rcpp::export]]
NumericVector RdistC(NumericMatrix Q1, NumericVector Q2){
	/*Compute the geodesic distance between quaternions Q1 and Q2*/
	/* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single quaternion*/
	
	int n = Q1.nrow(), i=0; 
	double cp;
	NumericVector rs(n);
	
	for(i=0;i<n;i++){
		
		cp = sum(Q1(i,_)*Q2);
		rs[i] = acos(2*cp*cp-1);
		
	}
	
	return rs;
}

// [[Rcpp::export]]
double oneRdistC(NumericMatrix Q1, NumericVector Q2){
	/*Compute the geodesic distance between quaternions Q1 and Q2*/
	/* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single quaternion*/
	
	double cp=0.0;
	double rs=0.0;
	cp = sum(Q1*Q2);
	rs = acos(2*cp*cp-1);
		
	return rs;
}


// [[Rcpp::export]]
NumericVector cdfunsC(NumericMatrix Qs, NumericVector Qhat){
	
	//Compute the values c and d to form the pivotal test statistic
	//This estimates c=2E(1-cos(r^2))/3 and d=E(1+2cos(r))/3 from a sample
	int n = Qs.nrow(), i;
	double crs;
	
	NumericVector cds(2);
	cds[0]=0.0;
	cds[1]=0.0;
	
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
NumericVector zhangQ4(NumericMatrix Q, int m){
	
	
	int n=Q.nrow(), i=0, j=0;
	NumericVector cdstar;
	IntegerVector samp(n);
	NumericVector unSamp;
	int numUn=0, maxSamp=0;
	
	NumericVector testStat(m);
	double sqrth=0.0;
	arma::mat Qstar(n,4);
	NumericVector QhatStar;
  NumericMatrix QhatStarMat(1,4);
	
	arma::mat QSamp = as<arma::mat>(Q); //Convert the sample into armadillo mode
	
	NumericMatrix QstarRcpp;
	
	NumericVector Qhat = as<NumericVector>(wrap(rotations2::meanQ4C(QSamp)));
	//NumericVector Qhat = as<NumericVector>(wrap(meanQ4C(QSamp)));
	
	for(j=0;j<m;j++){
		
		samp = floor(runif(n,0,n));			//Bootstrap sample of size n, with replacement
	  unSamp = unique(samp);
    numUn = unSamp.size();
    maxSamp = max(samp);
    
    while(numUn<4 || maxSamp>n-1){
      samp = floor(runif(n,0,n));	 //If bootstrap samp is less than 4 obs then							
	    unSamp = unique(samp);       //draw a new sample
      numUn = unSamp.size();
      maxSamp = max(samp);
    }
    
		for(i=0;i<n;i++){
			
			Qstar.row(i) = QSamp.row(samp[i]);		//Copying a matrix row by row produces a bunch of junk messages
		}																				//so I do it with arma instead of standard Rcpp
	
		QhatStar = as<NumericVector>(wrap(rotations2::meanQ4C(Qstar))); //Both of these functinos return arma variables so
		//QhatStar = as<NumericVector>(wrap(meanQ4C(Qstar)));
		QstarRcpp = as<NumericMatrix>(wrap(Qstar));					//They need to be converted to Rcpp type
		
		cdstar = cdfunsC(QstarRcpp,QhatStar);
		
		QhatStarMat = as<NumericMatrix>(wrap(QhatStar)); /*QhatStar needs to be a matrix to be used in RdistC*/
		sqrth = oneRdistC(QhatStarMat,Qhat);
		
		if(cdstar[0]<0.00001){
			printf("c is too small");
		}
		
		testStat[j] = 2*n*pow(cdstar[1],2)*pow(sqrth,2)/cdstar[0];
		
	}
	
	return testStat;
	
}

