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
	
	//Compute the values c and d to form the pivotal test statistic
	
	int n = Qs.nrow();
	double crs;
	
	NumericVector cds(2);
	NumericVector rs(n);

	
	for(int i=0; i<n; i++){
		
		rs[i] = RdistC(Qs(i,_),Qhat); //Get misorientation angle of Qs[i,] and hat(Q)
		
		crs = cos(rs[i]);
		
		cds[0] += 1-pow(crs,2);				//c=2E[1-cos(r)^2]/3
		cds[1] += 1+2*crs;						//d=E[1+2cos(r)]/3
	}
	
	cds[0] = (2*cds[0])/(3*n);
	cds[1] = (cds[1]/(3*n));
	
	return cds;
}

// [[Rcpp::export]]
NumericVector bootQhat(NumericMatrix Q, int m){
	
	int n=Q.nrow(), i=0, j=0, numUn;
	NumericVector samp, cdstar;
	
	NumericVector testStat(m);
	
	arma::mat Qstar(n,4);
	NumericVector QhatStar;
  NumericVector unSamp;
	
	arma::mat QSamp = as<arma::mat>(Q); //Convert the sample into armadillo mode
	
	NumericMatrix QstarRcpp;
	
	NumericVector Qhat = as<Rcpp::NumericVector>(wrap(meanQ4C(QSamp)));
	
	for(j=0;j<m;j++){
		
		samp = floor(runif(n,0,n));			//Bootstrap sample of size n, with replacement
	  unSamp = unique(samp);
    numUn = unSamp.size();
    
    while(numUn<4){
      samp = floor(runif(n,0,n));	 //If bootstrap samp is less than 3 obs then							
	    unSamp = unique(samp);       //draw a new sample
      numUn = unSamp.size();
    }
    
    
		for(i=0;i<n;i++){
			Qstar.row(i) = QSamp.row(samp[i]);		//Copying a matrix row by row produces a bunch of junk messages
		}																				//so I do it with arma instead of standard Rcpp
	
		QhatStar = as<Rcpp::NumericVector>(wrap(meanQ4C(Qstar)));
		QstarRcpp = as<Rcpp::NumericMatrix>(wrap(Qstar));
		
		cdstar = cdfunsC(QstarRcpp,QhatStar);
		
		testStat[j] = 2*n*pow(cdstar[1],2)*pow(RdistC(QhatStar,Qhat),2)/cdstar[0];
		
	}
	
	return testStat;
	
}

/*** R

#library(microbenchmark)
#library(rotations)
#source("U:/Thesis/Intervals/Code/IntervalFuns.R")
#n<-10
#rs<-rcayley(n)
#Qs<-genR(rs,space='Q4')
#Rs<-SO3(Qs)

#cTest<-bootQhat(Qs,100)
#cTest
#xs<-seq(0,max(cTest),length=1000)
#hist(cTest,breaks=100,prob=T)
#lines(xs,dchisq(xs,3))

#tim<-microbenchmark(
#	bootQhat(Qs,300),
#	ZhangCI(Rs,300,.95)
#)
#print(tim)
#plot(tim)

*/

