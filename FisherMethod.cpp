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
double fisherAxisC(arma::mat Qs, arma::rowvec Qhat){

	NumericMatrix Qss = as<NumericMatrix>(wrap(Qs));
	int cq4 = checkQ4(Qss);
	if(cq4){
		throw Rcpp::exception("The data are not in Q4.");
	}

  int n=Qs.n_rows;
  
	arma::mat Qsq=(Qs.t()*Qs)/n;
	arma::mat eigvec, Mhat(4,3);
	arma::vec eigval, eta(3);
  
  arma::eig_sym(eigval,eigvec,Qsq);   

	arma::mat G(3,3);
	G.zeros();

  int i, j, k;
  double Tm, denom;

  
  for(i=0;i<3;i++){
    Mhat.col(i)=eigvec.col(i);
    eta(i)=eigval(i);
  }
  
  for(j=0;j<3;j++){
		for(k=j;k<3;k++){
      
			denom = pow((n*(eigval[3]-eta[j])*(eigval[3]-eta[k])),-1);
			
			for(i=0;i<n;i++){
				G(j,k) = G(j,k) + arma::as_scalar(Qs.row(i)*Mhat.col(j)*Qs.row(i)*Mhat.col(k)*pow(Qs.row(i)*eigvec.col(3),2));
			}
      
      G(j,k) = G(j,k)*denom;
      G(k,j)=G(j,k);
		}
	}
  
  arma::mat Ginv = G.i();  
  
  Tm = arma::as_scalar(n*Qhat*Mhat*Ginv*Mhat.t()*Qhat.t());
  
  return Tm;
}


// [[Rcpp::export]]   
arma::rowvec meanQ4C(arma::mat Q) { 
	/*Compute the projected mean of the sample Q*/
	
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


//[[Rcpp::export]]
arma::vec fisherBootC(arma::mat Qs, int m){

  int n = Qs.n_rows;
  int i , j , numUn;
  
  arma::rowvec qhat = meanQ4C(Qs);

  arma::vec Tm(m);
  NumericVector samp, unSamp;
  arma::mat Qstar(n,4);
  Qstar.zeros();
  
  for(i=0;i<m;i++){
    
    samp = floor(runif(n,0,n));			//Bootstrap sample of size n, with replacement
	  unSamp = unique(samp);
    numUn = unSamp.size();
    
    while(numUn<4){
      samp = floor(runif(n,0,n));	 //If bootstrap samp is less than 4 obs then							
	    unSamp = unique(samp);       //draw a new sample
      numUn = unSamp.size();
    }
    
    for(j=0;j<n;j++){
      Qstar.row(j) = Qs.row(samp[j]);
    }
    
    Tm[i]=fisherAxisC(Qstar,qhat);
    
  }
  return Tm;
}

/*** R

#Rcpp::sourceCpp("ZhangMethod.cpp")

library(rotations)
#library(microbenchmark)
Qs<-ruars(10,rcayley,space='Q4',kappa=5)
meanQ4C(Qs)

Fisher<-fisherBootC(Qs,300)
#hist(Fisher,breaks=100,prob=T)
#ss<-seq(0,max(Fisher),length=1000)
#lines(ss,dchisq(ss,3))

#tim<-microbenchmark(
#fisherAxisC(Qs,Qhat),
#fisherAxisCompute(Qs,Qhat))
#print(tim)
#plot(tim)
*/

