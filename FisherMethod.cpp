#include <RcppArmadillo.h>   
#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]  
double fisherAxisC(arma::mat Qs, arma::vec Qhat){
	
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
			
			//printf(" %lf ",denom);
			
			for(i=0;i<n;i++){
				G(j,k) = G(j,k) + arma::as_scalar(Qs.row(i)*Mhat.col(j)*Qs.row(i)*Mhat.col(k)*pow(Qs.row(i)*eigvec.col(3),2));
			}
      
      G(j,k) = G(j,k)*denom;
      G(k,j)=G(j,k);
		}
	}
  

  //G.print("G:");
  arma::mat Ginv = G.i();  
  
  Tm = arma::as_scalar(n*Qhat.t()*Mhat*Ginv*Mhat.t()*Qhat);
  
  return Tm;
}


// [[Rcpp::export]]   
arma::vec meanQ4C(arma::mat Q) { 
  
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


//[[Rcpp::export]]
arma::vec fisherBootC(arma::mat Qs, int m){
  
  int n = Qs.n_rows;
  int i , j , numUn;
  arma::vec qhat=meanQ4C(Qs);
  arma::vec Tm(m);
  NumericVector samp, unSamp;
  arma::mat Qstar(n,4);
  
  for(i=0;i<m;i++){
    
    samp = floor(runif(n,0,n));			//Bootstrap sample of size n, with replacement
	  unSamp = unique(samp);
    numUn = unSamp.size();
    
    while(numUn<4){
      samp = floor(runif(n,0,n));	 //If bootstrap samp is less than 3 obs then							
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

#library(rotations)
#library(microbenchmark)
#Qs<-ruars(100,rcayley,space='Q4',kappa=5)
#Qhat<-mean(Qs)

#Fisher<-fisherBootC(Qs,300)
#hist(Fisher,breaks=100,prob=T)
#ss<-seq(0,max(Fisher),length=1000)
#lines(ss,dchisq(ss,3))

#tim<-microbenchmark(
#fisherAxisC(Qs,Qhat),
#fisherAxisCompute(Qs,Qhat))
#print(tim)
#plot(tim)
*/

