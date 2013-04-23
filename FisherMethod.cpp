#include <RcppArmadillo.h>   
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]  
double fisherAxisC(arma::mat Qs, arma::vec Qhat){
	
	arma::mat Qsq=Qs.t()*Qs;
	arma::mat eigvec;
	arma::vec eigval;
  
  arma::eig_sym(eigval,eigvec,Qsq);   
  arma::vec qhat=eigvec.col(3);
	arma::mat G(3,3);
  
  int i, j, k, n=Qs.n_rows;
  double Tm, denom;
  arma::vec m=eigvec.col(3);
  
  /*if(m[0]<0){
  	m = -m;
  }*/
  
  arma::mat Mhat = eigvec.submat(0,0,3,2);
  
  /*for(i=0;i<3;i++){
    Mhat.col(i)=eigvec.col(i);
  }*/
  
  for(j=0;j<3;j++){
		for(k=0;k<3;k++){
      
			denom=arma::as_scalar(pow((n*(eigval(3)-eigval(j))*(eigval(3)-eigval(k))),-1));
				
			for(i=0;i<n;i++){
				G(j,k) = G(j,k) + pow(sum(Qs.row(i)*Mhat.col(j))*sum(Qs.row(i)*Mhat.col(k))*sum(Qs.row(i)*eigvec.col(3)),2)*denom;
			}
      G(k,j)=G(j,k);
		}
	}
  
  arma::mat Ginv = G.i();
  
  Tm=arma::as_scalar(n*Qhat.t()*Mhat*Ginv*Mhat.t()*Qhat);
  
  return Tm;
}

/*** R
library(rotations)
library(microbenchmark)
Qs<-ruars(20,rcayley,space='Q4')
Qhat<-mean(Qs)

fisherAxisCompute<-function(Qs,Shat){
  
	n<-nrow(Qs)
	
	svdQs<-svd(t(Qs)%*%Qs/n)
	mhat<-svdQs$v[,1]
	Mhat<-t(svdQs$v[,-1])
	etad<-svdQs$d[1]
	etas<-svdQs$d[-1]
	G<-matrix(0,3,3)
		
	for(j in 1:3){
		for(k in j:3){
			denom<-1/(n*(etad-etas[j])*(etad-etas[k]))
				
			for(i in 1:n){
				G[j,k]<-G[k,j]<-G[j,k]+(Mhat[j,]%*%Qs[i,])*(Mhat[k,]%*%Qs[i,])*(mhat%*%Qs[i,])^2*denom
			}
		}
	}
	
	Ginv<-solve(G)
	
	#Ginv<-try(solve(G),silent=T)
		
	#if(class(Ginv)!='matrix'){
	#	Ginv<-diag(1/diag(G))
	#}
		
	Tm<-n*Shat%*%t(Mhat)%*%Ginv%*%Mhat%*%t(Shat)
	
	return(Tm)
}
tim<-microbenchmark(
fisherAxisC(Qs,Qhat),
fisherAxisCompute(Qs,Qhat))
print(tim)
plot(tim)
*/

