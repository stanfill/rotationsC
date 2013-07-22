#include <RcppArmadillo.h>   
#include <Rcpp.h>
//#include "../inst/include/rotations2.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]   
arma::rowvec meanQ4C(arma::mat Q) { 
	//Compute the projected mean of the sample Q.
	//NumericMatrix Qss = as<NumericMatrix>(wrap(Q));
	//int cq4 = checkQ4(Qss);
	//if(cq4){
	//	throw Rcpp::exception("The data are not in Q4.");
//	}
	
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

//' This estimates c=2E(1-cos(r^2))/3 and d=E(1+2cos(r))/3 from a sample
// [[Rcpp::export]]
NumericVector cdfunsC(NumericMatrix Qs, NumericVector Qhat){
	
	//Compute the values c and d to form the pivotal test statistic
	
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

//'Zhang bootstrap method paramaterized for quaternions
// [[Rcpp::export]]
NumericVector zhangQ4(NumericMatrix Q, int m){
	
	RNGScope scope; // using runif requires this to be set...I think  								
									// This has been shown to cause problems in the past so consider using the next line in its place
  //GetRNGstate();PutRNGstate();
  
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
	
	//NumericVector Qhat = as<NumericVector>(wrap(rotations2::meanQ4C(QSamp)));
	NumericVector Qhat = as<NumericVector>(wrap(meanQ4C(QSamp)));
	
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
	
		//QhatStar = as<NumericVector>(wrap(rotations2::meanQ4C(Qstar))); //Both of these functinos return arma variables so
		QhatStar = as<NumericVector>(wrap(meanQ4C(Qstar)));
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


//' Below here are all the functions for the Projected median

//' compute the riemannian distance between R1 and R2
// [[Rcpp::export]]
arma::rowvec rdistSO3C(arma::mat Rs, arma::mat R2){
  
  int n = Rs.n_rows, m=Rs.n_cols , i,j;
  
  if(m==3){
    Rs = Rs * R2.t();
    arma::rowvec theta(1); 
    theta(0) = acos(0.5*trace(Rs)-0.5);
    return theta;
  }
  
  
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

//' This estimates c=2E(1-cos(r^2))/3 and d=E(1+2cos(r))/3 from a sample
// [[Rcpp::export]]
NumericVector cdfunsCSO3(arma::mat Rs, arma::mat Rhat){
  
	//Compute the values c and d to form the pivotal test statistic
	
	int n = Rs.n_rows, i;
	double crs, OnemCrs;
	
	NumericVector cds(2);
	cds[0]=0.0;
	cds[1]=0.0;
	
	NumericVector rs(n);
	
	rs = rdistSO3C(Rs,Rhat);

	for(i=0; i<n; i++){
		
		crs = cos(rs[i]);
		
		cds[0] += crs;				//c=E[1+cos(r)]/6
		
		OnemCrs = 1-crs;
		OnemCrs = std::max(OnemCrs,1e-5); //I think sqrt(1-crs) is close to zero and causing the crash, for now max sure
															//that doesn't happen by taking at least 1e-5
															
		cds[1] += (1+3*crs)*(pow(OnemCrs,-0.5));							//d=E([1+3cos(r)]/12*sqrt[1-cos(r)])
	}
	
	cds[0] = ((cds[0]/n)+1)/6;
	cds[1] = (cds[1]/n)/12;
	
	return cds;
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


// [[Rcpp::export]]
NumericVector zhangMedianC(arma::mat Rs, int m){
	
  RNGScope scope; // using runif requires this to be set...I think.  
  								// This has been shown to cause problems in the past so consider using the next line in its place
  
  //GetRNGstate();PutRNGstate();
  
  int n = Rs.n_rows, i,j;
  arma::mat Shat = medianSO3C(Rs,2000,1e-5);
  arma::mat Rstar(n,9);
  arma::mat Sstar(3,3);
  NumericVector cdstar(2);
	IntegerVector samp(n);
	NumericVector unSamp;
	int numUn=0, maxSamp=0;
  NumericVector hsqrtMedian;
  NumericVector hstar(m);
  double hsq=0.0;
  
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
			
			Rstar.row(i) = Rs.row(samp[i]);		//Copying a matrix row by row produces a bunch of junk messages
		}	
    
    Sstar = medianSO3C(Rstar,2000,1e-5);
    
    cdstar = cdfunsCSO3(Rstar,Sstar);
    hsqrtMedian = rdistSO3C(Shat,Sstar);
    hsq = hsqrtMedian[0];
    
    hstar[j]=2*n*pow(cdstar[1],2)*pow(hsq,2)/cdstar[0];
  
  }
  
  return hstar;
}
