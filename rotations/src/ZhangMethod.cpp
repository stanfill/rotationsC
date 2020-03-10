#include "ZhangMethod.h"
#include "estimators.h"

Rcpp::NumericVector RdistC(Rcpp::NumericMatrix Q1, Rcpp::NumericVector Q2)
{
	/*Compute the geodesic distance between quaternions Q1 and Q2*/
	/* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single quaternion*/

	int n = Q1.nrow(), i=0;
	double cp;
	Rcpp::NumericVector rs(n);

	for(i=0;i<n;i++){

		cp = sum(Q1(i, Rcpp::_)*Q2);
		rs[i] = acos(2*cp*cp-1);

	}

	return rs;
}

arma::rowvec rdistSO3C(arma::mat Rs, arma::mat R2)
{
  int n = Rs.n_rows, m=Rs.n_cols , i,j;
  double tri;
  arma::mat R2t = R2.t();

  if(m==3){

  	Rs = Rs * R2t;
    //Rs.print("Rs2:");

    arma::rowvec theta(1);
    tri = trace(Rs);

    if((3-tri)<10e-10){
      theta(0) = 0;
    }else{
      theta(0) = acos(0.5*tri-0.5);
    }

    return theta;
  }

  arma::rowvec theta(n);
  theta.zeros();
  arma::mat33 Rsi;

  for(i=0; i<n ; i++){

    for(j = 0; j<9 ;j++){
      Rsi(j)=Rs(i,j);
    }

    Rsi = Rsi * R2t;

    tri = trace(Rsi);
    //printf("Trace: %lf",tri);
    if((3-tri)<10e-10){
      theta(i) = 0;
    }else{
      theta(i) = acos(0.5*tri-0.5);
    }

  }
  return theta;
}

Rcpp::NumericVector EdistC(Rcpp::NumericMatrix Q1, Rcpp::NumericVector Q2)
{
	/*Compute the Euclidean distance between quaternions Q1 and Q2*/
	/* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single quaternion*/

	int n = Q1.nrow(), i=0;
	double cp, rsi;
	Rcpp::NumericVector rs(n);

	for(i=0;i<n;i++){

		cp = sum(Q1(i, Rcpp::_)*Q2);
		rsi = 8*(1-(cp*cp));
		rs[i] = pow(rsi,0.5);

	}

	return rs;
}

double oneRdistC(Rcpp::NumericMatrix Q1, Rcpp::NumericVector Q2)
{
	/*Compute the geodesic distance between quaternions Q1 and Q2*/
	/* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single quaternion*/

	double cp=0.0;
	double rs=0.0;
	cp = sum(Q1*Q2);
	rs = acos(2*cp*cp-1);

	return rs;
}

Rcpp::NumericVector cdfunsC(Rcpp::NumericMatrix Qs, Rcpp::NumericVector Qhat)
{
	//Compute the projected mean values of c and d to form the pivotal test statistic
	//This estimates c=2E(1-cos(r^2))/3 and d=E(1+2cos(r))/3 from a sample
	int n = Qs.nrow(), i;
	double crs;

	Rcpp::NumericVector cds(2);
	cds[0]=0.0;
	cds[1]=0.0;

	Rcpp::NumericVector rs(n);

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

Rcpp::NumericVector cdfunsCMedian(Rcpp::NumericMatrix Qs, Rcpp::NumericVector Qhat)
{
	//Compute the values c and d to form the pivotal test statistic for the median using quaternions

	int n = Qs.nrow(), i;
	double crs, OnemCrs;

	Rcpp::NumericVector cds(2);
	cds[0]=0.0;
	cds[1]=0.0;

	Rcpp::NumericVector rs(n);

	rs = RdistC(Qs,Qhat);

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

Rcpp::NumericVector zhangQ4(Rcpp::NumericMatrix Q, int m)
{
	int n=Q.nrow(), i=0, j=0;
  Rcpp::NumericVector cdstar;
  Rcpp::IntegerVector samp(n);
  Rcpp::NumericVector unSamp;
	int numUn=0, maxSamp=0;

	Rcpp::NumericVector testStat(m);
	double sqrth=0.0;
	arma::mat Qstar(n,4);
	Rcpp::NumericVector QhatStar;
	Rcpp::NumericMatrix QhatStarMat(1,4);

	arma::mat QSamp = Rcpp::as<arma::mat>(Q); //Convert the sample into armadillo mode

	Rcpp::NumericMatrix QstarRcpp;

	Rcpp::NumericVector Qhat = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(meanQ4C(QSamp)));

	for(j=0;j<m;j++){

		samp = floor(Rcpp::runif(n,0,n));			//Bootstrap sample of size n, with replacement
	  unSamp = unique(samp);
    numUn = unSamp.size();
    maxSamp = max(samp);

    while(numUn<4 || maxSamp>n-1){
      samp = floor(Rcpp::runif(n,0,n));	 //If bootstrap samp is less than 4 obs then
	    unSamp = unique(samp);       //draw a new sample
      numUn = unSamp.size();
      maxSamp = max(samp);
    }

		for(i=0;i<n;i++){

			Qstar.row(i) = QSamp.row(samp[i]);		//Copying a matrix row by row produces a bunch of junk messages
		}																				//so I do it with arma instead of standard Rcpp

		QhatStar = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(meanQ4C(Qstar))); //Both of these functinos return arma variables so
		QstarRcpp = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(Qstar));					//They need to be converted to Rcpp type

		cdstar = cdfunsC(QstarRcpp,QhatStar);

		QhatStarMat = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(QhatStar)); /*QhatStar needs to be a matrix to be used in RdistC*/
		sqrth = oneRdistC(QhatStarMat,Qhat);

		if(cdstar[0]<0.0000001){
			//printf("c is too small");
      cdstar[0]=0.0000001;
		}

		testStat[j] = 2*n*pow(cdstar[1],2)*pow(sqrth,2)/cdstar[0];

	}

	return testStat;
}

Rcpp::NumericVector cdfunsCSO3(arma::mat Rs, arma::mat Rhat)
{
	//Compute the projected median values of c and d to form the pivotal test statistic
  //for SO3 data, used by the zhangMedianC function

	int n = Rs.n_rows, i;
	double crs, OnemCrs;

	Rcpp::NumericVector cds(2);
	cds[0]=0.0;
	cds[1]=0.0;

	Rcpp::NumericVector rs(n);

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

Rcpp::NumericVector zhangMedianC(arma::mat Rs, int m)
{
	//Compute the bootstrap version of the chang regions for SO3 data because that is what the
	//median function is written for

	Rcpp::RNGScope scope; // using runif requires this to be set...I think.
  								// This has been shown to cause problems in the past so consider using the next line in its place

  //GetRNGstate();PutRNGstate();

  int n = Rs.n_rows, i,j;
  arma::mat Shat = medianSO3C(Rs,2000,1e-5);
  arma::mat Rstar(n,9);
  arma::mat Sstar(3,3);
  Rcpp::NumericVector cdstar(2);
  Rcpp::IntegerVector samp(n);
  Rcpp::NumericVector unSamp;
	int numUn=0, maxSamp=0;
	Rcpp::NumericVector hsqrtMedian;
	Rcpp::NumericVector hstar(m);
  double hsq=0.0;

  for(j=0;j<m;j++){

		samp = floor(Rcpp::runif(n,0,n));			//Bootstrap sample of size n, with replacement
	  unSamp = unique(samp);
    numUn = unSamp.size();
    maxSamp = max(samp);

    while(numUn<4 || maxSamp>n-1){
      samp = floor(Rcpp::runif(n,0,n));	 //If bootstrap samp is less than 4 obs then
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
