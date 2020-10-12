#include <Rcpp.h>
//#include <RcppEigen.h>
//#include <RcppParallel.h>
#include <cmath>
#include <algorithm>
//#include <vector>
using namespace std;
using namespace Rcpp;
//using namespace RcppParallel;
//using namespace Eigen;

/* 
 *** By Rong (Mar.1st, 2019)
 *** Function: a Gibbs sampler to generate random samples from TPGM
*/

/* ***** Sub-function: calculate the log value of factorials ***** */
//double log_factorial(double x){
//	if(x < 0){
//	    stop("The log value of factorials is not available. Please input a value >= 0.");
//	}
	
//	double sum = 0;
//	for(double k = 1; k <= x; k++){
//			sum = sum + log(k);
//	}
//	return sum;
//}

// [[Rcpp::export]]
NumericMatrix Gibbs_TPGM(NumericVector psi, NumericMatrix theta, NumericVector D, size_t nSample, Nullable<double> burn_in=R_NilValue, Nullable<double> thin=R_NilValue){
	size_t p = theta.nrow();  // the variable size p
	NumericMatrix sample(nSample, p);
	
	double rcpp_burn_in;
	if(burn_in.isNull()){
		rcpp_burn_in = 1000;
		Rcout << "Use default burn-in = 1000" << endl;
	}
	else{
		rcpp_burn_in = as<double>(burn_in);
		Rcout << "In this case, burn-in = " << rcpp_burn_in << endl;
	}
	
    double rcpp_thin;
	if(thin.isNull()){
		rcpp_thin = 100;
		Rcout << "Use default thin = 100" << endl;
	}
	else{
		rcpp_thin = as<double>(thin);
		Rcout << "In this case, thin = " << rcpp_thin << endl;
	}
	
	// initialize the sample with the original Poisson distribution with parameter = 1	
    NumericVector x = rpois(p,1);	
	
		// burn-in period
		for(size_t i = 0; i < rcpp_burn_in; i++){			
			// update for each coordinate x_[j] conditional on the other coordinates x_{-j}
			for(size_t j = 0; j < p; j++){
				NumericVector inner = theta(j,_)*x;
				double mu = psi[j] + sum(inner);
				NumericVector numerator(D[j]+1);
				double log_factorial = 0; 
				for(double m = 0; m <= D[j]; m++){
					if(m != 0 && m != 1){
						log_factorial = log_factorial + log(m);
					}
					numerator[m] = exp(m*mu - log_factorial);
				}
				double denominator = sum(numerator);
			    
				// define the proposal distribution for each X_[j], which is the conditional distribution of X_[j] given X_{-j}
				NumericVector prob = numerator / denominator;
				
				// cumulative probability
				NumericVector cum_prob = cumsum(prob);
				double tmp = R::runif(0,1);
				for(size_t l = 0; l < (D[j]+1); l++){
					if(tmp <= cum_prob[l]){
						x[j] = l;
						break;
					}
				}				
			}
		}
	
	
		// start take samples after burn-in period with a specified thin
	    for(size_t k = 0; k < nSample; k++){	
	    	for(size_t i = 0; i < rcpp_thin; i++){
	    		for(size_t j = 0; j < p; j++){
	    				NumericVector inner = theta(j,_)*x;
	    				double mu = psi[j] + sum(inner);
	    				NumericVector numerator(D[j]+1);
	    				double log_factorial = 0; 
	    				for(double m = 0; m <= D[j]; m++){
	    					if(m != 0 && m != 1){
	    						log_factorial = log_factorial + log(m);
	    					}
	    					numerator[m] = exp(m*mu - log_factorial);
	    				}
	    				double denominator = sum(numerator);
	    					    
	    				// define the proposal distribution for each X_[j], which is the conditional distribution of X_[j] given X_{-j}
	    				NumericVector prob = numerator / denominator;
	    						
	    				// cumulative probability
	    				NumericVector cum_prob = cumsum(prob);
	    				double tmp = R::runif(0,1);
	    				for(size_t l = 0; l < (D[j]+1); l++){
	    					if(tmp <= cum_prob[l]){
	    						x[j] = l;
	    						break;
	    					}
	    				}				
	    		}
	      }
		
	
              sample(k,_) = x;
	
		      if(k % 10 == 0){
			       Rcout <<"sample " << (k+1) << endl;
		      }
       }
	
   return sample;
}




