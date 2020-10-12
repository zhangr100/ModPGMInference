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
 *** By Rong (Mar.6th, 2019)
 *** Function: a Gibbs sampler to generate random samples from SqrtPGM
*/

/* ***** Sub-function: calculate log-factorial ***** */
//double log_factorial(double x){
//	double sum = 0;
//	if(x > 0){
//		for(double k = 1; k <= x; k++){
//				sum = sum + log(k);
//		}
//	}
//	return sum;
//}

// [[Rcpp::export]]
NumericMatrix Gibbs_SqrtPGM(NumericVector psi, NumericMatrix theta, size_t nSample, Nullable<double> burn_in=R_NilValue, Nullable<double> thin=R_NilValue){
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
	NumericVector Sqrt = sqrt(x);

	// burn-in period
	for(size_t i = 0; i < rcpp_burn_in; i++){
			
			// update for each coordinate x_[j] conditional on the other coordinates x_{-j}
			for(size_t j = 0; j < p; j++){
				NumericVector inner = theta(j,_)*Sqrt;
				double mu = psi[j] + sum(inner) - theta(j,j)*Sqrt[j];
				
				// calculate the sum in the denominator with error < 1e-3
				double sum = 0;
				double value;
				double m = 0;
				double log_factorial = 0;
				do{
					if(m != 0 && m != 1){
						log_factorial = log_factorial + log(m);
					}
					value = sum;
					sum = sum + exp(sqrt(m)*mu + theta(j,j)*m - log_factorial);
					double diff = sum - value;
					if(abs(diff) < 1e-3){
						break;
					}
					m = m + 1;
				}
				while(m <= 1e6);
				
				if(m == 1e6 + 1){
					m = m - 1;
				}
				
				double log_factorial_1 = 0;
				NumericVector numerator(m+1);
				for(double q = 0; q < (m+1); q++){
					if(q != 0 && q != 1){
						log_factorial_1 = log_factorial_1 + log(q);
					}
					numerator[q] = exp(sqrt(q)*mu + theta(j,j)*q - log_factorial_1);
				}
			    
				// define the proposal distribution for each X_[j], which is the conditional distribution of X_[j] given X_{-j}
				NumericVector prob = numerator / sum;
				
				// cumulative probability
				NumericVector cum_prob = cumsum(prob);
				double tmp = R::runif(0,max(cum_prob));
				for(size_t l = 0; l < (m+1); l++){
					if(tmp <= cum_prob[l]){
						x[j] = l;
						Sqrt[j] = sqrt(x[j]);
						break;
					}
				}				
			}
		}
	
	// start take samples after burn-in period with a specified thin
		for(size_t k = 0; k < nSample; k++){	
			 for(size_t i = 0; i < rcpp_thin; i++){
				 for(size_t j = 0; j < p; j++){
				 			NumericVector inner = theta(j,_)*Sqrt;
				 			double mu = psi[j] + sum(inner) - theta(j,j)*Sqrt[j];
				 				
				 			// calculate the sum in the denominator with error < 1e-3
				 			double sum = 0;
				 			double value;
				 			double m = 0;
				 			double log_factorial = 0;
				 			do{
				 					if(m != 0 && m != 1){
				 						log_factorial = log_factorial + log(m);
				 					}
				 					value = sum;
				 					sum = sum + exp(sqrt(m)*mu + theta(j,j)*m - log_factorial);
				 					double diff = sum - value;
				 					if(abs(diff) < 1e-3){
				 						break;
				 					}
				 					m = m + 1;
				 			}
				 			while(m <= 1e6);
				 				
				 			if(m == 1e6 + 1){
				 				m = m - 1;
				 			}
				 				
				 			double log_factorial_1 = 0;
				 			NumericVector numerator(m+1);
				 			for(double q = 0; q < (m+1); q++){
				 				if(q != 0 && q != 1){
				 					log_factorial_1 = log_factorial_1 + log(q);
				 				}
				 				numerator[q] = exp(sqrt(q)*mu + theta(j,j)*q - log_factorial_1);
				 			}
				 			    
				 			// define the proposal distribution for each X_[j], which is the conditional distribution of X_[j] given X_{-j}
				 			NumericVector prob = numerator / sum;
				 				
				 			// cumulative probability
				 			NumericVector cum_prob = cumsum(prob);
				 			double tmp = R::runif(0,max(cum_prob));
				 			for(size_t l = 0; l < (m+1); l++){
				 				if(tmp <= cum_prob[l]){
				 					x[j] = l;
				 					Sqrt[j] = sqrt(x[j]);
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




