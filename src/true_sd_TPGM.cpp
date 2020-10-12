#include <Rcpp.h>
//#include <RcppParallel.h>
#include <cmath>
#include <algorithm>
using namespace std;
using namespace Rcpp;
//using namespace RcppParallel;

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

/*
*** By Rong (Mar. 5th, 2019)
*** Function: true standard deviation of each edge parameter from TPGM
*/

/* ****** Sub-function, log factorial ****** */
double log_factorial_sd_TPGM(double m){
	double log_factorial_value = 0;
	if(m != 0 && m != 1){
		for(double k = 2; k <= m; k++){
			log_factorial_value = log_factorial_value + log(k);
		}
	}
	return log_factorial_value;
}

/* ****** Sub-function, P1 calculation ****** */
NumericVector P1_calculation_sd_TPGM(NumericVector predict, double threshold){
	size_t n = predict.size();
	NumericVector p1(n);
	NumericVector numerator(n);
	NumericVector denominator(n);

	//for(size_t i = 0; i < n; i++){
		for(double j = 0; j <= threshold; j++){
			numerator = numerator + j*exp(j*predict-log_factorial_sd_TPGM(j));
			denominator = denominator + exp(j*predict-log_factorial_sd_TPGM(j));
		}
		p1 = numerator / denominator;
	//}
	return p1;
}


/* ****** Sub-function, P2 calculation ****** */
NumericVector P2_calculation_sd_TPGM(NumericVector predict, double threshold){
	size_t n = predict.size();
	NumericVector p2(n);
	NumericVector numerator(n);
	NumericVector denominator(n);
	//for(size_t i = 0; i < n; i++){
		for(double j = 0; j <= threshold; j++){
			numerator = numerator + pow(j,2)*exp(j*predict-log_factorial_sd_TPGM(j));
			denominator = denominator + exp(j*predict-log_factorial_sd_TPGM(j));
		}
		p2 = numerator / denominator;
	//}
	return p2;
}



/* ***** Function 1 ******** */
// [[Rcpp::export]]
NumericMatrix true_sd_TPGM(NumericMatrix x, NumericVector psi, NumericMatrix theta, NumericVector D){

	  //List data_split = data_splitting(x);

	  //NumericMatrix x_step1 = data_split["x_step1"];
	  //NumericMatrix x_step2 = data_split["x_step2"];

	  //size_t n1 = x_step1.nrow();
	  //size_t p = x_step1.ncol();
	  size_t n = x.nrow();
	  size_t p = x.ncol();

	  NumericMatrix true_stand_deviation(p,p);
	  NumericVector phi(p);

	 for(size_t i = 0; i < p; i++){
		 phi[i] = psi[i];
	 }

	 //Symmetrization of initial estimates
	 //for(size_t i = 0; i < p; i++){
		 //for(size_t j = 0; j < p; j++){
			 //theta_initial(i,j) = (theta_initial(i,j)+theta_initial(j,i))/2;
			 //theta_initial(j,i) = theta_initial(i,j);
		 //}
	 //}

	 //Bias Correction on beta_initial and symmetrization
	 //size_t n2 = x_step2.nrow();
	 NumericMatrix m(n,p);
	 for (size_t i = 0; i < p; i++){
		 for (size_t j = 0; j < n; j++){
			 m(j,i) = inner_product(x(j,_).begin(), x(j,_).end(), theta(i,_).begin(), 0.0) + phi[i];
		 }
	  }

	 bool display_progress = true;
	 Progress d(p*(p-1)/2, display_progress);

	 for (size_t i = 0; i < p-1; i++){

		 	 	if (Progress::check_abort() ){
		 		 		return -1.0;
		 		}

		        NumericVector tau = m(_,i);
			    NumericVector psi1 = P1_calculation_sd_TPGM(tau, D[i]);
			    NumericVector psi2 = P2_calculation_sd_TPGM(tau, D[i]) - pow(psi1,2);
			    NumericVector v1 = x(_,i)-psi1;

	    	for(size_t j = i+1; j < p; j++){

	    		d.increment();

	            NumericVector m1(n);
	            NumericVector m2(n);

	            m1 = m(_,i)-theta(i,j)*x(_,j)-phi[i];
	            m2 = m(_,j)-theta(j,i)*x(_,i)-phi[j];

	            NumericVector g_numerator(n);
	            NumericVector g_denominator(n);
	            for(size_t k2 = 0; k2 <= D[j]; k2++){
		            NumericVector capital_m(n);
	            	for(size_t k1 = 0; k1 <= D[i]; k1++){
	            		capital_m = capital_m + exp(phi[i]*k1-log_factorial_sd_TPGM(k1)+phi[j]*k2-log_factorial_sd_TPGM(k2)+theta(i,j)*k1*k2+k1*m1+k2*m2);
	            	}
	            	NumericVector inter_step = (P2_calculation_sd_TPGM(k2*theta(i,j)+m1+phi[i],D[i])-pow(P1_calculation_sd_TPGM(k2*theta(i,j)+m1+phi[i],D[i]),2))*capital_m;
	            	g_numerator = g_numerator + k2*inter_step;
	            	g_denominator = g_denominator + inter_step;
	            }

	            NumericVector g = g_numerator / g_denominator;

		        NumericVector v = x(_,j) - g;


	            double var_new = 1/(sum(pow(v,2) * psi2));
	            double standard_new = sqrt(var_new);


	            // Estimated standard deviation
	            true_stand_deviation(i,j) = standard_new;
	            true_stand_deviation(j,i) = true_stand_deviation(i,j);
	    	}
	  }
	  return true_stand_deviation;
}


