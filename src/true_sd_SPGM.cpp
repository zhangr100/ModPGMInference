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
*** Function: true standard deviation of each edge parameter from SPGM
*/

/* ***** Sub-function: transformation of the vector of the original sufficient statistics ***** */
NumericVector sublinear_sd_SPGM(NumericVector x, NumericVector lower, NumericVector upper){
	size_t p = x.size();
	NumericVector y(p);
	for(size_t i = 0; i < p; i++){
		if(x[i] <= lower[i]){
			y[i] = x[i];
		}
		else if(x[i] >= upper[i]){
			y[i] = (lower[i] + upper[i])/2;
		}
		else{
			y[i] = -1/(2*(upper[i] - lower[i]))*pow(x[i],2) + (upper[i]/(upper[i] - lower[i]))*x[i] - pow(lower[i],2)/(2*(upper[i]-lower[i]));
		}
	}
	return y;
}


/* ***** Sub-function: transformation of the original sufficient statistic ***** */
double sublinear2_sd_SPGM(double x, double lower, double upper){
	double y;
	if(x <= lower){
		y = x;
	}
	else if(x >= upper){
		y = (lower + upper)/2;
	}
	else{
		y = -1/(2*(upper - lower))*pow(x,2) + (upper/(upper - lower))*x - pow(lower,2)/(2*(upper - lower));
	}
	return y;
}



/* ****** Sub-function, log factorial ****** */
double log_factorial_sd_SPGM(double m){
	double log_factorial_value = 0;
	if(m != 0 && m != 1){
		for(double k = 2; k <= m; k++){
			log_factorial_value = log_factorial_value + log(k);
		}
	}
	return log_factorial_value;
}


/* ****** Sub-function, P1 calculation ****** */
NumericVector P1_calculation_sd_SPGM(NumericVector predict, double threshold1, double threshold2){
	size_t n = predict.size();
	NumericVector p1(n);
	NumericVector numerator(n);
	NumericVector denominator(n);

	NumericVector numerator_tmp(n);
	NumericVector denominator_tmp(n);
	double j = 0;
	do{
		for(size_t i = 0; i < n; i++){
				numerator_tmp[i] = numerator[i];
				denominator_tmp[i] = denominator[i];
		}
		numerator = numerator + sublinear2_sd_SPGM(j,threshold1,threshold2)*exp(sublinear2_sd_SPGM(j,threshold1,threshold2)*predict-log_factorial_sd_SPGM(j));
		denominator = denominator + exp(sublinear2_sd_SPGM(j,threshold1,threshold2)*predict-log_factorial_sd_SPGM(j));
		NumericVector diff1 = numerator - numerator_tmp;
		NumericVector diff2 = denominator - denominator_tmp;
		if((max(abs(diff1)) <= 1e-2) && (max(abs(diff2)) <= 1e-2)){
			break;
		}
		j = j + 1;
	}
	while(j <= 1e4);

	p1 = numerator / denominator;

	return p1;
}


/* ****** Sub-function, P2 calculation ****** */
NumericVector P2_calculation_sd_SPGM(NumericVector predict, double threshold1, double threshold2){
	size_t n = predict.size();
	NumericVector p2(n);
	NumericVector numerator(n);
	NumericVector denominator(n);

	NumericVector numerator_tmp(n);
	NumericVector denominator_tmp(n);
	double j = 0;
	do{
		for(size_t i = 0; i < n; i++){
				numerator_tmp[i] = numerator[i];
				denominator_tmp[i] = denominator[i];
		}
		numerator = numerator + pow(sublinear2_sd_SPGM(j,threshold1,threshold2),2)*exp(sublinear2_sd_SPGM(j,threshold1,threshold2)*predict-log_factorial_sd_SPGM(j));
		denominator = denominator + exp(sublinear2_sd_SPGM(j,threshold1,threshold2)*predict-log_factorial_sd_SPGM(j));
		NumericVector diff1 = numerator - numerator_tmp;
		NumericVector diff2 = denominator - denominator_tmp;
		if((max(abs(diff1)) <= 1e-2) && (max(abs(diff2)) <= 1e-2)){
			break;
		}
		j = j + 1;
	}
	while(j <= 1e4);

	p2 = numerator / denominator;

	return p2;
}



/* ***** Function 1 ******** */
// [[Rcpp::export]]
NumericMatrix true_sd_SPGM(NumericMatrix x, NumericVector psi, NumericMatrix theta, NumericVector D_0, NumericVector D_1){

	  //List data_split = data_splitting(x);

	  //NumericMatrix x_step1 = data_split["x_step1"];
	  //NumericMatrix x_step2 = data_split["x_step2"];

	  //size_t n1 = x_step1.nrow();
	  //size_t p = x_step1.ncol();
	  size_t n = x.nrow();
	  size_t p = x.ncol();

	  // Sublinear transformation on the original data matrix x
	  NumericMatrix sub_x(n,p);
	  for(size_t i = 0; i < n; i++){
		  sub_x(i,_) = sublinear_sd_SPGM(x(i,_), D_0, D_1);
	  }


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
			 m(j,i) = inner_product(sub_x(j,_).begin(), sub_x(j,_).end(), theta(i,_).begin(), 0.0) + phi[i];
		 }
	 }

	 bool display_progress = true;
	 Progress d(p*(p-1)/2, display_progress);

	 for (size_t i = 0; i < p-1; i++){

		 	 	if (Progress::check_abort() ){
		 		 		return -1.0;
		 		}

		 	 	NumericVector tau = m(_,i);
		 	 	NumericVector psi1 = P1_calculation_sd_SPGM(tau, D_0[i], D_1[i]);
		 	 	NumericVector psi2 = P2_calculation_sd_SPGM(tau, D_0[i], D_1[i]) - pow(psi1,2);
		 	 	NumericVector v1 = sub_x(_,i)-psi1;

	    	for(size_t j = i+1; j < p; j++){

	    		d.increment();

	            NumericVector m1(n);
	            NumericVector m2(n);

	            m1 = m(_,i)-theta(i,j)*sub_x(_,j)-phi[i];
	           	m2 = m(_,j)-theta(j,i)*sub_x(_,i)-phi[j];

	           	NumericVector g_numerator(n);
	           	NumericVector g_denominator(n);

	           	double k2 = 0;
	           	NumericVector g_numerator_tmp(n);
	           	NumericVector g_denominator_tmp(n);
	           	do{
	           	    //for(size_t k2 = 0; k2 <= D[j]; k2++){
	           	     for(size_t index = 0; index < n; index++){
	           	            	g_numerator_tmp[index] = g_numerator[index];
	           	            	g_denominator_tmp[index] = g_denominator[index];
	           	     }
	           		 NumericVector capital_m(n);
	           	     double k1 = 0;
	           		 NumericVector capital_m_tmp(n);
	           	     do{
	           	    	 //for(size_t k1 = 0; k1 <= D[i]; k1++){
	           	    	 for(size_t index1 = 0; index1 < n; index1++){
	           	    		 	 capital_m_tmp[index1] = capital_m[index1];
	           	         }
	           	         capital_m = capital_m + exp(phi[i]*sublinear2_sd_SPGM(k1,D_0[i],D_1[i])-log_factorial_sd_SPGM(k1)+phi[j]*sublinear2_sd_SPGM(k2,D_0[j],D_1[j])-log_factorial_sd_SPGM(k2)+theta(i,j)*sublinear2_sd_SPGM(k1,D_0[i],D_1[i])*sublinear2_sd_SPGM(k2,D_0[j],D_1[j])+sublinear2_sd_SPGM(k1,D_0[i],D_1[i])*m1+sublinear2_sd_SPGM(k2,D_0[j],D_1[j])*m2);
	           	         NumericVector diff = capital_m - capital_m_tmp;
	           	         if(max(abs(diff)) <= 1e-2){
	           	        	 break;
	           	         }
	           	         k1 = k1 + 1;
	           	     }
	           		 while(k1 <= 1e4);

	           	     NumericVector inter_step = (P2_calculation_sd_SPGM(sublinear2_sd_SPGM(k2,D_0[j],D_1[j])*theta(i,j)+m1+phi[i],D_0[i],D_1[i])-pow(P1_calculation_sd_SPGM(sublinear2_sd_SPGM(k2,D_0[j],D_1[j])*theta(i,j)+m1+phi[i],D_0[i],D_1[i]),2))*capital_m;
	           	     g_numerator = g_numerator + sublinear2_sd_SPGM(k2,D_0[j],D_1[j])*inter_step;
	           	     g_denominator = g_denominator + inter_step;

	           	     NumericVector diff1 = g_numerator - g_numerator_tmp;
	           	     NumericVector diff2 = g_denominator - g_denominator_tmp;

	           	     if((max(abs(diff1)) <= 1e-2) && (max(abs(diff2)) <= 1e-2)){
	           	         break;
	           	     }
	           	     k2 = k2 + 1;
	           	}
	           	while(k2 <= 1e4);

	            NumericVector g = g_numerator / g_denominator;

		        NumericVector v = sub_x(_,j) - g;

	            double var_new = 1/(sum(pow(v,2) * psi2));
	            double standard_new = sqrt(var_new);


	            // Estimated standard deviation
	            true_stand_deviation(i,j) = standard_new;
	            true_stand_deviation(j,i) = true_stand_deviation(i,j);
	    	}
	  }
	  return true_stand_deviation;
}







