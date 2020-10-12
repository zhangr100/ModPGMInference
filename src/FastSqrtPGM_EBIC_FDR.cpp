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
*** Statistical inference of large SqrtPGM with EBIC
***
By Rong Zhang, 2020.09.13
*/


/* ****** Sub-function, log factorial ****** */
double log_factorial_SqrtPGM_EBIC(double m){
	double log_factorial_value = 0;
	if(m != 0 && m != 1){
		for(double k = 2; k <= m; k++){
			log_factorial_value = log_factorial_value + log(k);
		}
	}
	return log_factorial_value;
}

/* ****** Sub-function, log factorial vector ****** */
NumericVector log_factorial_vector_SqrtPGM_EBIC(NumericVector m){
	size_t n = m.size();
	NumericVector log_factorial_value(n);

	for(size_t i = 0; i < n; i++){
		if(m[i] != 0 && m[i] != 1){
			for(double k = 2; k <= m[i]; k++){
				log_factorial_value[i] = log_factorial_value[i] + log(k);
			}
		}
	}

	return log_factorial_value;
}

/* ****** Sub-function, P1 calculation ****** */
NumericVector P1_calculation_SqrtPGM_EBIC(NumericVector predict, double diag_interact){
	size_t n = predict.size();
	NumericVector numerator(n);
	NumericVector denominator(n);
	NumericVector p1(n);

	NumericVector numerator_tmp(n);
	NumericVector denominator_tmp(n);
	double j = 0;
	do{
		for(size_t i = 0; i < n; i++){
				numerator_tmp[i] = numerator[i];
				denominator_tmp[i] = denominator[i];
		}
		numerator = numerator + sqrt(j)*exp(sqrt(j)*predict + diag_interact*j - log_factorial_SqrtPGM_EBIC(j));
		denominator = denominator + exp(sqrt(j)*predict + diag_interact*j - log_factorial_SqrtPGM_EBIC(j));
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
NumericVector P2_calculation_SqrtPGM_EBIC(NumericVector predict, double diag_interact){
	size_t n = predict.size();
	NumericVector numerator(n);
	NumericVector denominator(n);
	NumericVector p2(n);

	NumericVector numerator_tmp(n);
	NumericVector denominator_tmp(n);
	double j = 0;
	do{
		for(size_t i = 0; i < n; i++){
				numerator_tmp[i] = numerator[i];
				denominator_tmp[i] = denominator[i];
		}
		numerator = numerator + j*exp(sqrt(j)*predict + diag_interact*j - log_factorial_SqrtPGM_EBIC(j));
		denominator = denominator + exp(sqrt(j)*predict + diag_interact*j - log_factorial_SqrtPGM_EBIC(j));
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

/* ****** Sub-function, log term ****** */
NumericVector log_term_SqrtPGM_EBIC(NumericVector predict, double diag_interact){
	size_t n = predict.size();
	NumericVector log_term_SqrtPGM_EBIC(n);
	NumericVector log_term_SqrtPGM_EBIC_tmp(n);

	double j = 0;
	do{
		for(size_t i = 0; i < n; i++){
			log_term_SqrtPGM_EBIC_tmp[i] = log_term_SqrtPGM_EBIC[i];
		}
		log_term_SqrtPGM_EBIC = log_term_SqrtPGM_EBIC + exp(sqrt(j)*predict + diag_interact*j - log_factorial_SqrtPGM_EBIC(j));
		NumericVector diff = log_term_SqrtPGM_EBIC - log_term_SqrtPGM_EBIC_tmp;
		if(max(abs(diff)) <= 1e-2){
				break;
		}
		j = j + 1;
	}
	while(j <= 1e4);

	return log_term_SqrtPGM_EBIC;
}

/* ****** Sub-function, fast Lasso regression for SqrtPGM with intercept ******* */
List FastSqrtPGMLasso_Rcpp(NumericVector beta, NumericVector predictor, NumericVector p1, NumericVector p2, NumericVector ipy, NumericVector y, NumericMatrix x, double lambda, double step_size, size_t N){
  double tol = 1e-3; //threshold for stopping
  size_t p = ipy.size();
  double dbm;  //maximum of beta difference for while loop
  //NumericVector beta(p+1); //initialize output beta vector with length p and filled with 0.0
  size_t iter = 0; //number of iterations

  //NumericVector predictor(N);
  //NumericVector p1 = P1_calculation_SqrtPGM_EBIC(predictor, beta[0]);
  //NumericVector p2 = P2_calculation_SqrtPGM_EBIC(predictor, beta[0]);
  //NumericVector weight_predictor = weight*predictor;

  NumericVector x_plus(p+1);
  for(size_t j = 0; j < p+1; j++){
	  x_plus[j] = beta[j];
  }

  //update beta vector with proximal gradient descent
  do{

	double current_g = -sum(y*predictor + beta[0]*pow(y,2) - log_factorial_vector_SqrtPGM_EBIC(pow(y,2)) - log(log_term_SqrtPGM_EBIC(predictor,beta[0])))/N;

	// compute the gradient
	NumericVector gradient(p+1);
    for(size_t j = 0; j < p+1; j++){
      if(j == 0){
    	  //gradient[j] = sum(p2 - pow(y,2))/N;  // gradient of theta_ii
    	  gradient[j] = 0;
      }
      else{
    	  gradient[j] = (sum(x(_,j-1)*p1) - ipy[j-1])/N; // gradient of phi_i and theta_ij
      }
    }

    // use backtracking line search to choose step size
    double g_new;
    double g_approximate;
    double alpha = 1;
    NumericVector predictor_update(N);
    size_t iter1 = 0;
    do{
    	for(size_t i = 0; i < N; i++){
    		predictor_update[i] = predictor[i];
    	}
    	for(size_t j = 0; j < p+1; j++){
    		double beta_tmp;
    		double diffBeta;
    		double z = beta[j] - alpha*gradient[j];
    		if(j == 0){
    			beta_tmp = max(0.0, z) - max(0.0, -z);
    			diffBeta = beta_tmp - beta[j];
    			if(abs(diffBeta) > 0){
    				x_plus[j] = beta_tmp;
    			}
    		}
    		if(j == 1){
    			beta_tmp = max(0.0, z) - max(0.0, -z);
    			diffBeta = beta_tmp - beta[j];
    			if(abs(diffBeta) > 0){
    				x_plus[j] = beta_tmp;
    				predictor_update = predictor_update + x(_,j-1)*diffBeta;
    			}
    		}
    		if(j > 1){
    			beta_tmp = max(0.0, z - alpha*lambda) - max(0.0, -z - alpha*lambda);
    			diffBeta = beta_tmp - beta[j];
    			if(abs(diffBeta) > 0){
    				x_plus[j] = beta_tmp;
    				predictor_update = predictor_update + x(_,j-1)*diffBeta;
    			}
    		}
    	}
    	g_new = -sum(y*predictor_update + x_plus[0]*pow(y,2) - log_factorial_vector_SqrtPGM_EBIC(pow(y,2)) - log(log_term_SqrtPGM_EBIC(predictor_update,x_plus[0])))/N;
    	g_approximate = current_g + inner_product(gradient.begin(), gradient.end(), (x_plus-beta).begin(), 0.0) + 1/(2*alpha)*(sum(pow(x_plus-beta,2)));
    	alpha = step_size*alpha;
    	iter1++;
    }
    while(g_new > g_approximate && iter1 <= 1e6);

    dbm = max(abs(x_plus - beta));

    for(size_t i = 0; i < N; i++){
    	predictor[i] = predictor_update[i];
    }
    for(size_t j = 0; j < p+1; j++){
    	beta[j] = x_plus[j];
    }

    p1 = P1_calculation_SqrtPGM_EBIC(predictor, beta[0]);
    p2 = P2_calculation_SqrtPGM_EBIC(predictor, beta[0]);

    iter++;
  }
  while(dbm >= tol && iter <= 1e6);

  return List::create(_["beta"] = beta, _["predictor"] = predictor, _["p1"] = p1, _["p2"] = p2);
}

/* ****** Sub-function, fast Lasso regression for SqrtPGM without intercept ****** */
List FastSqrtPGMLasso_Rcpp2(NumericVector beta, NumericVector predictor, NumericVector p1, NumericVector p2, NumericVector ipy, NumericVector y, NumericMatrix x, double lambda, double step_size, size_t N){
	  double tol = 1e-3; //threshold for stopping
	  size_t p = ipy.size();
	  double dbm;  //maximum of beta difference for while loop
	  //NumericVector beta(p+1); //initialize output beta vector with length p and filled with 0.0
	  size_t iter = 0; //number of iterations

	  //NumericVector predictor(N);
	  //NumericVector p1 = P1_calculation_SqrtPGM_EBIC(predictor, beta[0]);
	  //NumericVector p2 = P2_calculation_SqrtPGM_EBIC(predictor, beta[0]);
	  //NumericVector weight_predictor = weight*predictor;

	  NumericVector x_plus(p+1);
	  for(size_t j = 0; j < p+1; j++){
	  	  	  x_plus[j] = beta[j];
	  }

	  //update beta vector with proximal gradient descent
	  do{

		double current_g = -sum(y*predictor + beta[0]*pow(y,2) - log_factorial_vector_SqrtPGM_EBIC(pow(y,2)) - log(log_term_SqrtPGM_EBIC(predictor,beta[0])))/N;

		// compute the gradient
		NumericVector gradient(p+1);
	    for(size_t j = 0; j < p+1; j++){
	      if(j == 0){
	    	  //gradient[j] = sum(p2 - pow(y,2))/N;  // gradient of theta_ii
	    	  gradient[j] = 0;
	      }
	      else{
	    	  gradient[j] = (sum(x(_,j-1)*p1) - ipy[j-1])/N; // gradient of theta_ij
	      }
	    }

	    // use backtracking line search to choose step size
	    double g_new;
	    double g_approximate;
	    double alpha = 1;
	    NumericVector predictor_update(N);
        size_t iter1 = 0;
	    do{
	    	for(size_t i = 0; i < N; i++){
	    	    predictor_update[i] = predictor[i];
	    	}
	    	for(size_t j = 0; j < p+1; j++){
	    		double beta_tmp;
	    		double diffBeta;
	    		double z = beta[j] - alpha*gradient[j];
	    		if(j == 0){
	    			beta_tmp = max(0.0, z) - max(0.0, -z);
	    			diffBeta = beta_tmp - beta[j];
	    			if(abs(diffBeta) > 0){
	    				x_plus[j] = beta_tmp;
	    			}
	    		}
	    		else{
	    			beta_tmp = max(0.0, z - alpha*lambda) - max(0.0, -z - alpha*lambda);
	    			diffBeta = beta_tmp - beta[j];
	    			if(abs(diffBeta) > 0){
	    				x_plus[j] = beta_tmp;
	    				predictor_update = predictor_update + x(_,j-1)*diffBeta;
	    			}
	    		}
	    	}
	    	g_new = -sum(y*predictor_update + x_plus[0]*pow(y,2) - log_factorial_vector_SqrtPGM_EBIC(pow(y,2)) - log(log_term_SqrtPGM_EBIC(predictor_update,x_plus[0])))/N;
	    	g_approximate = current_g + inner_product(gradient.begin(), gradient.end(), (x_plus-beta).begin(), 0.0) + 1/(2*alpha)*(sum(pow(x_plus-beta,2)));
	    	alpha = step_size*alpha;
	    	iter1++;
	    }
	    while(g_new > g_approximate && iter1 <= 1e6);

	    dbm = max(abs(x_plus - beta));

	    for(size_t i = 0; i < N; i++){
	        	predictor[i] = predictor_update[i];
	    }
	    for(size_t j = 0; j < p+1; j++){
	        	beta[j] = x_plus[j];
	    }

	    p1 = P1_calculation_SqrtPGM_EBIC(predictor, beta[0]);
	    p2 = P2_calculation_SqrtPGM_EBIC(predictor, beta[0]);

	    iter++;
	  }
	  while(dbm >= tol && iter <= 1e6);

	  return List::create(_["beta"] = beta, _["predictor"] = predictor, _["p1"] = p1, _["p2"] = p2);
}


/* ***** n choose k ***** */
//size_t ChooseC(size_t n,size_t k){
//	return Rf_choose(n,k);
//}

/* ***** Sub-function, return index of minimum element ***** */
double min_position_SqrtPGM_EBIC(NumericVector x){
	 NumericVector::iterator it = min_element(x.begin(), x.end());
	 double position = it - x.begin();
	 return position;
}

/* ***** sort the vector in descending order ***** */
NumericVector sortDescending_SqrtPGM_EBIC(NumericVector x){
	sort(x.begin(), x.end(), greater<double>());
	return x;
}


/* ***** Function 1 ******** */
// [[Rcpp::export]]
List FastSqrtPGM_EBIC_FDR(NumericMatrix x, double gamma = NA_REAL, double nlambda = NA_REAL, double step_size = NA_REAL, bool intercept=true, bool global=false, Nullable<NumericVector> regularization=R_NilValue, Nullable<NumericVector> alpha=R_NilValue, double N=NA_REAL, double delta_upper=NA_REAL, Nullable<NumericMatrix> true_graph=R_NilValue){

	  size_t n = x.nrow();
	  size_t p = x.ncol();

	  if(NumericVector::is_na(step_size)){
	    step_size = 0.5;
	    Rcout << "Use default step size = 0.5" << endl;
	  }
	  else{
	    Rcout << "In this case, step size = " << step_size << endl;
	  }

if(regularization.isNull()){
	  if(global==false){
	  NumericMatrix theta_cor(p,p);
	  NumericMatrix CI_low_theta(p,p);
	  NumericMatrix CI_high_theta(p,p);
	  NumericMatrix p_thetaCor(p,p);
	  NumericMatrix est_stand_deviation(p,p);
	  NumericMatrix test_statistic(p,p);
	  NumericMatrix theta_initial(p,p);
	  NumericVector phi(p);

	  if(intercept == true){

		  // Take square root of each element in the matrix
		  NumericMatrix sqrt_x(n,p);
		  for(size_t i = 0; i < p; i++){
			  sqrt_x(_,i) = sqrt(x(_,i));
		  }

		  // Include intercept in the square root matrix
		  NumericMatrix sqrt_x_new(n,p+1);
		  NumericVector intercept(n,1);
		  for(size_t i = 0; i < (p+1); i++){
			  if(i == 0){
				  sqrt_x_new(_,i) = intercept;
			  }
			  else{
			  	  sqrt_x_new(_,i) = sqrt_x(_,i-1);
			  }
		  }

	  //penalty parameter for L1 norm of betas
		  if(NumericVector::is_na(gamma)){
		  		  gamma = 0.5;
		  		  Rcout << "Use default gamma = 0.5" << endl;
		  }
		  else{
		  		  Rcout << "In this case, gamma =" << gamma << endl;
		  }

		  if(NumericVector::is_na(nlambda)){
		  		  nlambda = 20;
		  		  Rcout << "Use default number of tuning parameters = 20" << endl;
		  }
		  else{
		  		  Rcout << "In this case, number of tuning parameters = " << nlambda << endl;
		  }

		  NumericMatrix sqrt_IP_YX_new(p+1,p+1);
		  for(size_t i = 0; i < p+1; i++){
		 //NumericVector response = x_step1(_,i)+1;
			  for(size_t j = 0; j < p+1; j++){
				  sqrt_IP_YX_new(i,j) = inner_product(sqrt_x_new(_,i).begin(), sqrt_x_new(_,i).end(), sqrt_x_new(_,j).begin(), 0.0);
			  }
		  }

		  NumericMatrix sqrt_IP_YX_origin(p,p);
		  for(size_t i = 0; i < p; i++){
		  	//NumericVector response = x_step1(_,i)+1;
		  	  for(size_t j = 0; j < p; j++){
		  		  sqrt_IP_YX_origin(i,j) = inner_product(sqrt_x(_,i).begin(), sqrt_x(_,i).end(), sqrt_x(_,j).begin(), 0.0);
		  	  }
		   }


	  //NumericMatrix supnorm(100,p);

		  for(size_t i = 1; i < p+1; i++){

			  //if(i % 100 == 0){
				  Rcout <<"Lasso of Square Root PGM for variable " << i << endl;
			  //}

		 // Extract IP_YX[i,-i] and x[,-i]
	     NumericVector sqrt_tmpyx(p);
	     NumericVector sqrt_tmpyx_origin(p-1);
		 NumericMatrix sqrt_tmp_x(n,p);
		 size_t t = 0;
		 for (size_t k = 0; k < p+1; k++){
		     if (k != i){
		        sqrt_tmpyx[t] = sqrt_IP_YX_new(i,k);
		        sqrt_tmp_x(_,t) = sqrt_x_new(_,k);
		        t++;
		      }
		 }

		 t = 0;
		 for(size_t k = 0; k < p; k++){
		 	 if(k != (i-1)){
		 		 sqrt_tmpyx_origin[t] = sqrt_IP_YX_origin(i-1,k);
		 		 t++;
		 	 }
		 }

		 double max_number = max(abs(sqrt_tmpyx_origin))/n;

		 // Take 100 tuning parameters
		 NumericVector lambda(nlambda);

		 for(double k = 0; k < nlambda; k++){
			 lambda[k] = max_number/(pow(1000,k/(nlambda-1)));
		 }

		 NumericMatrix beta_tmp(nlambda,p+1);
		 NumericVector BIC_value(nlambda);

		 NumericVector beta(p+1);
		 NumericVector predictor(n);
		 NumericVector p1 = P1_calculation_SqrtPGM_EBIC(predictor, beta[0]);
		 NumericVector p2 = P2_calculation_SqrtPGM_EBIC(predictor, beta[0]);


		 for (size_t k = 0; k < nlambda; k++){

			 List Lasso_fit = FastSqrtPGMLasso_Rcpp(beta, predictor, p1, p2, sqrt_tmpyx, sqrt_x_new(_,i), sqrt_tmp_x, lambda[k], step_size, n);
			 beta = Lasso_fit["beta"];
			 predictor = Lasso_fit["predictor"];
			 p1 = Lasso_fit["p1"];
			 p2 = Lasso_fit["p2"];

			 beta_tmp(k,_) = beta;

			 //Return the number of non-zero coefficients
			 size_t count = 0;
			 for (size_t j = 2; j < p+1; j++){
				 if(beta[j] != 0){
					 count++;
		 		 }
		 	 }

		 	 NumericVector log_likelihood = sqrt_x_new(_,i)*predictor + beta[0]*pow(sqrt_x_new(_,i),2) - log_factorial_vector_SqrtPGM_EBIC(pow(sqrt_x_new(_,i),2)) - log(log_term_SqrtPGM_EBIC(predictor,beta[0]));
		 	 BIC_value[k] = -2*sum(log_likelihood) + count*(log(n)+2*gamma*log(p-1));
		 }

		 // Return the position of the minimum BIC
		 double position = min_position_SqrtPGM_EBIC(BIC_value);

		 NumericVector beta_select = beta_tmp(position,_);

		 phi[i-1] = beta_select[1];

		 for (size_t k = 0; k < p; k++){
			 	 if(i-1 == k){
					 theta_initial(i-1,k) = beta_select[0];
				 }
				 if(i-1 < k){
					 theta_initial(i-1,k) = beta_select[k+1];
				 }
				 if(i-1 > k){
					 theta_initial(i-1,k) = beta_select[k+2];
				 }
		}

	 }

	 //Symmetrization of initial estimates
	 //for(size_t i = 0; i < p; i++){
		 //for(size_t j = 0; j < p; j++){
			 //theta_initial(i,j) = (theta_initial(i,j)+theta_initial(j,i))/2;
			 //theta_initial(j,i) = theta_initial(i,j);
		 //}
	 //}

	 Rcout << "Bias correction for initial estimators. " << endl;

	 //Bias Correction on beta_initial and symmetrization
	 //size_t n2 = x_step2.nrow();
	 NumericMatrix m(n,p);
	 for (size_t i = 0; i < p; i++){
		 for (size_t j = 0; j < n; j++){
			 m(j,i) = inner_product(sqrt_x(j,_).begin(), sqrt_x(j,_).end(), theta_initial(i,_).begin(), 0.0) + phi[i];
		 }
	  }

	 bool display_progress = true;
	 Progress d(p*(p-1)/2, display_progress);

	 for (size_t i = 0; i < p-1; i++){

		 	 	if (Progress::check_abort() ){
		 		 		return -1.0;
		 		}

		        NumericVector tau = m(_,i) - theta_initial(i,i)*sqrt_x(_,i);
			    NumericVector psi1 = P1_calculation_SqrtPGM_EBIC(tau, theta_initial(i,i));
			    NumericVector psi2 = P2_calculation_SqrtPGM_EBIC(tau, theta_initial(i,i)) - pow(psi1,2);
			    NumericVector v1 = sqrt_x(_,i) - psi1;

	    	for(size_t j = i+1; j < p; j++){

	    		d.increment();

	            NumericVector m1(n);
	            NumericVector m2(n);

	            m1 = m(_,i)-theta_initial(i,j)*sqrt_x(_,j)-theta_initial(i,i)*sqrt_x(_,i)-phi[i];
	            m2 = m(_,j)-theta_initial(j,i)*sqrt_x(_,i)-theta_initial(j,j)*sqrt_x(_,j)-phi[j];

	            NumericVector g_numerator(n);
	            NumericVector g_denominator(n);

	            double k2 = 0;
	            NumericVector g_numerator_tmp(n);
	            NumericVector g_denominator_tmp(n);
	            do{
	            	   for(size_t index = 0; index < n; index++){
	            	           g_numerator_tmp[index] = g_numerator[index];
	            	           g_denominator_tmp[index] = g_denominator[index];
	            	   }
	            	   NumericVector capital_m(n);
	            	   double k1 = 0;
	            	   NumericVector capital_m_tmp(n);
	            	   do{
	            	           for(size_t index1 = 0; index1 < n; index1++){
	            	            	 capital_m_tmp[index1] = capital_m[index1];
	            	            }
	            	            capital_m = capital_m + exp(phi[i]*sqrt(k1)+theta_initial(i,i)*k1-log_factorial_SqrtPGM_EBIC(k1)+phi[j]*sqrt(k2)+theta_initial(j,j)*k2-log_factorial_SqrtPGM_EBIC(k2)+theta_initial(i,j)*sqrt(k1)*sqrt(k2)+sqrt(k1)*m1+sqrt(k2)*m2);
	            	            NumericVector diff = capital_m - capital_m_tmp;
	            	            if(max(abs(diff)) <= 1e-2){
	            	            	break;
	            	            }
	            	            k1 = k1 + 1;
	            	   }
	            	   while(k1 <= 1e4);

	            	   NumericVector inter_step = (P2_calculation_SqrtPGM_EBIC(sqrt(k2)*theta_initial(i,j)+m1+phi[i],theta_initial(i,i))-pow(P1_calculation_SqrtPGM_EBIC(sqrt(k2)*theta_initial(i,j)+m1+phi[i],theta_initial(i,i)),2))*capital_m;
	            	   g_numerator = g_numerator + sqrt(k2)*inter_step;
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

		        NumericVector v = sqrt_x(_,j) - g;

	            NumericVector tau1 = v*v1;
	            NumericVector tau2 = v*psi2*sqrt_x(_,j);

	            theta_cor(i,j) = theta_initial(i,j) + sum(tau1)/sum(tau2);
	            theta_cor(j,i) = theta_cor(i,j);

	            double var_new = 1/(sum(pow(v,2) * psi2));
	            double standard_new = sqrt(var_new);

	            // 95% confidence interval of de-biased beta
	            double z_95CI = 1.96;
	            CI_low_theta(i,j) = theta_cor(i,j) - z_95CI * standard_new;
	            CI_low_theta(j,i) = CI_low_theta(i,j);
	            CI_high_theta(i,j) = theta_cor(i,j) + z_95CI * standard_new;
	            CI_high_theta(j,i) = CI_high_theta(i,j);

	            double zscore = theta_cor(i,j)/standard_new;

	            //test statistic
	            test_statistic(i,j) = zscore;
	            test_statistic(j,i) = test_statistic(i,j);

	            // P-value for de-biased beta
	            p_thetaCor(i,j) = 2 * Rf_pnorm5(-abs(zscore), 0.0, 1.0, 1, 0);
	            p_thetaCor(j,i) = p_thetaCor(i,j);

	            // Estimated standard deviation
	            est_stand_deviation(i,j) = standard_new;
	            est_stand_deviation(j,i) = est_stand_deviation(i,j);
	    	}
	  }
	}
	else{
			  //penalty parameter for L1 norm of betas
				if(NumericVector::is_na(gamma)){
						gamma = 0.5;
						Rcout << "Use default gamma = 0.5" << endl;
				}
				else{
						Rcout << "In this case, gamma =" << gamma << endl;
				}

				if(NumericVector::is_na(nlambda)){
						nlambda = 20;
						Rcout << "Use default number of tuning parameters = 20" << endl;
				}
				else{
						Rcout << "In this case, number of tuning parameters = " << nlambda << endl;
				}


				// Take square root of each element in the matrix
				  NumericMatrix sqrt_x(n,p);
				  for(size_t i = 0; i < p; i++){
				  	  sqrt_x(_,i) = sqrt(x(_,i));
				  }


				  NumericMatrix sqrt_IP_YX(p,p);
				  for(size_t i = 0; i < p; i++){
				 //NumericVector response = x_step1(_,i)+1;
					  for(size_t j = 0; j < p; j++){
						  sqrt_IP_YX(i,j) = inner_product(sqrt_x(_,i).begin(), sqrt_x(_,i).end(), sqrt_x(_,j).begin(), 0.0);
					  }
				  }


			  //NumericMatrix supnorm(100,p);

				  for(size_t i = 0; i < p; i++){

					  //if(i % 100 == 0){
						  Rcout <<"Lasso of Square Root PGM for variable " << (i+1) << endl;
					  //}

				 // Extract IP_YX[i,-i] and x[,-i]
			     NumericVector sqrt_tmpyx(p-1);
				 NumericMatrix sqrt_tmp_x(n,p-1);
				 size_t t = 0;
				 for (size_t k = 0; k < p; k++){
				     if (k != i){
				        sqrt_tmpyx[t] = sqrt_IP_YX(i,k);
				        sqrt_tmp_x(_,t) = sqrt_x(_,k);
				        t++;
				      }
				 }

				 double max_number = max(abs(sqrt_tmpyx))/n;

				 // Take the tuning parameters
				 NumericVector lambda(nlambda);

				 for (double k = 0; k < nlambda; k++){
				 	 lambda[k] = max_number/(pow(1000,k/(nlambda-1)));
				 }

				 NumericMatrix beta_tmp(nlambda,p);
				 NumericVector BIC_value(nlambda);

				 NumericVector beta(p);
				 NumericVector predictor(n);
				 NumericVector p1 = P1_calculation_SqrtPGM_EBIC(predictor, beta[0]);
				 NumericVector p2 = P2_calculation_SqrtPGM_EBIC(predictor, beta[0]);

				 for (size_t k = 0; k < nlambda; k++){

					 List Lasso_fit = FastSqrtPGMLasso_Rcpp2(beta, predictor, p1, p2, sqrt_tmpyx, sqrt_x(_,i), sqrt_tmp_x, lambda[k], step_size, n);
				 	 beta = Lasso_fit["beta"];
				 	 predictor = Lasso_fit["predictor"];
				 	 p1 = Lasso_fit["p1"];
				 	 p2 = Lasso_fit["p2"];

				 	 beta_tmp(k,_) = beta;

				 	 //Return the number of non-zero coefficients
				 	 size_t count = 0;
				 	 for (size_t j = 1; j < p; j++){
				 		 if(beta[j] != 0){
				 			 count++;
				 		 }
				 	 }

				 	 NumericVector log_likelihood = sqrt_x(_,i)*predictor + beta[0]*pow(sqrt_x(_,i),2) - log_factorial_vector_SqrtPGM_EBIC(pow(sqrt_x(_,i),2)) - log(log_term_SqrtPGM_EBIC(predictor,beta[0]));
				 	 BIC_value[k] = -2*sum(log_likelihood) + count*(log(n)+2*gamma*log(p-1));

				 }

				 // Return the position of the minimum BIC
				 double position = min_position_SqrtPGM_EBIC(BIC_value);

				 NumericVector beta_select = beta_tmp(position,_);

				 for (size_t k = 0; k < p; k++){
						 if(i == k){
							 theta_initial(i,k) = beta_select[0];
						 }
						 if(i < k){
							 theta_initial(i,k) = beta_select[k];
						 }
						 if(i > k){
							 theta_initial(i,k) = beta_select[k+1];
						 }
				 }

			 }

			 //Symmetrization of initial estimates
			// for(size_t i = 0; i < p; i++){
				// for(size_t j = 0; j < p; j++){
					// theta_initial(i,j) = (theta_initial(i,j)+theta_initial(j,i))/2;
					// theta_initial(j,i) = theta_initial(i,j);
				// }
			// }

			 Rcout << "Bias correction for initial estimators. " << endl;

			 //Bias Correction on beta_initial and symmetrization
			 //size_t n2 = x_step2.nrow();
			 NumericMatrix m(n,p);
			 for (size_t i = 0; i < p; i++){
				 for (size_t j = 0; j < n; j++){
					 m(j,i) = inner_product(sqrt_x(j,_).begin(), sqrt_x(j,_).end(), theta_initial(i,_).begin(), 0.0);
				 }
			  }

			 bool display_progress = true;
			 Progress d(p*(p-1)/2, display_progress);

			 for (size_t i = 0; i < p-1; i++){
				 	 	if (Progress::check_abort() ){
				 			 return -1.0;
				 		}

				        NumericVector tau = m(_,i) - theta_initial(i,i)*sqrt_x(_,i);
					    NumericVector psi1 = P1_calculation_SqrtPGM_EBIC(tau, theta_initial(i,i));
					    NumericVector psi2 = P2_calculation_SqrtPGM_EBIC(tau, theta_initial(i,i)) - pow(psi1,2);;
					    NumericVector v1 = sqrt_x(_,i)-psi1;

			    	for(size_t j = i+1; j < p; j++){

			            d.increment();

			            NumericVector m1(n);
			            NumericVector m2(n);

			            m1 = m(_,i)-theta_initial(i,j)*sqrt_x(_,j)-theta_initial(i,i)*sqrt_x(_,i);
			            m2 = m(_,j)-theta_initial(j,i)*sqrt_x(_,i)-theta_initial(j,j)*sqrt_x(_,j);

			            NumericVector g_numerator(n);
			            NumericVector g_denominator(n);

			            double k2 = 0;
			            NumericVector g_numerator_tmp(n);
			            NumericVector g_denominator_tmp(n);
			            do{
			            	 for(size_t index = 0; index < n; index++){
			            	        g_numerator_tmp[index] = g_numerator[index];
			            	        g_denominator_tmp[index] = g_denominator[index];
			            	 }
			            	 NumericVector capital_m(n);
			            	 double k1 = 0;
			            	 NumericVector capital_m_tmp(n);
			            	 do{
			            	        for(size_t index1 = 0; index1 < n; index1++){
			            	            	 capital_m_tmp[index1] = capital_m[index1];
			            	        }
			            	        capital_m = capital_m + exp(theta_initial(i,i)*k1-log_factorial_SqrtPGM_EBIC(k1)+theta_initial(j,j)*k2-log_factorial_SqrtPGM_EBIC(k2)+theta_initial(i,j)*sqrt(k1)*sqrt(k2)+sqrt(k1)*m1+sqrt(k2)*m2);
			            	        NumericVector diff = capital_m - capital_m_tmp;
			            	        if(max(abs(diff)) <= 1e-2){
			            	            	 break;
			            	        }
			            	        k1 = k1 + 1;
			            	 }
			            	 while(k1 <= 1e4);

			            	 NumericVector inter_step = (P2_calculation_SqrtPGM_EBIC(sqrt(k2)*theta_initial(i,j)+m1,theta_initial(i,i))-pow(P1_calculation_SqrtPGM_EBIC(sqrt(k2)*theta_initial(i,j)+m1,theta_initial(i,i)),2))*capital_m;
			            	 g_numerator = g_numerator + sqrt(k2)*inter_step;
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

				        NumericVector v = sqrt_x(_,j) - g;

			            NumericVector tau1 = v*v1;
			            NumericVector tau2 = v*psi2*sqrt_x(_,j);

			            theta_cor(i,j) = theta_initial(i,j) + sum(tau1)/sum(tau2);
			            theta_cor(j,i) = theta_cor(i,j);

			            double var_new = 1/(sum(pow(v,2) * psi2));
			            double standard_new = sqrt(var_new);

			            // 95% confidence interval of de-biased beta
			            double z_95CI = 1.96;
			            CI_low_theta(i,j) = theta_cor(i,j) - z_95CI * standard_new;
			            CI_low_theta(j,i) = CI_low_theta(i,j);
			            CI_high_theta(i,j) = theta_cor(i,j) + z_95CI * standard_new;
			            CI_high_theta(j,i) = CI_high_theta(i,j);

			            double zscore = theta_cor(i,j)/standard_new;

			            //test statistic
			            test_statistic(i,j) = zscore;
			            test_statistic(j,i) = test_statistic(i,j);

			            // P-value for de-biased beta
			            p_thetaCor(i,j) = 2 * Rf_pnorm5(-abs(zscore), 0.0, 1.0, 1, 0);
			            p_thetaCor(j,i) = p_thetaCor(i,j);

			            // Estimated standard deviation
			            est_stand_deviation(i,j) = standard_new;
			            est_stand_deviation(j,i) = est_stand_deviation(i,j);
			    	}
			  }

			  //return List::create(_["intercept"] = phi, _["theta_initial"] = theta_initial, _["theta_cor"] = theta_cor, _["CI_low_theta"] = CI_low_theta, _["CI_high_theta"] = CI_high_theta, _["z_score"] = test_statistic, _["p_thetaCor"] = p_thetaCor, _["est_sd"] = est_stand_deviation);

	}

	  return List::create(_["intercept"] = phi, _["theta_initial"] = theta_initial, _["theta_cor"] = theta_cor, _["CI_low_theta"] = CI_low_theta, _["CI_high_theta"] = CI_high_theta, _["z_score"] = test_statistic, _["p_thetaCor"] = p_thetaCor, _["est_sd"] = est_stand_deviation);
	}
	else{
		NumericVector rcpp_level = NumericVector::create();
		if(alpha.isNull()){
				rcpp_level = NumericVector::create(0.05, 0.1);
		}
		else{
				rcpp_level = as<NumericVector>(alpha);
		}

		size_t alpha_size = rcpp_level.size();

		NumericVector FDR(alpha_size);
		NumericVector Power(alpha_size);
		NumericVector threshold(alpha_size);
		List global_decision(alpha_size);

		if(intercept==true){

			// Take square root of each element in the matrix
					  NumericMatrix sqrt_x(n,p);
					  for(size_t i = 0; i < p; i++){
						  sqrt_x(_,i) = sqrt(x(_,i));
					  }

					  // Include intercept in the square root matrix
					  NumericMatrix sqrt_x_new(n,p+1);
					  NumericVector intercept(n,1);
					  for(size_t i = 0; i < (p+1); i++){
						  if(i == 0){
							  sqrt_x_new(_,i) = intercept;
						  }
						  else{
						  	  sqrt_x_new(_,i) = sqrt_x(_,i-1);
						  }
					  }

					 Rcout <<"Choose noise level for each node-wise regression." << endl;

				  //penalty parameter for L1 norm of betas
					  if(NumericVector::is_na(gamma)){
					  		  gamma = 0.5;
					  		  Rcout << "Use default gamma = 0.5" << endl;
					  }
					  else{
					  		  Rcout << "In this case, gamma =" << gamma << endl;
					  }

					  if(NumericVector::is_na(nlambda)){
					  		  nlambda = 20;
					  		  Rcout << "Use default number of tuning parameters = 20" << endl;
					  }
					  else{
					  		  Rcout << "In this case, number of tuning parameters = " << nlambda << endl;
					  }

					  NumericMatrix sqrt_IP_YX_new(p+1,p+1);
					  for(size_t i = 0; i < p+1; i++){
					 //NumericVector response = x_step1(_,i)+1;
						  for(size_t j = 0; j < p+1; j++){
							  sqrt_IP_YX_new(i,j) = inner_product(sqrt_x_new(_,i).begin(), sqrt_x_new(_,i).end(), sqrt_x_new(_,j).begin(), 0.0);
						  }
					  }

					  NumericMatrix sqrt_IP_YX_origin(p,p);
					  for(size_t i = 0; i < p; i++){
					  	//NumericVector response = x_step1(_,i)+1;
					  	  for(size_t j = 0; j < p; j++){
					  		  sqrt_IP_YX_origin(i,j) = inner_product(sqrt_x(_,i).begin(), sqrt_x(_,i).end(), sqrt_x(_,j).begin(), 0.0);
					  	  }
					   }


				  //NumericMatrix supnorm(100,p);
					  NumericVector noise(p);

					  for(size_t i = 1; i < p+1; i++){

						  //if(i % 100 == 0){
							  Rcout <<"Lasso of Square Root PGM for variable " << i << endl;
						  //}

					 // Extract IP_YX[i,-i] and x[,-i]
				     NumericVector sqrt_tmpyx(p);
				     NumericVector sqrt_tmpyx_origin(p-1);
					 NumericMatrix sqrt_tmp_x(n,p);
					 size_t t = 0;
					 for (size_t k = 0; k < p+1; k++){
					     if (k != i){
					        sqrt_tmpyx[t] = sqrt_IP_YX_new(i,k);
					        sqrt_tmp_x(_,t) = sqrt_x_new(_,k);
					        t++;
					      }
					 }

					 t = 0;
					 for(size_t k = 0; k < p; k++){
					 	 if(k != (i-1)){
					 		 sqrt_tmpyx_origin[t] = sqrt_IP_YX_origin(i-1,k);
					 		 t++;
					 	 }
					 }

					 double max_number = max(abs(sqrt_tmpyx_origin))/n;

					 // Take 100 tuning parameters
					 NumericVector lambda(nlambda);

					 for(double k = 0; k < nlambda; k++){
						 lambda[k] = max_number/(pow(1000,k/(nlambda-1)));
					 }

					 NumericMatrix beta_tmp(nlambda,p+1);
					 NumericVector BIC_value(nlambda);
					 NumericVector variance_record(nlambda);

					 NumericVector beta(p+1);
					 NumericVector predictor(n);
					 NumericVector p1 = P1_calculation_SqrtPGM_EBIC(predictor, beta[0]);
					 NumericVector p2 = P2_calculation_SqrtPGM_EBIC(predictor, beta[0]);


					 for (size_t k = 0; k < nlambda; k++){

						 List Lasso_fit = FastSqrtPGMLasso_Rcpp(beta, predictor, p1, p2, sqrt_tmpyx, sqrt_x_new(_,i), sqrt_tmp_x, lambda[k], step_size, n);
						 beta = Lasso_fit["beta"];
						 predictor = Lasso_fit["predictor"];
						 p1 = Lasso_fit["p1"];
						 p2 = Lasso_fit["p2"];

						 beta_tmp(k,_) = beta;
						 variance_record[k] = max(p2 - pow(p1,2));

						 //Return the number of non-zero coefficients
						 size_t count = 0;
						 for (size_t j = 2; j < p+1; j++){
							 if(beta[j] != 0){
								 count++;
					 		 }
					 	 }

					 	 NumericVector log_likelihood = sqrt_x_new(_,i)*predictor + beta[0]*pow(sqrt_x_new(_,i),2) - log_factorial_vector_SqrtPGM_EBIC(pow(sqrt_x_new(_,i),2)) - log(log_term_SqrtPGM_EBIC(predictor,beta[0]));
					 	 BIC_value[k] = -2*sum(log_likelihood) + count*(log(n)+2*gamma*log(p-1));
					 }

					 // Return the position of the minimum BIC
					 double position = min_position_SqrtPGM_EBIC(BIC_value);
					 noise[i-1] = variance_record[position];

				 }

			if(NumericVector::is_na(N)){
			  N = 10;
			  Rcout << "Use default N = 10" << endl;
			}
			else{
			  Rcout << "In this case, N = " << N << endl;
			}

			if(NumericVector::is_na(delta_upper)){
			  delta_upper = 2;
			  Rcout << "Use default L = 2" << endl;
			}
			else{
			  Rcout << "In this case, L = " << delta_upper << endl;
			}

			double ndelta = delta_upper*N;

			Rcout << "The number of delta = " << ndelta << endl;

			// Define a sequence of delta from 0 to delta_upper
			NumericVector delta(ndelta);
			for(double m = 0; m < ndelta; m++){
					delta[m] = (m+1)/(ndelta/delta_upper);
			}

			delta = sortDescending_SqrtPGM_EBIC(delta);

			NumericVector error_cum(ndelta);
			List intercept_value(ndelta);
			List theta_init(ndelta);
			List theta_correction(ndelta);
			List Statistic(ndelta);
			List upper_stat_tmp(ndelta);
			List CI_low(ndelta);
			List CI_high(ndelta);
			List p_value(ndelta);
			List standard_deviation(ndelta);

			Rcout <<"Perform global inference." << endl;
			Rcout << "Use pre-specified level(s): " << rcpp_level << endl;

			//Take out the upper triangular of true graph if available
			NumericMatrix rcpp_true_graph(p,p);
			NumericVector omega_true(p*(p-1)/2);
			size_t t = 0;
			if(true_graph.isNull()){
					Rcout << "True graph is not available." << endl;
			}
			else{
					Rcout << "True graph is available." << endl;
					rcpp_true_graph = as<NumericMatrix>(true_graph);
					for (size_t i = 1; i < p; i++){
						for (size_t j = 0; j < i; j++){
								omega_true[t]=rcpp_true_graph(j,i);
								t++;
						}
					}
		   }

		   NumericMatrix tmp_beta(p+1,p);

		   NumericMatrix tmp_predictor(n,p);
		   NumericMatrix tmp_p1(n,p);
		   for(size_t i = 0; i < p; i++){
					tmp_p1(_,i) = P1_calculation_SqrtPGM_EBIC(tmp_predictor(_,i), tmp_beta(0,i));
		   }
		   NumericMatrix tmp_p2(n,p);
		   for(size_t i = 0; i < p; i++){
				    tmp_p2(_,i) = P2_calculation_SqrtPGM_EBIC(tmp_predictor(_,i), tmp_beta(0,i));
		   }

		   size_t count = 0;
		   Rcout << "Calculate Lasso of each variable with tuning parameters under each delta." << endl;
		   Rcout << "Record test statistics under each delta." << endl;

		   for(size_t delta_count = 0; delta_count < ndelta; delta_count++){
			   NumericMatrix theta_cor(p,p);
			   NumericMatrix CI_low_theta(p,p);
			   NumericMatrix CI_high_theta(p,p);
			   NumericMatrix p_thetaCor(p,p);
			   NumericMatrix est_stand_deviation(p,p);
			   NumericMatrix test_statistic(p,p);
			   NumericMatrix theta_initial(p,p);
			   NumericVector phi(p);

			   count++;

			   Rcout <<"delta " << (count) << endl;
			   for(size_t i = 1; i < p+1; i++){
					double lambda = delta[delta_count]*sqrt(noise[i-1]*(log(p)/n));

					Rcout <<"Lasso of Square Root PGM for variable " << i << endl;

					// Extract IP_YX[i,-i] and x[,-i]
					NumericVector sqrt_tmpyx(p);
					NumericMatrix sqrt_tmp_x(n,p);
					size_t t = 0;
					for (size_t k = 0; k < p+1; k++){
						if (k != i){
								sqrt_tmpyx[t] = sqrt_IP_YX_new(i,k);
								sqrt_tmp_x(_,t) = sqrt_x_new(_,k);
								t++;
						}
					}

					// Lasso
					NumericVector beta_column = tmp_beta(_,i-1);
					NumericVector predictor_column = tmp_predictor(_,i-1);
					NumericVector p1_column = tmp_p1(_,i-1);
					NumericVector p2_column = tmp_p2(_,i-1);

				  	List beta_result = FastSqrtPGMLasso_Rcpp(beta_column, predictor_column, p1_column, p2_column, sqrt_tmpyx, sqrt_x_new(_,i), sqrt_tmp_x, lambda, step_size, n);

				  	NumericVector beta = beta_result["beta"];
				  	NumericVector predictor = beta_result["predictor"];
				  	NumericVector p1 = beta_result["p1"];
				  	NumericVector p2 = beta_result["p2"];

				  	tmp_beta(_,i-1) = beta;
				  	tmp_predictor(_,i-1) = predictor;
				  	tmp_p1(_,i-1) = p1;
				  	tmp_p2(_,i-1) = p2;

				  	phi[i-1] = beta[1];

				  	for (size_t k = 0; k < p; k++){
				  		if(i-1 == k){
				  				theta_initial(i-1,k) = beta[0];
				  		}
				  		if(i-1 < k){
				  				theta_initial(i-1,k) = beta[k+1];
				  		}
				  		if(i-1 > k){
				  				theta_initial(i-1,k) = beta[k+2];
				  		}
				    }

			   }

			   Rcout << "Bias correction for initial estimators. " << endl;

			   	 //Bias Correction on beta_initial and symmetrization
			   	 //size_t n2 = x_step2.nrow();
			   	 NumericMatrix m(n,p);
			   	 for (size_t i = 0; i < p; i++){
			   		 for (size_t j = 0; j < n; j++){
			   			 m(j,i) = inner_product(sqrt_x(j,_).begin(), sqrt_x(j,_).end(), theta_initial(i,_).begin(), 0.0) + phi[i];
			   		 }
			   	  }

			   	 bool display_progress = true;
			   	 Progress d(p*(p-1)/2, display_progress);

			   	 for (size_t i = 0; i < p-1; i++){

			   		 	 	if (Progress::check_abort() ){
			   		 		 		return -1.0;
			   		 		}

			   		        NumericVector tau = m(_,i) - theta_initial(i,i)*sqrt_x(_,i);
			   			    NumericVector psi1 = P1_calculation_SqrtPGM_EBIC(tau, theta_initial(i,i));
			   			    NumericVector psi2 = P2_calculation_SqrtPGM_EBIC(tau, theta_initial(i,i)) - pow(psi1,2);
			   			    NumericVector v1 = sqrt_x(_,i) - psi1;

			   	    	for(size_t j = i+1; j < p; j++){

			   	    		d.increment();

			   	            NumericVector m1(n);
			   	            NumericVector m2(n);

			   	            m1 = m(_,i)-theta_initial(i,j)*sqrt_x(_,j)-theta_initial(i,i)*sqrt_x(_,i)-phi[i];
			   	            m2 = m(_,j)-theta_initial(j,i)*sqrt_x(_,i)-theta_initial(j,j)*sqrt_x(_,j)-phi[j];

			   	            NumericVector g_numerator(n);
			   	            NumericVector g_denominator(n);

			   	            double k2 = 0;
			   	            NumericVector g_numerator_tmp(n);
			   	            NumericVector g_denominator_tmp(n);
			   	            do{
			   	            	   for(size_t index = 0; index < n; index++){
			   	            	           g_numerator_tmp[index] = g_numerator[index];
			   	            	           g_denominator_tmp[index] = g_denominator[index];
			   	            	   }
			   	            	   NumericVector capital_m(n);
			   	            	   double k1 = 0;
			   	            	   NumericVector capital_m_tmp(n);
			   	            	   do{
			   	            	           for(size_t index1 = 0; index1 < n; index1++){
			   	            	            	 capital_m_tmp[index1] = capital_m[index1];
			   	            	            }
			   	            	            capital_m = capital_m + exp(phi[i]*sqrt(k1)+theta_initial(i,i)*k1-log_factorial_SqrtPGM_EBIC(k1)+phi[j]*sqrt(k2)+theta_initial(j,j)*k2-log_factorial_SqrtPGM_EBIC(k2)+theta_initial(i,j)*sqrt(k1)*sqrt(k2)+sqrt(k1)*m1+sqrt(k2)*m2);
			   	            	            NumericVector diff = capital_m - capital_m_tmp;
			   	            	            if(max(abs(diff)) <= 1e-2){
			   	            	            	break;
			   	            	            }
			   	            	            k1 = k1 + 1;
			   	            	   }
			   	            	   while(k1 <= 1e4);

			   	            	   NumericVector inter_step = (P2_calculation_SqrtPGM_EBIC(sqrt(k2)*theta_initial(i,j)+m1+phi[i],theta_initial(i,i))-pow(P1_calculation_SqrtPGM_EBIC(sqrt(k2)*theta_initial(i,j)+m1+phi[i],theta_initial(i,i)),2))*capital_m;
			   	            	   g_numerator = g_numerator + sqrt(k2)*inter_step;
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

			   		        NumericVector v = sqrt_x(_,j) - g;

			   	            NumericVector tau1 = v*v1;
			   	            NumericVector tau2 = v*psi2*sqrt_x(_,j);

			   	            theta_cor(i,j) = theta_initial(i,j) + sum(tau1)/sum(tau2);
			   	            theta_cor(j,i) = theta_cor(i,j);

			   	            double var_new = 1/(sum(pow(v,2) * psi2));
			   	            double standard_new = sqrt(var_new);

			   	            // 95% confidence interval of de-biased beta
			   	            double z_95CI = 1.96;
			   	            CI_low_theta(i,j) = theta_cor(i,j) - z_95CI * standard_new;
			   	            CI_low_theta(j,i) = CI_low_theta(i,j);
			   	            CI_high_theta(i,j) = theta_cor(i,j) + z_95CI * standard_new;
			   	            CI_high_theta(j,i) = CI_high_theta(i,j);

			   	            double zscore = theta_cor(i,j)/standard_new;

			   	            //test statistic
			   	            test_statistic(i,j) = zscore;
			   	            test_statistic(j,i) = test_statistic(i,j);

			   	            // P-value for de-biased beta
			   	            p_thetaCor(i,j) = 2 * Rf_pnorm5(-abs(zscore), 0.0, 1.0, 1, 0);
			   	            p_thetaCor(j,i) = p_thetaCor(i,j);

			   	            // Estimated standard deviation
			   	            est_stand_deviation(i,j) = standard_new;
			   	            est_stand_deviation(j,i) = est_stand_deviation(i,j);
			   	    	}
			   	  }
			   	  intercept_value[delta_count] = phi;
			   	  theta_init[delta_count] = theta_initial;
			   	  theta_correction[delta_count] = theta_cor;
			   	  Statistic[delta_count] = test_statistic;
			   	  CI_low[delta_count] = CI_low_theta;
			   	  CI_high[delta_count] = CI_high_theta;
			   	  p_value[delta_count] = p_thetaCor;
			   	  standard_deviation[delta_count] = est_stand_deviation;

			   	  // Take out the upper triangular of test statistic
			   	  size_t t = 0;
			   	  NumericVector upper_stat(p*(p-1)/2);
			   	  for (size_t i = 1; i < p; i++){
			   		  for (size_t j = 0; j < i; j++){
			   			  	  upper_stat[t]=test_statistic(j,i);
			   			  	  t++;
			   		  }
			   	  }

			   	  // threshold for choosing delta
			   	  double error_1 = 0;
			   	  for (double k = 3; k < 10; k++){
			   		  	  double normal_tail = Rf_qnorm5(1-k/20, 0.0, 1.0, 1, 0);
			   		  	  double error = sum(abs(upper_stat)>=normal_tail)/(p*(p-1)*k/20)-1;
			   		  	  error_1 = error_1 + pow(error,2);
			   	  }

			   	  upper_stat_tmp[delta_count] = upper_stat;
			   	  error_cum[delta_count] = error_1;

		   }

		   	   	   	   //Reverse the error_cum
		   				reverse(error_cum.begin(), error_cum.end());


		   				//Return the minimum position for error

		   				Rcout << "Choose delta for FDR control." << endl;
		   				double position_select = min_position_SqrtPGM_EBIC(error_cum);


		   				NumericVector phi_select = intercept_value[ndelta-position_select-1];
		   				NumericMatrix theta_initial_select = theta_init[ndelta-position_select-1];
		   				NumericMatrix theta_cor_select = theta_correction[ndelta-position_select-1];
		   				NumericMatrix test_statistic_select = Statistic[ndelta-position_select-1];
		   				NumericMatrix CI_low_theta_select = CI_low[ndelta-position_select-1];
		   				NumericMatrix CI_high_theta_select = CI_high[ndelta-position_select-1];
		   				NumericMatrix p_thetaCor_select = p_value[ndelta-position_select-1];
		   				NumericMatrix est_stand_deviation_select = standard_deviation[ndelta-position_select-1];
		   				NumericVector upper_stat_select = upper_stat_tmp[ndelta-position_select-1];


		   				// Get threshold
		   				double upper = 2*sqrt(log(p));
		   				NumericVector z(2000);
		   				NumericVector numerator1(2000);
		   				NumericVector denominator(2000);
		   				for(double k = 0; k < 2000; k++){
		   						z[k] = ((k+1)/2000) * upper;
		   						numerator1[k] = 2 * (1-Rf_pnorm5(z[k], 0.0, 1.0, 1, 0)) *p*(p-1)/2;
		   						denominator[k] = sum(abs(upper_stat_select)>=z[k]);
		   						if(denominator[k] < 1){
		   								  denominator[k] = 1;
		   						}
		   				}
		   				NumericVector FDP = numerator1/denominator;


		   				NumericVector position(alpha_size);

		   				for(size_t i = 0; i < alpha_size; i++){
		   					for (double k = 0; k < 2000 ; k++){
		   							if(FDP[k] <= rcpp_level[i]){      // threshold for FDP <= pre-specified level
		   								  	threshold[i] = z[k];
		   								  	position[i] = k;
		   								  	break;
		   							}
		   							else{
		   								  if(k==1999){
		   								  		threshold[i] = z[k];
		   								  		position[i] = k;
		   								  }
		   							}
		   					}
		   				}

		   						if(true_graph.isNull()){
		   								  for(size_t i = 0; i < alpha_size; i++){
		   								  			FDR[i] = FDP[position[i]];
		   								  }
		   								  for(size_t k = 0; k < alpha_size; k++){
		   								  		  NumericMatrix decision(p,p);
		   								  		  for(size_t i = 0; i < p-1; i++){
		   								  			  for(size_t j = i+1; j < p; j++){
		   								  		   				   if(abs(test_statistic_select(i,j)) >= threshold[k]){
		   								  		   					   decision(i,j) = 1;
		   								  		   				   }
		   								  		   				   decision(j,i) = decision(i,j);
		   								  		   			   }
		   								  		   	  }
		   								  		   	  global_decision[k] = decision;
		   								  }
		   						}
		   						else{

		   								  for(size_t i = 0; i < alpha_size; i++){
		   								  	   	  double FDR_denominator = sum(abs(upper_stat_select)>= threshold[i]);
		   								  	   	  if(FDR_denominator < 1){
		   								  	   			FDR_denominator = 1;
		   								  	   	  }

		   								  	   	  FDR[i] = sum((abs(upper_stat_select) >= threshold[i]) & (omega_true == 0)) / FDR_denominator;

		   								  	   	  double Power_denominator = sum(omega_true!=0);
		   								  	   	  if(Power_denominator < 1){
		   								  	   			Power_denominator = 1;
		   								  	   	  }
		   								  	   	  Power[i] = sum((abs(upper_stat_select) >= threshold[i]) & (omega_true != 0)) / Power_denominator;
		   								  }
		   								  for(size_t k = 0; k < alpha_size; k++){
		   								  		  NumericMatrix decision(p,p);
		   								  		  for(size_t i = 0; i < p-1; i++){
		   								  		   		for(size_t j = i+1; j < p; j++){
		   								  		   				if(abs(test_statistic_select(i,j)) >= threshold[k]){
		   								  		   					   decision(i,j) = 1;
		   								  		   				}
		   								  		   				decision(j,i) = decision(i,j);
		   								  		   		}
		   								  		  }
		   								  		  global_decision[k] = decision;
		   								  }
		   						 }



		   						if(true_graph.isNull()){
		   								 return List::create(_["intercept"] = phi_select, _["theta_initial"] = theta_initial_select, _["theta_cor"] = theta_cor_select, _["CI_low_theta"] = CI_low_theta_select, _["CI_high_theta"] = CI_high_theta_select, _["z_score"] = test_statistic_select, _["p_thetaCor"] = p_thetaCor_select, _["est_sd"] = est_stand_deviation_select, _["threshold"] = threshold, _["FDR"] = FDR, _["global_decision"] = global_decision);
		   						}
		   						else{
		   								 return List::create(_["intercept"] = phi_select, _["theta_initial"] = theta_initial_select, _["theta_cor"] = theta_cor_select, _["CI_low_theta"] = CI_low_theta_select, _["CI_high_theta"] = CI_high_theta_select, _["z_score"] = test_statistic_select, _["p_thetaCor"] = p_thetaCor_select, _["est_sd"] = est_stand_deviation_select, _["threshold"] = threshold, _["FDR"] = FDR, _["power"] = Power, _["global_decision"] = global_decision);
		   						}

		}
		else{

			 	 	 	 	Rcout <<"Choose noise level for each node-wise regression." << endl;

							//penalty parameter for L1 norm of betas
							if(NumericVector::is_na(gamma)){
									gamma = 0.5;
									Rcout << "Use default gamma = 0.5" << endl;
							}
							else{
									Rcout << "In this case, gamma =" << gamma << endl;
							}

							if(NumericVector::is_na(nlambda)){
									nlambda = 20;
									Rcout << "Use default number of tuning parameters = 20" << endl;
							}
							else{
									Rcout << "In this case, number of tuning parameters = " << nlambda << endl;
							}


							// Take square root of each element in the matrix
							  NumericMatrix sqrt_x(n,p);
							  for(size_t i = 0; i < p; i++){
							  	  sqrt_x(_,i) = sqrt(x(_,i));
							  }


							  NumericMatrix sqrt_IP_YX(p,p);
							  for(size_t i = 0; i < p; i++){
							 //NumericVector response = x_step1(_,i)+1;
								  for(size_t j = 0; j < p; j++){
									  sqrt_IP_YX(i,j) = inner_product(sqrt_x(_,i).begin(), sqrt_x(_,i).end(), sqrt_x(_,j).begin(), 0.0);
								  }
							  }


						  //NumericMatrix supnorm(100,p);
							  NumericVector noise(p);

							  for(size_t i = 0; i < p; i++){

								  //if(i % 100 == 0){
									  Rcout <<"Lasso of Square Root PGM for variable " << (i+1) << endl;
								  //}

							 // Extract IP_YX[i,-i] and x[,-i]
						     NumericVector sqrt_tmpyx(p-1);
							 NumericMatrix sqrt_tmp_x(n,p-1);
							 size_t t = 0;
							 for (size_t k = 0; k < p; k++){
							     if (k != i){
							        sqrt_tmpyx[t] = sqrt_IP_YX(i,k);
							        sqrt_tmp_x(_,t) = sqrt_x(_,k);
							        t++;
							      }
							 }

							 double max_number = max(abs(sqrt_tmpyx))/n;

							 // Take the tuning parameters
							 NumericVector lambda(nlambda);

							 for (double k = 0; k < nlambda; k++){
							 	 lambda[k] = max_number/(pow(1000,k/(nlambda-1)));
							 }

							 NumericMatrix beta_tmp(nlambda,p);
							 NumericVector BIC_value(nlambda);
							 NumericVector variance_record(nlambda);

							 NumericVector beta(p);
							 NumericVector predictor(n);
							 NumericVector p1 = P1_calculation_SqrtPGM_EBIC(predictor, beta[0]);
							 NumericVector p2 = P2_calculation_SqrtPGM_EBIC(predictor, beta[0]);

							 for (size_t k = 0; k < nlambda; k++){

								 List Lasso_fit = FastSqrtPGMLasso_Rcpp2(beta, predictor, p1, p2, sqrt_tmpyx, sqrt_x(_,i), sqrt_tmp_x, lambda[k], step_size, n);
							 	 beta = Lasso_fit["beta"];
							 	 predictor = Lasso_fit["predictor"];
							 	 p1 = Lasso_fit["p1"];
							 	 p2 = Lasso_fit["p2"];

							 	 beta_tmp(k,_) = beta;
							 	 variance_record[k] = max(p2 - pow(p1,2));

							 	 //Return the number of non-zero coefficients
							 	 size_t count = 0;
							 	 for (size_t j = 1; j < p; j++){
							 		 if(beta[j] != 0){
							 			 count++;
							 		 }
							 	 }

							 	 NumericVector log_likelihood = sqrt_x(_,i)*predictor + beta[0]*pow(sqrt_x(_,i),2) - log_factorial_vector_SqrtPGM_EBIC(pow(sqrt_x(_,i),2)) - log(log_term_SqrtPGM_EBIC(predictor,beta[0]));
							 	 BIC_value[k] = -2*sum(log_likelihood) + count*(log(n)+2*gamma*log(p-1));

							 }

							 // Return the position of the minimum BIC
							 double position = min_position_SqrtPGM_EBIC(BIC_value);
							 noise[i] = variance_record[position];

						 }

			if(NumericVector::is_na(N)){
			  N = 10;
			  Rcout << "Use default N = 10" << endl;
			}
			else{
			  Rcout << "In this case, N = " << N << endl;
			}

			if(NumericVector::is_na(delta_upper)){
			  delta_upper = 2;
			  Rcout << "Use default L = 2" << endl;
			}
			else{
			  Rcout << "In this case, L = " << delta_upper << endl;
			}

			double ndelta = delta_upper*N;

			Rcout << "The number of delta = " << ndelta << endl;

			// Define a sequence of delta from 0 to delta_upper
			NumericVector delta(ndelta);
			for(double m = 0; m < ndelta; m++){
				    delta[m] = (m+1)/(ndelta/delta_upper);
		    }

			delta = sortDescending_SqrtPGM_EBIC(delta);

			NumericVector error_cum(ndelta);
			List intercept_value(ndelta);
			List theta_init(ndelta);
			List theta_correction(ndelta);
			List Statistic(ndelta);
			List upper_stat_tmp(ndelta);
			List CI_low(ndelta);
			List CI_high(ndelta);
			List p_value(ndelta);
			List standard_deviation(ndelta);

			Rcout <<"Perform global inference." << endl;
			Rcout << "Use pre-specified level(s): " << rcpp_level << endl;

			//Take out the upper triangular of true graph if available
			NumericMatrix rcpp_true_graph(p,p);
			NumericVector omega_true(p*(p-1)/2);
			size_t t = 0;
			if(true_graph.isNull()){
					Rcout << "True graph is not available." << endl;
			}
			else{
					Rcout << "True graph is available." << endl;
					rcpp_true_graph = as<NumericMatrix>(true_graph);
					for (size_t i = 1; i < p; i++){
							for (size_t j = 0; j < i; j++){
									omega_true[t]=rcpp_true_graph(j,i);
									t++;
							}
					}
			}

			NumericMatrix tmp_beta(p,p);
			NumericMatrix tmp_predictor(n,p);
			NumericMatrix tmp_p1(n,p);
			for(size_t i = 0; i < p; i++){
					tmp_p1(_,i) = P1_calculation_SqrtPGM_EBIC(tmp_predictor(_,i), tmp_beta(0,i));
			}
			NumericMatrix tmp_p2(n,p);
			for(size_t i = 0; i < p; i++){
					tmp_p2(_,i) = P2_calculation_SqrtPGM_EBIC(tmp_predictor(_,i), tmp_beta(0,i));
			}

			size_t count = 0;

			Rcout << "Calculate Lasso of each variable with tuning parameters under each delta." << endl;
			Rcout << "Record test statistics under each delta." << endl;
			for(size_t delta_count = 0; delta_count < ndelta; delta_count++){
					NumericMatrix theta_cor(p,p);
					NumericMatrix CI_low_theta(p,p);
					NumericMatrix CI_high_theta(p,p);
					NumericMatrix p_thetaCor(p,p);
					NumericMatrix est_stand_deviation(p,p);
					NumericMatrix test_statistic(p,p);
					NumericMatrix theta_initial(p,p);
					NumericVector phi(p);

					count++;

					Rcout <<"delta " << (count) << endl;

					for(size_t i = 0; i < p; i++){
						double lambda = delta[delta_count]*sqrt(noise[i]*(log(p)/n));

						Rcout <<"Lasso of Square Root PGM for variable " << (i+1) << endl;
						//}

						// Extract IP_YX[i,-i] and x[,-i]
						NumericVector sqrt_tmpyx(p-1);
						NumericMatrix sqrt_tmp_x(n,p-1);
						size_t t = 0;
						for (size_t k = 0; k < p; k++){
								if (k != i){
										sqrt_tmpyx[t] = sqrt_IP_YX(i,k);
										sqrt_tmp_x(_,t) = sqrt_x(_,k);
										t++;
								}
						}

						// Lasso
						NumericVector beta_column = tmp_beta(_,i);
						NumericVector predictor_column = tmp_predictor(_,i);
						NumericVector p1_column = tmp_p1(_,i);
						NumericVector p2_column = tmp_p2(_,i);

						List beta_result = FastSqrtPGMLasso_Rcpp2(beta_column, predictor_column, p1_column, p2_column, sqrt_tmpyx, sqrt_x(_,i), sqrt_tmp_x, lambda, step_size, n);

						NumericVector beta = beta_result["beta"];
						NumericVector predictor = beta_result["predictor"];
						NumericVector p1 = beta_result["p1"];
						NumericVector p2 = beta_result["p2"];

						tmp_beta(_,i) = beta;
						tmp_predictor(_,i) = predictor;
						tmp_p1(_,i) = p1;
						tmp_p2(_,i) = p2;

						//phi[i-1] = beta[1];

						for (size_t k = 0; k < p; k++){
								if(i == k){
									 theta_initial(i,k) = beta[0];
								}
								if(i < k){
									theta_initial(i,k) = beta[k];
								}
								if(i > k){
									theta_initial(i,k) = beta[k+1];
								}
						}
					}

								Rcout << "Bias correction for initial estimators. " << endl;

								 //Bias Correction on beta_initial and symmetrization
								 //size_t n2 = x_step2.nrow();
								 NumericMatrix m(n,p);
								 for (size_t i = 0; i < p; i++){
									 for (size_t j = 0; j < n; j++){
										 m(j,i) = inner_product(sqrt_x(j,_).begin(), sqrt_x(j,_).end(), theta_initial(i,_).begin(), 0.0);
									 }
								  }

								 bool display_progress = true;
								 Progress d(p*(p-1)/2, display_progress);

								 for (size_t i = 0; i < p-1; i++){
									 	 	if (Progress::check_abort() ){
									 			 return -1.0;
									 		}

									        NumericVector tau = m(_,i) - theta_initial(i,i)*sqrt_x(_,i);
										    NumericVector psi1 = P1_calculation_SqrtPGM_EBIC(tau, theta_initial(i,i));
										    NumericVector psi2 = P2_calculation_SqrtPGM_EBIC(tau, theta_initial(i,i)) - pow(psi1,2);;
										    NumericVector v1 = sqrt_x(_,i)-psi1;

								    	for(size_t j = i+1; j < p; j++){

								            d.increment();

								            NumericVector m1(n);
								            NumericVector m2(n);

								            m1 = m(_,i)-theta_initial(i,j)*sqrt_x(_,j)-theta_initial(i,i)*sqrt_x(_,i);
								            m2 = m(_,j)-theta_initial(j,i)*sqrt_x(_,i)-theta_initial(j,j)*sqrt_x(_,j);

								            NumericVector g_numerator(n);
								            NumericVector g_denominator(n);

								            double k2 = 0;
								            NumericVector g_numerator_tmp(n);
								            NumericVector g_denominator_tmp(n);
								            do{
								            	 for(size_t index = 0; index < n; index++){
								            	        g_numerator_tmp[index] = g_numerator[index];
								            	        g_denominator_tmp[index] = g_denominator[index];
								            	 }
								            	 NumericVector capital_m(n);
								            	 double k1 = 0;
								            	 NumericVector capital_m_tmp(n);
								            	 do{
								            	        for(size_t index1 = 0; index1 < n; index1++){
								            	            	 capital_m_tmp[index1] = capital_m[index1];
								            	        }
								            	        capital_m = capital_m + exp(theta_initial(i,i)*k1-log_factorial_SqrtPGM_EBIC(k1)+theta_initial(j,j)*k2-log_factorial_SqrtPGM_EBIC(k2)+theta_initial(i,j)*sqrt(k1)*sqrt(k2)+sqrt(k1)*m1+sqrt(k2)*m2);
								            	        NumericVector diff = capital_m - capital_m_tmp;
								            	        if(max(abs(diff)) <= 1e-2){
								            	            	 break;
								            	        }
								            	        k1 = k1 + 1;
								            	 }
								            	 while(k1 <= 1e4);

								            	 NumericVector inter_step = (P2_calculation_SqrtPGM_EBIC(sqrt(k2)*theta_initial(i,j)+m1,theta_initial(i,i))-pow(P1_calculation_SqrtPGM_EBIC(sqrt(k2)*theta_initial(i,j)+m1,theta_initial(i,i)),2))*capital_m;
								            	 g_numerator = g_numerator + sqrt(k2)*inter_step;
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

									        NumericVector v = sqrt_x(_,j) - g;

								            NumericVector tau1 = v*v1;
								            NumericVector tau2 = v*psi2*sqrt_x(_,j);

								            theta_cor(i,j) = theta_initial(i,j) + sum(tau1)/sum(tau2);
								            theta_cor(j,i) = theta_cor(i,j);

								            double var_new = 1/(sum(pow(v,2) * psi2));
								            double standard_new = sqrt(var_new);

								            // 95% confidence interval of de-biased beta
								            double z_95CI = 1.96;
								            CI_low_theta(i,j) = theta_cor(i,j) - z_95CI * standard_new;
								            CI_low_theta(j,i) = CI_low_theta(i,j);
								            CI_high_theta(i,j) = theta_cor(i,j) + z_95CI * standard_new;
								            CI_high_theta(j,i) = CI_high_theta(i,j);

								            double zscore = theta_cor(i,j)/standard_new;

								            //test statistic
								            test_statistic(i,j) = zscore;
								            test_statistic(j,i) = test_statistic(i,j);

								            // P-value for de-biased beta
								            p_thetaCor(i,j) = 2 * Rf_pnorm5(-abs(zscore), 0.0, 1.0, 1, 0);
								            p_thetaCor(j,i) = p_thetaCor(i,j);

								            // Estimated standard deviation
								            est_stand_deviation(i,j) = standard_new;
								            est_stand_deviation(j,i) = est_stand_deviation(i,j);
								    	}
								  }

								 intercept_value[delta_count] = phi;
								 theta_init[delta_count] = theta_initial;
								 theta_correction[delta_count] = theta_cor;
								 Statistic[delta_count] = test_statistic;
								 CI_low[delta_count] = CI_low_theta;
								 CI_high[delta_count] = CI_high_theta;
								 p_value[delta_count] = p_thetaCor;
								 standard_deviation[delta_count] = est_stand_deviation;

								 //Take out the upper triangular of test statistic
								 size_t t = 0;
								 NumericVector upper_stat(p*(p-1)/2);
								 for (size_t i = 1; i < p; i++){
									 for (size_t j = 0; j < i; j++){
										 upper_stat[t]=test_statistic(j,i);
								 		 t++;
								 	 }
								 }

								 //threshold for choosing delta
								 double error_1 = 0;
								 for (double k = 3; k < 10; k++){
								 		double normal_tail = Rf_qnorm5(1-k/20, 0.0, 1.0, 1, 0);
								 		double error = sum(abs(upper_stat)>=normal_tail)/(p*(p-1)*k/20)-1;
								 		error_1 = error_1 + pow(error,2);
								 }

								 upper_stat_tmp[delta_count] = upper_stat;
								 error_cum[delta_count] = error_1;


			}

			//Reverse the error_cum
										reverse(error_cum.begin(), error_cum.end());


										//Return the minimum position for error

										Rcout << "Choose delta for FDR control." << endl;
										double position_select = min_position_SqrtPGM_EBIC(error_cum);


										NumericVector phi_select = intercept_value[ndelta-position_select-1];
										NumericMatrix theta_initial_select = theta_init[ndelta-position_select-1];
										NumericMatrix theta_cor_select = theta_correction[ndelta-position_select-1];
										NumericMatrix test_statistic_select = Statistic[ndelta-position_select-1];
										NumericMatrix CI_low_theta_select = CI_low[ndelta-position_select-1];
										NumericMatrix CI_high_theta_select = CI_high[ndelta-position_select-1];
										NumericMatrix p_thetaCor_select = p_value[ndelta-position_select-1];
										NumericMatrix est_stand_deviation_select = standard_deviation[ndelta-position_select-1];
										NumericVector upper_stat_select = upper_stat_tmp[ndelta-position_select-1];


										// Get threshold
										double upper = 2*sqrt(log(p));
										NumericVector z(2000);
										NumericVector numerator1(2000);
										NumericVector denominator(2000);
										for(double k = 0; k < 2000; k++){
												z[k] = ((k+1)/2000) * upper;
												numerator1[k] = 2 * (1-Rf_pnorm5(z[k], 0.0, 1.0, 1, 0)) *p*(p-1)/2;
												denominator[k] = sum(abs(upper_stat_select)>=z[k]);
												if(denominator[k] < 1){
														  denominator[k] = 1;
												}
										}
										NumericVector FDP = numerator1/denominator;


										NumericVector position(alpha_size);

										for(size_t i = 0; i < alpha_size; i++){
											for (double k = 0; k < 2000 ; k++){
													if(FDP[k] <= rcpp_level[i]){      // threshold for FDP <= pre-specified level
														  	threshold[i] = z[k];
														  	position[i] = k;
														  	break;
													}
													else{
														  if(k==1999){
														  		threshold[i] = z[k];
														  		position[i] = k;
														  }
													}
											}
										}

												if(true_graph.isNull()){
														  for(size_t i = 0; i < alpha_size; i++){
														  			FDR[i] = FDP[position[i]];
														  }
														  for(size_t k = 0; k < alpha_size; k++){
														  		  NumericMatrix decision(p,p);
														  		  for(size_t i = 0; i < p-1; i++){
														  			  for(size_t j = i+1; j < p; j++){
														  		   				   if(abs(test_statistic_select(i,j)) >= threshold[k]){
														  		   					   decision(i,j) = 1;
														  		   				   }
														  		   				   decision(j,i) = decision(i,j);
														  		   			   }
														  		   	  }
														  		   	  global_decision[k] = decision;
														  }
												}
												else{

														  for(size_t i = 0; i < alpha_size; i++){
														  	   	  double FDR_denominator = sum(abs(upper_stat_select)>= threshold[i]);
														  	   	  if(FDR_denominator < 1){
														  	   			FDR_denominator = 1;
														  	   	  }

														  	   	  FDR[i] = sum((abs(upper_stat_select) >= threshold[i]) & (omega_true == 0)) / FDR_denominator;

														  	   	  double Power_denominator = sum(omega_true!=0);
														  	   	  if(Power_denominator < 1){
														  	   			Power_denominator = 1;
														  	   	  }
														  	   	  Power[i] = sum((abs(upper_stat_select) >= threshold[i]) & (omega_true != 0)) / Power_denominator;
														  }
														  for(size_t k = 0; k < alpha_size; k++){
														  		  NumericMatrix decision(p,p);
														  		  for(size_t i = 0; i < p-1; i++){
														  		   		for(size_t j = i+1; j < p; j++){
														  		   				if(abs(test_statistic_select(i,j)) >= threshold[k]){
														  		   					   decision(i,j) = 1;
														  		   				}
														  		   				decision(j,i) = decision(i,j);
														  		   		}
														  		  }
														  		  global_decision[k] = decision;
														  }
												 }



												if(true_graph.isNull()){
														 return List::create(_["intercept"] = phi_select, _["theta_initial"] = theta_initial_select, _["theta_cor"] = theta_cor_select, _["CI_low_theta"] = CI_low_theta_select, _["CI_high_theta"] = CI_high_theta_select, _["z_score"] = test_statistic_select, _["p_thetaCor"] = p_thetaCor_select, _["est_sd"] = est_stand_deviation_select, _["threshold"] = threshold, _["FDR"] = FDR, _["global_decision"] = global_decision);
												}
												else{
														 return List::create(_["intercept"] = phi_select, _["theta_initial"] = theta_initial_select, _["theta_cor"] = theta_cor_select, _["CI_low_theta"] = CI_low_theta_select, _["CI_high_theta"] = CI_high_theta_select, _["z_score"] = test_statistic_select, _["p_thetaCor"] = p_thetaCor_select, _["est_sd"] = est_stand_deviation_select, _["threshold"] = threshold, _["FDR"] = FDR, _["power"] = Power, _["global_decision"] = global_decision);
												}

		}

	}
}
else{
	NumericVector rcpp_regularization = NumericVector::create();
	rcpp_regularization = as<NumericVector>(regularization);
	size_t regularization_size = rcpp_regularization.size();
	rcpp_regularization = sortDescending_SqrtPGM_EBIC(rcpp_regularization);

	List intercept_value(regularization_size);
	List theta_init(regularization_size);
	List theta_correction(regularization_size);
	List Statistic(regularization_size);
	List CI_low(regularization_size);
    List CI_high(regularization_size);
	List p_value(regularization_size);
	List standard_deviation(regularization_size);

	if(intercept == true){
		NumericMatrix sqrt_x(n,p);
		for(size_t i = 0; i < p; i++){
			sqrt_x(_,i) = sqrt(x(_,i));
		}

		// Include intercept in the square root matrix
		NumericMatrix sqrt_x_new(n,p+1);
		NumericVector intercept(n,1);
		for(size_t i = 0; i < (p+1); i++){
			if(i == 0){
				sqrt_x_new(_,i) = intercept;
			}
		    else{
				sqrt_x_new(_,i) = sqrt_x(_,i-1);
			}
		}

		NumericMatrix sqrt_IP_YX_new(p+1,p+1);
		for(size_t i = 0; i < p+1; i++){
				//NumericVector response = x_step1(_,i)+1;
			for(size_t j = 0; j < p+1; j++){
				sqrt_IP_YX_new(i,j) = inner_product(sqrt_x_new(_,i).begin(), sqrt_x_new(_,i).end(), sqrt_x_new(_,j).begin(), 0.0);
			}
		}

		NumericMatrix tmp_beta(p+1,p);

		NumericMatrix tmp_predictor(n,p);
		NumericMatrix tmp_p1(n,p);
		for(size_t i = 0; i < p; i++){
			  tmp_p1(_,i) = P1_calculation_SqrtPGM_EBIC(tmp_predictor(_,i), tmp_beta(0,i));
		}
		NumericMatrix tmp_p2(n,p);
		for(size_t i = 0; i < p; i++){
			  tmp_p2(_,i) = P2_calculation_SqrtPGM_EBIC(tmp_predictor(_,i), tmp_beta(0,i));
		}

		size_t count = 0;
		Rcout << "Calculate Lasso of each variable with tuning parameters under each pre-specified regularization." << endl;
		Rcout << "Record initial estimators under each regularization." << endl;

		for(size_t delta_count = 0; delta_count < regularization_size; delta_count++){
					   NumericMatrix theta_cor(p,p);
					   NumericMatrix CI_low_theta(p,p);
					   NumericMatrix CI_high_theta(p,p);
					   NumericMatrix p_thetaCor(p,p);
					   NumericMatrix est_stand_deviation(p,p);
					   NumericMatrix test_statistic(p,p);
					   NumericMatrix theta_initial(p,p);
					   NumericVector phi(p);

					   count++;

					   Rcout <<"regularization " << (count) << endl;
					   for(size_t i = 1; i < p+1; i++){
							double lambda = rcpp_regularization[delta_count];

							Rcout <<"Lasso of Square Root PGM for variable " << i << endl;

							// Extract IP_YX[i,-i] and x[,-i]
							NumericVector sqrt_tmpyx(p);
							NumericMatrix sqrt_tmp_x(n,p);
							size_t t = 0;
							for (size_t k = 0; k < p+1; k++){
								if (k != i){
										sqrt_tmpyx[t] = sqrt_IP_YX_new(i,k);
										sqrt_tmp_x(_,t) = sqrt_x_new(_,k);
										t++;
								}
							}

							// Lasso
							NumericVector beta_column = tmp_beta(_,i-1);
							NumericVector predictor_column = tmp_predictor(_,i-1);
							NumericVector p1_column = tmp_p1(_,i-1);
							NumericVector p2_column = tmp_p2(_,i-1);

						  	List beta_result = FastSqrtPGMLasso_Rcpp(beta_column, predictor_column, p1_column, p2_column, sqrt_tmpyx, sqrt_x_new(_,i), sqrt_tmp_x, lambda, step_size, n);

						  	NumericVector beta = beta_result["beta"];
						  	NumericVector predictor = beta_result["predictor"];
						  	NumericVector p1 = beta_result["p1"];
						  	NumericVector p2 = beta_result["p2"];

						  	tmp_beta(_,i-1) = beta;
						  	tmp_predictor(_,i-1) = predictor;
						  	tmp_p1(_,i-1) = p1;
						  	tmp_p2(_,i-1) = p2;

						  	phi[i-1] = beta[1];

						  	for (size_t k = 0; k < p; k++){
						  		if(i-1 == k){
						  				theta_initial(i-1,k) = beta[0];
						  		}
						  		if(i-1 < k){
						  				theta_initial(i-1,k) = beta[k+1];
						  		}
						  		if(i-1 > k){
						  				theta_initial(i-1,k) = beta[k+2];
						  		}
						    }

					   }
					   intercept_value[delta_count] = phi;
					   theta_init[delta_count] = theta_initial;
					   theta_correction[delta_count] = theta_cor;
					   Statistic[delta_count] = test_statistic;
					   CI_low[delta_count] = CI_low_theta;
					   CI_high[delta_count] = CI_high_theta;
					   p_value[delta_count] = p_thetaCor;
					   standard_deviation[delta_count] = est_stand_deviation;
			}
		 	return List::create(_["intercept"] = intercept_value, _["theta_initial"] = theta_init);

	}
	else{
		// Take square root of each element in the matrix
		NumericMatrix sqrt_x(n,p);
		for(size_t i = 0; i < p; i++){
			sqrt_x(_,i) = sqrt(x(_,i));
		}

	    NumericMatrix sqrt_IP_YX(p,p);
		for(size_t i = 0; i < p; i++){
			//NumericVector response = x_step1(_,i)+1;
			for(size_t j = 0; j < p; j++){
				sqrt_IP_YX(i,j) = inner_product(sqrt_x(_,i).begin(), sqrt_x(_,i).end(), sqrt_x(_,j).begin(), 0.0);
			}
		}

		NumericMatrix tmp_beta(p,p);
		NumericMatrix tmp_predictor(n,p);
		NumericMatrix tmp_p1(n,p);
		for(size_t i = 0; i < p; i++){
				tmp_p1(_,i) = P1_calculation_SqrtPGM_EBIC(tmp_predictor(_,i), tmp_beta(0,i));
		}
		NumericMatrix tmp_p2(n,p);
		for(size_t i = 0; i < p; i++){
				tmp_p2(_,i) = P2_calculation_SqrtPGM_EBIC(tmp_predictor(_,i), tmp_beta(0,i));
		}

		size_t count = 0;

		Rcout << "Calculate Lasso of each variable with tuning parameters under each pre-specified regularization." << endl;
		Rcout << "Record initial estimators under each regularization." << endl;
		for(size_t delta_count = 0; delta_count < regularization_size; delta_count++){
				NumericMatrix theta_cor(p,p);
				NumericMatrix CI_low_theta(p,p);
				NumericMatrix CI_high_theta(p,p);
				NumericMatrix p_thetaCor(p,p);
				NumericMatrix est_stand_deviation(p,p);
				NumericMatrix test_statistic(p,p);
				NumericMatrix theta_initial(p,p);
				NumericVector phi(p);

				count++;

				Rcout <<"regularization " << (count) << endl;

				for(size_t i = 0; i < p; i++){
						double lambda = rcpp_regularization[delta_count];

					    Rcout <<"Lasso of Square Root PGM for variable " << (i+1) << endl;
								//}

						// Extract IP_YX[i,-i] and x[,-i]
						NumericVector sqrt_tmpyx(p-1);
					    NumericMatrix sqrt_tmp_x(n,p-1);
					    size_t t = 0;
						for (size_t k = 0; k < p; k++){
						    if (k != i){
								sqrt_tmpyx[t] = sqrt_IP_YX(i,k);
								sqrt_tmp_x(_,t) = sqrt_x(_,k);
								t++;
							}
						}

						// Lasso
						NumericVector beta_column = tmp_beta(_,i);
						NumericVector predictor_column = tmp_predictor(_,i);
						NumericVector p1_column = tmp_p1(_,i);
						NumericVector p2_column = tmp_p2(_,i);

						List beta_result = FastSqrtPGMLasso_Rcpp2(beta_column, predictor_column, p1_column, p2_column, sqrt_tmpyx, sqrt_x(_,i), sqrt_tmp_x, lambda, step_size, n);

						NumericVector beta = beta_result["beta"];
						NumericVector predictor = beta_result["predictor"];
						NumericVector p1 = beta_result["p1"];
						NumericVector p2 = beta_result["p2"];

						tmp_beta(_,i) = beta;
						tmp_predictor(_,i) = predictor;
						tmp_p1(_,i) = p1;
						tmp_p2(_,i) = p2;

						//phi[i-1] = beta[1];

						for (size_t k = 0; k < p; k++){
								if(i == k){
									theta_initial(i,k) = beta[0];
								}
								if(i < k){
									theta_initial(i,k) = beta[k];
								}
								if(i > k){
									theta_initial(i,k) = beta[k+1];
								}
						}
					}
					intercept_value[delta_count] = phi;
					theta_init[delta_count] = theta_initial;
					theta_correction[delta_count] = theta_cor;
					Statistic[delta_count] = test_statistic;
					CI_low[delta_count] = CI_low_theta;
					CI_high[delta_count] = CI_high_theta;
					p_value[delta_count] = p_thetaCor;
					standard_deviation[delta_count] = est_stand_deviation;
		}
		return List::create(_["intercept"] = intercept_value, _["theta_initial"] = theta_init);


	}

}

	return R_NilValue;
}







