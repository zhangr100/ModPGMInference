## Wrapper function:

ModPGMInference <- function(x, model = "SqrtPGM", D = NULL, D_0 = NULL, D_1 = NULL, tuning = "EBIC", gamma = NA_real_, kfold = NA_real_, nlambda = NA_real_, step_size = NA_real_, intercept = TRUE, global = FALSE, alpha = NULL, regularization = NULL, N = NA_real_, delta_upper = NA_real_, true_graph = NULL){
  
  dimension <- dim(x)[2] 
  
  if(model == "TPGM"){
    if(is.null(D)){
        D <- apply(x,2,max)  
    }
    if(tuning == "EBIC"){
      result_ModPGM <- FastTPGM_EBIC_FDR(x = x, D = D, gamma = gamma, nlambda = nlambda, step_size = step_size, intercept = intercept, global = global, alpha = alpha, regularization = regularization, N = N, delta_upper = delta_upper, true_graph = true_graph)
    }
    if(tuning == "CV"){
      result_ModPGM <- FastTPGM_CV_FDR(x = x, D = D, kfold = kfold, nlambda = nlambda, step_size = step_size, intercept = intercept, global = global, alpha = alpha, regularization = regularization, N = N, delta_upper = delta_upper, true_graph = true_graph)
    }
  }
  
  if(model == "SPGM"){
    if(is.null(D_0)){
        D_0 <- rep(0,dimension)
    }
    if(is.null(D_1)){
        D_1 <- apply(x,2,max)
    }
    if(tuning == "EBIC"){
      result_ModPGM <- FastSPGM_EBIC_FDR(x = x, D_0 = D_0, D_1 = D_1, gamma = gamma, nlambda = nlambda, step_size = step_size, intercept = intercept, global = global, regularization = regularization, alpha = alpha, N = N, delta_upper = delta_upper, true_graph = true_graph)
    }
    if(tuning == "CV"){
      result_ModPGM <- FastSPGM_CV_FDR(x = x, D_0 = D_0, D_1 = D_1, kfold = kfold, nlambda = nlambda, step_size = step_size, intercept = intercept, global = global, regularization = regularization, alpha = alpha, N = N, delta_upper = delta_upper, true_graph = true_graph)
    }
  }
  
  if(model == "SqrtPGM"){
    if(tuning == "EBIC"){
      result_ModPGM <- FastSqrtPGM_EBIC_FDR(x = x, gamma = gamma, nlambda = nlambda, step_size = step_size, intercept = intercept, global = global, regularization = regularization, alpha = alpha, N = N, delta_upper = delta_upper, true_graph = true_graph)
    }
    if(tuning == "CV"){
      result_ModPGM <- FastSqrtPGM_CV_FDR(x = x, kfold = kfold, nlambda = nlambda, step_size = step_size, intercept = intercept, global = global, regularization = regularization, alpha = alpha, N = N, delta_upper = delta_upper, true_graph = true_graph)
    }
  }
  
  return(result_ModPGM)
  
}