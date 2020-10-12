## Wrapper function

ModPGMSampler <- function(psi = NULL, true_graph, model = "SqrtPGM", D = NULL, D_0 = NULL, D_1 = NULL, nSample, burn_in = NULL, thin = NULL){
  
  p <- dim(true_graph)[2]
  
  if(is.null(psi)){
    psi <- rep(0,p)
  }
  
  if(model == "TPGM"){
    cat("Generate samples for TPGM \n")
    if(is.null(D)){
      D <- rep(3,p) 
    }
    X <- Gibbs_TPGM(psi = psi, theta = true_graph, D = D, nSample = nSample, burn_in = burn_in, thin = thin)
  }
  
  if(model == "SPGM"){
    cat("Generate samples for SPGM \n")
    if(is.null(D_0)){
      D_0 <- rep(3,p)
    }
    if(is.null(D_1)){
      D_1 <- rep(6,p)
    }
    X <- Gibbs_SPGM(psi = psi, theta = true_graph, D_0 = D_0, D_1 = D_1, nSample = nSample, burn_in = burn_in, thin = thin)
  }
  
  if(model == "SqrtPGM"){
    cat("Generate samples for SqrtPGM \n")
    X <- Gibbs_SqrtPGM(psi = psi, theta = true_graph, nSample = nSample, burn_in = burn_in, thin = thin)
  }
  
  return(X)
  
}