## Wrapper function

ModPGM_true_sd <- function(x, psi = NULL, model = "SqrtPGM", true_graph, D = NULL, D_0 = NULL, D_1 = NULL){

  p <- dim(x)[2]

  if(is.null(psi)){
    psi <- rep(0,p)
  }

  if(model == "TPGM"){
    if(is.null(D)){
      D <- apply(x,2,max)
    }
    true_sd <- true_sd_TPGM(x = x, psi = psi, theta = true_graph, D = D)
  }

  if(model == "SPGM"){
    if(is.null(D_0)){
      D_0 <- rep(0,p)
    }
    if(is.null(D_1)){
      D_1 <- apply(x,2,max)
    }
    true_sd <- true_sd_SPGM(x = x, psi = psi, theta = true_graph, D_0 = D_0, D_1 = D_1)
  }

  if(model == "SqrtPGM"){
    true_sd <- true_sd_SqrtPGM(x = x, psi = psi, theta = true_graph)
  }

  return(true_sd)
}
