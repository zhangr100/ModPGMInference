library(ModPGMInference)

######## Graph settings: Chain, Grid, E-R, Scale-free ########

## Chain graph
p <- 200 ## p <- 400

## TPGM and SPGM
Omega.tmp=matrix(0,p,p)
block1 <- matrix(0,p/2,p/2)
block2 <- matrix(0,p/2,p/2)
for(i in 1:(p/2-1)){
  j=i+1
  block1[i,j]=sample(c(-0.3,0.3),1)
  block1[j,i] <- block1[i,j]
}

for(i in 1:(p/2-1)){
  j=i+1
  block2[i,j]=sample(c(-0.4,0.4),1)
  block2[j,i]=block2[i,j]
}

Omega.tmp[1:(p/2),1:(p/2)] <- block1
Omega.tmp[(p/2+1):p,(p/2+1):p] <- block2
Omega=Omega.tmp

## SqrtPGM
#Omega.tmp=matrix(0,p,p)
#block1 <- matrix(0,p/2,p/2)
#block2 <- matrix(0,p/2,p/2)
#for(i in 1:(p/2-1)){
#  j=i+1
#  block1[i,j]=sample(c(-0.6,0.6),1)
#  block1[j,i] <- block1[i,j]
#}

#for(i in 1:(p/2-1)){
#  j=i+1
#  block2[i,j]=sample(c(-0.9,0.9),1)
#  block2[j,i]=block2[i,j]
#}

#Omega.tmp[1:(p/2),1:(p/2)] <- block1
#Omega.tmp[(p/2+1):p,(p/2+1):p] <- block2
#Omega=Omega.tmp

###########################################

## Grid graph
#p <- 200 ## p <- 400

## TPGM and SPGM
#Omega.tmp=matrix(0,p,p)
#block1 <- matrix(0,p/2,p/2)
#block2 <- matrix(0,p/2,p/2)
#Grid_Graph <- function(x, y, set){
#  p <- x*y
#  M <- matrix(0, p, p)
  
#  for(i in 1:p){
#    k1=sample(set,1)
#    k2=sample(set,1)
#    if(i+1 < p+1) {M[i,i+1 ] <- k1}
#    if(i+y < p+1) {M[i,i+y ] <- k2}
#  }
#  for(i in 1:(p-1)){
#    for(j in (i+1):p){
#      M[j,i]=M[i,j]
#    }
#  }
#  return(M)
#}
# Input:
#x <- 10 # Number of rows     ## x <- 10, y <- 20 for p <- 400
#y <- 10 # Number of columns
#set <- c(-0.3,0.3)
#block1 <- Grid_Graph(x,y,set)

#set <- c(-0.4,0.4)
#block2 <- Grid_Graph(x,y,set)

#Omega.tmp[1:(p/2),1:(p/2)] <- block1
#Omega.tmp[(p/2+1):p,(p/2+1):p] <- block2
#Omega=Omega.tmp

## SqrtPGM
#Omega.tmp=matrix(0,p,p)
#block1 <- matrix(0,p/2,p/2)
#block2 <- matrix(0,p/2,p/2)
#Grid_Graph <- function(x, y, set){
#  p <- x*y
#  M <- matrix(0, p, p)
  
#  for(i in 1:p){
#    k1=sample(set,1)
#    k2=sample(set,1)
#    if(i+1 < p+1) {M[i,i+1 ] <- k1}
#    if(i+y < p+1) {M[i,i+y ] <- k2}
#  }
#  for(i in 1:(p-1)){
#    for(j in (i+1):p){
#      M[j,i]=M[i,j]
#    }
#  }
#  return(M)
#}
# Input:
#x <- 10 # Number of rows     ## x <- 10, y <- 20 for p <- 400
#y <- 10 # Number of columns

#set <- c(-0.6,0.6)
#block1 <- Grid_Graph(x,y,set)
#set <- c(-0.9,0.9)
#block2 <- Grid_Graph(x,y,set)

#Omega.tmp[1:(p/2),1:(p/2)] <- block1
#Omega.tmp[(p/2+1):p,(p/2+1):p] <- block2
#Omega=Omega.tmp

###########################################

## E-R graph
#p <- 200 ## p <- 400

## TPGM and SPGM
#Omega.tmp=matrix(0,p,p)

#block1 <- matrix(0,p/2,p/2)
#block2 <- matrix(0,p/2,p/2)

#prop = 4/(p/2) # sparseness
#for(i in 1:(p/2-1)){
#  for(j in (i+1):(p/2)){
#    block1[i, j]=sample(c(-0.3,0.3),1)*rbinom(1,1,prop)
#    block1[j, i]=block1[i, j]
#  }
#}

#for(i in 1:(p/2-1)){
#  for(j in (i+1):(p/2)){
#    block2[i, j]=sample(c(-0.4,0.4),1)*rbinom(1,1,prop)
#    block2[j, i]=block2[i, j]
#  }
#}

#Omega.tmp[1:(p/2),1:(p/2)] <- block1
#Omega.tmp[(p/2+1):p,(p/2+1):p] <- block2
#Omega=Omega.tmp

## SqrtPGM
#Omega.tmp=matrix(0,p,p)

#block1 <- matrix(0,p/2,p/2)
#block2 <- matrix(0,p/2,p/2)

#prop = 4/(p/2) # sparseness
#for(i in 1:(p/2-1)){
#  for(j in (i+1):(p/2)){
#    block1[i, j]=sample(c(-0.6,0.6),1)*rbinom(1,1,prop)
#    block1[j, i]=block1[i, j]
#  }
#}

#for(i in 1:(p/2-1)){
#  for(j in (i+1):(p/2)){
#    block2[i, j]=sample(c(-0.9,0.9),1)*rbinom(1,1,prop)
#    block2[j, i]=block2[i, j]
#  }
#}

#Omega.tmp[1:(p/2),1:(p/2)] <- block1
#Omega.tmp[(p/2+1):p,(p/2+1):p] <- block2
#Omega=Omega.tmp

###########################################

## Scale-free graph
#p <- 200 ## p <- 400

## TPGM and SPGM
#Omega.tmp=matrix(0,p,p)

#library(igraph)
#g <- barabasi.game(n=p/2, power=0.01, directed=F)
#block1 <- as.matrix(get.adjacency(g, type="both"))
#for(i in 1:(p/2-1)){
#  for(j in (i+1):(p/2)){
#    if(block1[i,j]!=0){
#      block1[i,j] = sample(c(-0.3,0.3),1)
#      block1[j,i] = block1[i,j]
#    }
#  }
#}
#g <- barabasi.game(n=p/2, power=0.01, directed=F)
#block2 <- as.matrix(get.adjacency(g, type="both"))
#for(i in 1:(p/2-1)){
#  for(j in (i+1):(p/2)){
#    if(block2[i,j]!=0){
#      block2[i,j] = sample(c(-0.4,0.4),1)
#      block2[j,i] = block2[i,j]
#    }
#  }
#}

#Omega.tmp[1:(p/2),1:(p/2)] <- block1
#Omega.tmp[(p/2+1):p,(p/2+1):p] <- block2
#Omega=Omega.tmp

## SqrtPGM
#Omega.tmp=matrix(0,p,p)

#library(igraph)
#g <- barabasi.game(n=p/2, power=0.01, directed=F)
#block1 <- as.matrix(get.adjacency(g, type="both"))
#for(i in 1:(p/2-1)){
#  for(j in (i+1):(p/2)){
#    if(block1[i,j]!=0){
#      block1[i,j] = sample(c(-0.6,0.6),1)
#      block1[j,i] = block1[i,j]
#    }
#  }
#}
#g <- barabasi.game(n=p/2, power=0.01, directed=F)
#block2 <- as.matrix(get.adjacency(g, type="both"))
#for(i in 1:(p/2-1)){
#  for(j in (i+1):(p/2)){
#    if(block2[i,j]!=0){
#      block2[i,j] = sample(c(-0.9,0.9),1)
#      block2[j,i] = block2[i,j]
#    }
#  }
#}

#Omega.tmp[1:(p/2),1:(p/2)] <- block1
#Omega.tmp[(p/2+1):p,(p/2+1):p] <- block2
#Omega=Omega.tmp


######## Generate random samples for multiple testing ######## 

n <- 400 ## n <- 150
p <- 200 ## p <- 400
iteration <- 100

## TPGM 
X <- list()
for(i in 1:iteration){
  X[[i]] <- ModPGMSampler(psi = rep(-0.5,p), true_graph = Omega, model = "TPGM", D = rep(3,p), nSample = n, burn_in = 5000)
}

## SPGM for Chain and Scale-free graphs
X <- list()
for(i in 1:iteration){
  X[[i]] <- ModPGMSampler(psi = rep(-0.5,p), true_graph = Omega, model = "SPGM", D_0 = rep(2,p), D_1 = rep(5,p), nSample = n, burn_in = 5000)
}

## SPGM for Grid and E-R graphs
X <- list()
for(i in 1:iteration){
  X[[i]] <- ModPGMSampler(psi = rep(-1,p), true_graph = Omega, model = "SPGM", D_0 = rep(2,p), D_1 = rep(5,p), nSample = n, burn_in = 5000)
}

## SqrtPGM
X <- list()
for(i in 1:iteration){
  X[[i]] <- ModPGMSampler(psi = rep(0,p), true_graph = Omega, model = "SqrtPGM", nSample = n, burn_in = 5000)
}


######## Perform multiple testing on random samples ######## 

p <- 200 ## p <- 400
iteration <- 100

## Tuning parameter selection: EBIC
## TPGM
FDR_record_TPGM_EBIC <- list()
power_record_TPGM_EBIC <- list()
FPR_record_TPGM_EBIC <- list()
TPR_record_TPGM_EBIC <- list()

for(i in 1:iteration){
  fit <- ModPGMInference(x = X[[i]], model = "TPGM", tuning = "EBIC", D = rep(3,p), nlambda = 100, true_graph = Omega, global = TRUE, alpha = seq(0.001,1,by = 0.001))
  FPR <- NULL
  TPR <- NULL
  FDR_record_TPGM_EBIC[[i]] <- fit$FDR
  power_record_TPGM_EBIC[[i]] <- fit$power
  
  for(j in 1:length(fit$threshold)){
    
    est <- fit$global_decision[[j]]
    est <- est[upper.tri(est)]
    
    result <- est
    true <- Omega[upper.tri(Omega)]
    
    FP <- sum(result!=0 & true==0)
    TN <- sum(result==0 & true==0)
    FPR[j] <- FP/(FP+TN)
    TP <- sum(result!=0 & true!=0)
    FN <- sum(result==0 & true!=0)
    TPR[j] <- TP/(TP+FN)
  }
  
  FPR_record_TPGM_EBIC[[i]] <- FPR
  TPR_record_TPGM_EBIC[[i]] <- TPR
  
  print(i)
}

## SPGM
#FDR_record_SPGM_EBIC <- list()
#power_record_SPGM_EBIC <- list()
#FPR_record_SPGM_EBIC <- list()
#TPR_record_SPGM_EBIC <- list()

#for(i in 1:iteration){
#  fit <- ModPGMInference(x = X[[i]], model = "SPGM", tuning = "EBIC", D_0 = rep(2,p), D_1 = rep(5,p), nlambda = 100, true_graph = Omega, global = TRUE, alpha = seq(0.001,1,by = 0.001))
#  FPR <- NULL
#  TPR <- NULL
#  FDR_record_SPGM_EBIC[[i]] <- fit$FDR
#  power_record_SPGM_EBIC[[i]] <- fit$power
  
#  for(j in 1:length(fit$threshold)){
    
#    est <- fit$global_decision[[j]]
#    est <- est[upper.tri(est)]
    
#    result <- est
#    true <- Omega[upper.tri(Omega)]
    
#    FP <- sum(result!=0 & true==0)
#    TN <- sum(result==0 & true==0)
#    FPR[j] <- FP/(FP+TN)
#    TP <- sum(result!=0 & true!=0)
#    FN <- sum(result==0 & true!=0)
#    TPR[j] <- TP/(TP+FN)
#  }
  
#  FPR_record_SPGM_EBIC[[i]] <- FPR
#  TPR_record_SPGM_EBIC[[i]] <- TPR
  
#  print(i)
#}

## SqrtPGM
#FDR_record_SqrtPGM_EBIC <- list()
#power_record_SqrtPGM_EBIC <- list()
#FPR_record_SqrtPGM_EBIC <- list()
#TPR_record_SqrtPGM_EBIC <- list()

#for(i in 1:iteration){
#  fit <- ModPGMInference(x = X[[i]], model = "SqrtPGM", tuning = "EBIC", nlambda = 100, true_graph = Omega, global = TRUE, alpha = seq(0.001,1,by = 0.001))
#  FPR <- NULL
#  TPR <- NULL
#  FDR_record_SqrtPGM_EBIC[[i]] <- fit$FDR
#  power_record_SqrtPGM_EBIC[[i]] <- fit$power
  
#  for(j in 1:length(fit$threshold)){
    
#    est <- fit$global_decision[[j]]
#    est <- est[upper.tri(est)]
    
#    result <- est
#    true <- Omega[upper.tri(Omega)]
    
#    FP <- sum(result!=0 & true==0)
#    TN <- sum(result==0 & true==0)
#    FPR[j] <- FP/(FP+TN)
#    TP <- sum(result!=0 & true!=0)
#    FN <- sum(result==0 & true!=0)
#    TPR[j] <- TP/(TP+FN)
#  }
  
#  FPR_record_SqrtPGM_EBIC[[i]] <- FPR
#  TPR_record_SqrtPGM_EBIC[[i]] <- TPR
  
#  print(i)
#}

## Tuning parameter selection: Cross validation (CV)
## TPGM 
FDR_record_TPGM_CV <- list()
power_record_TPGM_CV <- list()
FPR_record_TPGM_CV <- list()
TPR_record_TPGM_CV <- list()

for(i in 1:iteration){
  fit <- ModPGMInference(x = X[[i]], model = "TPGM", tuning = "CV", kfold = 10, D = rep(3,p), nlambda = 100, true_graph = Omega, global = TRUE, alpha = seq(0.001,1,by = 0.001))
  FPR <- NULL
  TPR <- NULL
  FDR_record_TPGM_CV[[i]] <- fit$FDR
  power_record_TPGM_CV[[i]] <- fit$power
  
  for(j in 1:length(fit$threshold)){
    
    est <- fit$global_decision[[j]]
    est <- est[upper.tri(est)]
    
    result <- est
    true <- Omega[upper.tri(Omega)]
    
    FP <- sum(result!=0 & true==0)
    TN <- sum(result==0 & true==0)
    FPR[j] <- FP/(FP+TN)
    TP <- sum(result!=0 & true!=0)
    FN <- sum(result==0 & true!=0)
    TPR[j] <- TP/(TP+FN)
  }
  
  FPR_record_TPGM_CV[[i]] <- FPR
  TPR_record_TPGM_CV[[i]] <- TPR
  
  print(i)
}

## SPGM
#FDR_record_SPGM_CV <- list()
#power_record_SPGM_CV <- list()
#FPR_record_SPGM_CV <- list()
#TPR_record_SPGM_CV <- list()

#for(i in 1:iteration){
#  fit <- ModPGMInference(x = X[[i]], model = "SPGM", tuning = "CV", kfold = 10, D_0 = rep(2,p), D_1 = rep(5,p), nlambda = 100, true_graph = Omega, global = TRUE, alpha = seq(0.001,1,by = 0.001))
#  FPR <- NULL
#  TPR <- NULL
#  FDR_record_SPGM_CV[[i]] <- fit$FDR
#  power_record_SPGM_CV[[i]] <- fit$power

#  for(j in 1:length(fit$threshold)){

#    est <- fit$global_decision[[j]]
#    est <- est[upper.tri(est)]

#    result <- est
#    true <- Omega[upper.tri(Omega)]

#    FP <- sum(result!=0 & true==0)
#    TN <- sum(result==0 & true==0)
#    FPR[j] <- FP/(FP+TN)
#    TP <- sum(result!=0 & true!=0)
#    FN <- sum(result==0 & true!=0)
#    TPR[j] <- TP/(TP+FN)
#  }

#  FPR_record_SPGM_CV[[i]] <- FPR
#  TPR_record_SPGM_CV[[i]] <- TPR

#  print(i)
#}

## SqrtPGM
#FDR_record_SqrtPGM_CV <- list()
#power_record_SqrtPGM_CV <- list()
#FPR_record_SqrtPGM_CV <- list()
#TPR_record_SqrtPGM_CV <- list()

#for(i in 1:iteration){
#  fit <- ModPGMInference(x = X[[i]], model = "SqrtPGM", tuning = "CV", kfold = 10, nlambda = 100, true_graph = Omega, global = TRUE, alpha = seq(0.001,1,by = 0.001))
#  FPR <- NULL
#  TPR <- NULL
#  FDR_record_SqrtPGM_CV[[i]] <- fit$FDR
#  power_record_SqrtPGM_CV[[i]] <- fit$power

#  for(j in 1:length(fit$threshold)){

#    est <- fit$global_decision[[j]]
#    est <- est[upper.tri(est)]

#    result <- est
#    true <- Omega[upper.tri(Omega)]

#    FP <- sum(result!=0 & true==0)
#    TN <- sum(result==0 & true==0)
#    FPR[j] <- FP/(FP+TN)
#    TP <- sum(result!=0 & true!=0)
#    FN <- sum(result==0 & true!=0)
#    TPR[j] <- TP/(TP+FN)
#  }

#  FPR_record_SqrtPGM_CV[[i]] <- FPR
#  TPR_record_SqrtPGM_CV[[i]] <- TPR

#  print(i)
#}

## Note that "delta_upper" and "N" are adjustable in real implementations
## In our simulation, we use different "delta_upper" for several settings
## For SqrtPGM with Grid graph and p = 200
## For TPGM, SPGM and SqrtPGM with Grid ang E-R graphs and p = 400 
## We use "delta_upper = 4"


## Medians (SDs) of FDR values at the 0.1 and the 0.2 levels and their power values 
FDR_TPGM_EBIC_0.1 <- NULL
power_TPGM_EBIC_0.1 <- NULL
FDR_TPGM_EBIC_0.2 <- NULL
power_TPGM_EBIC_0.2 <- NULL
for(i in 1:iteration){
  FDR_TPGM_EBIC_0.1[i] <- FDR_record_TPGM_EBIC[[i]][100]
  power_TPGM_EBIC_0.1[i] <- power_record_TPGM_EBIC[[i]][100]
  FDR_TPGM_EBIC_0.2[i] <- FDR_record_TPGM_EBIC[[i]][200]
  power_TPGM_EBIC_0.2[i] <- power_record_TPGM_EBIC[[i]][200]
}

FDR_TPGM_EBIC_0.1_median <- median(FDR_TPGM_EBIC_0.1)
FDR_TPGM_EBIC_0.1_sd <- sd(FDR_TPGM_EBIC_0.1)
power_TPGM_EBIC_0.1_median <- median(power_TPGM_EBIC_0.1)
power_TPGM_EBIC_0.1_sd <- sd(power_TPGM_EBIC_0.1)
FDR_TPGM_EBIC_0.2_median <- median(FDR_TPGM_EBIC_0.2)
FDR_TPGM_EBIC_0.2_sd <- sd(FDR_TPGM_EBIC_0.2)
power_TPGM_EBIC_0.2_median <- median(power_TPGM_EBIC_0.2)
power_TPGM_EBIC_0.2_sd <- sd(power_TPGM_EBIC_0.2)

#FDR_SPGM_EBIC_0.1 <- NULL
#power_SPGM_EBIC_0.1 <- NULL
#FDR_SPGM_EBIC_0.2 <- NULL
#power_SPGM_EBIC_0.2 <- NULL
#for(i in 1:iteration){
#  FDR_SPGM_EBIC_0.1[i] <- FDR_record_SPGM_EBIC[[i]][100]
#  power_SPGM_EBIC_0.1[i] <- power_record_SPGM_EBIC[[i]][100]
#  FDR_SPGM_EBIC_0.2[i] <- FDR_record_SPGM_EBIC[[i]][200]
#  power_SPGM_EBIC_0.2[i] <- power_record_SPGM_EBIC[[i]][200]
#}

#FDR_SPGM_EBIC_0.1_median <- median(FDR_SPGM_EBIC_0.1)
#FDR_SPGM_EBIC_0.1_sd <- sd(FDR_SPGM_EBIC_0.1)
#power_SPGM_EBIC_0.1_median <- median(power_SPGM_EBIC_0.1)
#power_SPGM_EBIC_0.1_sd <- sd(power_SPGM_EBIC_0.1)
#FDR_SPGM_EBIC_0.2_median <- median(FDR_SPGM_EBIC_0.2)
#FDR_SPGM_EBIC_0.2_sd <- sd(FDR_SPGM_EBIC_0.2)
#power_SPGM_EBIC_0.2_median <- median(power_SPGM_EBIC_0.2)
#power_SPGM_EBIC_0.2_sd <- sd(power_SPGM_EBIC_0.2)

#FDR_SqrtPGM_EBIC_0.1 <- NULL
#power_SqrtPGM_EBIC_0.1 <- NULL
#FDR_SqrtPGM_EBIC_0.2 <- NULL
#power_SqrtPGM_EBIC_0.2 <- NULL
#for(i in 1:iteration){
#  FDR_SqrtPGM_EBIC_0.1[i] <- FDR_record_SqrtPGM_EBIC[[i]][100]
#  power_SqrtPGM_EBIC_0.1[i] <- power_record_SqrtPGM_EBIC[[i]][100]
#  FDR_SqrtPGM_EBIC_0.2[i] <- FDR_record_SqrtPGM_EBIC[[i]][200]
#  power_SqrtPGM_EBIC_0.2[i] <- power_record_SqrtPGM_EBIC[[i]][200]
#}

#FDR_SqrtPGM_EBIC_0.1_median <- median(FDR_SqrtPGM_EBIC_0.1)
#FDR_SqrtPGM_EBIC_0.1_sd <- sd(FDR_SqrtPGM_EBIC_0.1)
#power_SqrtPGM_EBIC_0.1_median <- median(power_SqrtPGM_EBIC_0.1)
#power_SqrtPGM_EBIC_0.1_sd <- sd(power_SqrtPGM_EBIC_0.1)
#FDR_SqrtPGM_EBIC_0.2_median <- median(FDR_SqrtPGM_EBIC_0.2)
#FDR_SqrtPGM_EBIC_0.2_sd <- sd(FDR_SqrtPGM_EBIC_0.2)
#power_SqrtPGM_EBIC_0.2_median <- median(power_SqrtPGM_EBIC_0.2)
#power_SqrtPGM_EBIC_0.2_sd <- sd(power_SqrtPGM_EBIC_0.2)


FDR_TPGM_CV_0.1 <- NULL
power_TPGM_CV_0.1 <- NULL
FDR_TPGM_CV_0.2 <- NULL
power_TPGM_CV_0.2 <- NULL
for(i in 1:iteration){
  FDR_TPGM_CV_0.1[i] <- FDR_record_TPGM_CV[[i]][100]
  power_TPGM_CV_0.1[i] <- power_record_TPGM_CV[[i]][100]
  FDR_TPGM_CV_0.2[i] <- FDR_record_TPGM_CV[[i]][200]
  power_TPGM_CV_0.2[i] <- power_record_TPGM_CV[[i]][200]
}

FDR_TPGM_CV_0.1_median <- median(FDR_TPGM_CV_0.1)
FDR_TPGM_CV_0.1_sd <- sd(FDR_TPGM_CV_0.1)
power_TPGM_CV_0.1_median <- median(power_TPGM_CV_0.1)
power_TPGM_CV_0.1_sd <- sd(power_TPGM_CV_0.1)
FDR_TPGM_CV_0.2_median <- median(FDR_TPGM_CV_0.2)
FDR_TPGM_CV_0.2_sd <- sd(FDR_TPGM_CV_0.2)
power_TPGM_CV_0.2_median <- median(power_TPGM_CV_0.2)
power_TPGM_CV_0.2_sd <- sd(power_TPGM_CV_0.2)

#FDR_SPGM_CV_0.1 <- NULL
#power_SPGM_CV_0.1 <- NULL
#FDR_SPGM_CV_0.2 <- NULL
#power_SPGM_CV_0.2 <- NULL
#for(i in 1:iteration){
#  FDR_SPGM_CV_0.1[i] <- FDR_record_SPGM_CV[[i]][100]
#  power_SPGM_CV_0.1[i] <- power_record_SPGM_CV[[i]][100]
#  FDR_SPGM_CV_0.2[i] <- FDR_record_SPGM_CV[[i]][200]
#  power_SPGM_CV_0.2[i] <- power_record_SPGM_CV[[i]][200]
#}

#FDR_SPGM_CV_0.1_median <- median(FDR_SPGM_CV_0.1)
#FDR_SPGM_CV_0.1_sd <- sd(FDR_SPGM_CV_0.1)
#power_SPGM_CV_0.1_median <- median(power_SPGM_CV_0.1)
#power_SPGM_CV_0.1_sd <- sd(power_SPGM_CV_0.1)
#FDR_SPGM_CV_0.2_median <- median(FDR_SPGM_CV_0.2)
#FDR_SPGM_CV_0.2_sd <- sd(FDR_SPGM_CV_0.2)
#power_SPGM_CV_0.2_median <- median(power_SPGM_CV_0.2)
#power_SPGM_CV_0.2_sd <- sd(power_SPGM_CV_0.2)

#FDR_SqrtPGM_CV_0.1 <- NULL
#power_SqrtPGM_CV_0.1 <- NULL
#FDR_SqrtPGM_CV_0.2 <- NULL
#power_SqrtPGM_CV_0.2 <- NULL
#for(i in 1:iteration){
#  FDR_SqrtPGM_CV_0.1[i] <- FDR_record_SqrtPGM_CV[[i]][100]
#  power_SqrtPGM_CV_0.1[i] <- power_record_SqrtPGM_CV[[i]][100]
#  FDR_SqrtPGM_CV_0.2[i] <- FDR_record_SqrtPGM_CV[[i]][200]
#  power_SqrtPGM_CV_0.2[i] <- power_record_SqrtPGM_CV[[i]][200]
#}

#FDR_SqrtPGM_CV_0.1_median <- median(FDR_SqrtPGM_CV_0.1)
#FDR_SqrtPGM_CV_0.1_sd <- sd(FDR_SqrtPGM_CV_0.1)
#power_SqrtPGM_CV_0.1_median <- median(power_SqrtPGM_CV_0.1)
#power_SqrtPGM_CV_0.1_sd <- sd(power_SqrtPGM_CV_0.1)
#FDR_SqrtPGM_CV_0.2_median <- median(FDR_SqrtPGM_CV_0.2)
#FDR_SqrtPGM_CV_0.2_sd <- sd(FDR_SqrtPGM_CV_0.2)
#power_SqrtPGM_CV_0.2_median <- median(power_SqrtPGM_CV_0.2)
#power_SqrtPGM_CV_0.2_sd <- sd(power_SqrtPGM_CV_0.2)


######## Perform sole estimation on random samples ######## 

p <- 200 ## p <- 400
iteration <- 100

## Tuning parameter selection: EBIC
## TPGM
fit_TPGM_EBIC <- list()

for(i in 1:iteration){
  fit_TPGM_EBIC[[i]] <- ModPGMInference(x = X[[i]], model = "TPGM", tuning = "EBIC", D = rep(3,p), nlambda = 100)
  print(i)
}

FDR_TPGM_EBIC_est <- NULL
power_TPGM_EBIC_est <- NULL

for(i in 1:iteration){
  est <- fit_TPGM_EBIC[[i]]$theta_initial
  
  ## symmetrize
  for(k in 1:(p-1)){
    for(j in (k+1):p){
      est[k,j] = (est[k,j]+est[j,k])/2
      est[j,k] = est[k,j]
    }
  }
  est <- est[upper.tri(est)]
  true <- Omega[upper.tri(Omega)]
  
  FDR_TPGM_EBIC_est[i] <- sum(est!=0 & true==0)/max(sum(est!=0),1)
  power_TPGM_EBIC_est[i] <- sum(est!=0 & true!=0)/sum(true!=0)
}

## SPGM
#fit_SPGM_EBIC <- list()

#for(i in 1:iteration){
#  fit_SPGM_EBIC[[i]] <- ModPGMInference(x = X[[i]], model = "SPGM", tuning = "EBIC", D_0 = rep(2,p), D_1 = rep(5,p), nlambda = 100)
#  print(i)
#}

#FDR_SPGM_EBIC_est <- NULL
#power_SPGM_EBIC_est <- NULL

#for(i in 1:iteration){
#  est <- fit_SPGM_EBIC[[i]]$theta_initial
  
  ## symmetrize
#  for(k in 1:(p-1)){
#    for(j in (k+1):p){
#      est[k,j] = (est[k,j]+est[j,k])/2
#      est[j,k] = est[k,j]
#    }
#  }
#  est <- est[upper.tri(est)]
#  true <- Omega[upper.tri(Omega)]
  
#  FDR_SPGM_EBIC_est[i] <- sum(est!=0 & true==0)/max(sum(est!=0),1)
#  power_SPGM_EBIC_est[i] <- sum(est!=0 & true!=0)/sum(true!=0)
#}

## SqrtPGM
#fit_SqrtPGM_EBIC <- list()

#for(i in 1:iteration){
#  fit_SqrtPGM_EBIC[[i]] <- ModPGMInference(x = X[[i]], model = "SqrtPGM", tuning = "EBIC", nlambda = 100)
#  print(i)
#}

#FDR_SqrtPGM_EBIC_est <- NULL
#power_SqrtPGM_EBIC_est <- NULL

#for(i in 1:iteration){
#  est <- fit_SqrtPGM_EBIC[[i]]$theta_initial
  
  ## symmetrize
#  for(k in 1:(p-1)){
#    for(j in (k+1):p){
#      est[k,j] = (est[k,j]+est[j,k])/2
#      est[j,k] = est[k,j]
#    }
#  }
#  est <- est[upper.tri(est)]
#  true <- Omega[upper.tri(Omega)]
  
#  FDR_SqrtPGM_EBIC_est[i] <- sum(est!=0 & true==0)/max(sum(est!=0),1)
#  power_SqrtPGM_EBIC_est[i] <- sum(est!=0 & true!=0)/sum(true!=0)
#}


## Tuning parameter selection: CV
## TPGM
fit_TPGM_CV <- list()

for(i in 1:iteration){
  fit_TPGM_CV[[i]] <- ModPGMInference(x = X[[i]], model = "TPGM", tuning = "CV", kfold = 10, D = rep(3,p), nlambda = 100)
  print(i)
}

FDR_TPGM_CV_est <- NULL
power_TPGM_CV_est <- NULL

for(i in 1:iteration){
  est <- fit_TPGM_CV[[i]]$theta_initial
  
  ## symmetrize
  for(k in 1:(p-1)){
    for(j in (k+1):p){
      est[k,j] = (est[k,j]+est[j,k])/2
      est[j,k] = est[k,j]
    }
  }
  est <- est[upper.tri(est)]
  true <- Omega[upper.tri(Omega)]
  
  FDR_TPGM_CV_est[i] <- sum(est!=0 & true==0)/max(sum(est!=0),1)
  power_TPGM_CV_est[i] <- sum(est!=0 & true!=0)/sum(true!=0)
}

## SPGM
#fit_SPGM_CV <- list()

#for(i in 1:iteration){
#  fit_SPGM_CV[[i]] <- ModPGMInference(x = X[[i]], model = "SPGM", tuning = "CV", kfold = 10, D_0 = rep(2,p), D_1 = rep(5,p), nlambda = 100)
#  print(i)
#}

#FDR_SPGM_CV_est <- NULL
#power_SPGM_CV_est <- NULL

#for(i in 1:iteration){
#  est <- fit_SPGM_CV[[i]]$theta_initial

## symmetrize
#  for(k in 1:(p-1)){
#    for(j in (k+1):p){
#      est[k,j] = (est[k,j]+est[j,k])/2
#      est[j,k] = est[k,j]
#    }
#  }
#  est <- est[upper.tri(est)]
#  true <- Omega[upper.tri(Omega)]

#  FDR_SPGM_CV_est[i] <- sum(est!=0 & true==0)/max(sum(est!=0),1)
#  power_SPGM_CV_est[i] <- sum(est!=0 & true!=0)/sum(true!=0)
#}

## SqrtPGM
#fit_SqrtPGM_CV <- list()

#for(i in 1:iteration){
#  fit_SqrtPGM_CV[[i]] <- ModPGMInference(x = X[[i]], model = "SqrtPGM", tuning = "CV", kfold = 10, nlambda = 100)
#  print(i)
#}

#FDR_SqrtPGM_CV_est <- NULL
#power_SqrtPGM_CV_est <- NULL

#for(i in 1:iteration){
#  est <- fit_SqrtPGM_CV[[i]]$theta_initial

## symmetrize
#  for(k in 1:(p-1)){
#    for(j in (k+1):p){
#      est[k,j] = (est[k,j]+est[j,k])/2
#      est[j,k] = est[k,j]
#    }
#  }
#  est <- est[upper.tri(est)]
#  true <- Omega[upper.tri(Omega)]

#  FDR_SqrtPGM_CV_est[i] <- sum(est!=0 & true==0)/max(sum(est!=0),1)
#  power_SqrtPGM_CV_est[i] <- sum(est!=0 & true!=0)/sum(true!=0)
#}


## Medians (SDs) of FDR values and their power values from sole estimation
FDR_TPGM_EBIC_est_median <- median(FDR_TPGM_EBIC_est)
FDR_TPGM_EBIC_est_sd <- sd(FDR_TPGM_EBIC_est)
power_TPGM_EBIC_est_median <- median(power_TPGM_EBIC_est)
power_TPGM_EBIC_est_sd <- sd(power_TPGM_EBIC_est)

#FDR_SPGM_EBIC_est_median <- median(FDR_SPGM_EBIC_est)
#FDR_SPGM_EBIC_est_sd <- sd(FDR_SPGM_EBIC_est)
#power_SPGM_EBIC_est_median <- median(power_SPGM_EBIC_est)
#power_SPGM_EBIC_est_sd <- sd(power_SPGM_EBIC_est)

#FDR_SqrtPGM_EBIC_est_median <- median(FDR_SqrtPGM_EBIC_est)
#FDR_SqrtPGM_EBIC_est_sd <- sd(FDR_SqrtPGM_EBIC_est)
#power_SqrtPGM_EBIC_est_median <- median(power_SqrtPGM_EBIC_est)
#power_SqrtPGM_EBIC_est_sd <- sd(power_SqrtPGM_EBIC_est)

FDR_TPGM_CV_est_median <- median(FDR_TPGM_CV_est)
FDR_TPGM_CV_est_sd <- sd(FDR_TPGM_CV_est)
power_TPGM_CV_est_median <- median(power_TPGM_CV_est)
power_TPGM_CV_est_sd <- sd(power_TPGM_CV_est)

#FDR_SPGM_CV_est_median <- median(FDR_SPGM_CV_est)
#FDR_SPGM_CV_est_sd <- sd(FDR_SPGM_CV_est)
#power_SPGM_CV_est_median <- median(power_SPGM_CV_est)
#power_SPGM_CV_est_sd <- sd(power_SPGM_CV_est)

#FDR_SqrtPGM_CV_est_median <- median(FDR_SqrtPGM_CV_est)
#FDR_SqrtPGM_CV_est_sd <- sd(FDR_SqrtPGM_CV_est)
#power_SqrtPGM_CV_est_median <- median(power_SqrtPGM_CV_est)
#power_SqrtPGM_CV_est_sd <- sd(power_SqrtPGM_CV_est)

######## ROC curves ######## 

p <- 200 ## p <- 400
iteration <- 100

## TPGM
FPR_record_TPGM_est <- list()
TPR_record_TPGM_est <- list()

max_lambda <- 10
min_lambda <- 0.001*max_lambda
lambda_path <- NULL
for(i in 1:1000){
  lambda_path[i] <- max_lambda/(1000^((i-1)/(1000-1))) 
}

for(i in 1:iteration){
  TPGM_fit <- ModPGMInference(x = X[[i]], model = "TPGM", D = rep(3,p), regularization = lambda_path)
  
  FPR1 <- NULL
  TPR1 <- NULL
  
  for(m in 1:1000){
    est <- TPGM_fit$theta_initial[[m]]
    
    ## symmetrize
    for(k in 1:(p-1)){
      for(j in (k+1):p){
        est[k,j] = (est[k,j]+est[j,k])/2
        est[j,k] = est[k,j]
      }
    }
    est <- est[upper.tri(est)]
    true <- Omega[upper.tri(Omega)]
    FP <- sum(est!=0 & true==0)
    TN <- sum(est==0 & true==0)
    FPR1[m] <- FP/(FP+TN)
    TP <- sum(est!=0 & true!=0)
    FN <- sum(est==0 & true!=0)
    TPR1[m] <- TP/(TP+FN)
    
  }	
  
  
  FPR_record_TPGM_est[[i]] <- FPR1
  TPR_record_TPGM_est[[i]] <- TPR1
  
  print(i)
}

## SPGM
#FPR_record_SPGM_est <- list()
#TPR_record_SPGM_est <- list()

#max_lambda <- 10
#min_lambda <- 0.001*max_lambda
#lambda_path <- NULL
#for(i in 1:1000){
#  lambda_path[i] <- max_lambda/(1000^((i-1)/(1000-1))) 
#}

#for(i in 1:iteration){
#  SPGM_fit <- ModPGMInference(x = X[[i]], model = "SPGM", D_0 = rep(2,p), D_1 = rep(5,p), nlambda = 100, regularization = lambda_path)	
  
#  FPR1 <- NULL
#  TPR1 <- NULL
  
#  for(m in 1:1000){
#    est <- SPGM_fit$theta_initial[[m]]
    
    ## symmetrize
#    for(k in 1:(p-1)){
#      for(j in (k+1):p){
#        est[k,j] = (est[k,j]+est[j,k])/2
#        est[j,k] = est[k,j]
#      }
#    }
    
#    est <- est[upper.tri(est)]
#    true <- Omega[upper.tri(Omega)]
#    FP <- sum(est!=0 & true==0)
#    TN <- sum(est==0 & true==0)
#    FPR1[m] <- FP/(FP+TN)
#    TP <- sum(est!=0 & true!=0)
#    FN <- sum(est==0 & true!=0)
#    TPR1[m] <- TP/(TP+FN)
    
#  }	
  
#  FPR_record_SPGM_est[[i]] <- FPR1
#  TPR_record_SPGM_est[[i]] <- TPR1
  
#  print(i)
#}

## SqrtPGM
#FPR_record_SqrtPGM_est <- list()
#TPR_record_SqrtPGM_est <- list()

#max_lambda <- 10
#min_lambda <- 0.001*max_lambda
#lambda_path <- NULL
#for(i in 1:1000){
#  lambda_path[i] <- max_lambda/(1000^((i-1)/(1000-1))) 
#}

#for(i in 1:iteration){
#  SqrtPGM_fit <- ModPGMInference(x = X[[i]], model = "SqrtPGM", nlambda = 100, regularization = lambda_path)
  
#  FPR1 <- NULL
#  TPR1 <- NULL
  
#  for(m in 1:1000){
#    est <- SqrtPGM_fit$theta_initial[[m]]
    
    ## symmetrize
#    for(k in 1:(p-1)){
#      for(j in (k+1):p){
#        est[k,j] = (est[k,j]+est[j,k])/2
#        est[j,k] = est[k,j]
#      }
#    }
    
#    est <- est[upper.tri(est)]
#    true <- Omega[upper.tri(Omega)]
#    FP <- sum(est!=0 & true==0)
#    TN <- sum(est==0 & true==0)
#    FPR1[m] <- FP/(FP+TN)
#    TP <- sum(est!=0 & true!=0)
#    FN <- sum(est==0 & true!=0)
#    TPR1[m] <- TP/(TP+FN)
    
#  }	
  
#  FPR_record_SqrtPGM_est[[i]] <- FPR1
#  TPR_record_SqrtPGM_est[[i]] <- TPR1
  
#  print(i)
#}

## ROC curve for our method with EBIC and the sole estimation
iteration <- 100

## TPGM
A <- matrix(0,iteration,1000)
for(i in 1:iteration){
  A[i,] <- FPR_record_TPGM_EBIC[[i]]
}

FPR_TPGM_EBIC <- apply(A, 2, median)

A <- matrix(0,iteration,1000)
for(i in 1:iteration){
  A[i,] <- TPR_record_TPGM_EBIC[[i]]
}

TPR_TPGM_EBIC <- apply(A, 2, median)

A <- matrix(0,iteration,1000)
for(i in 1:iteration){
  A[i,] <- FPR_record_TPGM_est[[i]]
}

FPR_TPGM_est <- apply(A, 2, median)

A <- matrix(0,iteration,1000)
for(i in 1:iteration){
  A[i,] <- TPR_record_TPGM_est[[i]]
}

TPR_TPGM_est <- apply(A, 2, median)

## SPGM
#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- FPR_record_SPGM_EBIC[[i]]
#}

#FPR_SPGM_EBIC <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- TPR_record_SPGM_EBIC[[i]]
#}

#TPR_SPGM_EBIC <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- FPR_record_SPGM_est[[i]]
#}

#FPR_SPGM_est <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- TPR_record_SPGM_est[[i]]
#}

#TPR_SPGM_est <- apply(A, 2, median)

## SqrtPGM
#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- FPR_record_SqrtPGM_EBIC[[i]]
#}

#FPR_SqrtPGM_EBIC <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- TPR_record[[i]]
#}

#TPR_SqrtPGM_EBIC <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- FPR_record_SqrtPGM_est[[i]]
#}

#FPR_SqrtPGM_est <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- TPR_record_SqrtPGM_est[[i]]
#}

#TPR_SqrtPGM_est <- apply(A, 2, median)

plot(FPR_TPGM_EBIC,TPR_TPGM_EBIC,type="l",col="red", main = "Chain (TPGM)", xlim = c(0,0.2),lwd = 4, xlab = "FPR", ylab = "TPR", cex.lab = 1.5, cex.main = 2.5)
lines(FPR_TPGM_est,TPR_TPGM_est,col="blue",lty = 2, lwd = 4)
legend("bottomright", legend=c("Proposed", "Estimation"),
       col=c("red", "blue"), lty = 1:2, lwd = 2:2, cex=1.5)
#plot(FPR_SPGM_EBIC,TPR_SPGM_EBIC,type="l",col="red", xlim = c(0,0.2),main = "Chain (SPGM)", lwd = 4,xlab = "FPR", ylab = "TPR",cex.lab = 1.5, cex.main = 2.5)
#lines(FPR_SPGM_est,TPR_SPGM_est,col="blue",lty = 2, lwd = 4)
#legend("bottomright", legend=c("Proposed", "Estimation"),
#       col=c("red", "blue"), lty = 1:2, lwd = 2:2, cex=1.5)
#plot(FPR_SqrtPGM_EBIC,TPR_SqrtPGM_EBIC,type="l",col="red", xlim = c(0,0.2),main = "Chain (SqrtPGM)",lwd = 4,xlab = "FPR", ylab = "TPR",cex.lab = 1.5, cex.main = 2.5)
#lines(FPR_SqrtPGM_est,TPR_SqrtPGM_est,col="blue",lty = 2, lwd = 4)
#legend("bottomright", legend=c("Proposed", "Estimation"),
#       col=c("red", "blue"), lty = 1:2, lwd = 2:2, cex=1.5)


## ROC curve for our method with EBIC and CV
iteration <- 100

## TPGM
A <- matrix(0,iteration,1000)
for(i in 1:iteration){
  A[i,] <- FPR_record_TPGM_EBIC[[i]]
}

FPR_TPGM_EBIC <- apply(A, 2, median)

A <- matrix(0,iteration,1000)
for(i in 1:iteration){
  A[i,] <- TPR_record_TPGM_EBIC[[i]]
}

TPR_TPGM_EBIC <- apply(A, 2, median)

A <- matrix(0,iteration,1000)
for(i in 1:iteration){
  A[i,] <- FPR_record_TPGM_CV[[i]]
}

FPR_TPGM_CV <- apply(A, 2, median)

A <- matrix(0,iteration,1000)
for(i in 1:iteration){
  A[i,] <- TPR_record_TPGM_CV[[i]]
}

TPR_TPGM_CV <- apply(A, 2, median)

## SPGM
#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- FPR_record_SPGM_EBIC[[i]]
#}

#FPR_SPGM_EBIC <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- TPR_record_SPGM_EBIC[[i]]
#}

#TPR_SPGM_EBIC <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- FPR_record_SPGM_CV[[i]]
#}

#FPR_SPGM_CV <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- TPR_record_SPGM_CV[[i]]
#}

#TPR_SPGM_CV <- apply(A, 2, median)

## SqrtPGM
#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- FPR_record_SqrtPGM_EBIC[[i]]
#}

#FPR_SqrtPGM_EBIC <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- TPR_record_SqrtPGM_EBIC[[i]]
#}

#TPR_SqrtPGM_EBIC <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- FPR_record_SqrtPGM_CV[[i]]
#}

#FPR_SqrtPGM_CV <- apply(A, 2, median)

#A <- matrix(0,iteration,1000)
#for(i in 1:iteration){
#  A[i,] <- TPR_record_SqrtPGM_CV[[i]]
#}

#TPR_SqrtPGM_CV <- apply(A, 2, median)

plot(FPR_TPGM_EBIC,TPR_TPGM_EBIC,type="l",col="red", main = "Chain (TPGM)", xlim = c(0,0.2),lwd = 4, xlab = "FPR", ylab = "TPR", cex.lab = 1.5, cex.main = 2.5)
lines(FPR_TPGM_CV,TPR_TPGM_CV,col="black",lty = 2, lwd = 4)
legend("bottomright", legend=c("Proposed (EBIC)", "Proposed (Cross validation)"),
       col=c("red", "black"), lty = 1:2, lwd = 2:2, cex=1.5)
#plot(FPR_SPGM_EBIC,TPR_SPGM_EBIC,type="l",col="red", xlim = c(0,0.2),main = "Chain (SPGM)", lwd = 4,xlab = "FPR", ylab = "TPR",cex.lab = 1.5, cex.main = 2.5)
#lines(FPR_SPGM_CV,TPR_SPGM_CV,col="black",lty = 2, lwd = 4)
#legend("bottomright", legend=c("Proposed (EBIC)", "Proposed (Cross validation)"),
#       col=c("red", "black"), lty = 1:2, lwd = 2:2, cex=1.5)
#plot(FPR_SqrtPGM_EBIC,TPR_SqrtPGM_EBIC,type="l",col="red", xlim = c(0,0.2),main = "Chain (SqrtPGM)",lwd = 4,xlab = "FPR", ylab = "TPR",cex.lab = 1.5, cex.main = 2.5)
#lines(FPR_SqrtPGM_CV,TPR_SqrtPGM_CV,col="black",lty = 2, lwd = 4)
#legend("bottomright", legend=c("Proposed (EBIC)", "Proposed (Cross validation)"),
#       col=c("red", "black"), lty = 1:2, lwd = 2:2, cex=1.5)
