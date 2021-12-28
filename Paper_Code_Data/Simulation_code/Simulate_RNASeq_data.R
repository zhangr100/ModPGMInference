library(ModPGMInference)
library(SILGGM)
library(huge)

######## Evaluation on simulated RNA-seq data ########
## Scale-free
p <- 400
Omega.tmp=matrix(0,p,p)

library(igraph)
g <- barabasi.game(n=p/2, power=0.01, directed=F)
block1 <- as.matrix(get.adjacency(g, type="both"))
for(i in 1:(p/2-1)){
  for(j in (i+1):(p/2)){
    if(block1[i,j]!=0){
      block1[i,j] = sample(c(-0.6,0.6),1)
      block1[j,i] = block1[i,j]
    }
  }
}
g <- barabasi.game(n=p/2, power=0.01, directed=F)
block2 <- as.matrix(get.adjacency(g, type="both"))
for(i in 1:(p/2-1)){
  for(j in (i+1):(p/2)){
    if(block2[i,j]!=0){
      block2[i,j] = sample(c(-0.9,0.9),1)
      block2[j,i] = block2[i,j]
    }
  }
}

Omega.tmp[1:(p/2),1:(p/2)] <- block1
Omega.tmp[(p/2+1):p,(p/2+1):p] <- block2
Omega=Omega.tmp

## Generate count matrix of RNA-seq data
n <- 300  ## n <- 150
p <- 400
iteration <- 100

count <- list() 

for(i in 1:iteration){
  ## Step 1: Simulate normalized RNA-seq data
  X <- ModPGMSampler(psi = rep(0,p), model = "SqrtPGM", true_graph = Omega, nSample = n, burn_in = 5000)
  
  ## Step 2: Pseudo-random number addition
  uniform <- matrix(0,n,p)
  for(k in 1:n){
    uniform[k,] <- runif(p,0,1)
  }
  
  X_new <- X + uniform
  
  ## Step 3: Inverse power transform
  count_value <- exp(log(X_new)/0.2517)
  
  ## Step 4: Final count generation
  count[[i]] <- floor(count_value)
}

## Implement our method on the simulated RNA-seq data
## SqrtPGM
FDR_record_SqrtPGM <- list()
power_record_SqrtPGM <- list()
FPR_record_SqrtPGM <- list()
TPR_record_SqrtPGM <- list()

for(i in 1:iteration){
  power_transformation <- count[[i]]^(0.2517)
  power_transformation <- floor(power_transformation)
  
  fit <- ModPGMInference(x = power_transformation, model = "SqrtPGM", nlambda = 100, global = TRUE, true_graph = Omega, alpha = seq(0.001,1,by = 0.001))
  FPR <- NULL
  TPR <- NULL
  FDR_record_SqrtPGM[[i]] <- fit$FDR
  power_record_SqrtPGM[[i]] <- fit$power
  
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
  
  FPR_record_SqrtPGM[[i]] <- FPR
  TPR_record_SqrtPGM[[i]] <- TPR
  
  print(i)
}

## SPGM
FDR_record_SPGM <- list()
power_record_SPGM <- list()
FPR_record_SPGM <- list()
TPR_record_SPGM <- list()

for(i in 1:iteration){
  power_transformation <- count[[i]]^(0.2517)
  power_transformation <- floor(power_transformation)
  
  fit <- ModPGMInference(x = power_transformation, model = "SPGM", nlambda = 100, D_0 = rep(0,dim(power_transformation)[2]), D_1 = apply(power_transformation,2,max), global = TRUE, true_graph = Omega, alpha = seq(0.001,1,by = 0.001))
  FPR <- NULL
  TPR <- NULL
  FDR_record_SPGM[[i]] <- fit$FDR
  power_record_SPGM[[i]] <- fit$power
  
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
  
  FPR_record_SPGM[[i]] <- FPR
  TPR_record_SPGM[[i]] <- TPR
  
  print(i)
}

## Implement GFC_L on the simulated RNA-seq data
FDR_record_GFC_L <- list()
power_record_GFC_L <- list()
FPR_record_GFC_L <- list()
TPR_record_GFC_L <- list()

for(i in 1:iteration){
  
  fit <- SILGGM(x = huge.npn(log(count[[i]]+1)),method = "GFC_L", ndelta = 20, true_graph = Omega, alpha = seq(0.001,1,by = 0.001))
  FPR <- NULL
  TPR <- NULL
  FDR_record_GFC_L[[i]] <- fit$FDR
  power_record_GFC_L[[i]] <- fit$power
  
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
  
  FPR_record_GFC_L[[i]] <- FPR
  TPR_record_GFC_L[[i]] <- TPR
  
  print(i)
}

## Medians (SDs) of FDR values at the 0.1 and the 0.2 levels and their power values
FDR_SqrtPGM_0.1 <- NULL
power_SqrtPGM_0.1 <- NULL
FDR_SqrtPGM_0.2 <- NULL
power_SqrtPGM_0.2 <- NULL

FDR_SPGM_0.1 <- NULL
power_SPGM_0.1 <- NULL
FDR_SPGM_0.2 <- NULL
power_SPGM_0.2 <- NULL

FDR_GFC_L_0.1 <- NULL
power_GFC_L_0.1 <- NULL
FDR_GFC_L_0.2 <- NULL
power_GFC_L_0.2 <- NULL

for(i in 1:iteration){
  FDR_SqrtPGM_0.1[i] <- FDR_record_SqrtPGM[[i]][100]
  power_SqrtPGM_0.1[i] <- power_record_SqrtPGM[[i]][100]
  FDR_SqrtPGM_0.2[i] <- FDR_record_SqrtPGM[[i]][200]
  power_SqrtPGM_0.2[i] <- power_record_SqrtPGM[[i]][200]
  
  FDR_SPGM_0.1[i] <- FDR_record_SPGM[[i]][100]
  power_SPGM_0.1[i] <- power_record_SPGM[[i]][100]
  FDR_SPGM_0.2[i] <- FDR_record_SPGM[[i]][200]
  power_SPGM_0.2[i] <- power_record_SPGM[[i]][200]
  
  FDR_GFC_L_0.1[i] <- FDR_record_GFC_L[[i]][100]
  power_GFC_L_0.1[i] <- power_record_GFC_L[[i]][100]
  FDR_GFC_L_0.2[i] <- FDR_record_GFC_L[[i]][200]
  power_GFC_L_0.2[i] <- power_record_GFC_L[[i]][200]
}

FDR_SqrtPGM_0.1_median <- median(FDR_SqrtPGM_0.1)
FDR_SqrtPGM_0.1_sd <- sd(FDR_SqrtPGM_0.1)
FDR_SqrtPGM_0.2_median <- median(FDR_SqrtPGM_0.2)
FDR_SqrtPGM_0.2_sd <- sd(FDR_SqrtPGM_0.2)
FDR_SPGM_0.1_median <- median(FDR_SPGM_0.1)
FDR_SPGM_0.1_sd <- sd(FDR_SPGM_0.1)
FDR_SPGM_0.2_median <- median(FDR_SPGM_0.2)
FDR_SPGM_0.2_sd <- sd(FDR_SPGM_0.2)
FDR_GFC_L_0.1_median <- median(FDR_GFC_L_0.1)
FDR_GFC_L_0.1_sd <- sd(FDR_GFC_L_0.1)
FDR_GFC_L_0.2_median <- median(FDR_GFC_L_0.2)
FDR_GFC_L_0.2_sd <- sd(FDR_GFC_L_0.2)

## ROC-type curve (FDR and power) for our method (SqrtPGM and SPGM) and GFC_L
iteration <- 100

A <- matrix(0,iteration,1000)
B <- matrix(0,iteration,1000)

for(i in 1:iteration){
  A[i,] <- FDR_record_SqrtPGM[[i]]
  B[i,] <- power_record_SqrtPGM[[i]]
}

FDR_SqrtPGM <- apply(A,2,median)
power_SqrtPGM <- apply(B,2,median)

A <- matrix(0,iteration,1000)
B <- matrix(0,iteration,1000)

for(i in 1:iteration){
  A[i,] <- FDR_record_SPGM[[i]]
  B[i,] <- power_record_SPGM[[i]]
}

FDR_SPGM <- apply(A,2,median)
power_SPGM <- apply(B,2,median)

A <- matrix(0,iteration,1000)
B <- matrix(0,iteration,1000)

for(i in 1:iteration){
  A[i,] <- FDR_record_GFC_L[[i]]
  B[i,] <- power_record_GFC_L[[i]]
}

FDR_GFC_L <- apply(A,2,median)
power_GFC_L <- apply(B,2,median)

plot(c(0,FDR_SqrtPGM),c(0,power_SqrtPGM),type="l",col="red", main = NULL, xlim = c(0,0.5),lwd = 4, xlab = "FDR", ylab = "TPR", cex.lab = 1.5, cex.main = 2.5)
lines(c(0,FDR_SPGM),c(0,power_SPGM),col="green",lty = 4, lwd = 4)
lines(c(0,FDR_GFC_L),c(0,power_GFC_L),col="black",lty = 3, lwd = 4)
legend("bottomright", legend=c("Proposed (SqrtPGM)", "Proposed (SPGM)", "GFC_L"),
       col=c("red", "green", "black"), lty = c(1,4,3), lwd = c(4,4,4), cex=1.5)


## Histograms of RNA-seq data for childhood atopic asthma and simulated RNA-seq data
load("ge_atopic_asthma_original.RData")

par(mfrow = c(1,2))
A <- hist(ge_atopic_asthma_original,xlab = "Gene count", main = "Real RNA-seq data of childhood atopic asthma", xlim = c(0,500), breaks = 5000, probability = TRUE, cex.lab = 1.5, cex.main = 1.5)
hist(count[[1]],xlab = "Gene count", main = "Simulated RNA-seq data", xlim = c(0,500), ylim = c(0, max(A$density)), breaks = 725, probability = TRUE, cex.lab = 1.5, cex.main = 1.5)
