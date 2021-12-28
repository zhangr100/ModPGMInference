library(ModPGMInference)

######## Graph settings: Chain, Grid, E-R, Scale-free ########

## Chain graph
set <- c(-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4)
p <- 100 ## p <- 400

Omega.tmp=matrix(0,p,p)
for(i in 1:(p-1)){
  j=i+1
  Omega.tmp[i,j]=sample(set,1)
}
for(i in 1:(p-1)){
  for(j in (i+1):p){
    Omega.tmp[j,i]=Omega.tmp[i,j]
  }
}
Omega=Omega.tmp

###########################################

## Grid graph
#set <- c(-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4)
#p <- 100 ## p <- 400

#Grid_Graph <- function(x, y){
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
#x <- 10 # Number of rows    ## x <- 20, y <- 20 for p <- 400
#y <- 10 # Number of columns 
#Omega <- Grid_Graph(x,y)

###########################################

## E-R graph
#set <- c(-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4)
#p <- 100 ## p <- 400

#Omega.tmp <- matrix(0,p,p)
#prop = 4/p # sparseness
#for(i in 1:(p-1)){
#  for(j in (i+1):p){
#    Omega.tmp[i, j]=sample(set,1)*rbinom(1,1,prop)
#    Omega.tmp[j, i]=Omega.tmp[i, j]
#  }
#}
#Omega <- Omega.tmp

###########################################

## Scale-free graph
#set <- c(-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4)
#p <- 100 ## p <- 400

#library(igraph)
#g <- barabasi.game(n=p, power=0.01, directed=F)
#theta <- as.matrix(get.adjacency(g, type="both"))
#diag(theta) <- 0
#for(i in 1:(p-1)){
#  for(j in (i+1):p){
#    if(theta[i,j]!=0){
#      theta[i,j] = sample(set,1)
#      theta[j,i] = theta[i,j]
#    }
#  }
#}
#Omega <- theta


######## Generate random samples for individual inference ######## 

n <- 300  ## n <- 100, 150
p <- 100  ## p <- 400
iteration <- 100

## TPGM
X <- list()
for(i in 1:iteration){
  X[[i]] <- ModPGMSampler(psi = rep(0,p), true_graph = Omega, model = "TPGM", D = rep(3,p), nSample = n, burn_in = 5000)
}

## SPGM
#X <- list()
#for(i in 1:iteration){
#  X[[i]] <- ModPGMSampler(psi = rep(-0.5,p), true_graph = Omega, model = "SPGM", D_0 = rep(3,p), D_1 = rep(6,p), nSample = n, burn_in = 5000)
#}

## SqrtPGM
#X <- list()
#for(i in 1:iteration){
#  X[[i]] <- ModPGMSampler(psi = rep(0,p), true_graph = Omega, model = "SqrtPGM", nSample = n, burn_in = 5000)
#}


######## Perform individual inference on random samples ######## 

p <- 100  ## p <- 400
iteration <- 100

## Tuning parameter selection: EBIC
## TPGM
fit_TPGM_EBIC <- list()

for(i in 1:iteration){
  fit_TPGM_EBIC[[i]] <- ModPGMInference(x = X[[i]], model = "TPGM", tuning = "EBIC", D = rep(3,p), nlambda = 100)
  print(i)
}

## SPGM
#fit_SPGM_EBIC <- list()

#for(i in 1:iteration){
#  fit_SPGM_EBIC[[i]] <- ModPGMInference(x = X[[i]], model = "SPGM", tuning = "EBIC", D_0 = rep(3,p), D_1 = rep(6,p), nlambda = 100)
#  print(i)
#}

## SqrtPGM
#fit_SqrtPGM_EBIC <- list()

#for(i in 1:iteration){
#  fit_SqrtPGM_EBIC[[i]] <- ModPGMInference(x = X[[i]], model = "SqrtPGM", tuning = "EBIC", nlambda = 100)
#  print(i)
#}

## Tuning parameter selection: Cross validation (CV)
## TPGM 
fit_TPGM_CV <- list()

for(i in 1:iteration){
  fit_TPGM_CV[[i]] <- ModPGMInference(x = X[[i]], model = "TPGM", tuning = "CV", kfold = 10, D = rep(3,p), nlambda = 100)
  print(i)
}

## SPGM
#fit_SPGM_CV <- list()

#for(i in 1:iteration){
#  fit_SPGM_CV[[i]] <- ModPGMInference(x = X[[i]], model = "SPGM", tuning = "CV", kfold = 10, D_0 = rep(3,p), D_1 = rep(6,p), nlambda = 100)
#  print(i)
#}

## SqrtPGM
#fit_SqrtPGM_CV <- list()

#for(i in 1:iteration){
#  fit_SqrtPGM_CV[[i]] <- ModPGMInference(x = X[[i]], model = "SqrtPGM", tuning = "CV", kfold = 10, nlambda = 100)
#  print(i)
#}

## Evaluation of empirical coverage of 95% confidence intervals
CI_nonzero_TPGM_EBIC <- NULL   ## confidence intervals for non-zero edge parameters 
CI_zero_TPGM_EBIC <- NULL      ## confidence intervals for zero edge parameters
#CI_nonzero_SPGM_EBIC <- NULL
#CI_zero_SPGM_EBIC <- NULL
#CI_nonzero_SqrtPGM_EBIC <- NULL
#CI_zero_SqrtPGM_EBIC <- NULL

for(i in 1:iteration){
  
  ## TPGM
  true <- Omega[upper.tri(Omega)]
  CI_low <- fit_TPGM_EBIC[[i]]$CI_low_theta[upper.tri(fit_TPGM_EBIC[[i]]$CI_low_theta)]
  CI_high <- fit_TPGM_EBIC[[i]]$CI_high_theta[upper.tri(fit_TPGM_EBIC[[i]]$CI_high_theta)]
  CI_low_nonzero <- CI_low[true!=0]
  CI_high_nonzero <- CI_high[true!=0]
  true_non_zero <- true[true!=0]
  
  CI_nonzero_TPGM_EBIC[i] <- sum(CI_low_nonzero <= true_non_zero & true_non_zero <= CI_high_nonzero)/sum(true!=0)
  
  CI_low_zero <- CI_low[true==0]
  CI_high_zero <- CI_high[true==0]
  true_zero <- true[true==0]
  
  CI_zero_TPGM_EBIC[i] <- sum(CI_low_zero <= true_zero & true_zero <= CI_high_zero)/sum(true==0)

  ## SPGM
  #true <- Omega[upper.tri(Omega)]
  #CI_low <- fit_SPGM_EBIC[[i]]$CI_low_theta[upper.tri(fit_SPGM_EBIC[[i]]$CI_low_theta)]
  #CI_high <- fit_SPGM_EBIC[[i]]$CI_high_theta[upper.tri(fit_SPGM_EBIC[[i]]$CI_high_theta)]
  #CI_low_nonzero <- CI_low[true!=0]
  #CI_high_nonzero <- CI_high[true!=0]
  #true_non_zero <- true[true!=0]
  
  #CI_nonzero_SPGM_EBIC[i] <- sum(CI_low_nonzero <= true_non_zero & true_non_zero <= CI_high_nonzero)/sum(true!=0)
  
  #CI_low_zero <- CI_low[true==0]
  #CI_high_zero <- CI_high[true==0]
  #true_zero <- true[true==0]
  
  #CI_zero_SPGM_EBIC[i] <- sum(CI_low_zero <= true_zero & true_zero <= CI_high_zero)/sum(true==0)
  
  ## SqrtPGM
  #true <- Omega[upper.tri(Omega)]
  #CI_low <- fit_SqrtPGM_EBIC[[i]]$CI_low_theta[upper.tri(fit_SqrtPGM_EBIC[[i]]$CI_low_theta)]
  #CI_high <- fit_SqrtPGM_EBIC[[i]]$CI_high_theta[upper.tri(fit_SqrtPGM_EBIC[[i]]$CI_high_theta)]
  #CI_low_nonzero <- CI_low[true!=0]
  #CI_high_nonzero <- CI_high[true!=0]
  #true_non_zero <- true[true!=0]
  
  #CI_nonzero_SqrtPGM_EBIC[i] <- sum(CI_low_nonzero <= true_non_zero & true_non_zero <= CI_high_nonzero)/sum(true!=0)
  
  #CI_low_zero <- CI_low[true==0]
  #CI_high_zero <- CI_high[true==0]
  #true_zero <- true[true==0]
  
  #CI_zero_SqrtPGM_EBIC[i] <- sum(CI_low_zero <= true_zero & true_zero <= CI_high_zero)/sum(true==0)

}

CI_nonzero_TPGM_CV <- NULL   ## confidence intervals for non-zero edge parameters 
CI_zero_TPGM_CV <- NULL      ## confidence intervals for zero edge parameters
#CI_nonzero_SPGM_CV <- NULL
#CI_zero_SPGM_CV <- NULL
#CI_nonzero_SqrtPGM_CV <- NULL
#CI_zero_SqrtPGM_CV <- NULL

for(i in 1:iteration){
  
  ## TPGM
  true <- Omega[upper.tri(Omega)]
  CI_low <- fit_TPGM_CV[[i]]$CI_low_theta[upper.tri(fit_TPGM_CV[[i]]$CI_low_theta)]
  CI_high <- fit_TPGM_CV[[i]]$CI_high_theta[upper.tri(fit_TPGM_CV[[i]]$CI_high_theta)]
  CI_low_nonzero <- CI_low[true!=0]
  CI_high_nonzero <- CI_high[true!=0]
  true_non_zero <- true[true!=0]
  
  CI_nonzero_TPGM_CV[i] <- sum(CI_low_nonzero <= true_non_zero & true_non_zero <= CI_high_nonzero)/sum(true!=0)
  
  CI_low_zero <- CI_low[true==0]
  CI_high_zero <- CI_high[true==0]
  true_zero <- true[true==0]
  
  CI_zero_TPGM_CV[i] <- sum(CI_low_zero <= true_zero & true_zero <= CI_high_zero)/sum(true==0)
  
  ## SPGM
  #true <- Omega[upper.tri(Omega)]
  #CI_low <- fit_SPGM_CV[[i]]$CI_low_theta[upper.tri(fit_SPGM_CV[[i]]$CI_low_theta)]
  #CI_high <- fit_SPGM_CV[[i]]$CI_high_theta[upper.tri(fit_SPGM_CV[[i]]$CI_high_theta)]
  #CI_low_nonzero <- CI_low[true!=0]
  #CI_high_nonzero <- CI_high[true!=0]
  #true_non_zero <- true[true!=0]
  
  #CI_nonzero_SPGM_CV[i] <- sum(CI_low_nonzero <= true_non_zero & true_non_zero <= CI_high_nonzero)/sum(true!=0)
  
  #CI_low_zero <- CI_low[true==0]
  #CI_high_zero <- CI_high[true==0]
  #true_zero <- true[true==0]
  
  #CI_zero_SPGM_CV[i] <- sum(CI_low_zero <= true_zero & true_zero <= CI_high_zero)/sum(true==0)
  
  ## SqrtPGM
  #true <- Omega[upper.tri(Omega)]
  #CI_low <- fit_SqrtPGM_CV[[i]]$CI_low_theta[upper.tri(fit_SqrtPGM_CV[[i]]$CI_low_theta)]
  #CI_high <- fit_SqrtPGM_CV[[i]]$CI_high_theta[upper.tri(fit_SqrtPGM_CV[[i]]$CI_high_theta)]
  #CI_low_nonzero <- CI_low[true!=0]
  #CI_high_nonzero <- CI_high[true!=0]
  #true_non_zero <- true[true!=0]
  
  #CI_nonzero_SqrtPGM_CV[i] <- sum(CI_low_nonzero <= true_non_zero & true_non_zero <= CI_high_nonzero)/sum(true!=0)
  
  #CI_low_zero <- CI_low[true==0]
  #CI_high_zero <- CI_high[true==0]
  #true_zero <- true[true==0]
  
  #CI_zero_SqrtPGM_CV[i] <- sum(CI_low_zero <= true_zero & true_zero <= CI_high_zero)/sum(true==0)
  
}

## Median values (SDs) of empirical coverage of 95% confidence intervals from 100 replications
CI_nonzero_TPGM_EBIC_median <- median(CI_nonzero_TPGM_EBIC)
CI_nonzero_TPGM_EBIC_sd <- sd(CI_nonzero_TPGM_EBIC)
#CI_nonzero_SPGM_EBIC_median <- median(CI_nonzero_SPGM_EBIC)
#CI_nonzero_SPGM_EBIC_sd <- sd(CI_nonzero_SPGM_EBIC)
#CI_nonzero_SqrtPGM_EBIC_median <- median(CI_nonzero_SqrtPGM_EBIC)
#CI_nonzero_SqrtPGM_EBIC_sd <- sd(CI_nonzero_SqrtPGM_EBIC)

CI_zero_TPGM_EBIC_median <- median(CI_zero_TPGM_EBIC)
CI_zero_TPGM_EBIC_sd <- sd(CI_zero_TPGM_EBIC)
#CI_zero_SPGM_EBIC_median <- median(CI_zero_SPGM_EBIC)
#CI_zero_SPGM_EBIC_sd <- sd(CI_zero_SPGM_EBIC)
#CI_zero_SqrtPGM_EBIC_median <- median(CI_zero_SqrtPGM_EBIC)
#CI_zero_SqrtPGM_EBIC_sd <- sd(CI_zero_SqrtPGM_EBIC)


CI_nonzero_TPGM_CV_median <- median(CI_nonzero_TPGM_CV)
CI_nonzero_TPGM_CV_sd <- sd(CI_nonzero_TPGM_CV)
#CI_nonzero_SPGM_CV_median <- median(CI_nonzero_SPGM_CV)
#CI_nonzero_SPGM_CV_sd <- sd(CI_nonzero_SPGM_CV)
#CI_nonzero_SqrtPGM_CV_median <- median(CI_nonzero_SqrtPGM_CV)
#CI_nonzero_SqrtPGM_CV_sd <- sd(CI_nonzero_SqrtPGM_CV)

CI_zero_TPGM_CV_median <- median(CI_zero_TPGM_CV)
CI_zero_TPGM_CV_sd <- sd(CI_zero_TPGM_CV)
#CI_zero_SPGM_CV_median <- median(CI_zero_SPGM_CV)
#CI_zero_SPGM_CV_sd <- sd(CI_zero_SPGM_CV)
#CI_zero_SqrtPGM_CV_median <- median(CI_zero_SqrtPGM_CV)
#CI_zero_SqrtPGM_CV_sd <- sd(CI_zero_SqrtPGM_CV)


######## Evaluation of asymptotic normality ######## 
## Note: due to a different true graph, the true values of corresponding edge parameter might not be exactly ones shown in the manuscript
## Therefore, the true values of edge parameter in titles of histograms may need to be changed to the corrsponding elements from the generated true graph settings

p <- 100  ## p <- 400
iteration <- 100

## Obtain true standard deviation of each edge parameter
## TPGM
sd_TPGM <- list()

for(i in 1:iteration){
  TPGM_sd <- ModPGM_true_sd(x=X[[i]], psi=rep(0,p), model = "TPGM", true_graph = Omega, D = rep(3,p)) 
  sd_TPGM[[i]] <- TPGM_sd
  print(i)
}

## SPGM
#sd_SPGM <- list()

#for(i in 1:iteration){
#  SPGM_sd <- ModPGM_true_sd(x=X[[i]], psi=rep(-0.5,p), model = "SPGM", true_graph = Omega, D_0 = rep(3,p), D_1 = rep(6,p)) 
#  sd_SPGM[[i]] <- SPGM_sd
#  print(i)
#}

## SqrtPGM
#sd_SqrtPGM <- list()

#for(i in 1:iteration){
#  SqrtPGM_sd <- ModPGM_true_sd(x=X[[i]], psi = rep(0,p), true_graph = Omega) 
#  sd_SqrtPGM[[i]] <- SqrtPGM_sd
#  print(i)
#}

## TPGM
theta12_TPGM <- NULL
sd12_TPGM_true <- NULL
theta23_TPGM <- NULL
sd23_TPGM_true <- NULL
theta78_TPGM <- NULL
sd78_TPGM_true <- NULL
theta47_TPGM <- NULL
sd47_TPGM_true <- NULL
for(i in 1:iteration){
  theta12_TPGM[i] <- fit_TPGM_EBIC[[i]]$theta_cor[1,2]
  theta23_TPGM[i] <- fit_TPGM_EBIC[[i]]$theta_cor[2,3]
  theta78_TPGM[i] <- fit_TPGM_EBIC[[i]]$theta_cor[7,8]
  theta47_TPGM[i] <- fit_TPGM_EBIC[[i]]$theta_cor[4,7]
  sd12_TPGM_true[i] <- sd_TPGM[[i]][1,2]
  sd23_TPGM_true[i] <- sd_TPGM[[i]][2,3]
  sd78_TPGM_true[i] <- sd_TPGM[[i]][7,8]
  sd47_TPGM_true[i] <- sd_TPGM[[i]][4,7]
}

## SPGM
#theta12_SPGM <- NULL
#sd12_SPGM_true <- NULL
#theta34_SPGM <- NULL
#sd34_SPGM_true <- NULL
#theta1516_SPGM <- NULL
#sd1516_SPGM_true <- NULL
#theta48_SPGM <- NULL
#sd48_SPGM_true <- NULL
#for(i in 1:iteration){
#  theta12_SPGM[i] <- fit_SPGM_EBIC[[i]]$theta_cor[1,2]
#  theta34_SPGM[i] <- fit_SPGM_EBIC[[i]]$theta_cor[3,4]
#  theta1516_SPGM[i] <- fit_SPGM_EBIC[[i]]$theta_cor[15,16]
#  theta48_SPGM[i] <- fit_SPGM_EBIC[[i]]$theta_cor[4,8]
#  sd12_SPGM_true[i] <- sd_SPGM[[i]][1,2]
#  sd34_SPGM_true[i] <- sd_SPGM[[i]][3,4]
#  sd1516_SPGM_true[i] <- sd_SPGM[[i]][15,16]
#  sd48_SPGM_true[i] <- sd_SPGM[[i]][4,8]
#}

## SqrtPGM
#theta67_SqrtPGM <- NULL
#sd67_SqrtPGM_true <- NULL
#theta89_SqrtPGM <- NULL
#sd89_SqrtPGM_true <- NULL
#theta34_SqrtPGM <- NULL
#sd34_SqrtPGM_true <- NULL
#theta59_SqrtPGM <- NULL
#sd59_SqrtPGM_true <- NULL
#for(i in 1:iteration){
#  theta67_SqrtPGM[i] <- fit_SqrtPGM_EBIC[[i]]$theta_cor[6,7]
#  theta89_SqrtPGM[i] <- fit_SqrtPGM_EBIC[[i]]$theta_cor[8,9]
#  theta34_SqrtPGM[i] <- fit_SqrtPGM_EBIC[[i]]$theta_cor[3,4]
#  theta59_SqrtPGM[i] <- fit_SqrtPGM_EBIC[[i]]$theta_cor[5,9]
#  sd67_SqrtPGM_true[i] <- sd_SqrtPGM[[i]][6,7]
#  sd89_SqrtPGM_true[i] <- sd_SqrtPGM[[i]][8,9]
#  sd34_SqrtPGM_true[i] <- sd_SqrtPGM[[i]][3,4]
#  sd59_SqrtPGM_true[i] <- sd_SqrtPGM[[i]][5,9]
#}

sd12_TPGM_true <- median(sd12_TPGM_true)
sd23_TPGM_true <- median(sd23_TPGM_true)
sd78_TPGM_true <- median(sd78_TPGM_true)
sd47_TPGM_true <- median(sd47_TPGM_true)
#sd12_SPGM_true <- median(sd12_SPGM_true)
#sd34_SPGM_true <- median(sd34_SPGM_true)
#sd1516_SPGM_true <- median(sd1516_SPGM_true)
#sd48_SPGM_true <- median(sd48_SPGM_true)
#sd67_SqrtPGM_true <- median(sd67_SqrtPGM_true)
#sd89_SqrtPGM_true <- median(sd89_SqrtPGM_true)
#sd34_SqrtPGM_true <- median(sd34_SqrtPGM_true)
#sd59_SqrtPGM_true <- median(sd59_SqrtPGM_true)

true_Omega_TPGM <- Omega
#true_Omega_SPGM <- Omega
#true_Omega_SqrtPGM <- Omega


par(mfrow=c(3,4))
hist(theta12_TPGM, prob=TRUE, 
     xlab=expression(tilde(theta)[12]),xlim = c(true_Omega_TPGM[1,2]-4*sd12_TPGM_true,true_Omega_TPGM[1,2]+4*sd12_TPGM_true), main=bquote(paste("TPGM"," (",theta[12]," = ",.(true_Omega_TPGM[1,2]),")")),cex.lab = 1.5, cex.main = 2.5)
curve(dnorm(x, mean=true_Omega_TPGM[1,2], sd=sd12_TPGM_true), 
      col="red", lwd=4, add=TRUE, yaxt="n")
hist(theta23_TPGM, prob=TRUE,
     xlab=expression(tilde(theta)[23]),xlim = c(true_Omega_TPGM[2,3]-4*sd23_TPGM_true,true_Omega_TPGM[2,3]+4*sd23_TPGM_true),main=bquote(paste("TPGM"," (",theta[23]," = ",.(true_Omega_TPGM[2,3]),")")),cex.lab = 1.5, cex.main = 2.5)
curve(dnorm(x, mean=true_Omega_TPGM[2,3], sd=sd23_TPGM_true),
      col="red", lwd=4, add=TRUE, yaxt="n")
hist(theta78_TPGM, prob=TRUE,
     xlab=expression(tilde(theta)[78]),xlim = c(true_Omega_TPGM[7,8]-4*sd78_TPGM_true,true_Omega_TPGM[7,8]+4*sd78_TPGM_true),main=bquote(paste("TPGM"," (",theta[78]," = ",.(true_Omega_TPGM[7,8]),")")),cex.lab = 1.5, cex.main = 2.5)
curve(dnorm(x, mean=true_Omega_TPGM[7,8], sd=sd78_TPGM_true),
      col="red", lwd=4, add=TRUE, yaxt="n")
hist(theta47_TPGM, prob=TRUE,
     xlab=expression(tilde(theta)[47]),xlim = c(true_Omega_TPGM[4,7]-4*sd47_TPGM_true,true_Omega_TPGM[4,7]+4*sd47_TPGM_true),main=bquote(paste("TPGM"," (",theta[47]," = ",.(true_Omega_TPGM[4,7]),")")),cex.lab = 1.5, cex.main = 2.5)
curve(dnorm(x, mean=true_Omega_TPGM[4,7], sd=sd47_TPGM_true),
      col="red", lwd=4, add=TRUE, yaxt="n")

#hist(theta12_SPGM, prob=TRUE,
#     xlab=expression(tilde(theta)[12]),xlim = c(true_Omega_SPGM[1,2]-4*sd12_SPGM_true,true_Omega_SPGM[1,2]+4*sd12_SPGM_true),main=bquote(paste("SPGM"," (",theta[12]," = ",.(true_Omega_SPGM[1,2]),")")),cex.lab = 1.5, cex.main = 2.5)
#curve(dnorm(x, mean=true_Omega_SPGM[1,2], sd=sd12_SPGM_true),
#      col="red", lwd=4, add=TRUE, yaxt="n")
#hist(theta34_SPGM, prob=TRUE,
#     xlab=expression(tilde(theta)[34]),xlim = c(true_Omega_SPGM[3,4]-4*sd34_SPGM_true,true_Omega_SPGM[3,4]+4*sd34_SPGM_true),main=bquote(paste("SPGM"," (",theta[34]," = ",.(true_Omega_SPGM[3,4]),")")),cex.lab = 1.5, cex.main = 2.5)
#curve(dnorm(x, mean=true_Omega_SPGM[3,4], sd=sd34_SPGM_true),
#      col="red", lwd=4, add=TRUE, yaxt="n")
#hist(theta1516_SPGM, prob=TRUE,
#     xlab=expression(tilde(theta)[list(15,16)]),xlim = c(true_Omega_SPGM[15,16]-4*sd1516_SPGM_true,true_Omega_SPGM[15,16]+4*sd1516_SPGM_true),main=bquote(paste("SPGM"," (",theta[list(15,16)]," = ",.(true_Omega_SPGM[15,16]),")")),cex.lab = 1.5, cex.main = 2.5)
#curve(dnorm(x, mean=true_Omega_SPGM[15,16], sd=sd1516_SPGM_true),
#      col="red", lwd=4, add=TRUE, yaxt="n")
#hist(theta48_SPGM, prob=TRUE,
#     xlab=expression(tilde(theta)[48]),xlim = c(true_Omega_SPGM[4,8]-4*sd48_SPGM_true,true_Omega_SPGM[4,8]+4*sd48_SPGM_true),main=bquote(paste("SPGM"," (",theta[48]," = ",.(true_Omega_SPGM[4,8]),")")),cex.lab = 1.5, cex.main = 2.5)
#curve(dnorm(x, mean=true_Omega_SPGM[4,8], sd=sd48_SPGM_true),
#      col="red", lwd=4, add=TRUE, yaxt="n")

#hist(theta67_SqrtPGM, prob=TRUE,
#     xlab=expression(tilde(theta)[67]),xlim = c(true_Omega_SqrtPGM[6,7]-4*sd67_SqrtPGM_true,true_Omega_SqrtPGM[6,7]+4*sd67_SqrtPGM_true),main=bquote(paste("SqrtPGM"," (",theta[67]," = ",.(true_Omega_SqrtPGM[6,7]),")")),cex.lab = 1.5, cex.main = 2.5)
#curve(dnorm(x, mean=true_Omega_SqrtPGM[6,7], sd=sd67_SqrtPGM_true),
#      col="red", lwd=4, add=TRUE, yaxt="n")
#hist(theta89_SqrtPGM, prob=TRUE,
#     xlab=expression(tilde(theta)[89]),xlim = c(true_Omega_SqrtPGM[8,9]-4*sd89_SqrtPGM_true,true_Omega_SqrtPGM[8,9]+4*sd89_SqrtPGM_true),main=bquote(paste("SqrtPGM"," (",theta[89]," = ",.(true_Omega_SqrtPGM[8,9]),")")),cex.lab = 1.5, cex.main = 2.5)
#curve(dnorm(x, mean=true_Omega_SqrtPGM[8,9], sd=sd89_SqrtPGM_true),
#      col="red", lwd=4, add=TRUE, yaxt="n")
#hist(theta34_SqrtPGM, prob=TRUE,
#     xlab=expression(tilde(theta)[34]),xlim = c(true_Omega_SqrtPGM[3,4]-4*sd34_SqrtPGM_true,true_Omega_SqrtPGM[3,4]+4*sd34_SqrtPGM_true),main=bquote(paste("SqrtPGM"," (",theta[34]," = ",.(true_Omega_SqrtPGM[3,4]),")")),cex.lab = 1.5, cex.main = 2.5)
#curve(dnorm(x, mean=true_Omega_SqrtPGM[3,4], sd=sd34_SqrtPGM_true),
#      col="red", lwd=4, add=TRUE, yaxt="n")
#hist(theta59_SqrtPGM, prob=TRUE,
#     xlab=expression(tilde(theta)[59]),xlim = c(true_Omega_SqrtPGM[5,9]-4*sd59_SqrtPGM_true,true_Omega_SqrtPGM[5,9]+4*sd59_SqrtPGM_true),ylim=c(0,3),main=bquote(paste("SqrtPGM"," (",theta[59]," = ",.(true_Omega_SqrtPGM[5,9]),")")),cex.lab = 1.5, cex.main = 2.5)
#curve(dnorm(x, mean=true_Omega_SqrtPGM[5,9], sd=sd59_SqrtPGM_true),
#      col="red", lwd=4, add=TRUE, yaxt="n")
