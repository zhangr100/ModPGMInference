library(ModPGMInference)

######## Real data analysis for childhood atopic asthma ########

## Preprocessing
load("ge_atopic_asthma_original.RData")
ge_atopic_asthma_process <- Preprocess(X = ge_atopic_asthma_original, log_power_trans_only = TRUE)
ge_atopic_asthma_process <- ge_atopic_asthma_process[, colnames(ge_atopic_asthma_original)]

## Implement our method using EBIC with FDR control at 0.001 level
## TPGM
D = apply(ge_atopic_asthma_process,2,max)
D = as.numeric(D)

fit_TPGM_realdata <- ModPGMInference(x = ge_atopic_asthma_process, model = "TPGM", D = D, global = TRUE, alpha = 0.001)

## SPGM
D1 = apply(ge_atopic_asthma_process,2,max)
D1 = as.numeric(D1)
D0 = rep(0,dim(ge_atopic_asthma_process)[2])

fit_SPGM_realdata <- ModPGMInference(x = ge_atopic_asthma_process, model = "SPGM", D_0 = D0, D_1 = D1, global = TRUE, alpha = 0.001)

## SqrtPGM
fit_SqrtPGM_realdata <- ModPGMInference(x = ge_atopic_asthma_process, model = "SqrtPGM", global = TRUE, alpha = 0.001)


## Implement our method using CV with FDR control at 0.001 level
## TPGM
D = apply(ge_atopic_asthma_process,2,max)
D = as.numeric(D)

fit_TPGM_realdata_CV <- ModPGMInference(x = ge_atopic_asthma_process, model = "TPGM", tuning = "CV", kfold = 10, D = D, global = TRUE, alpha = 0.001)

## SPGM
D1 = apply(ge_atopic_asthma_process,2,max)
D1 = as.numeric(D1)
D0 = rep(0,dim(ge_atopic_asthma_process)[2])

fit_SPGM_realdata_CV <- ModPGMInference(x = ge_atopic_asthma_process, model = "SPGM", tuning = "CV", kfold = 10, D_0 = D0, D_1 = D1,  global = TRUE, alpha = 0.001)

## SqrtPGM
fit_SqrtPGM_realdata_CV <- ModPGMInference(x = ge_atopic_asthma_process, model = "SqrtPGM", tuning = "CV", kfold = 10, global = TRUE, alpha = 0.001)


## Implement sole estimation with EBIC
## TPGM
D = apply(ge_atopic_asthma_process,2,max)
D = as.numeric(D)

fit_TPGM_realdata_EBIC_estimation <- ModPGMInference(x = ge_atopic_asthma_process, model = "TPGM", D = D)

## SPGM
D1 = apply(ge_atopic_asthma_process,2,max)
D1 = as.numeric(D1)
D0 = rep(0,dim(ge_atopic_asthma_process)[2])

fit_SPGM_realdata_EBIC_estimation <- ModPGMInference(x = ge_atopic_asthma_process, model = "SPGM", D_0 = D0, D_1 = D1)

## SqrtPGM
fit_SqrtPGM_realdata_EBIC_estimation <- ModPGMInference(x = ge_atopic_asthma_process, model = "SqrtPGM")


## Implement sole estimation with CV
## TPGM
D = apply(ge_atopic_asthma_process,2,max)
D = as.numeric(D)

fit_TPGM_realdata_CV_estimation <- ModPGMInference(x = ge_atopic_asthma_process, model = "TPGM", tuning = "CV", kfold = 10, D = D)

## SPGM
D1 = apply(ge_atopic_asthma_process,2,max)
D1 = as.numeric(D1)
D0 = rep(0,dim(ge_atopic_asthma_process)[2])

fit_SPGM_realdata_CV_estimation <- ModPGMInference(x = ge_atopic_asthma_process, model = "SPGM", tuning = "CV", kfold = 10, D_0 = D0, D_1 = D1)

## SqrtPGM
fit_SqrtPGM_realdata_CV_estimation <- ModPGMInference(x = ge_atopic_asthma_process, model = "SqrtPGM", tuning = "CV", kfold = 10)


p <- dim(fit_TPGM_realdata_EBIC_estimation$theta_initial)[2]
A <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    A[i,j] = (fit_TPGM_realdata_EBIC_estimation$theta_initial[i,j]+fit_TPGM_realdata_EBIC_estimation$theta_initial[j,i])/2
    A[j,i] = A[i,j]
  }
}
A[A!=0] = 1

node_degree_TPGM_EBIC_estimation <- apply(A,1,sum)

A <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    A[i,j] = (fit_SPGM_realdata_EBIC_estimation$theta_initial[i,j]+fit_SPGM_realdata_EBIC_estimation$theta_initial[j,i])/2
    A[j,i] = A[i,j]
  }
}
A[A!=0] = 1

node_degree_SPGM_EBIC_estimation <- apply(A,1,sum)

A <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    A[i,j] = (fit_SqrtPGM_realdata_EBIC_estimation$theta_initial[i,j]+fit_SqrtPGM_realdata_EBIC_estimation$theta_initial[j,i])/2
    A[j,i] = A[i,j]
  }
}
A[A!=0] = 1

node_degree_SqrtPGM_EBIC_estimation <- apply(A,1,sum)

node_degree_TPGM <- apply(fit_TPGM_realdata$global_decision[[1]],1,sum)
node_degree_SPGM <- apply(fit_SPGM_realdata$global_decision[[1]],1,sum)
node_degree_SqrtPGM <- apply(fit_SqrtPGM_realdata$global_decision[[1]],1,sum)

p <- dim(fit_TPGM_realdata_CV_estimation$theta_initial)[2]
A <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    A[i,j] = (fit_TPGM_realdata_CV_estimation$theta_initial[i,j]+fit_TPGM_realdata_CV_estimation$theta_initial[j,i])/2
    A[j,i] = A[i,j]
  }
}
A[A!=0] = 1

node_degree_TPGM_CV_estimation <- apply(A,1,sum)

p <- dim(fit_SPGM_realdata_CV_estimation$theta_initial)[2]
A <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    A[i,j] = (fit_SPGM_realdata_CV_estimation$theta_initial[i,j]+fit_SPGM_realdata_CV_estimation$theta_initial[j,i])/2
    A[j,i] = A[i,j]
  }
}
A[A!=0] = 1

node_degree_SPGM_CV_estimation <- apply(A,1,sum)

p <- dim(fit_SqrtPGM_realdata_CV_estimation$theta_initial)[2]
A <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    A[i,j] = (fit_SqrtPGM_realdata_CV_estimation$theta_initial[i,j]+fit_SqrtPGM_realdata_CV_estimation$theta_initial[j,i])/2
    A[j,i] = A[i,j]
  }
}
A[A!=0] = 1

node_degree_SqrtPGM_CV_estimation <- apply(A,1,sum)

node_degree_TPGM_CV <- apply(fit_TPGM_realdata_CV$global_decision[[1]],1,sum)
node_degree_SPGM_CV <- apply(fit_SPGM_realdata_CV$global_decision[[1]],1,sum)
node_degree_SqrtPGM_CV <- apply(fit_SqrtPGM_realdata_CV$global_decision[[1]],1,sum)

######################################################################

## Use power law to evaluate the overall network structure
## EBIC
table_TPGM <- table(node_degree_TPGM)
Prob_TPGM <- table_TPGM/sum(table_TPGM)
Prob_TPGM <- Prob_TPGM[-1]
degree_number_TPGM <- as.numeric(names(Prob_TPGM))
Corr_TPGM <- cor(log2(degree_number_TPGM),log2(Prob_TPGM))

table_SPGM <- table(node_degree_SPGM)
Prob_SPGM <- table_SPGM/sum(table_SPGM)
Prob_SPGM <- Prob_SPGM[-1]
degree_number_SPGM <- as.numeric(names(Prob_SPGM))
Corr_SPGM <- cor(log2(degree_number_SPGM),log2(Prob_SPGM))

table_SqrtPGM <- table(node_degree_SqrtPGM)
Prob_SqrtPGM <- table_SqrtPGM/sum(table_SqrtPGM)
Prob_SqrtPGM <- Prob_SqrtPGM[-1]
degree_number_SqrtPGM <- as.numeric(names(Prob_SqrtPGM))
Corr_SqrtPGM <- cor(log2(degree_number_SqrtPGM),log2(Prob_SqrtPGM))

Prob_TPGM <- as.numeric(Prob_TPGM)
Prob_SPGM <- as.numeric(Prob_SPGM)
Prob_SqrtPGM <- as.numeric(Prob_SqrtPGM)

table_TPGM_EBIC_estimation <- table(node_degree_TPGM_EBIC_estimation)
Prob_TPGM_EBIC_estimation <- table_TPGM_EBIC_estimation/sum(table_TPGM_EBIC_estimation)
Prob_TPGM_EBIC_estimation <- Prob_TPGM_EBIC_estimation[-1]
degree_number_TPGM_EBIC_estimation <- as.numeric(names(Prob_TPGM_EBIC_estimation))
Corr_TPGM_EBIC_estimation <- cor(log2(degree_number_TPGM_EBIC_estimation),log2(Prob_TPGM_EBIC_estimation))

table_SPGM_EBIC_estimation <- table(node_degree_SPGM_EBIC_estimation)
Prob_SPGM_EBIC_estimation <- table_SPGM_EBIC_estimation/sum(table_SPGM_EBIC_estimation)
Prob_SPGM_EBIC_estimation <- Prob_SPGM_EBIC_estimation[-1]
degree_number_SPGM_EBIC_estimation <- as.numeric(names(Prob_SPGM_EBIC_estimation))
Corr_SPGM_EBIC_estimation <- cor(log2(degree_number_SPGM_EBIC_estimation),log2(Prob_SPGM_EBIC_estimation))

table_SqrtPGM_EBIC_estimation <- table(node_degree_SqrtPGM_EBIC_estimation)
Prob_SqrtPGM_EBIC_estimation <- table_SqrtPGM_EBIC_estimation/sum(table_SqrtPGM_EBIC_estimation)
Prob_SqrtPGM_EBIC_estimation <- Prob_SqrtPGM_EBIC_estimation[-1]
degree_number_SqrtPGM_EBIC_estimation <- as.numeric(names(Prob_SqrtPGM_EBIC_estimation))
Corr_SqrtPGM_EBIC_estimation <- cor(log2(degree_number_SqrtPGM_EBIC_estimation),log2(Prob_SqrtPGM_EBIC_estimation))

Prob_TPGM_EBIC_estimation <- as.numeric(Prob_TPGM_EBIC_estimation)
Prob_SPGM_EBIC_estimation <- as.numeric(Prob_SPGM_EBIC_estimation)
Prob_SqrtPGM_EBIC_estimation <- as.numeric(Prob_SqrtPGM_EBIC_estimation)


## CV
table_TPGM_CV <- table(node_degree_TPGM_CV)
Prob_TPGM_CV <- table_TPGM_CV/sum(table_TPGM_CV)
Prob_TPGM_CV <- Prob_TPGM_CV[-1]
degree_number_TPGM_CV <- as.numeric(names(Prob_TPGM_CV))
Corr_TPGM_CV <- cor(log2(degree_number_TPGM_CV),log2(Prob_TPGM_CV))

table_SPGM_CV <- table(node_degree_SPGM_CV)
Prob_SPGM_CV <- table_SPGM_CV/sum(table_SPGM_CV)
Prob_SPGM_CV <- Prob_SPGM_CV[-1]
degree_number_SPGM_CV <- as.numeric(names(Prob_SPGM_CV))
Corr_SPGM_CV <- cor(log2(degree_number_SPGM_CV),log2(Prob_SPGM_CV))

table_SqrtPGM_CV <- table(node_degree_SqrtPGM_CV)
Prob_SqrtPGM_CV <- table_SqrtPGM_CV/sum(table_SqrtPGM_CV)
Prob_SqrtPGM_CV <- Prob_SqrtPGM_CV[-1]
degree_number_SqrtPGM_CV <- as.numeric(names(Prob_SqrtPGM_CV))
Corr_SqrtPGM_CV <- cor(log2(degree_number_SqrtPGM_CV),log2(Prob_SqrtPGM_CV))

Prob_TPGM_CV <- as.numeric(Prob_TPGM_CV)
Prob_SPGM_CV <- as.numeric(Prob_SPGM_CV)
Prob_SqrtPGM_CV <- as.numeric(Prob_SqrtPGM_CV)

table_TPGM_CV_estimation <- table(node_degree_TPGM_CV_estimation)
Prob_TPGM_CV_estimation <- table_TPGM_CV_estimation/sum(table_TPGM_CV_estimation)
degree_number_TPGM_CV_estimation <- as.numeric(names(Prob_TPGM_CV_estimation))
Corr_TPGM_CV_estimation <- cor(log2(degree_number_TPGM_CV_estimation),log2(Prob_TPGM_CV_estimation))

table_SPGM_CV_estimation <- table(node_degree_SPGM_CV_estimation)
Prob_SPGM_CV_estimation <- table_SPGM_CV_estimation/sum(table_SPGM_CV_estimation)
degree_number_SPGM_CV_estimation <- as.numeric(names(Prob_SPGM_CV_estimation))
Corr_SPGM_CV_estimation <- cor(log2(degree_number_SPGM_CV_estimation),log2(Prob_SPGM_CV_estimation))

table_SqrtPGM_CV_estimation <- table(node_degree_SqrtPGM_CV_estimation)
Prob_SqrtPGM_CV_estimation <- table_SqrtPGM_CV_estimation/sum(table_SqrtPGM_CV_estimation)
degree_number_SqrtPGM_CV_estimation <- as.numeric(names(Prob_SqrtPGM_CV_estimation))
Corr_SqrtPGM_CV_estimation <- cor(log2(degree_number_SqrtPGM_CV_estimation),log2(Prob_SqrtPGM_CV_estimation))

Prob_TPGM_CV_estimation <- as.numeric(Prob_TPGM_CV_estimation)
Prob_SPGM_CV_estimation <- as.numeric(Prob_SPGM_CV_estimation)
Prob_SqrtPGM_CV_estimation <- as.numeric(Prob_SqrtPGM_CV_estimation)


## Log2-log2 plots of node degree distributions for EBIC 
par(mfrow=c(2,3))
plot(log2(degree_number_TPGM),log2(Prob_TPGM),xlab="log2-degree",ylab="log2-probability",main="Proposed (TPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_TPGM),log2(Prob_TPGM)), col = "red", lwd = 4)
legend("topright", bty="n", legend=paste("Corr = ", 
                                         format(Corr_TPGM, digits=2)), cex = 2)
plot(log2(degree_number_SPGM),log2(Prob_SPGM),xlab="log2-degree",ylab="log2-probability",main="Proposed (SPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_SPGM),log2(Prob_SPGM)), col = "red", lwd = 4)
legend("topright", bty="n", legend=paste("Corr = ", 
                                         format(Corr_SPGM, digits=2)), cex = 2)
plot(log2(degree_number_SqrtPGM),log2(Prob_SqrtPGM),xlab="log2-degree",ylab="log2-probability",main="Proposed (SqrtPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_SqrtPGM),log2(Prob_SqrtPGM)), col = "red", lwd = 4)
legend("topright", bty="n", legend=paste("Corr = ", 
                                         format(Corr_SqrtPGM, digits=2)), cex = 2)
plot(log2(degree_number_TPGM_EBIC_estimation),log2(Prob_TPGM_EBIC_estimation),xlab="log2-degree",ylab="log2-probability",main="Estimation (TPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_TPGM_EBIC_estimation),log2(Prob_TPGM_EBIC_estimation)), col = "red", lwd = 4)
legend("topright", bty="n", legend=paste("Corr = ", 
                                         format(Corr_TPGM_EBIC_estimation, digits=2)), cex = 2)
plot(log2(degree_number_SPGM_EBIC_estimation),log2(Prob_SPGM_EBIC_estimation),xlab="log2-degree",ylab="log2-probability",main="Estimation (SPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_SPGM_EBIC_estimation),log2(Prob_SPGM_EBIC_estimation)), col = "red", lwd = 4)
legend("topright", bty="n", legend=paste("Corr = ", 
                                         format(Corr_SPGM_EBIC_estimation, digits=2)), cex = 2)
plot(log2(degree_number_SqrtPGM_EBIC_estimation),log2(Prob_SqrtPGM_EBIC_estimation),xlab="log2-degree",ylab="log2-probability",main="Estimation (SqrtPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_SqrtPGM_EBIC_estimation),log2(Prob_SqrtPGM_EBIC_estimation)), col = "red", lwd = 4)
legend("topright", bty="n", legend=paste("Corr = ", 
                                         format(Corr_SqrtPGM_EBIC_estimation, digits=2)), cex = 2)


## Log2-log2 plots of node degree distributions for CV
par(mfrow=c(2,3))
plot(log2(degree_number_TPGM_CV),log2(Prob_TPGM_CV),xlab="log2-degree",ylab="log2-probability",main="Proposed (TPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_TPGM_CV),log2(Prob_TPGM_CV)), col = "red", lwd = 4)
legend("topright", bty="n", legend=paste("Corr = ", 
                                         format(Corr_TPGM_CV, digits=2)), cex = 2)
plot(log2(degree_number_SPGM_CV),log2(Prob_SPGM_CV),xlab="log2-degree",ylab="log2-probability",main="Proposed (SPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_SPGM_CV),log2(Prob_SPGM_CV)), col = "red", lwd = 4)
legend("topright", bty="n", legend=paste("Corr = ", 
                                         format(Corr_SPGM_CV, digits=2)), cex = 2)
plot(log2(degree_number_SqrtPGM_CV),log2(Prob_SqrtPGM_CV),xlab="log2-degree",ylab="log2-probability",main="Proposed (SqrtPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_SqrtPGM_CV),log2(Prob_SqrtPGM_CV)), col = "red", lwd = 4)
legend("topright", bty="n", legend=paste("Corr = ", 
                                         format(Corr_SqrtPGM_CV, digits=2)), cex = 2)
plot(log2(degree_number_TPGM_CV_estimation),log2(Prob_TPGM_CV_estimation),xlab="log2-degree",ylab="log2-probability",main="Estimation (TPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_TPGM_CV_estimation),log2(Prob_TPGM_CV_estimation)), col = "red", lwd = 4)
legend("topleft", bty="n", legend=paste("Corr = ", 
                                        format(Corr_TPGM_CV_estimation, digits=2)), cex = 2)
plot(log2(degree_number_SPGM_CV_estimation),log2(Prob_SPGM_CV_estimation),xlab="log2-degree",ylab="log2-probability",main="Estimation (SPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_SPGM_CV_estimation),log2(Prob_SPGM_CV_estimation)), col = "red", lwd = 4)
legend("topleft", bty="n", legend=paste("Corr = ", 
                                        format(Corr_SPGM_CV_estimation, digits=2)), cex = 2)
plot(log2(degree_number_SqrtPGM_CV_estimation),log2(Prob_SqrtPGM_CV_estimation),xlab="log2-degree",ylab="log2-probability",main="Estimation (SqrtPGM)",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_SqrtPGM_CV_estimation),log2(Prob_SqrtPGM_CV_estimation)), col = "red", lwd = 4)
legend("topleft", bty="n", legend=paste("Corr = ", 
                                        format(Corr_SqrtPGM_CV_estimation, digits=2)), cex = 2)


######################################################################

## Gene modularity
## Our method with EBIC
TPGM_module <- fit_TPGM_realdata$global_decision[[1]]
SPGM_module <- fit_SPGM_realdata$global_decision[[1]]
SqrtPGM_module <- fit_SqrtPGM_realdata$global_decision[[1]]

library(igraph)

net_TPGM = graph.adjacency(adjmatrix = TPGM_module, mode= "undirected", weighted=NULL)  
net_SPGM = graph.adjacency(adjmatrix = SPGM_module, mode= "undirected", weighted=NULL)
net_SqrtPGM = graph.adjacency(adjmatrix = SqrtPGM_module, mode= "undirected", weighted=NULL)

sgc.net_TPGM = leading.eigenvector.community(net_TPGM)
sgc.net_SPGM = leading.eigenvector.community(net_SPGM)
sgc.net_SqrtPGM = leading.eigenvector.community(net_SqrtPGM)

## membership of network contains community assignments of the network nodes
membership_TPGM <- membership(sgc.net_TPGM)
membership_SPGM <- membership(sgc.net_SPGM)
membership_SqrtPGM <- membership(sgc.net_SqrtPGM)


p <- dim(fit_TPGM_realdata_EBIC_estimation$theta_initial)[2]
A <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    A[i,j] = (fit_TPGM_realdata_EBIC_estimation$theta_initial[i,j]+fit_TPGM_realdata_EBIC_estimation$theta_initial[j,i])/2
    A[j,i] = A[i,j]
  }
}
A[A!=0] = 1
TPGM_estimation_module <- A

A <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    A[i,j] = (fit_SPGM_realdata_EBIC_estimation$theta_initial[i,j]+fit_SPGM_realdata_EBIC_estimation$theta_initial[j,i])/2
    A[j,i] = A[i,j]
  }
}
A[A!=0] = 1
SPGM_estimation_module <- A

A <- matrix(0,p,p)
for(i in 1:(p-1)){
  for(j in (i+1):p){
    A[i,j] = (fit_SqrtPGM_realdata_EBIC_estimation$theta_initial[i,j]+fit_SqrtPGM_realdata_EBIC_estimation$theta_initial[j,i])/2
    A[j,i] = A[i,j]
  }
}
A[A!=0] = 1
SqrtPGM_estimation_module <- A


library(igraph)

net_TPGM_estimation = graph.adjacency(adjmatrix = TPGM_estimation_module, mode= "undirected", weighted=NULL)  
net_SPGM_estimation = graph.adjacency(adjmatrix = SPGM_estimation_module, mode= "undirected", weighted=NULL)
net_SqrtPGM_estimation = graph.adjacency(adjmatrix = SqrtPGM_estimation_module, mode= "undirected", weighted=NULL)

sgc.net_TPGM_estimation = leading.eigenvector.community(net_TPGM_estimation)
sgc.net_SPGM_estimation = leading.eigenvector.community(net_SPGM_estimation)
sgc.net_SqrtPGM_estimation = leading.eigenvector.community(net_SqrtPGM_estimation)

## membership of network contains community assignments of the network nodes
membership_TPGM_estimation <- membership(sgc.net_TPGM_estimation)
membership_SPGM_estimation <- membership(sgc.net_SPGM_estimation)
membership_SqrtPGM_estimation <- membership(sgc.net_SqrtPGM_estimation)

## GFC_L and nonparanormal SKEPTIC
library(huge)
ge_atopic_asthma_GGM <- huge::huge.npn(log(ge_atopic_asthma_original+1))

library(SILGGM)

fit_GFC_L_realdata <- SILGGM::SILGGM(x = ge_atopic_asthma_GGM, method = "GFC_L", alpha = 0.001, ndelta = 20)
GFC_L_module <- fit_GFC_L_realdata$global_decision[[1]]
net_GFC_L = graph.adjacency(adjmatrix = GFC_L_module, mode= "undirected", weighted=NULL)  

sgc.net_GFC_L = leading.eigenvector.community(net_GFC_L)

## membership of network contains community assignments of the network nodes
membership_GFC_L <- membership(sgc.net_GFC_L)


skeptic <- huge::huge.npn(ge_atopic_asthma_original, npn.func = "skeptic")
library(Matrix)
skeptic <- as.matrix(nearPD(skeptic)$mat)

p <- dim(skeptic)[1]
n <- dim(ge_atopic_asthma_original)[1]

fit_SKEPTIC_realdata <- huge::huge(skeptic, method = "glasso")
nlamda <- length(fit_SKEPTIC_realdata$lambda)

## EBIC criterion (minimize)
EBIC.score <- NULL
EBIC.gamma <- 0.5

for(i in 1:nlamda){
  EBIC.score[i] = -n*fit_SKEPTIC_realdata$loglik[i] + log(n)*fit_SKEPTIC_realdata$df[i] + 4*EBIC.gamma*log(p)*fit_SKEPTIC_realdata$df[i]
}

index <- which(EBIC.score == min(EBIC.score))

## Selected graph structure
SKEPTIC_select <- fit_SKEPTIC_realdata$icov[[index]]
diag(SKEPTIC_select) <- 0
SKEPTIC_select[SKEPTIC_select!=0] <- 1

SKEPTIC_module <- SKEPTIC_select
net_SKEPTIC = graph.adjacency(adjmatrix = SKEPTIC_module, mode= "undirected", weighted=NULL)  

sgc.net_SKEPTIC = leading.eigenvector.community(net_SKEPTIC)

## membership of network contains community assignments of the network nodes
membership_SKEPTIC <- membership(sgc.net_SKEPTIC)


######################################################################

## Gene interactions within the JAK-STAT signaling pathway

JAK_STAT_gene <- c("IL6R","IL6","LIFR","IL5RA","IL7","IL20","IL6ST","IL3RA","SOCS1","CSF3R","CSF3","CSF2RA")

fit_TPGM_matrix <- fit_TPGM_realdata$global_decision[[1]]
fit_SPGM_matrix <- fit_SPGM_realdata$global_decision[[1]]
fit_SqrtPGM_matrix <- fit_SqrtPGM_realdata$global_decision[[1]]
colnames(fit_TPGM_matrix) <- colnames(ge_atopic_asthma_original)
rownames(fit_TPGM_matrix) <- colnames(ge_atopic_asthma_original)
colnames(fit_SPGM_matrix) <- colnames(ge_atopic_asthma_original)
rownames(fit_SPGM_matrix) <- colnames(ge_atopic_asthma_original)
colnames(fit_SqrtPGM_matrix) <- colnames(ge_atopic_asthma_original)
rownames(fit_SqrtPGM_matrix) <- colnames(ge_atopic_asthma_original)

fit_TPGM_estimation_matrix <- TPGM_estimation_module
fit_SPGM_estimation_matrix <- SPGM_estimation_module
fit_SqrtPGM_estimation_matrix <- SqrtPGM_estimation_module
colnames(fit_TPGM_estimation_matrix) <- colnames(ge_atopic_asthma_original)
rownames(fit_TPGM_estimation_matrix) <- colnames(ge_atopic_asthma_original)
colnames(fit_SPGM_estimation_matrix) <- colnames(ge_atopic_asthma_original)
rownames(fit_SPGM_estimation_matrix) <- colnames(ge_atopic_asthma_original)
colnames(fit_SqrtPGM_estimation_matrix) <- colnames(ge_atopic_asthma_original)
rownames(fit_SqrtPGM_estimation_matrix) <- colnames(ge_atopic_asthma_original)

fit_GFC_L_matrix <- fit_GFC_L_realdata$global_decision[[1]]
fit_SKEPTIC_matrix <- SKEPTIC_module
colnames(fit_GFC_L_matrix) <- colnames(ge_atopic_asthma_original)
rownames(fit_GFC_L_matrix) <- colnames(ge_atopic_asthma_original)
colnames(fit_SKEPTIC_matrix) <- colnames(ge_atopic_asthma_original)
rownames(fit_SKEPTIC_matrix) <- colnames(ge_atopic_asthma_original)

TPGM_JAK_STAT_matrix <- fit_TPGM_matrix[JAK_STAT_gene, JAK_STAT_gene]
SPGM_JAK_STAT_matrix <- fit_SPGM_matrix[JAK_STAT_gene, JAK_STAT_gene]
SqrtPGM_JAK_STAT_matrix <- fit_SqrtPGM_matrix[JAK_STAT_gene, JAK_STAT_gene]
TPGM_estimation_JAK_STAT_matrix <- fit_TPGM_estimation_matrix[JAK_STAT_gene, JAK_STAT_gene]
SPGM_estimation_JAK_STAT_matrix <- fit_SPGM_estimation_matrix[JAK_STAT_gene, JAK_STAT_gene]
SqrtPGM_estimation_JAK_STAT_matrix <- fit_SqrtPGM_estimation_matrix[JAK_STAT_gene, JAK_STAT_gene]
GFC_L_JAK_STAT_matrix <- fit_GFC_L_matrix[JAK_STAT_gene, JAK_STAT_gene]
SKEPTIC_JAK_STAT_matrix <- fit_SKEPTIC_matrix[JAK_STAT_gene, JAK_STAT_gene]

g1 <- igraph::graph_from_adjacency_matrix(adjmatrix = TPGM_JAK_STAT_matrix,mode = "undirected")
g2 <- igraph::graph_from_adjacency_matrix(adjmatrix = SPGM_JAK_STAT_matrix,mode = "undirected")
g3 <- igraph::graph_from_adjacency_matrix(adjmatrix = SqrtPGM_JAK_STAT_matrix,mode = "undirected")
g4 <- igraph::graph_from_adjacency_matrix(adjmatrix = TPGM_estimation_JAK_STAT_matrix,mode = "undirected")
g5 <- igraph::graph_from_adjacency_matrix(adjmatrix = SPGM_estimation_JAK_STAT_matrix,mode = "undirected")
g6 <- igraph::graph_from_adjacency_matrix(adjmatrix = SqrtPGM_estimation_JAK_STAT_matrix,mode = "undirected")
g7 <- igraph::graph_from_adjacency_matrix(adjmatrix = GFC_L_JAK_STAT_matrix,mode = "undirected")
g8 <- igraph::graph_from_adjacency_matrix(adjmatrix = SKEPTIC_JAK_STAT_matrix,mode = "undirected")
coords <- igraph::layout.davidson.harel(g1)

par(mfrow=c(2,4))
igraph::plot.igraph(g1, vertex.color= NA, vertex.label.cex = 1.5, vertex.frame.color= NA, edge.color="red", edge.width = 4, layout = coords)
title("Proposed (TPGM)", cex.main = 2)
igraph::plot.igraph(g2, vertex.color= NA, vertex.label.cex = 1.5, vertex.frame.color= NA, edge.color="red", edge.width = 4, layout = coords)
title("Proposed (SPGM)", cex.main = 2)
igraph::plot.igraph(g3, vertex.color= NA, vertex.label.cex = 1.5, vertex.frame.color= NA, edge.color="red", edge.width = 4, layout = coords)
title("Proposed (SqrtPGM)", cex.main = 2)
igraph::plot.igraph(g7, vertex.color= NA, vertex.label.cex = 1.5, vertex.frame.color= NA, edge.color="red", edge.width = 4, layout = coords)
title("GFC_L", cex.main = 2)
igraph::plot.igraph(g4, vertex.color= NA, vertex.label.cex = 1.5, vertex.frame.color= NA, edge.color="red", edge.width = 4, layout = coords)
title("Estimation (TPGM)", cex.main = 2)
igraph::plot.igraph(g5, vertex.color= NA, vertex.label.cex = 1.5, vertex.frame.color= NA, edge.color="red", edge.width = 4, layout = coords)
title("Estimation (SPGM)", cex.main = 2)
igraph::plot.igraph(g6, vertex.color= NA, vertex.label.cex = 1.5, vertex.frame.color= NA, edge.color="red", edge.width = 4, layout = coords)
title("Estimation (SqrtPGM)", cex.main = 2)
igraph::plot.igraph(g8, vertex.color= NA, vertex.label.cex = 1.5, vertex.frame.color= NA, edge.color="red", edge.width = 4, layout = coords)
title("Nonparanormal SKEPTIC", cex.main = 2)


######################################################################

## Additional real data results

## Other two methods under normal assumption
X <- ge_atopic_asthma_GGM

## Additional method 1: c-level partial correlation
library(scalreg)

n = dim(X)[1]; p = dim(X)[2]
t0 = 2; tau = seq(0, 2, 0.01); smax = n / 2; lentau = length(tau); 
IndMatrix = matrix(1, p, p) - diag(rep(1, p))
Eresidual = matrix(0, n, p)
CoefMatrix = matrix(0, p, p - 1)
meanX = colMeans(X)
X = t(t(X) - meanX)
XS = matrix(0, n, p)
for (i in 1 : p){
  XS[, i] = X[, i] / sd(X[, i])
}

for (i in 1 : p){
  out = scalreg::scalreg(X = XS[, -i], y = X[, i], lam0 = sqrt(2 * 2.01 * log(p * (log(p))^(1.5) / sqrt(n)) / n))
  Eresidual[, i] = out$residuals
  CoefMatrix[i, ] = out$coefficients / apply(X[, -i], 2, sd)
  print(i)
}

CovRes = t(Eresidual) %*% Eresidual / n
Est = matrix(1, p, p)

for (i in 1 : (p - 1)){
  for (j in (i + 1) : p){
    temp = Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1]
    Est[i, j] = mean(temp) / sqrt(diag(CovRes)[i] * diag(CovRes)[j])
    Est[j, i] = Est[i, j]
  }
}

EstThresh = Est * ( abs(Est) >= (t0 * sqrt(log(p) / n) * IndMatrix) )
kappa = (n / 3) * mean( colSums(Eresidual^4) / (colSums(Eresidual^2))^2 )

resprop = list(); rejectprop = c();
for (i in 1 : lentau){
  Threshold = tau[i] * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
  SRec = 1 * (abs(Est) > Threshold); NoNSRec = 1 * (SRec == 0)
  resprop[[i]] = which(SRec == 1, arr.ind = TRUE)
  rejectprop = c(rejectprop, max(1, (sum(SRec) - p)))
}

FDPprop = 2 * p * (p - 1) * ( 1 - pnorm( tau * sqrt(log(p)) ) ) / rejectprop

alpha = 0.001
lenalpha = length(alpha)

FDPresprop = list();

for (i in 1 : lenalpha){
  if (sum(FDPprop <= alpha[i]) > 0) tauprop = min(c(2, tau[FDPprop <= alpha[i]]))
  if (sum(FDPprop <= alpha[i]) == 0) tauprop = 2
  Threshold = tauprop * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
  SRec = 1 * (abs(Est) > Threshold); NoNSRec = 1 * (SRec == 0)
  FDPresprop[[i]] = which(SRec == 1, arr.ind = TRUE)
}

diag(SRec) <- 0
c_level_select <- SRec
colnames(c_level_select) <- colnames(X)
row.names(c_level_select) <- colnames(X)

## Additional method 2: SPACE (Sparse PArtial Correlation Estimation)
alpha=1
l1=1/sqrt(n)*qnorm(1-alpha/(2*p^2))
result_SPACE <- space::space.joint(Y.m = X, lam1 = l1*n*1.61, lam2 = 0, weight = 2)

SPACE_select <- result_SPACE$ParCor
diag(SPACE_select) <- 0
SPACE_select[SPACE_select!=0] <- 1
colnames(SPACE_select) <- colnames(X)
row.names(SPACE_select) <- colnames(X)

## Overall network structure
node_degree_c_level <- apply(c_level_select,1,sum)
node_degree_SPACE <- apply(SPACE_select,1,sum)

table_c_level <- table(node_degree_c_level)
table_SPACE <- table(node_degree_SPACE)
Prob_c_level <- table_c_level/sum(table_c_level)
Prob_SPACE <- table_SPACE/sum(table_SPACE)
Prob_c_level <- Prob_c_level[-1]
Prob_SPACE <- Prob_SPACE[-1]
degree_number_c_level <- as.numeric(names(Prob_c_level))
degree_number_SPACE <- as.numeric(names(Prob_SPACE))
Prob_c_level <- as.numeric(Prob_c_level)
Prob_SPACE <- as.numeric(Prob_SPACE)

Corr_c_level <- cor(log2(degree_number_c_level),log2(Prob_c_level))
Corr_SPACE <- cor(log2(degree_number_SPACE),log2(Prob_SPACE))

par(mfrow = c(1,2))
plot(log2(degree_number_c_level),log2(Prob_c_level),xlab="log2-degree",ylab="log2-probability",main="c-level PC",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_c_level),log2(Prob_c_level)), col = "red", lwd = 4)
legend("topleft", bty="n", legend=paste("Corr = ", 
                                        format(Corr_c_level, digits=2)), cex = 2)
plot(log2(degree_number_SPACE),log2(Prob_SPACE),xlab="log2-degree",ylab="log2-probability",main="SPACE",cex.lab=1.5,cex.axis=1.5,cex.main=2, pch = 19)
lines(lowess(log2(degree_number_SPACE),log2(Prob_SPACE)), col = "red", lwd = 4)
legend("topleft", bty="n", legend=paste("Corr = ", 
                                        format(Corr_SPACE, digits=2)), cex = 2)


## Results using TPM values
load("ge_atopic_asthma_tpm.RData")

GFC_L_tpm <- SILGGM::SILGGM(x = log(ge_atopic_asthma_tpm+1), method = "GFC_L",  alpha = 0.001, ndelta = 20)
GFC_L_module_tpm <- GFC_L_tpm$global_decision[[1]]

net_GFC_L_tpm <- graph.adjacency(adjmatrix = GFC_L_module_tpm, mode= "undirected", weighted=NULL)
sgc.net_GFC_L_tpm = leading.eigenvector.community(net_GFC_L_tpm)
membership_GFC_L_tpm <- membership(sgc.net_GFC_L_tpm)

JAK_STAT_gene <- c("IL6R","IL6","LIFR","IL5RA","IL7","IL20","IL6ST","IL3RA","SOCS1","CSF3R","CSF3","CSF2RA")
colnames(GFC_L_module_tpm) <- colnames(ge_atopic_asthma_tpm)
row.names(GFC_L_module_tpm) <- colnames(ge_atopic_asthma_tpm)

GFC_L_JAK_STAT_matrix_tpm <- GFC_L_module_tpm[JAK_STAT_gene, JAK_STAT_gene] 

g <- igraph::graph_from_adjacency_matrix(adjmatrix = GFC_L_JAK_STAT_matrix_tpm, mode = "undirected")

igraph::plot.igraph(g, vertex.color= NA, vertex.label.cex = 1.5, vertex.frame.color= NA, edge.color="red", edge.width = 4)
title("GFC_L", cex.main = 2)