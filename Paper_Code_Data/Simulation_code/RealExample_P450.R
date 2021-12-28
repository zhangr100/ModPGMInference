library(ModPGMInference)

Correlation_Scale_free <- function(X){
  degree <- apply(X,1,sum)
  prob <- table(degree)/sum(table(degree))
  node_degree <- as.numeric(names(prob))
  prob_degree <- as.numeric(prob)
  if(node_degree[1] == 0){
    node_degree <- node_degree[-1]
    prob_degree <- prob_degree[-1]
  } 
  correlation <- cor(log2(node_degree),log2(prob_degree))
  return(correlation)
}


######## Real example for P450 data ########

## Preprocessing
P450 <- read.table("RawCounts_subgene.txt", header = TRUE)
P450_gene <- P450[4:dim(P450)[1],]
P450_gene <- t(P450_gene)
P450_gene <- apply(P450_gene,2,as.numeric)
X <- Preprocess(X = P450_gene)

## "Ground truth" structure
Omega <- matrix(0,dim(P450_gene)[2],dim(P450_gene)[2])
colnames(Omega) <- colnames(P450_gene)
row.names(Omega) <- colnames(P450_gene)

## connection 1
Omega["EHHADH","ACSM3"] <- 1
Omega["ACSM3","EHHADH"] <- 1

## connection 2
Omega["ACSM3","CYP2C19"] <- 1
Omega["CYP2C19","ACSM3"] <- 1

## connection 3
Omega["CYP2C19","CYP2C9"] <- 1
Omega["CYP2C9","CYP2C19"] <- 1

## connection 4
Omega["CYP2C9","CYP3A7"] <- 1
Omega["CYP3A7","CYP2C9"] <- 1

## Connection 5
Omega["CYP3A4","CYP3A43"] <- 1
Omega["CYP3A43","CYP3A4"] <- 1

## Connection 6
Omega["CYP2B6","CYP2B7P1"] <- 1
Omega["CYP2B7P1","CYP2B6"] <- 1

## Connection 7
Omega["CYP2C9","CYP2C8"] <- 1
Omega["CYP2C8","CYP2C9"] <- 1

## Connection 8
Omega["CYP2A7","CYP2A13"] <- 1
Omega["CYP2A13","CYP2A7"] <- 1

## Connection 9
Omega["CYP2A7","CYP2A6"] <- 1
Omega["CYP2A6","CYP2A7"] <- 1

## Connection 10
Omega["CYP2C8","ZGPAT"] <- 1
Omega["ZGPAT","CYP2C8"] <- 1

## Connection 11
Omega["ZGPAT","GLYAT"] <- 1
Omega["GLYAT","ZGPAT"] <- 1

## Connection 12
Omega["GLYAT","NTHL1"] <- 1
Omega["NTHL1","GLYAT"] <- 1

## Connection 13
Omega["NTHL1","NR1I3"] <- 1
Omega["NR1I3","NTHL1"] <- 1

## Connection 14
Omega["AKR1D1","GLYAT"] <- 1
Omega["GLYAT","AKR1D1"] <- 1

## Connection 15
Omega["AKR1D1","CYP2C9"] <- 1
Omega["CYP2C9","AKR1D1"] <- 1

## Connection 16
Omega["SLC10A1","AKR1D1"] <- 1
Omega["AKR1D1","SLC10A1"] <- 1

## Connection 17
Omega["SLC10A1","FMO3"] <- 1
Omega["FMO3","SLC10A1"] <- 1

## Connection 18
Omega["FMO3","PGRMC1"] <- 1
Omega["PGRMC1","FMO3"] <- 1

## Connection 19
Omega["AKR1D1","HAAO"] <- 1
Omega["HAAO","AKR1D1"] <- 1

## Connection 20
Omega["HAAO","CYP4A11"] <- 1
Omega["CYP4A11","HAAO"] <- 1

## Connection 21
Omega["HAAO","CYP27A1"] <- 1
Omega["CYP27A1","HAAO"] <- 1

## Connection 22
Omega["AKR1D1","SLC16A2"] <- 1
Omega["SLC16A2","AKR1D1"] <- 1

## Connection 23
Omega["SLC16A2","CYP4F2"] <- 1
Omega["CYP4F2","SLC16A2"] <- 1

## Connection 24
Omega["SLC16A2","SLC27A5"] <- 1
Omega["SLC27A5","SLC16A2"] <- 1

## Connection 25
Omega["CLU","SLC27A5"] <- 1
Omega["SLC27A5","CLU"] <- 1

## Connection 26
Omega["CLU","NCOR1"] <- 1
Omega["NCOR1","CLU"] <- 1

## Connection 27
Omega["GLYAT","ETNK2"] <- 1
Omega["ETNK2","GLYAT"] <- 1

## Connection 28
Omega["ETNK2","HNF4A"] <- 1
Omega["HNF4A","ETNK2"] <- 1

## Connection 29
Omega["ETNK2","SELENBP1"] <- 1
Omega["SELENBP1","ETNK2"] <- 1

## Connection 30
Omega["ETNK2","NR1I2"] <- 1
Omega["NR1I2","ETNK2"] <- 1

## Connection 31
Omega["NCOA7","HNF4A"] <- 1
Omega["HNF4A","NCOA7"] <- 1

## Connection 32
Omega["ALDH1L1","HNF4A"] <- 1
Omega["HNF4A","ALDH1L1"] <- 1

## Connection 33
Omega["CYP8B1","SELENBP1"] <- 1
Omega["SELENBP1","CYP8B1"] <- 1

## Connection 34
Omega["NR1I2","PRDX6"] <- 1
Omega["PRDX6","NR1I2"] <- 1

## Connection 35
Omega["PRDX6","PPARG"] <- 1
Omega["PPARG","PRDX6"] <- 1

## Connection 36
Omega["NCOA7","BCL6"] <- 1
Omega["BCL6","NCOA7"] <- 1

## Connection 37
Omega["CEBPD","BCL6"] <- 1
Omega["BCL6","CEBPD"] <- 1

## Connection 38
Omega["ALDH1L1","FOXA2"] <- 1
Omega["FOXA2","ALDH1L1"] <- 1

## Connection 39
Omega["ALDH1L1","MAT1A"] <- 1
Omega["MAT1A","ALDH1L1"] <- 1

## Connection 40
Omega["MAT1A","CYP1A1"] <- 1
Omega["CYP1A1","MAT1A"] <- 1


## Implement our approach with TPGM, SPGM and SqrtPGM and GFC_L under normal assumption
result_TPGM <- ModPGMInference(x = X, model = "TPGM", D = apply(X,2,max), global = TRUE,  N = 10, delta_upper = 1, alpha = c(0.001,0.005,0.01,0.05,0.1,0.15), true_graph = Omega)
result_SPGM <- ModPGMInference(x = X, model = "SPGM", D_0 = rep(0,dim(X)[2]), D_1 = apply(X,2,max), global = TRUE, delta_upper = 10, N = 1, alpha = c(0.001,0.005,0.01,0.05,0.1,0.15), true_graph = Omega)
result_SqrtPGM <- ModPGMInference(x = X, model = "SqrtPGM", global = TRUE, delta_upper = 0.5, N = 20, alpha = c(0.001,0.005,0.01,0.05,0.1,0.15), true_graph = Omega)
result_GFC_L <- SILGGM::SILGGM(huge::huge.npn(log(P450_gene+1)), method = "GFC_L", global = TRUE, ndelta = 10, alpha = c(0.001,0.005,0.01,0.05,0.1,0.15), true_graph = Omega)


## Part I: Evaluation of scale-free topology for the overall inferred network structure

Correlation_TPGM <- sapply(result_TPGM$global_decision,Correlation_Scale_free)
Correlation_SPGM <- sapply(result_SPGM$global_decision,Correlation_Scale_free)
Correlation_SqrtPGM <- sapply(result_SqrtPGM$global_decision,Correlation_Scale_free)
Correlation_GFC_L <- sapply(result_GFC_L$global_decision,Correlation_Scale_free)


## Part II: Performance of identifying "ground truth" interactions with FDR control at level 0.001 
result_SPGM_FDR_0.001 <- result_SPGM$global_decision[[1]]  
result_GFC_L_FDR_0.001 <- result_GFC_L$global_decision[[1]]

pair <- reshape::melt(upper.tri(Omega))
pair <- pair[pair[,3] %in% TRUE,c(1,2)]
row.names(pair) <- NULL

global_decision_SPGM <- result_SPGM_FDR_0.001[upper.tri(result_SPGM_FDR_0.001)]
global_decision_GFC_L <- result_GFC_L_FDR_0.001[upper.tri(result_GFC_L_FDR_0.001)]

result_SPGM_FDR_0.001_table <- cbind(pair, global_decision_SPGM)
result_GFC_L_FDR_0.001_table <- cbind(pair, global_decision_GFC_L)

gene_symbol <- colnames(Omega)
result_SPGM_FDR_0.001_table[,1] <- gene_symbol[as.numeric(result_SPGM_FDR_0.001_table[,1])]
result_SPGM_FDR_0.001_table[,2] <- gene_symbol[as.numeric(result_SPGM_FDR_0.001_table[,2])]
result_GFC_L_FDR_0.001_table[,1] <- gene_symbol[as.numeric(result_GFC_L_FDR_0.001_table[,1])]
result_GFC_L_FDR_0.001_table[,2] <- gene_symbol[as.numeric(result_GFC_L_FDR_0.001_table[,2])]

colnames(result_SPGM_FDR_0.001_table)[c(1,2)] <- c("gene1","gene2")
colnames(result_GFC_L_FDR_0.001_table)[c(1,2)] <- c("gene1","gene2")

true_interaction <- Omega[upper.tri(Omega)]

result_SPGM_FDR_0.001_table <- cbind(result_SPGM_FDR_0.001_table,true_interaction)
result_GFC_L_FDR_0.001_table <- cbind(result_GFC_L_FDR_0.001_table,true_interaction)

## Overlapped gene interactions with the "ground truth"
result_SPGM_FDR_0.001_table_overlap <- result_SPGM_FDR_0.001_table[result_SPGM_FDR_0.001_table$global_decision_SPGM==1 & result_SPGM_FDR_0.001_table$true_interaction==1,]
result_GFC_L_FDR_0.001_table_overlap <- result_GFC_L_FDR_0.001_table[result_GFC_L_FDR_0.001_table$global_decision_GFC_L==1 & result_GFC_L_FDR_0.001_table$true_interaction==1,]