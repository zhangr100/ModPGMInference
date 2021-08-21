# ModPGMInference
Statistical Inference of Modified Poisson-Type Graphical Models

## Description
Provides a novel method for both edge-wise and global statistical inference
with FDR control based on three modified Poisson-type graphical models: the truncated
Poisson graphical model (TPGM), the sub-linear Poisson graphical model (SPGM) and the
square-root Poisson graphical model (SqrtPGM). We focus on the high-dimensional settings
where dimension p is allowed to be far larger than sample size n and implement the method
using efficient proximal gradient descent. The method will be potentially useful in analysis
of large count-valued or discrete omics data (e.g. RNA-seq). Other functions in the package
including sample generation and data preprocessing for count-valued data. Windows users
should install 'Rtools' before the installation of this package.

## Installation
Windows users should install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first. `ModPGMInference` can be installed into `R` environment using the following steps:

* **Step 1: Install the `devtools` package**

```r
install.packages("devtools")
```

* **Step 2: Install the package from this repository**

```r
library(devtools)
install_github("zhangr100/ModPGMInference")
```

## Functions

There are 4 main functions in this package:

* `ModPGMInference`: Statistical Inference of Modified Poisson-Type Graphical Models.
* `ModPGMSampler`: Sample Generator for the Modified Poisson-Type Graphical Models.
* `ModPGM_true_sd`: True Standard Deviation of Each Edge for Modified Poisson-Type Graphical Models.
* `Preprocess`: Preprocessing the Count-Valued RNA-Seq Data.

## Examples

Generate random samples from Chain graph using TPGM.

```r
library(ModPGMInference)
set <- c(-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4)
p <- 20

Omega.tmp <- matrix(0,p,p)
for(i in 1:(p-1)){
  j <- i+1
  Omega.tmp[i,j] <- sample(set,1)
}
for(i in 1:(p-1)){
  for(j in (i+1):p){
    Omega.tmp[j,i] <- Omega.tmp[i,j]
  }
}
Omega <- Omega.tmp

n <- 100
X <- ModPGMSampler(psi = rep(0,p), true_graph = Omega, model = "TPGM", D = rep(3,p), nSample = n, burn_in = 5000)
```

Obtain true standard deviation of each edge.

```r
true_sd <- ModPGM_true_sd(x = X, psi = rep(0,p), model = "TPGM", true_graph = Omega, D = rep(3,p))
```

Perform statistical inference on random samples.

```r
result <- ModPGMInference(x = X, model = "TPGM", tuning = "EBIC", D = rep(3,p), nlambda = 100)
```

Simulate raw count-valued data and perform preprocessing.

```r
uniform <- matrix(0,n,p)
for(k in 1:n){
uniform[k,] <- runif(p,0,1)
}
X_new <- X + uniform

count_value <- exp(log(X_new)/0.2517)
count_value <- floor(count_value)

X_process <- Preprocess(X = count_value)
```

## Author(s)

Rong Zhang (<roz16@pitt.edu>), Zhao Ren (<zren@pitt.edu>), Juan C. Celed\'on (<juan.celedon@chp.edu>) and Wei Chen (<wei.chen@chp.edu>)

## Maintainer

Rong Zhang (<roz16@pitt.edu>)

## References

Zhang, R., Ren, Z., Celed\'on, J. C. and Chen, W. (2021). Inference of Large Modified Poisson-Type Graphical Models: Application to RNA-Seq Data in Childhood Atopic Asthma Studies. *The Annals of Applied Statistics.* **15**(2): 831-855. [Link](https://doi.org/10.1214/20-AOAS1413)
