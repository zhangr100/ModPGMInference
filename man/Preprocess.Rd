\name{Preprocess}
\alias{Preprocess}
\title{Preprocessing the Count-Valued RNA-Seq Data}
\usage{
Preprocess <- function(X, quanNorm = 0.75, nLowCount = 20, percentLowCount = 0.95,
NumGenes = 500, log_power_trans_only = FALSE)
}
\arguments{
  \item{X}{
An \code{n} by \code{p} raw count matrix.
}
  \item{quanNorm}{
A parameter to control the quantile normalization for the data. The default value is \code{0.75}.
}
  \item{nLowCount}{
The minimum count value across samples to filter out a gene. The default value is \code{20}.
}
  \item{percentLowCount}{
The percent of samples for a gene to be less than \code{nLowCount}. The default value is \code{0.95}.
}
  \item{NumGenes}{
Number of genes to keep in the data after preprocessing. The default value is \code{500}.
}
  \item{log_power_trans_only}{
Only log or power transform is performed and retains \code{NumGenes} genes with largest inter-sample variances if \code{log_power_trans_only = TRUE}. The default is \code{FALSE}.
}
}
\value{
An \code{n} by \code{NumGenes} or \code{p} or number of retained genes data matrix after preprocessing.
}
\details{
The preprocessing steps for the raw count-valued data include: 1. Quantile normalization for samples; 2. Filter out genes with low counts across all samples; 3. Retain genes with large inter-sample variances; 4. Use log or power transform which is selected based on Kolmogorov-Smirnov goodness of fit test to make the data closer to the Poisson distribution.
}
\description{
Performs preprocessing on raw count-valued data (e.g. RNA-seq).
}

\references{
Zhang, R., Ren, Z., Celedon, J. C. and Chen, W. (2021). Inference of Large Modified Poisson-Type Graphical Models: Application to RNA-Seq Data in Childhood Atopic Asthma Studies. \emph{The Annals of Applied Statistics}, \bold{15}(2), 831-855.
}

\examples{
## Chain graph
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

## Simulate a count-valued data set
n <- 100
X <- ModPGMSampler(psi = rep(0,p), true_graph = Omega, model = "TPGM", D = rep(3,p),
nSample = n, burn_in = 5000)

uniform <- matrix(0,n,p)
for(k in 1:n){
  uniform[k,] <- runif(p,0,1)
}

X_new <- X + uniform

count_value <- exp(log(X_new)/0.2517)

count_value <- floor(count_value)

## Pre-processing
X_process <- Preprocess(X = count_value)
}
