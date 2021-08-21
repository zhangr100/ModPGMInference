\name{ModPGM_true_sd}
\alias{ModPGM_true_sd}
\title{True Standard Deviation of Each Edge for Modified Poisson-Type Graphical Models}
\usage{
ModPGM_true_sd <- function(x, psi = NULL, model = "SqrtPGM", true_graph, D = NULL,
D_0 = NULL, D_1 = NULL)
}
\arguments{
  \item{x}{
An \code{n} by \code{p} data matrix (\code{n} is the sample size and \code{p} is the dimension).
}
  \item{psi}{
A \code{p}-length vector of user-supplied values of intercepts. The default is a vector of all 0s if no sequence is provided.
}
  \item{model}{
Specification of modified Poisson-type graphical models: \code{"TPGM"} for truncated Poisson, \code{"SPGM"} for sub-linear Poisson and \code{"SqrtPGM"} for square-root Poisson. The default value is \code{"SqrtPGM"}.
}
  \item{true_graph}{
A known true graph structure.
}
  \item{D}{
A pre-specified \code{p}-length vector of truncation levels for each data column for \code{"TPGM"}. All the values need to be positive. The default is a vector of maximum values for each column of the input data matrix.
}
  \item{D_0}{
A pre-specified \code{p}-length vector of lower-bound truncation levels for each data column for \code{"SPGM"}. The default is a vector with all 0s.
}
  \item{D_1}{
A pre-specified \code{p}-length vector of upper-bound truncation levels for each data column for \code{"SPGM"}. The default is a vector of maximum values for each column of the input data matrix.
}
}
\value{
A \code{p} by \code{p} matrix which encodes true standard deviations from all node pairs.
}
\description{
Returns true standard deviation of each \eqn{\theta_{ij}} for three modified Poisson-type graphical models if a true graph structure is known.
}

\references{
Zhang, R., Ren, Z., Celedon, J. C. and Chen, W. (2021). Inference of Large Modified Poisson-Type Graphical Models: Application to RNA-Seq Data in Childhood Atopic Asthma Studies. \emph{The Annals of Applied Statistics}, \bold{15}(2), 831-855.
}

\examples{
## Chain graph
set <- c(-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4)
p <- 10

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

## Generate samples
n <- 100
X <- ModPGMSampler(psi = rep(0,p), true_graph = Omega, model = "TPGM", D = rep(3,p),
nSample = n, burn_in = 5000)

## Obtain true standard deviation of each edge
true_sd <- ModPGM_true_sd(x = X, psi = rep(0,p), model = "TPGM", true_graph = Omega,
D = rep(3,p))
}
