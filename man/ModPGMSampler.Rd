\name{ModPGMSampler}
\alias{ModPGMSampler}
\title{Sample Generator for the Modified Poisson-Type Graphical Models}
\usage{
ModPGMSampler <- function(psi = NULL, true_graph, model = "SqrtPGM", D = NULL,
D_0 = NULL, D_1 = NULL, nSample, burn_in = NULL, thin = NULL)
}
\arguments{
  \item{psi}{
A \code{p}-length vector of user-supplied values of intercepts. The default is a vector of all 0s if no sequence is provided.
}
  \item{true_graph}{
A known true graph structure.
}
  \item{model}{
Specification of modified Poisson-type graphical models: \code{"TPGM"} for truncated Poisson, \code{"SPGM"} for sub-linear Poisson and \code{"SqrtPGM"} for square-root Poisson. The default value is \code{"SqrtPGM"}.
}
  \item{D}{
A pre-specified \code{p}-length vector of truncation levels for each data column for \code{"TPGM"}. All the values need to be positive. The default is a vector of all 3s when there is no specification.
}
  \item{D_0}{
A pre-specified \code{p}-length vector of lower-bound truncation levels for each data column for \code{"SPGM"}. The default is a vector with all 3s when there is no specification.
}
  \item{D_1}{
A pre-specified \code{p}-length vector of upper-bound truncation levels for each data column for \code{"SPGM"}. The default is a vector with all 6s when there is no specification.
}
  \item{nSample}{
The sample size \code{n}, or number of samples to be generated.
}
  \item{burn_in}{
The burn-in period of Gibbs sampler (or the number of samples to be discarded at the beginning). The default value is \code{1000}.
}
  \item{thin}{
The thinning degree of Gibbs sampler. By default \code{thin = 100} is used, which would result in keeping every 100th generated sample and discard all other samples.
}
}
\value{
An \code{n} by \code{p} data matrix.
}
\description{
Implements sample generation for three modified Poisson-type graphical models based on the Gibbs sampler.
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
}
