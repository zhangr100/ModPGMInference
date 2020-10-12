\name{ModPGMInference}
\alias{ModPGMInference}
\title{Statistical Inference of Modified Poisson-Type Graphical Models}
\usage{
ModPGMInference <- function(x, model = "SqrtPGM", D = NULL, D_0 = NULL, D_1 = NULL,
tuning = "EBIC", gamma = NA_real_, kfold = NA_real_, nlambda = NA_real_,
step_size = NA_real_, intercept = TRUE, global = FALSE, alpha = NULL,
regularization = NULL, N = NA_real_, delta_upper = NA_real_, true_graph = NULL)
}
\arguments{
  \item{x}{
An \code{n} by \code{p} data matrix (\code{n} is the sample size and \code{p} is the dimension, where \code{p} is allowed to be far larger than \code{n}).
}
  \item{model}{
Specification of modified Poisson-type graphical models: \code{"TPGM"} for truncated Poisson, \code{"SPGM"} for sub-linear Poisson and \code{"SqrtPGM"} for square-root Poisson. The default value is \code{"SqrtPGM"}.
}
  \item{D}{
A \code{p}-length vector of truncation levels for each column of the input data matrix for \code{"TPGM"}. All the values need to be positive. The default is a vector of maximum values for each column of the input data matrix.
}
  \item{D_0}{
A \code{p}-length vector of lower-bound truncation levels for each column of the input data matrix for \code{"SPGM"}. The default is a vector with all 0s.
}
  \item{D_1}{
A \code{p}-length vector of upper-bound truncation levels for each column of the input data matrix for \code{"SPGM"}. The default is a vector of maximum values for each column of the input data matrix.
}
  \item{tuning}{
Specification of methods for selection of tuning parameters when implementing the initial step of our approach: \code{"EBIC"} for the EBIC criterion and \code{"CV"} for the cross validation. The default value is \code{"EBIC"}.
}
  \item{gamma}{
The parameter in the EBIC criterion. The default value is \code{0.5}. ONLY applicable when \code{tuning = "EBIC"}.
}
  \item{kfold}{
Specification of fold numbers in the cross validation. The default value is \code{5}. ONLY applicable when \code{tuning = "CV"}.
}
  \item{nlambda}{
Number of tuning parameters for each node-wise regression when implementing the initial step of our approach. The default value is \code{20}.
}
  \item{step_size}{
A multiplicative parameter to decrease the step size during backtracking line search in the proximal gradient descent. Has to satisfy: \code{0 < step_size < 1}. The default value is \code{0.5}.
}
  \item{intercept}{
Should intercepts be estimated (\code{intercept = TRUE}) or set to 0s (\code{intercept = FALSE}) when implementing the initial step of our approach. The default is \code{TRUE}.
}
  \item{global}{
Should edge-wise (\code{global = FALSE}) or global statistical inference (\code{global = TRUE}) be performed based on the refined step of our approach. The default is \code{FALSE}.
}
  \item{alpha}{
A user-supplied sequence of pre-sepecified \code{alpha} levels for FDR control when \code{global = TRUE}. The default is \code{alpha = 0.05,0.1} if no sequence is provided.
}
  \item{regularization}{
A user-supplied sequence of tuning parameters when only sole estimation is performed.
}
  \item{N}{
A pre-specified value related to the number (e.g. \code{delta_upper*N}) of \eqn{\delta} values when selecting tuning parameters for global inference. The default value is \code{10}.
}
  \item{delta_upper}{
A pre-specified value for the upper-bound level of \eqn{\delta} (e.g. \eqn{0 \le \delta \le} \code{delta_upper}) when selecting tuning parameters for global inference. The default value is \code{2}.
}
  \item{true_graph}{
The true graph structure in a study if available. The default value is \code{NULL}. This argument is particularly for global inference. If a true graph is available, both FDR(s) and the corresponding power(s) will be provided in the outputs. Otherwise, only FDR(s) and the associated threshold(s) for all absolute values of test statistics will be provided.
}
}
\value{
A list is returned including:
  \item{intercept}{
A sequence of estimated intercepts from the initial step of our approach.
}
  \item{theta_initial}{
A matrix of estimated \eqn{\theta_{ij}}'s which depict conditional dependence of each \eqn{(i,j)^{th}} pair of nodes after only the initial step of our approach.
}
  \item{theta_cor}{
A matrix of estimated \eqn{\theta_{ij}}'s which depict conditional dependence of each \eqn{(i,j)^{th}} pair of nodes after the refined step of our approach for bias correction.
}
  \item{CI_low_theta}{
A matrix of lower values of 95\% confidence interval for each \eqn{\theta_{ij}} based on the estimated values in \code{theta_cor} after bias correction.
}
  \item{CI_high_theta}{
A matrix of higher values of 95\% confidence interval for each \eqn{\theta_{ij}} based on the estimated values in \code{theta_cor}.
}
  \item{z_score}{
A matrix of z-scores for each \eqn{\theta_{ij}} based on the estimated values in \code{theta_cor}.
}
  \item{p_thetaCor}{
A matrix of p-values for each \eqn{\theta_{ij}} based on the estimated values in \code{theta_cor}.
}
  \item{est_sd}{
A matrix of estimated standard deviations for each \eqn{\theta_{ij}} based on the estimated values in \code{theta_cor}.
}
  \item{threshold}{
The threshold sequence for absolute values of test statistics associated with the estimated FDR sequence.
}
  \item{FDR}{
The estimated FDR sequence for global inference of all pairs of \eqn{\theta_{ij}}'s  based on the pre-specified \code{alpha} level(s).
}
  \item{power}{
The estimated power sequence for global inference of all pairs of \eqn{\theta_{ij}}'s associated with the estimated FDR sequence. ONLY applicable if \code{true_graph} is available.
}
  \item{global_decision}{
A list of \code{p} by \code{p} adjacency matrices of inferred graphs under the global inference corrsponding to the sequence of pre-sepecified \code{alpha} levels. A value of \code{1} in the matrix means that there is conditional dependence (or an edge) between the node pair, while a value of \code{0} means condtional independence (or no edge).
} \cr

If \code{regularization} is available, only sole estimation is performed. Two lists named with \code{intercept} and \code{theta_initial} are returned, where each element records the estimated values corresponding to each user-defined tuning parameter with a descending order in the \code{regularization} sequence automatically.
}
\description{
The main function for both edge-wise statistical inference and multiple testing
with FDR control of three modified Poisson-type graphical models.
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

## Perform inference on random samples
result <- ModPGMInference(x = X, model = "TPGM", tuning = "EBIC", D = rep(3,p), nlambda = 100)
}
