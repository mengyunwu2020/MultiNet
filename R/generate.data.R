#' Generating the simulated data.
#'
#' @param N @ N: K0 * 1 vector, the sample sizes of subgroups.
#' @param Mu0.list a list including K0 mean vectors (p * 1).
#' @param Theta0.list a list including K0 precision matrices (p * p).
#' @param Sigma0.list a list including K0 correlation matrices (p * p).
#'
#' @return A list "whole.data" including:
#' \item{L0}{n * 1 vector, the subgroup labels to which each sample belongs.}
#' \item{Mu0}{K0 * p matrix, K0 mean vectors.}
#' \item{Theta0}{K0 * p * p array, K0 precision matrices.}
#' \item{data}{n * p matrix, the design matrix.}
#' \item{n_all}{int, the total sample size.}
#' \item{K0}{int, the true number of subgroups.}
#' @export
#'
generate.data = function(N,Mu0.list,Theta0.list,Sigma0.list){
  K0 = length(Mu0.list)
  p = length(Mu0.list[[1]])
  Mu0=matrix(0,K0,p);L0=NULL;Theta0=array(0, dim = c(p, p, K0));data=NULL
  for (k in 1:K0) {
    Mu0[k,] <- Mu0.list[[k]]
    L0 <- c(L0,rep(k,N[k]))
  }
  for (k in 1:K0) {
    Theta0[,,k] <- as.matrix(Theta0.list[[k]])
  }
  for (k in 1:K0) {
    data <- rbind(data,mvrnorm(N[k],Mu0[k,],Sigma0.list[[k]]))
  }
  n_all = dim(data)[1]
  whole.data=list()
  whole.data$L0=L0
  whole.data$Mu0=Mu0
  whole.data$Theta0=Theta0
  whole.data$data=data
  whole.data$n_all=n_all
  whole.data$K0=K0
  return(whole.data)
}
