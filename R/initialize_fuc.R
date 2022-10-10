############################# Functions for generating the initial values ############################
initialize_fuc = function(data, K, n.start = 100){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: initialize_fuc
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Generating the initial values using K-means.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  Mu <- matrix(0, K, p)
  kmeans.clust <- kmeans(data, K, nstart = n.start)
  memb <- kmeans.clust$cluster
  prob <- kmeans.clust$size/n
  Theta <- array(0, dim = c(p, p, K))
  S <- array(0, dim = c(p, p, K))
  for(k in 1:K)
  {
    Mu[k,] <- t(colMeans(data[memb == k, , drop = FALSE]) )
    S[,,k]  <- cov(data[memb == k, , drop = FALSE])
    
    if(det(S[,,k])<1e-4) S[,,k]=S[,,k]+diag(0.01,p,p)
    Theta[,,k] <- solve(S[,,k])
  }
  
  int.res <- list()
  int.res$prob <-  prob
  int.res$Mu <-  Mu
  int.res$Theta <- Theta
  int.res$S <- S
  int.res$memb <- memb
  return(int.res)
}