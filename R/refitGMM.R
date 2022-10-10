#' Refit function
#'
#' @param X An n by p data matrix.
#' @param Theta A p-by-p-by-K array of K precision matrices.
#' @param mu A K by p matrix of K cluster mean.
#' @param L.mat A n by K matrix. Element (i, j) represents whether the i-th sample belongs to the j-th category.
#' @param A.hat A p-by-p-by-K array of K edge matrices.
#'
#' @return A list containing the precision matrices array and mean.
#' @export
#'

refitGMM<-function(X,Theta,mu,L.mat,A.hat){
  K=dim(A.hat)[3]
  members=apply(L.mat,1,which.max)
  p=dim(X)[2]
  for (k.ind in 1:K) {
    data1=matrix(X[which(members==k.ind),],ncol=p)

      if(dim(data1)[1]==1){
        mu[k.ind,]=data1

    }
    L_ikx = t(t(data1) - mu[k.ind,])
    Stemp= t(L_ikx) %*% L_ikx / dim(data1)[1]
    A.hat[,,k.ind]=(abs(Theta[,,k.ind])>0)*1
    tmp=A.hat[,,k.ind]
    tmp[upper.tri(tmp)]=1
    zero=which(tmp==0,arr.ind = T)

    if (dim(zero)[1]== p*(p-1)/2){
      warning("One or more matrices are constrained to be zero")
    }else if(dim(zero)[1]!=0){
        set.seed(1)
        fit <- glasso(Stemp, rho = 0, zero = zero, penalize.diagonal=FALSE, maxit = 50)
        Theta[,,k.ind]= (fit$wi + t(fit$wi))/2


    }
    # while (kappa(Stemp) > 1e+2){
    #   Stemp = Stemp + 0.05 * diag(p)
    # }

  }
  return(list(Theta=Theta,mu=mu))
}
