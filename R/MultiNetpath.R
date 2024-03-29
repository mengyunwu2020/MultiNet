#' Estimation of multiple networks with common structures in heterogeneous subgroups (MultiNet) along path.
#'
#' @param data An n x p data matrix.
#' @param Kseq Vector of number of clusters.
#' @param lambda_vec A vector of non-negative tuning parameters for the composite MCP. It should be decreasing from large to small.
#' @param gamma A vector of non-negative regularization parameters for the composite MCP. Default is 3.
#' @param initial.selection The method of initialization. Can be "k-means". Otherwise, it must be initialized with the initialize parameter.
#' @param initialize A list, which includes the initialization values of Theta, Mu and prob.
#' @param eps Tolerance for the EM algorithm. The default value is 1e-3.
#' @param maxiter The maximum number of iterations.
#' @param trace.inter Whether to track progress in the main function.
#' @param refit Whether to perform the refit process after the EM algorithm converges.
#' @param tunoption One character. If tunoption=='pathwise', it means we used a path optimization scheme.
#' @param seed The seed used in the initialization step.
#' @param numthreads Number of threads to use.
#'
#' @return A list "resultall" with each element including:
#' \item{K}{Subgroups' number.}
#' \item{Mu_hat.list}{List of estimated subgroups' means.}
#' \item{Theta_hat.list}{List of estimated subgroups' precision matrices.}
#' \item{prob.list}{List of estimated subgroups' prior probabilities.}
#' \item{member.list}{List of estimated memberships of each subject.}
#' \item{L.mat.list}{List of estimated probabilities that each subject belongs to each subgroup.}
#'
#'
#' @import glasso
#' @import MASS
#' @importFrom igraph ba.game
#' @importFrom igraph get.edgelist
#' @export
#'
#' @examples
#' n <- 200              # The sample size of each subgroup
#' p <- 100              # The dimension of the precision matrix
#' K0 <- 3
#' set.seed(1)
#' mue <- .8
#' nonnum <- 2
#' mu01 <- c(rep(mue,nonnum),rep(-mue,nonnum),rep(-0,p-2*nonnum))
#' mu02 <- c(rep(mue,2*nonnum),rep(-0,p-2*nonnum))
#' mu03 <- c(rep(-mue,2*nonnum),rep(-0,p-2*nonnum))
#' num.differ=8
#' N <- rep(n,K0)
#' A.list <- Power.law.network(p,s=10,umin=0.3,umax=0.6,num.differ=num.differ)
#' Theta01 <- A.list$A1
#' Theta02 <- A.list$A2
#' Theta03 <- A.list$A3
#' sigma01 <- solve(Theta01)
#' sigma02 <- solve(Theta02)
#' sigma03 <- solve(Theta03)
#' Mu0.list <- list(mu01,mu02,mu03)
#' Sigma0.list <- list(sigma01,sigma02,sigma03)
#' Theta0.list <- list(Theta01,Theta02,Theta03)
#' set.seed(123)
#' whole.data <- generate.data(N,Mu0.list,Theta0.list,Sigma0.list)
#' lambda <- list()
#' opt2=sqrt(log((p-1)*K0))/sqrt(nrow(whole.data$data)*K0)
#' lambda$gamma=3
#' lambda$lambda_vec=seq(1*opt2,0.8*opt2,length.out=3)
#' res<-MultiNetpath(whole.data$data,K0,lambda$lambda_vec,lambda$gamma,trace.inter= FALSE,numthreads=1)




MultiNetpath = function(data, Kseq,lambda_vec,gamma,initial.selection="K-means",
                     initialize=NULL, eps = 1e-3, maxiter=100, trace.inter, refit=TRUE,tunoption='pathwise',
                     seed=1,numthreads=1){
  L1 = length(lambda_vec)
  L2 = length(gamma)
  L3 = length(Kseq)

  L = L1+L2+L3

  resultall=vector('list',L3)

  for(ppp in 1:L3){

    K=Kseq[ppp]
    cat('K=',K,'\n')
    n.start=K
    n_all = dim(data)[1]
    p=dim(data)[2]
    if(initial.selection=="K-means"){
      set.seed(seed)
      out.initial = initialize_fuc(data,K,n.start=n.start)
      memb = out.initial$memb
      L.mat = matrix(0,n_all,K)
      for(jj in 1:n_all) L.mat[jj, memb[jj]]=1
      out.initial$L.mat = L.mat
      p=dim(data)[2]
      S=list()

      for(k in 1:K){
        S[[k]]  <- diag(1,p,p)

        out.initial$Theta[,,k] <-diag(1,p,p)
      }

    } else if(initial.selection=='random'){
      if(!(is.null(initialize))){
        memship_ini=initialize$memb
      }
      set.seed(seed)
      memship_ini = sample(seq(1,K),n_all, replace =TRUE)
      clust_size = as.numeric(table(memship_ini))
      center_ini = matrix(0,nrow =K, ncol=dim(data)[2])
      prob =rep(0,K)
      for(l in 1:K){
        center_ini[l,] = apply(data[memship_ini==l,], 2, mean)
        prob[l]=sum(memship_ini==l)/n_all
      }

      p=dim(data)[2]
      Mu= center_ini # initialize mu
      Theta <- array(0, dim = c(p, p, K))
      S <- array(0, dim = c(p, p, K))
      for(k in 1:K)
      {
        S[,,k]  <- diag(1,p,p)

        Theta[,,k] <- diag(1,p,p)
      }

      out.initial<- list()
      out.initial$Mu <-  Mu
      out.initial$Theta <- Theta
      out.initial$S <- S
      out.initial$prob<-prob
      out.initial$memb <-memship_ini
      L.mat = matrix(0,n_all,K)
      for(jj in 1:n_all) L.mat[jj, memship_ini[jj]]=1
      out.initial$L.mat = L.mat

    }else if(initial.selection=="input_memb"){
      memb = initialize$memb
      n=length(memb)
      Mu = matrix(0,nrow =K, ncol=p)
      prob =rep(0,K)
      for(l in 1:K){
        Mu[l,] = apply(data[memb==l,], 2, mean)
        prob[l]=sum(memb==l)/n
      }


      Theta = array(NA,dim=c(p,p,K))
      S=list()

      for(k in 1:K){
        S[[k]]  <- diag(1,p,p)

        if(det(S[[k]])<1e-4) S[[k]]=S[[k]]+diag(0.01,p,p)
        Theta[,,k] <- solve(S[[k]])
      }
      L.mat = matrix(0,n,K)
      for(jj in 1:n) L.mat[jj, memb[jj]]=1

      out.initial=list()
      out.initial$L.mat=L.mat
      out.initial$Theta=Theta
      out.initial$prob=prob

      out.initial$Mu = Mu
    }else {
      out.initial = initialize
    }

    Mu_hat.list = list()
    Theta_hat.list = list()
    prob.list = list()
    L.mat.list = list()
    member.list = list()

    l=0
    lam=matrix(0,(L1*L2),2)
    for (ll in 1:L1) {
      for(jj in 1:L2){
          l=l+1
          lam1 = lambda_vec[ll]
          lam2=gamma[jj]

          cat('tuning parameter:',round(lam1,3),'is running\n')

          if(l==1) beta=NULL
          PP = MultiNet(data, K,lam1, lam2, eps = eps, niter=maxiter, initialization='input',
                     initialize=out.initial,beta=beta,traces=trace.inter, refit=refit,
                     n.start=n.start,re=seed,numthreads=numthreads)

          beta=PP$beta
          mu_hat=PP$mu;Theta_hat=PP$Theta;
          L.mat = PP$L.mat;prob = PP$prob;
          member = PP$member

          if(tunoption=='pathwise'){

          L.mat = matrix(0,n_all,K)
          for(pj in 1:n_all) L.mat[pj,  member[pj]]=1
          out.initial$L.mat = L.mat
          out.initial$prob <-  prob
          out.initial$Mu <-  mu_hat
          out.initial$Theta <- Theta_hat
          out.initial$memb <-member
          }


          Mu_hat.list[[l]]=mu_hat; Theta_hat.list[[l]]=Theta_hat; prob.list[[l]]=prob; L.mat.list[[l]] = L.mat; member.list[[l]]=member

          lam[l,]=c(lam1,lam2)

          resultall[[ppp]]=list(K=K, Mu_hat.list=Mu_hat.list,Theta_hat.list=Theta_hat.list,prob.list=prob.list,
                          member.list=member.list,L.mat.list=L.mat.list)


      }
    }


  }



  return(list(resultall=resultall))
}





