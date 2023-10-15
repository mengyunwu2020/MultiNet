#' Estimation of multiple networks with common structures in heterogeneous subgroups (MultiNet).
#'@details
#' This function estimates multiple networks based on SQRT-sparse group lasso.
#'
#' @param data An n by p data matrix.
#' @param K Number of clusters.
#' @param lambda_mcp A non-negative tuning parameters for the composite MCP.
#' @param gamma A non-negative regularization parameter for the composite MCP. Default is 3.
#' @param eps Tolerance for the EM algorithm. The default value is 1e-3.
#' @param niter Maximum number of iterations.
#' @param initialization The method of initialization. Can be "k-means". Otherwise, it must be initialized with the initialize parameter.
#' @param initialize A list, which includes the initialization values of Theta, Mu and prob.
#' @param beta (p-1)*K by p Matrix of initialized regression coefficients. Default is NULL.
#' @param traces Whether to trace intermediate results.
#' @param refit Whether to perform the refit process after the EM algorithm converges.
#' @param n.start The parameter in k-means initialization, indicating how many random sets should be selected.
#' @param re The seed used in the initialization step.
#' @param numthreads Number of threads to use.
#'
#'
#'
#' @return a list with entries
#' \item{mu}{Estimated mean matrix. }
#' \item{Theta}{Array of estimated precision matrices. }
#' \item{prob}{A vector of estimated probabilities. }
#' \item{L.mat}{n by K Matrix of estimated probability that each sample belongs to each subgroup. }
#' \item{member}{A vector representing the estimated membership of each sample. }
#' \item{beta}{(p-1)*K by p Matrix of estimated regression coefficients. }
#'
#' @import glasso
#' @import MASS
#' @import parallel
#' @import doParallel
#' @import grpreg
#' @import foreach
#' @importFrom igraph ba.game
#' @importFrom igraph get.edgelist
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 200              # The sample size of each subgroup
#' p <- 100              # The dimension of the precision matrix
#' K0 <- 3
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
#' set.seed(1)
#' whole.data <- generate.data(N,Mu0.list,Theta0.list,Sigma0.list)
#' lambda <- list()
#' opt2=sqrt(log((p-1)*K0))/sqrt(nrow(whole.data$data)*K0)
#' lambda$gamma=3
#' lambda$lambda_mcp=opt2
#' res<-MultiNet(whole.data$data,K0,lambda$lambda_mcp,lambda$gamma,numthreads=1)


MultiNet =function(data, K,lambda_mcp,gamma=3,
                   eps = 5e-2, niter = 20, initialization="K-means", initialize,
                   beta=NULL,traces=FALSE, refit=TRUE,n.start=NULL,
                   re=1,numthreads=1)

{

  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])


  initial.selection=initialization



  if(is.null(beta))  beta=matrix(0,(p-1)*K,p)
  if(is.null(n.start)) n.start=K
  index=rep(1:(p-1),K)

  if(initial.selection=="K-means"){
    set.seed(re)
    out.initial = initialize_fuc(data,K,n.start)
    prob = out.initial$prob
    mu = out.initial$Mu
    Theta = out.initial$Theta
    memb = out.initial$memb
    L.mat = matrix(0,n,K)
    for(jj in 1:n) L.mat[jj, memb[jj]]=1


    S=list()

    for(k in 1:K){
      S[[k]]  <- diag(1,p,p)

      Theta[,,k] <- diag(1,p,p)
    }
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
      Theta[,,k] <- diag(1,p,p)
    }
    L.mat = matrix(0,n,K)
    for(jj in 1:n) L.mat[jj, memb[jj]]=1

    mu = Mu
  }else {
    Theta = initialize$Theta
    mu = initialize$Mu
    prob = initialize$prob
    L.mat = initialize$L.mat
    if(is.null(L.mat)){
      memb=initialize$memship_ini
    }else{
      memb = apply(L.mat,1,which.max)
    }

    if(min(table(memb))<10){
      set.seed(re)
      out.initial = initialize_fuc(data,K,n.start)
      prob = out.initial$prob
      mu = out.initial$Mu
      Theta = out.initial$Theta
      memb = out.initial$memb
      L.mat = matrix(0,n,K)
      for(jj in 1:n) L.mat[jj, memb[jj]]=1


      S=list()

      for(k in 1:K){
        S[[k]]  <- diag(1,p,p)

        Theta[,,k] <- diag(1,p,p)
      }


    }
  }

  # EM algorithm
  t = 0
  diff_mu = 10
  diff_theta = 10
  A.hat=array(1,dim=c(p,p,K))
  nK= apply(L.mat,2,sum)
  tmp=c()
  for(kkk in 1:K){
    tmp=c(tmp,1/diag(Theta[,,kkk]))
  }
  stand_sigma=residual=matrix(tmp,K,p,byrow = T)
  inverse_intercept=matrix(0,K,p)
  det_vec=rep(0,K)
  for(kk in 1:K){
    det_vec[kk]=det(Theta[,,kk])
    while(det_vec[kk]<1e-4)
    {
      Theta[,,kk]=Theta[,,kk]+diag(0.01,p,p)
      det_vec[kk]=det(Theta[,,kk])
    }


    if(is.infinite(det_vec[kk])) det_vec[kk]=2^1023
    for(j in 1:p){

      stand_sigma[kk,j] <-1/sqrt(residual[kk,j])
      temp= -Theta[j,-j,kk]*residual[kk,j]
      inverse_intercept[kk,j]=(mu[kk,j]-sum(temp*mu[kk,-j]))*stand_sigma[kk,j]
    }
  }


  cl <- makeCluster(numthreads)
  registerDoParallel(cl)
  while(t < niter)
  {
    prob.old = prob
    mu.old = mu
    Theta.old = Theta
    L.mat.old = L.mat
    nK.old = nK
    A.hat.old=A.hat
    ptmp = c()

    # print(prob)
    tmpp=apply(mu.old, 1,function(x) (t(data) -x))
    Theta.old2=indtemp=list()
    for(kind in 1:K){
      temp=abs(Theta.old[,,kind])>0
      diag(temp)=0
      indtemp[[kind]]=which(colSums(temp)!=0)
      Theta.old2[[kind]]=Theta.old[indtemp[[kind]],indtemp[[kind]],kind]
    }

    ptmp = c()
    for(j in 1:K){
      probtmp=rep(0,n)
      for(i in 1:n){
        probtmp[i]=0
        for(l in 1:K){
          tm2=tmpp[((i-1)*p+1):(i*p),j]
          tm1=tmpp[((i-1)*p+1):(i*p),l]
          tmpind1=setdiff(1:p,indtemp[[j]])
          tmpind2=setdiff(1:p,indtemp[[l]])
          con=matrix(tm2[indtemp[[j]]],nrow=1)%*%Theta.old2[[j]]%*%matrix(tm2[indtemp[[j]]],ncol=1)+sum(diag(Theta.old[,,j])[tmpind1]*(tm2[tmpind1])^2)-matrix(tm1[indtemp[[l]]],nrow=1)%*%Theta.old2[[l]]%*%matrix(tm1[indtemp[[l]]],ncol=1)-sum(diag(Theta.old[,,l])[tmpind2]*(tm1[tmpind2])^2)
          tmp=exp(con*0.5)
          if(is.infinite(tmp)) tmp=1e100
          if(j==l){
            mml=1}else{mml=det_vec[l]/det_vec[j]}
          if(is.infinite(mml)) mml=1e100
          probtmp[i]=probtmp[i]+prob[l]*(mml)^0.5*tmp
        }
        probtmp[i]=prob[j]/probtmp[i]
      }
      ptmp = cbind(ptmp,probtmp)
    }
    L.mat = ptmp

    prob= colMeans(L.mat)
    if(sum(is.na(prob))>0){
      warning('something wrong with probability calculation!')
      break;
    }

    if(sum(prob)!=1) prob=prob/sum(prob)


    nK= apply(L.mat,2,sum)

    sqrt_L=sqrt(L.mat)


    set.seed(1)
    noder=nodewise_cmcp(data,K,nrow=n,p,stand_sigma,sqrt_L,index,lambdas=lambda_mcp,gamma=gamma,inverse_intercept,nK,numthreads)


    coef.m=array(0,dim=c(p,p,K))


 

    for(i in 1:p){
      tmp_ind=setdiff(1:p,i)
      thresbeta=as.numeric(noder[1:(K*(p-1)),i])
      beta[,i]=thresbeta


      a<-rep(1,K)
      b<-rep(1,K)
      c2<-rep(1,K)


      for(kk in 1:K){

        besm=noder[(K*(p-1)+1):(K*(p-1+n)),i]
        tmp2=sum(L.mat[,kk]*stand_sigma[kk,i]*data[,i])-sum(sqrt_L[,kk]*besm[(n*(kk-1)+1):(n*kk)])


        inverse_intercept[kk,i]=tmp2/nK[kk]

        tmp=(besm[(n*(kk-1)+1):(n*kk)]+sqrt(L.mat[,kk])*inverse_intercept[kk,i])*sqrt(L.mat[,kk])*data[,i]
        b[kk]=sum(tmp)
        a[kk]=sum(L.mat[,kk]*data[,i]^2)

        stand_sigma[kk,i]=(b[kk])/(2*a[kk])+(sqrt(b[kk]^2+4*a[kk]*nK[kk]))/(2*a[kk])

        Theta[i,i,kk] =stand_sigma[kk,i]^2

        mu[kk,i]=(inverse_intercept[kk,i]+sum(thresbeta[((p-1)*(kk-1)+1):((p-1)*kk)]*mu[kk,-i]))/stand_sigma[kk,i]


      }

      coef.m[tmp_ind,i,]=thresbeta/rep(stand_sigma[,i],each=p-1)

      Theta[tmp_ind,i,]= -coef.m[tmp_ind,i,]*(stand_sigma[,i]^2)

    }




    A.hat=array(NA,dim=c(p,p,K))


    for(kk in 1:K){
      Theta[,,kk]=(Theta[,,kk] + t(Theta[,,kk]))/2
      A.hat[,,kk]=(abs(Theta[,,kk])>0)*1
      diag(A.hat[,,kk])=0
    } 

    for(kk in 1:K){
      det_vec[kk]=det(Theta[,,kk])
      while(det_vec[kk]<1e-4)
      {
         
        Theta[,,kk]=Theta[,,kk]+diag(0.01,p,p)
        det_vec[kk]=det(Theta[,,kk])
      }
      
      
      if(is.infinite(det_vec[kk])) det_vec[kk]=2^1023
    }



    t = t + 1
    diff_mu = norm(mu.old-mu,type="2")/(norm(mu,type="2")+0.001)
    diff_theta = norm(Theta.old-Theta,type="2")/(norm(Theta,type="2")+0.001)
    diff_Ahat=sum(abs(A.hat.old-A.hat))
    if(traces){
      cat('the ', t,'th: diff_Ahat:',diff_Ahat,'diff_mu',diff_mu,'diff_theta',diff_theta,'\n')
    }
    if(max(diff_mu,diff_theta)<eps&diff_Ahat==0){
      break;
    }
  }
  stopCluster(cl)
  if(refit){
    restmp<-refitGMM(data,Theta,mu,L.mat,A.hat)
    Theta=restmp$Theta
    mu=restmp$mu
  }


  member = apply(L.mat,1,which.max)

  return(list(mu=mu, Theta= Theta,prob=prob, L.mat= L.mat,member=member,beta=beta))

}








