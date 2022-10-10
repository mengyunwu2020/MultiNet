#' Estimation of multiple networks with common structures in heterogeneous subgroups (MultiNet).
#'@details
#' This function estimates multiple networks based on SQRT-sparse group lasso.
#'
#' @param data An n by p data matrix.
#' @param K Number of clusters.
#' @param lambda.lasso Penalty parameter for lasso term.
#' @param lambda.similar Penalty parameter for group item.
#' @param eps Tolerance for the EM algorithm. The default value is 1e-3.
#' @param niter Maximum number of iterations.
#' @param initialization The method of initialization. Can be "k-means". Otherwise, it must be initialized with the initialize parameter.
#' @param initialize A list, which includes the initialization values of Theta, Mu and prob.
#' @param beta (p-1)*K by p Matrix of initialized regression coefficients. Default is NULL.
#' @param traces Whether to trace intermediate results.
#' @param tau.theta Threshold parameter for precision matrix.
#' @param refit Whether to perform the refit process after the EM algorithm converges.
#' @param n.start The parameter in k-means initialization, indicating how many random sets should be selected.
#' @param re The seed used in the initialization step.
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
#' @importFrom igraph ba.game
#' @importFrom igraph get.edgelist
#' @export
#'
#' @examples
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
#' opt=sqrt(log(p-1))/sqrt(nrow(whole.data$data)*K0)
#' opt2=sqrt(log((p-1)*K0))/sqrt(nrow(whole.data$data)*K0)
#' lambda$lambda.similar=0.05*opt
#' lambda$lambda.lasso=0.1*opt2
#' res<-MultiNet(whole.data$data,K0,lambda$lambda.lasso,lambda$lambda.similar)


MultiNet = function(data, K, lambda.lasso, lambda.similar,
                eps = 1e-3, niter = 100, initialization="K-means", initialize,
                beta=NULL,traces=FALSE,tau.theta=1e-3,refit=TRUE,
                n.start=NULL,re=1){

  n <- as.integer(dim(data)[1])
  p <- as.integer(dim(data)[2])

  initial.selection=initialization



  if(is.null(beta))  beta=matrix(0,(p-1)*K,p)
  if(is.null(n.start)) n.start=K


    index=rep(1:(p-1),K)

    ord <- order(index)
    unOrd <- match(1:length(ord),ord)
    groupLen=rep(K,p-1)
    rangeGroupInd=c(0:(p-1))*K


    if(initial.selection=="K-means"){
      set.seed(re)
      out.initial = initialize_fuc(data,K,n.start)
      prob = out.initial$prob
      mu = out.initial$Mu
      Theta = out.initial$Theta
      memb = out.initial$memb
      L.mat = matrix(0,n,K)
      for(jj in 1:n) L.mat[jj, memb[jj]]=1
    } else {
      Theta = initialize$Theta
      mu = initialize$Mu
      prob = initialize$prob
      L.mat = initialize$L.mat
      if(is.null(L.mat)){
        memb=initialize$memship_ini
      }else{
        memb = apply(L.mat,1,which.max)
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
    residual=matrix(tmp,K,p,byrow = TRUE)

    cc_hat=matrix(0,K,p)
    det_vec=rep(0,K)
    for(kk in 1:K){
      det_vec[kk]=det(Theta[,,kk])
      for(j in 1:p){
        temp= -Theta[j,-j,kk]*residual[kk,j]
        cc_hat[kk,j]=mu[kk,j]-sum(temp*mu[kk,-j])
      }
    }



    while(t < niter)
    {
      prob.old = prob
      mu.old = mu
      Theta.old = Theta
      L.mat.old = L.mat
      nK.old = nK
      A.hat.old=A.hat
      ptmp = c()

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

      tmp=sqrt(L.mat)


      mis=matrix(1,K,p)
      mi=sqrt(mis)

      noder=nodewise(data,K,nrow=n,p,num_relaxation_round=5,mi,max_iter=5,innerIter=5,tmp,
                     index,
                     lambda.lasso,
                     lambda.similar,
                     nlam=1,
                     groupLen,
                     rangeGroupInd,
                     beta=beta,reorderfun,tau=0,
                     standard=0,
                     cc_hat)


      nodewiseres=noder$coeflist
      nodewiseresidual=noder$residual


      coef.m=array(0,dim=c(p,p,K))
      for(i in 1:p){
        tmp=setdiff(1:p,i)
        thresbeta=as.numeric(nodewiseres[[i]])
        beta[,i]=thresbeta
        coef.m[tmp,i,]=thresbeta[unOrd]

        for(kk in 1:K){
          wtm=nodewiseresidual[[i]][(n*(kk-1)+1):(n*kk)]*sqrt(mis[kk,i])
          residual[kk,i]= sum((wtm)^2)/nK[kk]
          Theta[i,i,kk] =1/residual[kk,i]
          Theta[tmp,i,kk]= -coef.m[tmp,i,kk]/residual[kk,i]
          cc_hat[kk,i]=sum(sqrt(L.mat[,kk])*wtm)/nK[kk]+cc_hat[kk,i]
          atmp=thresbeta[unOrd]
          mu[kk,i]=cc_hat[kk,i]+sum(atmp[((p-1)*(kk-1)+1):((p-1)*kk)]*mu[kk,-i])
        }


      }

      for(kk in 1:K){
          Theta[,,kk]=(Theta[,,kk] + t(Theta[,,kk]))/2
          if(det(Theta[,,kk])<0)
            warning('the current precision matrix is not definite!')
      }

      A.hat=apply(coef.m,3,turn.symmetric)
      A.hat=(abs(A.hat)>0)*1
      A.hat=array(A.hat,dim=c(p,p,K))
      for(k in 1:K){
        diag(A.hat[,,k])=1
      }
      det_vec=rep(0,K)
      for (k.ind in 1:K) {
        L_ikx = sqrt(L.mat[,k.ind])*t(t(data) - mu[k.ind,])
        Stemp= t(L_ikx) %*% L_ikx / nK[k.ind]
        tmp=A.hat[,,k.ind]
        tmp[upper.tri(tmp)]=1
        zero=which(tmp==0,arr.ind = TRUE)

        if (dim(zero)[1]== p*(p-1)/2){
          warning("One or more matrices are constrained to be zero")
        }

        set.seed(1)

        fit <- glasso(Stemp, rho = 0, zero = zero, penalize.diagonal=FALSE, maxit = 30)

        Theta[,,k.ind]= (fit$wi + t(fit$wi))/2


        det_vec[k.ind]=det(Theta[,,k.ind])


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


    if(refit){
      restmp<-refitGMM(data,Theta,mu,L.mat,A.hat)
      Theta=restmp$Theta
      mu=restmp$mu
    }

    Theta[abs(Theta) < tau.theta] <- 0



  member = apply(L.mat,1,which.max)

  return(list(mu=mu, Theta= Theta,prob=prob, L.mat= L.mat,member=member,beta=beta))
}








