
nodewise_cmcp<-function(data, K, nrow,p,stand_sigma,sqrt_L,index,lambdas,gamma,inverse_intercept,nK,numthreads){



  weightX=matrix(0,nrow,(K*p))

  XX=matrix(0,(nrow*K),(K*(p-1)))


  Y=matrix(0,(nrow*K),1)


  l=1;
  for(k in 1:K){
    for(j in 1:p){
      weightX[,l]=sqrt_L[,k]*data[,j]
      l=l+1
    }
  } 


  a=foreach(i=1:p,.packages = c('grpreg'),.inorder=TRUE,.combine = cbind,  .export = c("CMCP","K","p","weightX","nrow","stand_sigma","sqrt_L","index","lambdas","gamma","inverse_intercept","nK"))%dopar% {


    XX=matrix(0,(nrow*K),(K*(p-1)))


    Y=matrix(0,(nrow*K),1)

    for(k in 1:K){
      tmp=weightX[1:nrow,((k-1)*p+1):(k*p)]
      XX[((k-1)*nrow+1):(k*nrow), ((k-1)*(p-1)+1):((k)*(p-1))]=tmp[1:nrow,-i]

      Y[((k-1)*nrow+1):(k*nrow),1]=tmp[1:nrow,i]*stand_sigma[(k),i]-inverse_intercept[k,i]*sqrt_L[,k]
    }

    out= CMCP(XX,Y,lambdas,gamma,index,stand_sigma,nK)


    temp= XX%*%matrix(out,ncol = 1)
    res=c(out,as.vector(temp))
    res
  } 
  return(a)
}


