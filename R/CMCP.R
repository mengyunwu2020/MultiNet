
CMCP<-function(X,y,lambda_input,gamma,index,standsigma,nK){
  X=as.matrix(X)
  y=as.vector(y)
  res=grpreg(X, y, group=index, penalty="cMCP",
               family="gaussian",nlambda = 2,lambda.min = lambda_input,lambda =lambda_input,
               gamma=gamma)

  bb=res$beta[-1]
  return(bb)
}
