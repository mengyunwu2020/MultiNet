reorderfun<-function(index,X){
  ord <- order(index)
  X <- X[,ord]
  return(X)#
}