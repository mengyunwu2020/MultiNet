
turn.symmetric<-function(matrix,diag=F){
  temp=(matrix+t(matrix))/2
  if(diag) diag(temp)=1
  return(temp)
}
