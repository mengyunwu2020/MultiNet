
f.den.vec = function(data, mu, Theta){
  
  p = length(mu)
  if(is.infinite(det(Theta))){
    tmp=eigen(Theta)$values
  }else{
    tmp=det(Theta)
  }
  fden = as.numeric(exp((-p/2)*log(2*pi)+1/2*sum(log(tmp))-1/2*diag(t(t(data) - as.numeric(mu)) %*% Theta %*% (t(data) - as.numeric(mu)))))
  return(fden)
}

