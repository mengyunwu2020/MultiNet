f.den.vec = function(data, mu, Theta){
  p = length(mu)
  # if(p>300){
    ps=eigen(Theta)
    eigenval=ps$values

    fden = as.numeric( prod((2*pi/eigenval)^(-1/2)) * exp(-1/2*diag(t(t(data) - as.numeric(mu)) %*% Theta %*% (t(data) - as.numeric(mu)))) )


  # }else{
  #   fden = as.numeric( (2*pi)^(-p/2) * (det(Theta))^(1/2) * exp(-1/2*diag(t(t(data) - as.numeric(mu)) %*% Theta %*% (t(data) - as.numeric(mu)))) )
  #
  # }
  return(fden)
}




