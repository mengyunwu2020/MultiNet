#'Generating the s-block power-law precision matrices.
#' @param p Dimensions of the precision matrix.
#' @param s The number of sub-networks.
#' @param umin The lower bound of non-zero elements on non-diagonal elements.
#' @param umax The upper bound of non-zero elements on non-diagonal elements.
#' @param num.differ number of sub-networks which are different across subgroups.
#'
#' @return A list including The precision matrices of three subgroups.
#'
#' @importFrom igraph ba.game
#' @importFrom igraph get.edgelist
#' @import MASS
#' @import stats
#' @importFrom Matrix bdiag
#' @export


Power.law.network = function(p,s=10,umin=0.4,umax=0.7,num.differ){
  m=2
  Dnum=3
  pp=p/s
  if(p%%s != 0){
    print("warning! Matrix dimensions cannot be rounded by sub-matrix dimensions.")
  }
  submatrix=submatrix2=submatrix3=list()

  DD=diag(1,pp,pp)
  for(im in (pp/2+1):pp){
    DD[im,im]=Dnum
  }

  for (ss in 1:num.differ) {
    g = ba.game(pp, m=m,directed = F)
    Eg= as.data.frame(get.edgelist(g))
    subi = diag(0,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
    # for (i in 1:pp) {

    temppp=subi+diag((abs(min(eigen(subi)$values))+0.2),pp,pp)
    subi = DD%*%temppp%*%DD
    # }
    submatrix[[ss]]=subi


    g = ba.game(pp, m=m,directed = F)
    Eg = as.data.frame(get.edgelist(g))
    subi = diag(0,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
    # for (i in 1:pp) {
    temppp=subi+diag((abs(min(eigen(subi)$values))+0.2),pp,pp)
    subi = DD%*%temppp%*%DD
    # }
    submatrix2[[ss]] = subi



    g = ba.game(pp, m=m,directed = F)
    Eg = as.data.frame(get.edgelist(g))
    subi = diag(0,pp)
    for (q in 1:dim(Eg)[1]) {
      i=Eg[q,1];j=Eg[q,2]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
    # for (i in 1:pp) {
    temppp=subi+diag((abs(min(eigen(subi)$values))+0.2),pp,pp)
    subi = DD%*%temppp%*%DD
    # }
    submatrix3[[ss]] = subi


  }
  if(num.differ!=s){
      for (ss in (num.differ+1):s){
        g = ba.game(pp, m=m,directed = F)
        Eg= as.data.frame(get.edgelist(g))
        subi = subi2= subi3= diag(0,pp)
        for (q in 1:dim(Eg)[1]) {
          i=Eg[q,1];j=Eg[q,2]
          ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
          subi[i,j]=ij;subi[j,i]=ij

          ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
          subi2[i,j]=ij;subi2[j,i]=ij

          ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
          subi3[i,j]=ij;subi3[j,i]=ij
        }
        # for (i in 1:pp) {
        temppp=subi+diag((abs(min(eigen(subi)$values))+0.2),pp,pp)
        subi = DD%*%temppp%*%DD
        temppp=subi2+diag((abs(min(eigen(subi2)$values))+0.2),pp,pp)
        subi2 = DD%*%temppp%*%DD
        temppp=subi3+diag((abs(min(eigen(subi3)$values))+0.2),pp,pp)
        subi3 = DD%*%temppp%*%DD
        # }
        submatrix[[ss]]=subi
        submatrix2[[ss]]=subi2
        submatrix3[[ss]]=subi3
      }


  }

  if(ss==1){
    A=submatrix[[1]]
  }else{
    A=submatrix[[1]]
    for (ss in 2:s) {
      A=bdiag(A,submatrix[[ss]])
    }
  }
  A = as.matrix(A)

  if(ss==1){
    A2=submatrix2[[1]]
  }else{
    A2=submatrix2[[1]]
    for (ss in 2:s) {
      A2=bdiag(A2,submatrix2[[ss]])
    }
  }
  A2 = as.matrix(A2)

  if(ss==1){A3=submatrix3[[1]]}
  else{A3=submatrix3[[1]]
  for (ss in 2:s) {

    A3=bdiag(A3,submatrix3[[ss]])
  }
  }
  A3 = as.matrix(A3)

  return(list(A1=A,A2=A2,A3=A3))
}

