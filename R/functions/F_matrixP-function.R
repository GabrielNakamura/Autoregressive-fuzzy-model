
matrix.p <- function(L, phylo){
  match<-picante::match.phylo.comm(phylo,L)
  matrix.w <- as.matrix(match$comm)
  similar.phy<-ape::vcv(phy=match$phy,corr=TRUE)
  matrix.phy <- 1/colSums(similar.phy)
  matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
  if(diag==FALSE){
    diag(matrix.q)<-0
  } else { matrix.q=matrix.q
  }
  matrix.P <- matrix.w %*% matrix.q
  return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
}

