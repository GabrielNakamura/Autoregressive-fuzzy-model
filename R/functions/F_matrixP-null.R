
p.n.taxa <- function(samp, L, phylo){
  L.null<-L
  colnames(L.null)<-col.names[samp]
  match.null<-picante::match.phylo.comm(phylo,L.null)
  MP.null <- matrix.p(match.null$comm,match.null$phy)$matrix.P
  return(MP.null)
}
