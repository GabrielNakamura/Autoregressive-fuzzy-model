matrix.p<-function(comm,phylodist)
{
  matrix.w <- as.matrix(comm)
  phylodist <- as.matrix(phylodist)
  similar.phy <- 1 - (phylodist/max(phylodist))
  matrix.phy <- 1/colSums(similar.phy)
  matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
  matrix.P <- matrix.w %*% matrix.q
  return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
}

matrix.double.center <- function(mat){
  mean_row_partial <- apply(mat, MARGIN = 2, 
                            function(x){
                              x - rowMeans(mat)
                            }
  )
  mean_col_partial <- t(apply(mean_row_partial, MARGIN = 1, 
                              function(x){
                                x - colMeans(mat)
                              }))
  double_center_matrix <- mean_col_partial + mean(as.matrix(mat))
  return(double_center_matrix)
}

#' Metacommunity simulation with different relationship between phylogenetic signal and environmental variable
#'
#' @param tree A phylogenetic tree in \code{phylo} format
#' @param Nspp Scalar. The number of species that will be simulated in the metacommunities
#' @param Ncomm Scalar. The number of communities that will be simulated in the metacommunities
#' @param perm Scalar. The number of permutations to be used to calculate power and type I error rates
#' @param power Scalar. Power parameter to compute the degree of phylogenetic signal in the simulated trait. The same parameter to be passed to \code{\link[ape]{rTraitCont}} function
#' @param u Scalar. Mean value to be used in rnorm function
#' @param envir Logical. If TRUE an environmental vector with normal distribution will be generated and used in simulations
#' @param binary Logical. If TRUE (default) the community matrix will be transformed to presence and absence
#' @param test Logical. If TRUE (default) it will be computed power and type I error rates for beta parameter
#'
#' @return a list with all parameters returned from the simulation model 
#' 
#' @export
#'
#' @examples
simul_metacomm_autoregress<- function(tree,
                                      Nspp, 
                                      Ncomm, 
                                      perm, 
                                      power, 
                                      u, 
                                      envir= TRUE,
                                      binary= TRUE, 
                                      test= TRUE,
                                      parallel = ){
 
  a<- rnorm(Nspp,u,10)
  h<- runif(Nspp,0,30)
  L<- matrix(NA, Ncomm, Nspp)
  L.cent<- matrix(NA, Ncomm, Nspp)
  P.cent<- matrix(NA, nrow= Ncomm, ncol= Nspp)
  resid.L.cent<- matrix(NA, nrow= Ncomm, ncol= Nspp)
  N<- matrix(NA ,nrow=Nspp, ncol= 1, dimnames=list(1:Nspp, "niche_pos"))
  E<- matrix(NA, nrow= Ncomm, ncol= 1, dimnames=list(1:Ncomm, "env"))
  rownames(N)<- tree$tip.label
  niche<- ape::rTraitCont(ape::compute.brlen(tree,power=power),model="BM")
  N<- scales::rescale(niche,c(-1,101))
  K.N<- picante::Kcalc(N, tree)
  Res<- matrix(NA, nrow= 1, ncol=11, 
               dimnames= list(1, 
                              paste(c("R2.W","Coef.a.W","Coef.b.W","Beta.W",
                                      "Value.PCPS.1","Value.PCPS.2","R2.Res","Coef.a.Res",
                                      "Coef.b.Res","Beta.Res","p.value")
                              )
               )
  )
  if(envir == TRUE){
    E <- runif(Ncomm,0,100)
  } else {E = N
  }
  for(i in 1:Ncomm){
    for (j in 1:Nspp){
      L[i,j]<- rpois(1,h[j]*exp((-((E[i]-N[j])^2))/(2*(a[j]^2))))
    }
  }
  colnames(L)<- tree$tip.label
  rownames(L)<- rownames(as.matrix(L), FALSE, prefix="comm")
  if(binary == TRUE){
    L<- vegan::decostand(L,"pa")
  } else{L=L}
  for (c in 1:nrow(L)){
    for (s in 1:ncol(L)){
      L.cent[c,s]<- L[c,s]-mean(L[c,])-mean(L[,s])+mean(L)
    }
  }
  pcps.L<- PCPS::pcps(L, cophenetic(tree))
  P<- matrix.p(L, cophenetic(tree))$matrix.P
  for (c in 1:nrow(P)){
    for (s in 1:ncol(P)){
      P.cent[c,s]<-P[c,s]-mean(P[c,])-mean(P[,s])+mean(P)
    }
  }
  mod.L <- lm(as.numeric(L.cent) ~ as.numeric(P.cent)) # linear model
  
  #### power/type I error test ####     
  if(test == TRUE){
    P.null.cent <- matrix(NA, nrow= Ncomm, ncol=Nspp)
    mod.beta <- QuantPsyc::lm.beta(mod.L)
    seqpermutation.taxa <- SYNCSA::permut.vector(ncol(L), nset = perm)
    seqpermutation.taxa <- lapply(seq_len(nrow(seqpermutation.taxa)), function(i) seqpermutation.taxa[i,])
    phylodist<- cophenetic(tree)
    p.n.taxa <- function(samp, comm, phylodist){
      MP.null <- matrix.p(comm, phylodist[samp, samp])$matrix.P
      return(MP.null)
    }
    P.null <- lapply(seqpermutation.taxa, p.n.taxa, comm = L, phylodist = phylodist)
    
    P.null.cent_list<- lapply(P.null, function(x){
      for (c in 1:nrow(x)){
        for (s in 1:ncol(x)){
          P.null.cent[c,s]<-x[c,s]-mean(x[c,])-mean(x[,s])+mean(x)
        }
      }
      P.null.cent
    })
    mod.L.null_list<- lapply(P.null.cent_list, function(x){
      mod.L.null_list<- lm(as.numeric(L.cent) ~ as.numeric(x))
      mod.L.null_list
    })
    mod.beta.null<- matrix(unlist(lapply(mod.L.null_list, function(x){
      QuantPsyc::lm.beta(x)
    })), nrow = perm, ncol= 1)
    
    p.value<- (sum(ifelse(mod.beta.null >= mod.beta, 1, 0)) +1)/(perm + 1)
  } else {p.value=NA
  }
  pred.L<-matrix(predict(mod.L), nrow(L), ncol(L), byrow=TRUE, dimnames=list(rownames(L),colnames(L)))
  resid.L<- L.cent-pred.L
  for (c in 1:nrow(resid.L)){
    for (s in 1:ncol(resid.L)){
      resid.L.cent[c,s]<-resid.L[c,s]-mean(resid.L[c,])-mean(resid.L[,s])+mean(resid.L)
    }
  }
  mod.resid.L<- lm(as.numeric(resid.L.cent) ~ as.numeric(P.cent))
  Res[ ,1]<- summary(mod.L)$r.squared
  Res[ ,2]<- summary(mod.L)$coefficients[1]
  Res[ ,3]<- summary(mod.L)$coefficients[2]
  Res[ ,4]<- QuantPsyc::lm.beta(mod.L)
  Res[ ,5]<- pcps.L$values[1,2]
  Res[ ,6]<- pcps.L$values[2,2]
  Res[ ,7]<- summary(mod.resid.L)$r.squared
  Res[ ,8]<- summary(mod.resid.L)$coefficients[1]
  Res[ ,9]<- summary(mod.resid.L)$coefficients[2]
  Res[ ,10]<- QuantPsyc::lm.beta(mod.resid.L)
  Res[ ,11]<- p.value
  list_res<- list(Tree= tree, Niche= N, K.niche= K.N, Environment= E, W.matrices= L,
                  P.matrices= P, Predicted.W.matrices= pred.L, 
                  Residual.W.matrices= resid.L,
                  Model.Results= Res)
  return(list_res)
}
