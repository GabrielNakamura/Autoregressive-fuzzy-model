matrix.p <- function(L, phylo, diag){
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

p.n.taxa <- function(samp, L, phylo, diag, col.names){
  L.null<-L
  colnames(L.null)<-col.names[samp]
  match.null<-picante::match.phylo.comm(phylo,L.null)
  MP.null <- matrix.p(match.null$comm,match.null$phy, diag)$matrix.P
  return(MP.null)
}

#' Metacommunity autoregressive model based on phylogenetic fuzzy sets
#'
#' @param comm Matrix. Rows are communities and columns species
#' @param phylo A phylo object containing the phylogeny of species in comm
#' @param binary Logical. If TRUE (default) the community matrix will be transformed to presence and absence
#' @param diag Logical. If TRUE the diagonal in matrix q is included, if FALSE the diagonal is removed
#' @param test Logical. If TRUE (default) it will be computed power and type I error rates for beta parameter
#' @param nperm Scalar. An integer indicating the number of permutations to be used to calculate power and type I error rates
#' @param parallel Scalar. An integer indicating the number of cores to be used in parallel computation. If NULL parallel computing won't be used
#'
#' @return A list with results from autoregressive model containing the following itens:
#' \itemize{
#' Matrix.P=P,Prediction=pred.L, Residuals= resid.L,Test.Results=Res
#' \item `Matrix.P` phylogenetic fuzzy matrix
#' \item `Prediction` predicted value for occurrence matrix
#' \item `Residuals` residuals obtained from original occurrence matrix and predicted occurrence matrix
#' \item `Test.Results` a matrix with parameters from ARM: R2, slope and p-value
#' }
#' @export
#'
#' @examples
Fuzzy_ARM <- function(comm,
                      phylo, 
                      binary=TRUE,
                      diag=FALSE,
                      test = TRUE, 
                      nperm = 1000, 
                      parallel = NULL){
  match<-picante::match.phylo.comm(phylo,comm)
  phylo<-match$phy
  comm<-match$comm
  
  Res <- matrix(NA,
                nrow = 1, 
                ncol = 3,
                dimnames = list(1, paste(
                  c("R2.W", "Beta.W", "P.value")
                )
                )
  )
  if(binary==TRUE){
    comm <- vegan::decostand(comm, "pa")
  } else{comm = comm}
  
  
  L.cent <- matrix.double.center(comm)
  
  P <- matrix.p(comm, phylo=phylo, diag = diag)$matrix.P
  P.cent<-matrix.double.center(P)
  mod.L<-lm(as.numeric(L.cent)~as.numeric(P.cent))
  pred.L <- matrix(predict(mod.L), nrow(comm), ncol(comm), byrow = TRUE, dimnames = list(rownames(comm), colnames(comm)))
  resid.L<- L.cent-pred.L
  Res[,1]<-summary(mod.L)$r.squared
  Res[,2]<- QuantPsyc::lm.beta(mod.L)
  #### Test == TRUE  
  if(test == TRUE){
    seqpermutation.taxa_perm <- SYNCSA::permut.vector(ncol(comm), nset = nperm)
    seqpermutation.taxa <- lapply(seq_len(nrow(seqpermutation.taxa_perm)), function(i) seqpermutation.taxa_perm[i,])
    col.names<-colnames(comm)
    
    P.null <- lapply(seqpermutation.taxa, p.n.taxa, L = comm, phylo = phylo, diag = diag, col.names = col.names)
    P.null.cent_list<- lapply(P.null, FUN=function(x){
      P.null.cent<-matrix.double.center(x)
    })
    newClusters <- FALSE
    # function to be used in parallel computation - lm model
    FUN <- function(y, i){
      res <- lm(as.numeric(y) ~ as.numeric(i))
      return(res)
    } # function to be called in parallel process
    
    if (is.numeric(parallel)) {
      parallel <- parallel::makeCluster(parallel, type = "PSOCK")
      newClusters <- TRUE
    }
    if (!inherits(parallel, "cluster")) {
      mod.L.null_list<- lapply(P.null.cent_list, function(x){
        mod.L.null_list<- lm(as.numeric(L.cent) ~ as.numeric(x))
      }
      )
    } else{
      # parallel version
      mod.L.null_list<- parallel::parLapply(cl = parallel,
                                            X = P.null.cent_list, 
                                            fun = FUN,
                                            y = L.cent)
    }
    if (newClusters){
      parallel::stopCluster(parallel)
    }
    
    mod.beta.null <- matrix(unlist(lapply(mod.L.null_list, 
                                          function(x){
                                            QuantPsyc::lm.beta(x)})),
                            nrow = nperm, 
                            ncol= 1)
    Res[,3]<- (sum(ifelse(mod.beta.null > Res[,2], 1, 0)) +1)/(nperm + 1)
    return(list(Matrix.P=P, Prediction=pred.L, Residuals= resid.L,Test.Results=Res))
  } else{
    Res[,3]<-"Nops"
    return(list(Matrix.P=P,Prediction=pred.L, Residuals= resid.L,Test.Results=Res))
  }
}

