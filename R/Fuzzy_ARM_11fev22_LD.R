matrix.p <- function(comm, phylo){
  matrix.w <- as.matrix(comm)
  #matrix.w <- as.matrix(sweep(comm, 1, rowSums(comm,na.rm = TRUE), "/"))
  #w.NA <- apply(matrix.w, 2, is.na)
  #matrix.w[w.NA] <- 0
  similar.phy<-ape::vcv(phylo,corr=TRUE)
  matrix.phy <- 1/rowSums(similar.phy)
  matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
  if(diag==FALSE){
    diag(matrix.q)<-0
  } else { matrix.q=matrix.q}
  matrix.P <- matrix.w %*% matrix.q
  return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
}


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
  L.cent <- matrix.double.center(comm)
  
  P <- matrix.p(comm, phylo=phylo)$matrix.P
  P.cent<-matrix.double.center(P)
  #P.cent <- matrix(NA, nrow=nrow(P), ncol = ncol(P))
 #for (c in 1:nrow(P)){
    #for (s in 1:ncol(P)){
      #P.cent[c,s] <- P[c, s] - mean(P[c, ]) - mean(P[, s]) + mean(P)
    #}
  #}
  mod.L<-lm(as.numeric(L.cent)~as.numeric(P.cent))
  pred.L<-predict(mod.L)
  resid.L<-L.cent-pred.L
  Res[,1]<-summary(mod.L)$r.squared
  Res[,2]<- QuantPsyc::lm.beta(mod.L)
  #### Test == TRUE  
  if(test == TRUE){
    seqpermutation.taxa_perm <- SYNCSA::permut.vector(ncol(comm), nset = nperm)
    seqpermutation.taxa <- lapply(seq_len(nrow(seqpermutation.taxa_perm)), function(i) seqpermutation.taxa_perm[i,])
    col.names<-colnames(comm)
    p.n.taxa <- function(samp, comm, phylo){
      comm.null<-comm
      colnames(comm.null)<-col.names[samp]
      match.null<-picante::match.phylo.comm(phylo,comm.null)
      MP.null <- matrix.p(match.null$comm,match.null$phy)$matrix.P
      return(MP.null)
    }
    P.null <- lapply(seqpermutation.taxa, p.n.taxa, comm = comm, phylo = phylo)
    P.null.cent_list<- lapply(P.null, FUN=function(x){
      P.null.cent<-matrix.double.center(x)
    })
    newClusters <- FALSE
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
      mod.L.null_list<- parallel::parLapply(cl = parallel,
                                            X = P.null.cent_list, 
                                            fun = function(i){
                                              lm(as.numeric(L.cent) ~ as.numeric(i))
                                            }
      )
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

