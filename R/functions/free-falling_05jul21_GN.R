free.falling <- function(comm,
                         phylo, 
                         binary, 
                         test = TRUE, 
                         nperm = 1000, 
                         parallel = 8){
  matrix.p <- function(comm, phylodist)
  {
    matrix.w <- as.matrix(comm)
    phylodist <- as.matrix(phylodist)
    similar.phy <- 1 - (phylodist/max(phylodist))
    matrix.phy <- 1/colSums(similar.phy)
    matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
    matrix.P <- matrix.w %*% matrix.q
    return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
  }
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
  
  L.cent <- matrix(NA, 
                   nrow = nrow(comm), 
                   ncol = ncol(comm)
  )
  
  for (c in 1:nrow(comm)){
    for (s in 1:ncol(comm)){
      L.cent[c,s] <- comm[c, s] - mean(comm[c, ]) - mean(comm[, s]) + mean(comm)
    }
  }
  
  P <- matrix.p(comm, phylodist = cophenetic(phylo))$matrix.P
  P.cent <- matrix(NA, nrow=nrow(P), ncol = ncol(P))
  for (c in 1:nrow(P)){
    for (s in 1:ncol(P)){
      P.cent[c,s] <- P[c, s] - mean(P[c, ]) - mean(P[, s]) + mean(P)
    }
  }
  mod.L<-lm(as.numeric(L.cent)~as.numeric(P.cent))
  pred.L<-predict(mod.L)
  resid.L<-L.cent-pred.L
  Res[,1]<-summary(mod.L)$r.squared
  Res[,2]<- QuantPsyc::lm.beta(mod.L)
  #### Test == TRUE  
  if(test == TRUE){
    P.null.cent<- matrix(NA, nrow= nrow(P.cent), ncol=ncol(P.cent))
    seqpermutation.taxa <- SYNCSA::permut.vector(ncol(comm), nset = nperm)
    seqpermutation.taxa <- lapply(seq_len(nrow(seqpermutation.taxa)), function(i) seqpermutation.taxa[i,])
    phylodist<- cophenetic(phylo)
    p.n.taxa <- function(samp, comm, phylodist){
      MP.null <- matrix.p(comm, phylodist[samp, samp])$matrix.P
      return(MP.null)
    }
    P.null <- lapply(seqpermutation.taxa, p.n.taxa, comm = comm, phylodist = phylodist)
    P.null.cent_list<- lapply(P.null, FUN=function(x){
      for (c in 1:nrow(x)){
        for (s in 1:ncol(x)){
          P.null.cent[c,s]<-x[c,s]-mean(x[c,])-mean(x[,s])+mean(x)
        }
      }
      P.null.cent
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
    Res[,3]<- (sum(ifelse(mod.beta.null >= Res[,2], 1, 0)) +1)/(nperm + 1)
    return(Res)
  } else{
    return(Res)
  }
}

