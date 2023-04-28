
source(here::here("R", "functions", "F_fuzzy-ARM.R"))

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


#' Metacommunity simulation with different relationship between phylogenetic signal and environmental variable
#'
#' @param Ncomm Scalar. The number of communities that will be simulated in the metacommunities
#' @param Nspp Scalar. The number of species that will be simulated in the metacommunities
#' @param envir Logical. If TRUE an environmental vector with normal distribution will be generated and used in simulations
#' @param diag Logical. If TRUE the diagonal in matrix q is included, if FALSE the diagonal is removed
#' @param u Scalar. Mean value to be used in rnorm function
#' @param power Scalar. Power parameter to compute the degree of phylogenetic signal in the simulated trait. The same parameter to be passed to \code{\link[ape]{rTraitCont}} function
#' @param binary Logical. If TRUE (default) the community matrix will be transformed to presence and absence
#' @param runs Scalar. An integer informing the number of times to run the simulation model
#' @param test Logical. If TRUE (default) it will be computed power and type I error rates for beta parameter
#' @param nperm Scalar. An integer indicating the number of permutations to be used to calculate power and type I error rates
#' @param parallel Scalar. An integer indicating the number of cores to be used in parallel computation. If NULL parallel computing won't be used
#'
#' @return \itemize{
#' \item `Trees` list of simulated trees
#' \item `Niche` simulated niche values
#' \item `K.N` phylogenetic signal calculated for niche values accordingly to K statistic
#' \item `Environment` if environment is TRUE the environmental vector for each simulation will be returned
#' \item `W.matrices` simulated occurrence matrices
#' \item `P.matrices` simulated P matrices
#' \item `Predicted.W.matrices` predicted W matrices for each run
#' \item `Residual.W.matrices` residual matrices from the model between P and W matrices
#' \item `Model.Results` coefficient values from model between P and W matrices
#' 
#' }
#' a list with all parameters returned from the simulation model 
#' @export
#'
#' @examples
simul.comm <- function(Ncomm, 
                        Nspp,
                        envir=FALSE,
                        diag=FALSE,
                        u=5,
                        power=5,
                        binary=FALSE,
                        runs=30,
                        test=FALSE,
                        nperm=999,
                        parallel=NULL){
  # Compute matrix P:
  # Double center matrix transformation: 
  
  Res <- matrix(NA, 
                nrow = runs,
                ncol = 9,
                dimnames = list(1:runs, paste(c("R2.W","Coef.a.W","Coef.b.W","Beta.W","p.value","R2.Res","Coef.a.Res","Coef.b.Res","Beta.Res")
                )
                )
  )
  tree.list <- vector(mode = "list", length = runs) 
  N <- matrix(NA, nrow = Nspp, ncol = runs, dimnames = list(1:Nspp, 1:runs))
  K.N <- matrix(NA, nrow = runs, ncol = 1, dimnames = list(1:runs,"K"))
  E <- matrix(NA, nrow = Ncomm, ncol = runs, dimnames = list(1:Ncomm, 1:runs))
  L.list <- vector(mode = "list", length = runs)
  P.list <- list()
  pred.L.list <- list()
  resid.L.list <- list()
  
  # setting a progress bar 
  
  pb <- txtProgressBar(min = 0,   
                       max = runs, 
                       style = 3, 
                       width = 50,
                       char = "=")
  
  
  for(k in 1:runs){
    #k = 4
    phylo <- geiger::sim.bdtree(b=0.1, d=0, stop = "taxa", n = Nspp, extinct = FALSE)
    tree.list[[k]] <- phylo
    a <- rnorm(Nspp, u, 10)
    h <- runif(Nspp, 0, 30)
    L <- matrix(NA, Ncomm, Nspp, dimnames=list(paste("C", 1:Ncomm, sep=""), phylo$tip.label))
    rownames(N) <- phylo$tip.label
    
    # Simulate niche N:    
    niche <- ape::rTraitCont(ape::compute.brlen(phylo, power = power), model="BM")
    N[ ,k] <- scales::rescale(niche, c(-1,101)) # niche matrix
    K.N[k,] <- picante::Kcalc(N[,k], phylo)
    # Simulate environment E:
    if(envir==TRUE){
      E[,k]<-runif(Ncomm,0,100)
    } else {E[,k]=N[,k]
    }
    # Simulate L:	    
    for(i in 1:Ncomm){
      for (j in 1:Nspp){
        L[i,j]<-rpois(1,h[j]*exp((-((E[,k][i]-N[,k][j])^2))/(2*(a[j]^2))))
      }
    }
    colnames(L) <- phylo$tip.label
    rownames(L) <- rownames(as.matrix(L), FALSE, prefix = "comm")
    if(binary==TRUE){
      L<-vegan::decostand(L,"pa")
    } else{L=L}
    L.list[[k]]<-L
    L.cent<-matrix.double.center(L)
    
    # running fuzzy model
    test_res <- Fuzzy_ARM(comm = L, phylo = phylo, binary = binary, diag = diag, test = test, nperm = nperm, parallel = parallel)
    
    
    Res[k,5] <- test_res$Test.Results[, "P.value"]
    # Compute predicted and residual W matrices:
    P.list[[k]] <- test_res$Matrix.P
    pred.L.list[[k]] <- test_res$Prediction
    resid.L.list[[k]]<- test_res$Residuals
    resid.L.cent <- matrix.double.center(test_res$Residuals)
    P.cent <- matrix.double.center(test_res$Matrix.P)
    mod.resid.L <- lm(as.numeric(resid.L.cent) ~ as.numeric(P.cent))
    Res[k,6]<-summary(mod.resid.L)$r.squared
    Res[k,7]<-summary(mod.resid.L)$coefficients[1]
    Res[k,8]<-summary(mod.resid.L)$coefficients[2]
    Res[k,9]<-QuantPsyc::lm.beta(mod.resid.L)
    setTxtProgressBar(pb, k)
  }
  return(list(Trees = tree.list,
              Niche = N,
              K.niche = K.N, 
              Environment = E, 
              W.matrices = L.list, 
              P.matrices = P.list, 
              Predicted.W.matrices = pred.L.list, 
              Residual.W.matrices = resid.L.list, 
              Model.Results = Res)
  )
}