simul.comm<-function(Ncomm,Nspp,envir=FALSE,u=5,power=5,binary=FALSE,runs=30,test=FALSE,perm=999){
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
  Res<-matrix(NA,nrow=runs,ncol=11,dimnames=list(1:runs,paste(c("R2.W","Coef.a.W","Coef.b.W","Beta.W","Value.PCPS.1","Value.PCPS.2","R2.Res","Coef.a.Res","Coef.b.Res","Beta.Res","p.value"))))
  tree.list<-list()
  N<-matrix(NA,nrow=Nspp,ncol=runs,dimnames=list(1:Nspp,1:runs))
  K.N<-matrix(NA,nrow=runs,ncol=1,dimnames=list(1:runs,"K"))
  E<-matrix(NA,nrow=Ncomm,ncol=runs,dimnames=list(1:Ncomm,1:runs))
  L.cent<-matrix(NA,nrow=Ncomm,ncol=Nspp)
  P.cent<-matrix(NA,nrow=Ncomm,ncol=Nspp)
  P.null.cent<-matrix(NA,nrow=Ncomm,ncol=Nspp)
  resid.L.cent<-matrix(NA,nrow=Ncomm,ncol=Nspp)
  L.list<-list()
  P.list<-list()
  pred.L.list<-list()
  resid.L.list<-list()
    for(k in 1:runs){
      tree<-geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
      tree.list[[k]]<-tree
      a<-rnorm(Nspp,u,10)
      h<-runif(Nspp,0,30)
      L<-matrix(NA,Ncomm,Nspp)
      rownames(N)<-tree$tip.label
 	    niche<-ape::rTraitCont(ape::compute.brlen(tree,power=power),model="BM")
      N[,k]<-scales::rescale(niche,c(-1,101))
 	    K.N[k,]<-picante::Kcalc(N[,k],tree)
	       if(envir==TRUE){
	          E[,k]<-runif(Ncomm,0,100)
	       } else {E[,k]=N[,k]
	       }
		    for(i in 1:Ncomm){
			    for (j in 1:Nspp){
				    L[i,j]<-rpois(1,h[j]*exp((-((E[,k][i]-N[,k][j])^2))/(2*(a[j]^2))))
			    }
		    }
	    colnames(L)<-tree$tip.label
	    rownames(L)<-rownames(as.matrix(L),FALSE,prefix="comm")
	      if(binary==TRUE){
	        L<-vegan::decostand(L,"pa")
	      } else{L=L}
	    L.list[[k]]<-L
	      for (c in 1:nrow(L)){
	        for (s in 1:ncol(L)){
	          L.cent[c,s]<-L[c,s]-mean(L[c,])-mean(L[,s])+mean(L)
	        }
	      }
	    pcps.L<-PCPS::pcps(L, cophenetic(tree))
	    P<-matrix.p(L,cophenetic(tree))$matrix.P
	    P.list[[k]]<-P
	    L.list[[k]]<-L
	      for (c in 1:nrow(P)){
	        for (s in 1:ncol(P)){
	          P.cent[c,s]<-P[c,s]-mean(P[c,])-mean(P[,s])+mean(P)
	        }
	      }
	    mod.L<-lm(as.numeric(L.cent)~as.numeric(P.cent))
	    #### power/type I error test####     
	      if(test==TRUE){
	        mod.beta<-QuantPsyc::lm.beta(mod.L)
	        mod.beta.null<-matrix(NA,perm,1)
	          for(p in 1:perm){
	            dist_null<-picante::taxaShuffle(cophenetic(tree))
	            match.names<-match(colnames(L),colnames(dist_null))
	            P.null<-SYNCSA::matrix.p(L,as.matrix(dist_null[match.names, match.names]))$matrix.P
	              for (c in 1:nrow(P.null)){
	                for (s in 1:ncol(P.null)){
	                  P.null.cent[c,s]<-P.null[c,s]-mean(P.null[c,])-mean(P.null[,s])+mean(P.null)
	                }
	              }
	            mod.L.null<-lm(as.numeric(L.cent)~as.numeric(P.null.cent))
	            mod.beta.null[p,]<-QuantPsyc::lm.beta(mod.L.null)
	            print(p)
	          }
	        p.value<-(sum(ifelse(mod.beta.null>=mod.beta,1,0))+1)/(perm+1)
	      } else {p.value=NA
	        }
	    pred.L<-matrix(predict(mod.L),nrow(L),ncol(L),byrow=TRUE,dimnames=list(rownames(L),colnames(L)))
	    pred.L.list[[k]]<-pred.L
	    resid.L<-L.cent-pred.L
	    resid.L.list[[k]]<-resid.L
	      for (c in 1:nrow(resid.L)){
	        for (s in 1:ncol(resid.L)){
	          resid.L.cent[c,s]<-resid.L[c,s]-mean(resid.L[c,])-mean(resid.L[,s])+mean(resid.L)
	        }
	      }
	    mod.resid.L<-lm(as.numeric(resid.L.cent)~as.numeric(P.cent))
	    Res[k,1]<-summary(mod.L)$r.squared
	    Res[k,2]<-summary(mod.L)$coefficients[1]
      Res[k,3]<-summary(mod.L)$coefficients[2]
      Res[k,4]<-QuantPsyc::lm.beta(mod.L)
      Res[k,5]<-pcps.L$values[1,2]
      Res[k,6]<-pcps.L$values[2,2]
      Res[k,7]<-summary(mod.resid.L)$r.squared
      Res[k,8]<-summary(mod.resid.L)$coefficients[1]
      Res[k,9]<-summary(mod.resid.L)$coefficients[2]
      Res[k,10]<-QuantPsyc::lm.beta(mod.resid.L)
      Res[k,11]<-p.value
    }
  return(list(Trees=tree.list,Niche=N,K.niche=K.N,Environment=E,W.matrices=L.list,P.matrices=P.list,Predicted.W.matrices=pred.L.list,Residual.W.matrices=resid.L.list,Model.Results=Res))
  }