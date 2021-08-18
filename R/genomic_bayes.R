####################################################################################################################
#    Module 6: Bayesian models
####################################################################################################################
#'
#' Genomic prediction models implemented using Bayesian Methods (small data)
#'
#' @description
#' Genomic prediction models implemented using Bayesian Methods (small data).
#' The models are implemented using empirical Bayesian methods. The hyperparameters of the dispersion parameters of the Bayesian model can
#' be obtained from prior information or estimated by maximum likelihood, and conditional on these, the model is fitted using
#' Markov chain Monte Carlo. These functions are currently under development and future release will be able to handle large data sets.
#'
#'
#' @param y is a matrix of phenotypes
#' @param W is a matrix of centered and scaled genotypes
#' @param nsamp is the number of samples after burnin
#' @param sets is a list of markers defining a group
#' @param nsets is a list of number of marker groups
#' @param phi is the proportion of markers in each marker variance class (phi=c(0.999,0.001),used if method="ssvs")
#' @param h2 is the trait heritability
#' @param method specifies the methods used (method="mixed","bayesC",blasso","blr")
#' @param nburn is the number of burnin samples
#' @param nsave is the number of samples to save
#' @param tol is the tolerance
#'

#' @author Peter Sørensen


#' @examples
#'
#'
#' # Simulate data and test functions
#'
#' W <- matrix(rnorm(100000),nrow=1000)
#' set1 <- sample(1:ncol(W),5)
#' set2 <- sample(1:ncol(W),5)
#' sets <- list(set1,set2)
#' g <- rowSums(W[,c(set1,set2)])
#' e <- rnorm(nrow(W),mean=0,sd=1)
#' y <- g + e
#'
#'
#' fitM <- gbayes(y=y, W=W, method="mixed")
#' fitA <- gbayes(y=y, W=W, method="bayesA")
#' fitL <- gbayes(y=y, W=W, method="blasso")
#' fitC <- gbayes(y=y, W=W, method="bayesC")
#'
#' fitM <- gbayes(y=y, W=W, method="mixed", algorithm="fastbayes")
#' fitA <- gbayes(y=y, W=W, method="bayesA", algorithm="fastbayes")
#' fitL <- gbayes(y=y, W=W, method="blasso", algorithm="fastbayes")
#' fitC <- gbayes(y=y, W=W, method="bayesC", algorithm="fastbayes")


#'
#' @export
#'

gbayes <- function(y=NULL, X=NULL, W=NULL, Glist=NULL, chr=NULL, rsids=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                   vara=NULL, varb=NULL, vare=NULL, ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=TRUE,
                   h2=NULL, pi=0.001, updateB=TRUE, updateE=TRUE, updatePi=TRUE, models=NULL,
                   nub=4, nue=4, nit=100, nit_local=NULL,nit_global=NULL,
                   method="mixed", algorithm="default") {
     
     methods <- c("blup","mixed","bayesA","blasso","bayesC")
     method <- match(method, methods) - 1
     if( !sum(method%in%c(0:4))== 1 ) stop("Method specified not valid") 
     
     nt <- 1
     if(is.list(y)) nt <- length(y)
     
     if(nt==1 && !algorithm=="sbayes") fit <- bayes(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
                                                     vara=vara, varb=varb, vare=vare, 
                                                     ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
                                                     h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
                                                     nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm)  
     
     if(nt==1 && algorithm=="sbayes" && is.null(Glist)) fit <- sbayes(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
                                                                      vara=vara, varb=varb, vare=vare, 
                                                                      ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
                                                                      h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
                                                                      nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm)  
     
     if(nt>1 && is.null(Glist)) fit <- mtbayes(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
                                               vara=vara, varb=varb, vare=vare, 
                                               ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
                                               h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
                                               nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm) 
     
     # if(nt==1 && algorithm=="sbayes" && !is.null(Glist)) {
     #   
     #   if(is.matrix(y)) ids <- rownames(y)
     #   if(is.vector(y)) ids <- names(y)
     #   rws <- match(ids,Glist$ids)
     #   if(any(is.na(rws))) stop("some elements in names(y) does not match elements in Glist$ids ")       
     # 
     #   if(is.null(chr)) chromosomes <- 1:Glist$nchr
     #   if(!is.null(chr)) chromosomes <- chr
     #   fit <- vector(length=length(chromosomes),mode="list")
     #   names(fit) <- chromosomes
     #   
     # 
     #   e <- y
     #   g <- rep(0,Glist$n)
     #   
     #   for (chr in chromosomes) {
     #     rsidsCVS <- Glist$rsids[[chr]]
     #     if(!is.null(rsids)) rsidsCVS <- Glist$rsids[[chr]][Glist$rsids[[chr]]%in%rsids]
     #     clsCVS <- match(rsidsCVS,Glist$rsids[[chr]])
     #     clsCVS <- clsCVS[!is.na(clsCVS)]
     #     covs <- cvs(y=e,Glist=Glist,chr=chr, cls=clsCVS)
     #     mlogp <- -log10(covs$p[[1]])
     #     if(any(is.na(mlogp))) {
     #       print(paste("Number of marker removed:",sum(is.na(mlogp))))
     #       mlogp <- mlogp[!is.na(mlogp)]
     #     }
     #     if(!is.null(rsids)) {
     #       print(paste("Number of markers used:",sum(names(mlogp)%in%rsids),"from chromosome:",chr))
     #       mlogp <- mlogp[names(mlogp)%in%rsids]
     #     }
     #     cls <- match(names(mlogp),Glist$rsids[[chr]])
     #     m <- length(mlogp) 
     #     o <- order(mlogp,decreasing=TRUE)
     #     cls <- splitWithOverlap(cls[o],1000,0)
     #     sets <- splitWithOverlap((1:m)[o],1000,0)
     #     dm <- bm <- rep(0,m)
     #     names(dm) <- names(bm) <- names(mlogp)
     #     varem <- varbm <- pim <- vector(length=length(sets),mode="list")
     #     for (i in 1:length(sets)) {
     #       W <- getG(Glist, chr=chr, scale=TRUE, cls=cls[[i]])
     #       LD <- crossprod(W[rws,])
     #       fitset <- sbayes(y=e, X=X, W=W[rws,], b=b, badj=badj, seb=seb, LD=LD, n=n,
     #                        vara=vara, varb=varb, vare=vare, 
     #                        ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=FALSE,
     #                        h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
     #                        nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm)
     #       gset <- crossprod(t(W),fitset$bm)
     #       g <- g + gset
     #       e <- e - gset[rws,]
     #       dm[sets[[i]]] <- fitset$dm
     #       bm[sets[[i]]] <- fitset$bm
     #       varem[[i]] <- fitset$E
     #       varbm[[i]] <- fitset$B
     #       pim[[i]] <- fitset$Pi
     #       print(paste("Finished segment:",i,"out of",length(sets),"segments on chromosome:",chr))
     #     }
     #     fit[[chr]] <- list(bm=bm,dm=dm,E=varem,B=varbm,Pi=pim)
     #   }
     #   fit$g <- g[rws,]
     #   fit$gtrain <- g[rws,]
     #   fit$gtest <- g[-rws,]
     #   fit$e <- e
     # }
     if(nt==1 && algorithm=="sbayes" && !is.null(Glist)) {
       fit <- NULL
       if(is.matrix(y)) ids <- rownames(y)
       if(is.vector(y)) ids <- names(y)
       rws <- match(ids,Glist$ids)
       if(any(is.na(rws))) stop("some elements in names(y) does not match elements in Glist$ids ")       
       
       if(is.null(chr)) chromosomes <- 1:Glist$nchr
       if(!is.null(chr)) chromosomes <- chr
       
       bm <- dm <- covs <- fit <- stat <- vector(length=Glist$nchr,mode="list")
       names(covs) <- names(bm) <- names(dm) <- 1:Glist$nchr
       
       e <- y
       
       if(is.null(nit_local)) nit_local <- 500
       if(is.null(nit_global)) nit_global <- as.integer(nit/nit_local)
       for (chr in chromosomes){
         yy <- sum((e-mean(e))**2)
         n <- length(e)
         print(paste("Computing summary statistics for chromosome:",chr))
         LD <- getSparseLD(Glist = Glist, chr = chr)
         LD$indices <- lapply(LD$indices,function(x){x-1})
         LD$values <- lapply(LD$values,function(x){x*n})
         rsidsLD <- names(LD$values)
         clsLD <- match(rsidsLD,Glist$rsids[[chr]])
         covs[[chr]] <- cvs(y=e,Glist=Glist,chr=chr, cls=clsLD)
         wy <- covs[[chr]]$Xy[[1]]
         b <- rep(0,length(wy))
         fit[[chr]] <- sbayes_sparse(yy=yy, 
                              wy=wy,
                              b=b, 
                              LDvalues=LD$values, 
                              LDindices=LD$indices, 
                              method=method, 
                              nit=nit, 
                              n=n, 
                              pi=pi,
                              nue=nue, 
                              nub=nub, 
                              updateB=updateB, 
                              updateE=updateE, 
                              updatePi=updatePi)
         #b <- fit[[chr]]$b
         #bm[[chr]] <- fit[[chr]]$bm
         #dm[[chr]] <- fit[[chr]]$dm
         
         # for (iter in 1:nit_global) {
         #   b <- bmchr <- dmchr <- rep(0,length(wy))
         #   print(paste("Fit", methods[method+1],"for chromosome:",chr))
         #   for (i in 1:length(Glist$clsLD[[chr]])) {
         #     fitset <- sbayes_dense(yy=yy, wy=wy[ Glist$clsLD[[chr]][[i]] ],
         #                           b=b[ Glist$clsLD[[chr]][[i]] ], 
         #                           LD=Glist$LD[[chr]][[i]]*n, 
         #                           method=method, 
         #                           nit=nit_local, 
         #                           n=n, 
         #                           pi=pi,
         #                           nue=nue, 
         #                           nub=nub, 
         #                           updateB=updateB, 
         #                           updateE=updateE, 
         #                           updatePi=updatePi)
         #     b[Glist$clsLD[[chr]][[i]]] <- fitset$b
         #     bmchr[Glist$clsLD[[chr]][[i]]] <- fitset$bm
         #     dmchr[Glist$clsLD[[chr]][[i]]] <- fitset$dm
         #   }
         # }
         #bm[[chr]] <- cbind(bm[[chr]], bmchr)
         #dm[[chr]] <- cbind(dm[[chr]], dmchr)
         stat[[chr]] <- data.frame(chr=rep(chr,length(rsidsLD)),rsids=rsidsLD,alleles=Glist$a2[[chr]][clsLD], af=Glist$af[[chr]][clsLD],bm=fit$bm)
         rownames(stat[[chr]]) <- rsidsLD
       }
       #stat <- do.call(rbind, stat)
       fit$stat <- stat
       fit$covs <- covs
     }
     

     #fit$acc <- acc(yobs=y,ypred=fit$g)
     fit$method <- methods[method+1]
          
     return(fit)
     
          # ids <- NULL
          # if(is.matrix(y)) ids <- rownames(y)
          # if(is.vector(y)) ids <- names(y)
          # 
          # if(scaleY) y <- as.vector(scale(y)) 
          # 
          # if(is.null(pi)) pi <- 0.001
          # 
          # if(is.null(h2)) h2 <- 0.5
          # 
          # n <- nrow(W)
          # m <- ncol(W)
          # 
          # if(is.null(ids)) warning("No names/rownames provided for y")
          # if(is.null(rownames(W))) warning("No names/rownames provided for W")
          # if(!is.null(ids) & !is.null(rownames(W))) {
          #   if(any(is.na(match(ids,rownames(W))))) stop("Names/rownames for y does match rownames for W")
          #                         
          # }
          # 
          # 
          # 
          # if(is.null(b)) b <- rep(0,m)
          # e=y-mean(y)
          # if(is.null(vare)) vare <- var(e)
          # if(is.null(varb)) varb <- (vare/m)*h2
          # if(is.null(lambda)) lambda <- rep(vare/varb,m)
          # if(is.null(vara)) vara <- vare*h2
          # 
          # #if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/(pi*m*0.5))
          # if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/m)
          # if(is.null(sse_prior)) sse_prior <- nue*vare
          # 
          # 
          # if(algorithm=="default") {
          #      fit <- .Call("_qgg_bayes",
          #                   y=y, 
          #                   W=split(W, rep(1:ncol(W), each = nrow(W))), 
          #                   b=b,
          #                   lambda = lambda,
          #                   pi = pi,
          #                   vara = vara,
          #                   varb = varb,
          #                   vare = vare,
          #                   ssb_prior=ssb_prior,
          #                   sse_prior=sse_prior,
          #                   nub=nub,
          #                   nue=nue,
          #                   updateB = updateB,
          #                   updateE = updateE,
          #                   updatePi = updatePi,
          #                   nit=nit,
          #                   method=as.integer(method)) 
          #      names(fit[[1]]) <- colnames(W)
          #      fit[[7]] <- crossprod(t(W),fit[[10]])[,1]
          #      names(fit[[7]]) <- names(fit[[8]]) <- ids
          #      names(fit) <- c("bm","dm","mu","B","E","Pi","g","e","param","b")
          # } 
          # 
          # 
          # if(algorithm=="fastbayes") {
          #      LD <- crossprod(W)
          #      
          #      fit <- .Call("_qgg_fbayes",
          #                     y=y, 
          #                    W=split(W, rep(1:ncol(W), each = nrow(W))), 
          #                    LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
          #                    b=b,
          #                    lambda = lambda,
          #                    pi = pi,
          #                    vara = vara,
          #                    varb = varb,
          #                    vare = vare,
          #                    nub=nub,
          #                    nue=nue,
          #                    updateB = updateB,
          #                    updateE = updateE,
          #                    updatePi = updatePi,
          #                    nit=nit,
          #                    method=as.integer(method)) 
          #      names(fit[[1]]) <- colnames(W)
          #      names(fit) <- c("b","p","mu","B","E","Pi")
          # }

     
     
     # if(algorithm=="sbayes") {
     #   
     #   
     #   ids <- NULL
     #   if(is.matrix(y)) ids <- rownames(y)
     #   if(is.vector(y)) ids <- names(y)
     #   
     #   if(scaleY) y <- as.vector(scale(y)) 
     #   
     #   wy <- as.vector(crossprod(W,y))           
     #   yy <- sum((y-mean(y))**2)
     #   
     #   n <-nrow(W)       
     #   
     #   if(!is.null(W) && is.null(LD)) {
     #     n <- nrow(W)
     #     LD <- crossprod(W)
     #     #LD <- crossprod(W)/(n-1)
     #   }    
     #   
     #   
     #   m <- ncol(LD)
     #   
     #   if(is.null(pi)) pi <- 0.01
     #   
     #   if(is.null(h2)) h2 <- 0.5
     #   
     #   if(is.null(vare)) vare <- 1
     #   if(is.null(varb)) varb <- (vare/m)*h2
     #   if(is.null(lambda)) lambda <- rep(vare/varb,m)
     #   if(is.null(vara)) vara <- vare*h2
     #   
     #   if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/(pi*m*0.5))
     #   if(is.null(sse_prior)) sse_prior <- nue*vare
     #   
     #   if(is.null(b)) b <- rep(0,m)
     #   
     # 
     #      fit <- .Call("_qgg_sbayes",
     #                   wy=wy, 
     #                   LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
     #                   b = b,
     #                   lambda = lambda,
     #                   yy = yy,
     #                   pi = pi,
     #                   vara = vara,
     #                   varb = varb,
     #                   vare = vare,
     #                   ssb_prior=ssb_prior,
     #                   sse_prior=sse_prior,
     #                   nub=nub,
     #                   nue=nue,
     #                   updateB = updateB,
     #                   updateE = updateE,
     #                   updatePi = updatePi,
     #                   n=n,
     #                   nit=nit,
     #                   method=as.integer(method))
     #      names(fit[[1]]) <- 
     #      #names(fit) <- c("b","p","mu","B","E","Pi")
     #      names(fit) <- c("bm","dm","mu","B","E","Pi","g","e","param","b")
     # }
     
     # if(nt>1) {
     #      
     #      if(method==2) stop("Multiple trait not yet implemented for Bayes A") 
     #      if(method==3) stop("Multiple trait not yet implemented for Bayesian Lasso")
     #      
     #      n <- nrow(W)
     #      m <- ncol(W)
     #      
     #      if(scaleY) y <- lapply(y,function(x){as.vector(scale(x))})
     #      if(!scaleY) y <- lapply(y,function(x){x-mean(x) })
     #      
     #      
     #      if(is.null(b)) b <- lapply(1:nt,function(x){rep(0,m)}) 
     #      
     #      if(is.null(models)) {
     #           models <- rep(list(0:1), nt)
     #           models <- t(do.call(expand.grid, models))
     #           models <- split(models, rep(1:ncol(models), each = nrow(models)))
     #           pi <- c(0.999,rep(0.001,length(models)-1)) 
     #      } 
     #      if(is.character(models)) {
     #           if(models=="restrictive") {
     #                models <- list(rep(0,nt),rep(1,nt))
     #                pi <- c(0.999,0.001)
     #           }
     #      }
     #      
     #      #if(is.null(pi)) {
     #      #     pi <- c(0.999,rep(0.001,length(models)-1)) 
     #      #}
     #      
     #      if(is.null(h2)) h2 <- 0.5
     #      if(is.null(vare)) {
     #           vare <- diag(sapply(y,var))
     #      }
     #      if(is.null(varb)) varb <- diag(sapply(y,var)/(m*pi[length(models)]))*h2
     #      #if(is.null(varb)) varb <- (vare*h2)/(m*pi[length(models)])
     #      if(is.null(vara)) vara <- diag(diag(vare))*h2
     #      
     #      #if(is.null(ssb_prior)) ssb_prior <-  diag((nub-2.0)/nub * (vara/(m*pi[length(models)])))
     #      #if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/(m*pi[length(models)]))
     #      if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/m)
     #      #if(is.null(sse_prior)) sse_prior <- diag(nue*vare)
     #      if(is.null(sse_prior)) sse_prior <- nue*diag(diag(vare))
     #      
     #              
     #      fit <- .Call("_qgg_mtbayes",
     #                   y=y, 
     #                   W=split(W, rep(1:ncol(W), each = nrow(W))), 
     #                   b=b,
     #                   B = varb,
     #                   E = vare,
     #                   ssb_prior=split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
     #                   sse_prior=split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
     #                   models=models,
     #                   pi=pi,
     #                   nub=nub,
     #                   nue=nue,
     #                   updateB=updateB,
     #                   updateE=updateE,
     #                   updatePi=updatePi,
     #                   nit=nit,
     #                   method=as.integer(method))
     #      fit[[6]] <- matrix(unlist(fit[[6]]), ncol = nt, byrow = TRUE)
     #      fit[[7]] <- matrix(unlist(fit[[7]]), ncol = nt, byrow = TRUE)
     #      colnames(fit[[6]]) <- rownames(fit[[6]]) <- paste0("T",1:nt)
     #      colnames(fit[[7]]) <- rownames(fit[[7]]) <- paste0("T",1:nt)
     #      fit[[11]] <- matrix(unlist(fit[[11]]), ncol = nt, byrow = TRUE)
     #      fit[[12]] <- matrix(unlist(fit[[12]]), ncol = nt, byrow = TRUE)
     #      colnames(fit[[11]]) <- rownames(fit[[11]]) <- paste0("T",1:nt)
     #      colnames(fit[[12]]) <- rownames(fit[[12]]) <- paste0("T",1:nt)
     #      fit[[13]] <- fit[[13]][[1]]
     #      fit[[14]] <- fit[[14]][[1]]
     #      fit[[15]] <- fit[[15]][[1]]
     #      names(fit) <- c("bm","dm","mu","Bm","Em","rg","re","g","e","b","B","E","pi","pim","order")
     # 
     # }
     

}

bayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                   vara=NULL, varb=NULL, vare=NULL, 
                   ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                   h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                   nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL) {
  ids <- NULL
  if(is.matrix(y)) ids <- rownames(y)
  if(is.vector(y)) ids <- names(y)
  
  if(scaleY) y <- as.vector(scale(y)) 
  
  if(is.null(pi)) pi <- 0.001
  
  if(is.null(h2)) h2 <- 0.5
  
  n <- nrow(W)
  m <- ncol(W)
  
  if(is.null(ids)) warning("No names/rownames provided for y")
  if(is.null(rownames(W))) warning("No names/rownames provided for W")
  if(!is.null(ids) & !is.null(rownames(W))) {
    if(any(is.na(match(ids,rownames(W))))) stop("Names/rownames for y does match rownames for W")
    
  }
  
  
  
  if(is.null(b)) b <- rep(0,m)
  e=y-mean(y)
  if(is.null(vare)) vare <- var(e)
  if(method<4 && is.null(varb)) varb <- (vare*h2)/m
  if(method==4 && is.null(varb)) varb <- (vare*h2)/(m*pi)
  #if(is.null(varb)) varb <- (vare/m)*h2
  if(is.null(lambda)) lambda <- rep(vare/varb,m)
  if(is.null(vara)) vara <- vare*h2
  
  #if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/(pi*m*0.5))
  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/m)
  if(is.null(sse_prior)) sse_prior <- nue*vare
  
  
  if(algorithm=="default") {
    fit <- .Call("_qgg_bayes",
                 y=y, 
                 W=split(W, rep(1:ncol(W), each = nrow(W))), 
                 b=b,
                 lambda = lambda,
                 pi = pi,
                 vara = vara,
                 varb = varb,
                 vare = vare,
                 ssb_prior=ssb_prior,
                 sse_prior=sse_prior,
                 nub=nub,
                 nue=nue,
                 updateB = updateB,
                 updateE = updateE,
                 updatePi = updatePi,
                 nit=nit,
                 method=as.integer(method)) 
    names(fit[[1]]) <- colnames(W)
    fit[[7]] <- crossprod(t(W),fit[[10]])[,1]
     names(fit[[7]]) <- names(fit[[8]]) <- ids
    names(fit) <- c("bm","dm","mu","B","E","Pi","g","e","param","b")
  } 
  
  
  if(algorithm=="fastbayes") {
    LD <- crossprod(W)
    
    fit <- .Call("_qgg_fbayes",
                 y=y, 
                 W=split(W, rep(1:ncol(W), each = nrow(W))), 
                 LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
                 b=b,
                 lambda = lambda,
                 pi = pi,
                 vara = vara,
                 varb = varb,
                 vare = vare,
                 nub=nub,
                 nue=nue,
                 updateB = updateB,
                 updateE = updateE,
                 updatePi = updatePi,
                 nit=nit,
                 method=as.integer(method)) 
    names(fit[[1]]) <- colnames(W)
    names(fit) <- c("b","p","mu","B","E","Pi")
  }
  
  return(fit)
  
  
}


#sbayes <- function(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
#                    vara=vara, varb=varb, vare=vare, 
#                    ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
#                    h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
#                    nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm) {

sbayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                  vara=NULL, varb=NULL, vare=NULL, 
                  ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                  h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                  nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL) {
  
  ids <- NULL
  if(is.matrix(y)) ids <- rownames(y)
  if(is.vector(y)) ids <- names(y)
  
  if(scaleY) y <- as.vector(scale(y)) 
  
  wy <- as.vector(crossprod(W,y))           
  yy <- sum((y-mean(y))**2)
  
  n <-nrow(W)       
  
  if(!is.null(W) && is.null(LD)) {
    n <- nrow(W)
    LD <- crossprod(W)
    #LD <- crossprod(W)/(n-1)
  }    
  
  
  m <- ncol(LD)
  
  if(is.null(pi)) pi <- 0.01
  
  if(is.null(h2)) h2 <- 0.5
  
  if(is.null(vare)) vare <- 1
  #if(is.null(varb)) varb <- (vare/m)*h2
  if(method<4 && is.null(varb)) varb <- (vare*h2)/m
  if(method==4 && is.null(varb)) varb <- (vare*h2)/(m*pi)
  if(is.null(lambda)) lambda <- rep(vare/varb,m)
  if(is.null(vara)) vara <- vare*h2
  
  #if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/(pi*m*0.5))
  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/m)
  if(is.null(sse_prior)) sse_prior <- nue*vare
  
  if(is.null(b)) b <- rep(0,m)
  
  
  fit <- .Call("_qgg_sbayes",
               wy=wy, 
               LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
               b = b,
               lambda = lambda,
               yy = yy,
               pi = pi,
               vara = vara,
               varb = varb,
               vare = vare,
               ssb_prior=ssb_prior,
               sse_prior=sse_prior,
               nub=nub,
               nue=nue,
               updateB = updateB,
               updateE = updateE,
               updatePi = updatePi,
               n=n,
               nit=nit,
               method=as.integer(method))
  names(fit[[1]]) <- rownames(LD)
  if(!is.null(W)) fit[[7]] <- crossprod(t(W),fit[[10]])[,1]
  names(fit[[7]]) <- ids
  names(fit) <- c("bm","dm","mu","B","E","Pi","g","e","param","b")

  return(fit)
  
}


sbayes_dense <- function(yy=NULL, wy=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                        vara=NULL, varb=NULL, vare=NULL, 
                        ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                        h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                        nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL) {
  
  
  m <- ncol(LD)
  if(is.null(pi)) pi <- 0.01
  if(is.null(h2)) h2 <- 0.5
  if(is.null(vare)) vare <- 1
  if(method<4 && is.null(varb)) varb <- (vare*h2)/m
  if(method==4 && is.null(varb)) varb <- (vare*h2)/(m*pi)
  if(is.null(lambda)) lambda <- rep(vare/varb,m)
  if(is.null(vara)) vara <- vare*h2
  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/m)
  if(is.null(sse_prior)) sse_prior <- nue*vare
  if(is.null(b)) b <- rep(0,m)
  
  
  fit <- .Call("_qgg_sbayes",
               wy=wy, 
               LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
               b = b,
               lambda = lambda,
               yy = yy,
               pi = pi,
               vara = vara,
               varb = varb,
               vare = vare,
               ssb_prior=ssb_prior,
               sse_prior=sse_prior,
               nub=nub,
               nue=nue,
               updateB = updateB,
               updateE = updateE,
               updatePi = updatePi,
               n=n,
               nit=nit,
               method=as.integer(method))
  names(fit[[1]]) <- rownames(LD)
  names(fit) <- c("bm","dm","mu","B","E","Pi","g","e","param","b")
  
  return(fit)
  
}

sbayes_sparse <- function(yy=NULL, wy=NULL, b=NULL, badj=NULL, seb=NULL, LDvalues=NULL,LDindices=NULL, n=NULL,
                        vara=NULL, varb=NULL, vare=NULL, 
                        ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                        h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                        nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL) {
  
  
  m <- length(LDvalues)
  if(is.null(pi)) pi <- 0.01
  if(is.null(h2)) h2 <- 0.5
  if(is.null(vare)) vare <- 1
  if(method<4 && is.null(varb)) varb <- (vare*h2)/m
  if(method==4 && is.null(varb)) varb <- (vare*h2)/(m*pi)
  if(is.null(lambda)) lambda <- rep(vare/varb,m)
  if(is.null(vara)) vara <- vare*h2
  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/m)
  if(is.null(sse_prior)) sse_prior <- nue*vare
  if(is.null(b)) b <- rep(0,m)
  
  
  fit <- .Call("_qgg_sbayes_spa",
               wy=wy, 
               LDvalues=LDvalues, 
               LDindices=LDindices, 
               b = b,
               lambda = lambda,
               yy = yy,
               pi = pi,
               vara = vara,
               varb = varb,
               vare = vare,
               ssb_prior=ssb_prior,
               sse_prior=sse_prior,
               nub=nub,
               nue=nue,
               updateB = updateB,
               updateE = updateE,
               updatePi = updatePi,
               n=n,
               nit=nit,
               method=as.integer(method))
  names(fit[[1]]) <- names(LDvalues)
  names(fit) <- c("bm","dm","mu","B","E","Pi","g","e","param","b")
  
  return(fit)
  
}

#mtbayes <- function(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
#                    vara=vara, varb=varb, vare=vare, 
#                    ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
#                    h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
#                    nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm) {
mtbayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                  vara=NULL, varb=NULL, vare=NULL, 
                  ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                  h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                  nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL) {
  
  if(is.list(y)) nt <- length(y)
  if(!is.list(y)) stop("This is not a multiple trait analysis")
  
  if(method==2) stop("Multiple trait not yet implemented for Bayes A") 
  if(method==3) stop("Multiple trait not yet implemented for Bayesian Lasso")
  
  n <- nrow(W)
  m <- ncol(W)
  
  if(scaleY) y <- lapply(y,function(x){as.vector(scale(x))})
  if(!scaleY) y <- lapply(y,function(x){x-mean(x) })
  
  
  if(is.null(b)) b <- lapply(1:nt,function(x){rep(0,m)}) 
  
  if(is.null(models)) {
    models <- rep(list(0:1), nt)
    models <- t(do.call(expand.grid, models))
    models <- split(models, rep(1:ncol(models), each = nrow(models)))
    pi <- c(0.999,rep(0.001,length(models)-1)) 
  } 
  if(is.character(models)) {
    if(models=="restrictive") {
      models <- list(rep(0,nt),rep(1,nt))
      pi <- c(0.999,0.001)
    }
  }
  
  #if(is.null(pi)) {
  #     pi <- c(0.999,rep(0.001,length(models)-1)) 
  #}
  
  if(is.null(h2)) h2 <- 0.5
  if(is.null(vare)) {
    vare <- diag(sapply(y,var))
  }
  if(is.null(varb)) varb <- diag(sapply(y,var)/(m*pi[length(models)]))*h2
  #if(is.null(varb)) varb <- (vare*h2)/(m*pi[length(models)])
  if(is.null(vara)) vara <- diag(diag(vare))*h2
  
  #if(is.null(ssb_prior)) ssb_prior <-  diag((nub-2.0)/nub * (vara/(m*pi[length(models)])))
  #if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/(m*pi[length(models)]))
  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (vara/m)
  #if(is.null(sse_prior)) sse_prior <- diag(nue*vare)
  if(is.null(sse_prior)) sse_prior <- nue*diag(diag(vare))
  
  
  fit <- .Call("_qgg_mtbayes",
               y=y, 
               W=split(W, rep(1:ncol(W), each = nrow(W))), 
               b=b,
               B = varb,
               E = vare,
               ssb_prior=split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
               sse_prior=split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
               models=models,
               pi=pi,
               nub=nub,
               nue=nue,
               updateB=updateB,
               updateE=updateE,
               updatePi=updatePi,
               nit=nit,
               method=as.integer(method))
  fit[[6]] <- matrix(unlist(fit[[6]]), ncol = nt, byrow = TRUE)
  fit[[7]] <- matrix(unlist(fit[[7]]), ncol = nt, byrow = TRUE)
  colnames(fit[[6]]) <- rownames(fit[[6]]) <- paste0("T",1:nt)
  colnames(fit[[7]]) <- rownames(fit[[7]]) <- paste0("T",1:nt)
  fit[[11]] <- matrix(unlist(fit[[11]]), ncol = nt, byrow = TRUE)
  fit[[12]] <- matrix(unlist(fit[[12]]), ncol = nt, byrow = TRUE)
  colnames(fit[[11]]) <- rownames(fit[[11]]) <- paste0("T",1:nt)
  colnames(fit[[12]]) <- rownames(fit[[12]]) <- paste0("T",1:nt)
  fit[[13]] <- fit[[13]][[1]]
  fit[[14]] <- fit[[14]][[1]]
  fit[[15]] <- fit[[15]][[1]]
  names(fit) <- c("bm","dm","mu","Bm","Em","rg","re","g","e","b","B","E","pi","pim","order")
  
  return(fit)
  
}




#' Plot fit from gbayes
#'
#' @description
#' Summary plots from gbayes fit

#' @param fit object from gbayes

#'
#' @export
#'

plotBayes <- function(fit=NULL, causal=NULL) {
     if(!is.list(fit[[1]])) {
          layout(matrix(1:6,nrow=3,ncol=2))
          plot(fit[[1]],ylab="Posterior", xlab="Marker", main="Marker effect", frame.plot=FALSE)  
          if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch=4, cex=2, lwd=3 )
          #plot(fit[[2]],ylab="Posterior mean", xlab="Marker", pch="✈",,main="Marker indicator", frame.plot=FALSE)  
          plot(fit[[2]],ylab="Posterior mean", xlab="Marker", main="Marker indicator", ylim=c(0,1), frame.plot=FALSE)  
          if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch=4, cex=2, lwd=3 )
          #if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch="✈", cex=2, lwd=3 )
          hist(fit[[6]],xlab="Posterior", main="Pi")  
          hist(fit[[4]],xlab="Posterior", main="Marker variance")  
          hist(fit[[5]],xlab="Posterior", main="Residual variance")  
          plot(fit[[6]],xlab="Sample", ylab="Pi", frame.plot=FALSE)  
     } 
     if(is.list(fit[[1]])) {
          layout(matrix(1:4,2,2))
          matplot(as.data.frame(fit[[1]]),ylab="Marker effect", frame.plot=FALSE)  
          if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="green", pch=4, cex=2, lwd=3 )
          matplot(as.data.frame(fit[[2]]),ylab="Marker indicator", frame.plot=FALSE)  
          if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="green", pch=4, cex=2, lwd=3 )
          matplot(as.data.frame(fit[[4]]),ylab="Marker variance", frame.plot=FALSE)  
          matplot(as.data.frame(fit[[5]]),ylab="Residual variance", frame.plot=FALSE)  
     } 
     
}

#'
#' @export
#'

plotCvs <- function(fit=NULL, causal=NULL) {
  
  layout(matrix(1:length(fit$p),ncol=1))
  for (i in 1:length(fit$p)) {
    plot(-log10(fit$p[[i]]),ylab="mlogP", xlab="Marker", main="Marker effect", frame.plot=FALSE)  
    if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch=4, cex=2, lwd=3 )
  }
}


gsim <- function(nt=1,W=NULL,n=1000,m=1000) {
  if(is.null(W)) {
    W <- matrix(runif(n),ncol=1)
    for (i in 2:m) {
      W <- cbind(W,scale(W[,i-1]) + runif(n))  
    }
  }
  y <- e <- vector(length=nt,mode="list")
  names(y) <- paste0("D",1:nt)
  set0 <- sample(1:ncol(W),2)
  set1 <- b1 <- g1 <- vector(length=nt,mode="list")
  g <- NULL
  for (i in 1:nt){
    b0 <- sample(c(0.5,-0.5,1,-1),2)
    g0 <- W[,set0]%*%b0
    set1[[i]] <- sample(1:ncol(W),2)
    b1[[i]] <- sample(c(0.5,-0.5,1,-1),length(set1[[i]]))
    g1[[i]] <- W[,set1[[i]]]%*%b1[[i]]
    e[[i]] <- rnorm(nrow(W),mean=0,sd=1)
    y[[i]] <- g0+g1[[i]]+e[[i]]
    g <- cbind(g,g0+g1[[i]])
  }
  colnames(g) <- paste0("D",1:nt) 
  return( list( y=y,W=W, e=e,g=g,b0=b0,b1=b1,set0=set0,set1=set1,causal=c(set0,unlist(set1))))
}


splitWithOverlap <- function(vec, seg.length, overlap) {
     starts = seq(1, length(vec), by=seg.length-overlap)
     ends   = starts + seg.length - 1
     ends[ends > length(vec)] = length(vec)
     lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

neff <- function(seb=NULL,af=NULL,Vy=1) {
     seb2 <- seb**2
     vaf <- 2*af*(1-af)
     neff <- round(median(Vy/(vaf*seb2)))
     return(neff)
}

sortedSets <- function(o = NULL, msize = 500) {
  m <- length(o)
  sets <- vector(length=m,mode="list") 
  for (i in 1:m){
    sets[[i]] <- c((i-msize):(i-1),i:(i+msize))
    sets[[i]] <- sets[[i]][sets[[i]]>0]
    sets[[i]] <- sets[[i]][sets[[i]]<(m+1)]
  }
  
  indices <- 1:m
  keep <- rep(TRUE,m)
  osets <- NULL
  for (i in 1:m) {
    if(keep[o[i]]){
      keep[sets[[ o[i] ]]] <- FALSE
      osets <- c(osets,o[i])
    }
  }
  if(any(!indices%in%unique(unlist(sets[osets])))) stop("Some markers not in sets")
  return(sets[osets])
}
