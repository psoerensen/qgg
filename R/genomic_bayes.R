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

gbayes <- function(y=NULL, X=NULL, W=NULL, stat=NULL, covs=NULL, trait=NULL, fit=NULL, Glist=NULL, chr=NULL, rsids=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                   varg=NULL, varb=NULL, vare=NULL, ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=TRUE,
                   h2=NULL, pi=0.001, updateB=TRUE, updateE=TRUE, updatePi=TRUE, models=NULL,
                   nub=4, nue=4, nit=100, nit_local=NULL,nit_global=NULL,
                   method="mixed", algorithm="default") {
     
     methods <- c("blup","mixed","bayesA","blasso","bayesC","ssvs")
     method <- match(method, methods) - 1
     if( !sum(method%in%c(0:5))== 1 ) stop("Method specified not valid") 
     
     nt <- 1
     if(is.list(y)) nt <- length(y)

     if(nt==1 && algorithm=="default" && !is.null(W)) fit <- bayes(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
                                                     varg=varg, varb=varb, vare=vare, 
                                                     ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
                                                     h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
                                                     nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm)  
     
     if(nt==1 && algorithm=="sbayes" && is.null(Glist)) fit <- sbayes(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
                                                                      varg=varg, varb=varb, vare=vare, 
                                                                      ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
                                                                      h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
                                                                      nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm)  
     
     if(nt>1 && is.null(Glist)) fit <- mtbayes(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
                                               varg=varg, varb=varb, vare=vare, 
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
     #                        varg=varg, varb=varb, vare=vare, 
     #                        ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=FALSE,
     #                        h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
     #                        nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm)
     #       gset <- crossprod(t(W),fitset$bm)
     #       g <- g + gset
     #       e <- e - gset[rws,]
     #       dm[sets[[i]]] <- fitset$dm
     #       bm[sets[[i]]] <- fitset$bm
     #       varem[[i]] <- fitset$vare
     #       varbm[[i]] <- fitset$varb
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
     
     if( !is.null(y) && nt==1 && algorithm=="sbayes" && !is.null(Glist)) {
       fit <- NULL
       if(is.matrix(y)) ids <- rownames(y)
       if(is.vector(y)) ids <- names(y)
       rws <- match(ids,Glist$ids)
       if(any(is.na(rws))) stop("some elements in names(y) does not match elements in Glist$ids ")       
       n <- length(y)
       
       if(is.null(chr)) chromosomes <- 1:Glist$nchr
       if(!is.null(chr)) chromosomes <- chr
       
       bm <- dm <- fit <- stat <- vector(length=Glist$nchr,mode="list")
       names(bm) <- names(dm) <- names(fit) <- names(stat) <- 1:Glist$nchr
       
       yy <- sum((y-mean(y))**2)
       
       if(is.null(covs)) {
         covs <- vector(length=Glist$nchr,mode="list")
         names(covs) <- 1:Glist$nchr
         for (chr in chromosomes){
           print(paste("Computing summary statistics for chromosome:",chr))
           covs[[chr]] <- cvs(y=y,Glist=Glist,chr=chr)
         }
       } 

       for (chr in chromosomes){
         print(paste("Extract sparse LD matrix for chromosome:",chr))
         #LD <- getSparseLD(Glist = Glist, chr = chr)
         LD <- getSparseLD(Glist = Glist, chr = chr, onebased=FALSE)
         
         #LD$indices <- lapply(LD$indices,function(x){x-1})
         LD$values <- lapply(LD$values,function(x){x*n})
         rsidsLD <- names(LD$values)
         clsLD <- match(rsidsLD,Glist$rsids[[chr]])
         wy <- covs[[chr]][[1]][rsidsLD,"wy"]
         b <- rep(0,length(wy))
         print( paste("Fit",methods[method+1] ,"on chromosome:",chr))
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
                              h2=h2, 
                              lambda=lambda, 
                              varb=varb, 
                              vare=vare, 
                              updateB=updateB, 
                              updateE=updateE, 
                              updatePi=updatePi)
         stat[[chr]] <- data.frame(rsids=rsidsLD,chr=rep(chr,length(rsidsLD)),
                                   pos=Glist$pos[[chr]][clsLD], a1=Glist$a1[[chr]][clsLD],
                                   a2=Glist$a2[[chr]][clsLD], af=Glist$af[[chr]][clsLD],bm=fit[[chr]]$bm)
         rownames(stat[[chr]]) <- rsidsLD
       }
       stat <- do.call(rbind, stat)
       rownames(stat) <- stat$rsids
       fit$stat <- stat
       fit$covs <- covs
     }
  
     if( !is.null(stat) && nt==1 && !is.null(Glist)) {
       fit <- NULL
       
       if(is.null(chr)) chromosomes <- 1:Glist$nchr
       if(!is.null(chr)) chromosomes <- chr
       if(is.null(LD)) LD <- vector(length=Glist$nchr,mode="list")
       
       if(is.data.frame(stat)) {
         # this is for external summary statistics
         nt <- 1
         rsidsLD <- unlist(Glist$rsidsLD)
         b <- wy <- matrix(0,nrow=length(rsidsLD),ncol=nt)
         rownames(b) <- rownames(wy) <- rsidsLD
         trait_names <- "badj"       
         stat <- stat[rownames(stat)%in%rsidsLD,]
         b2 <- stat$b^2
         seb2 <- stat$seb^2
         if(!is.null(stat$dfe)) n <- as.integer(mean(stat$dfe)+2)
         if(!is.null(stat$n)) n <- as.integer(mean(stat$n))
         if(!is.null(stat$ww)) ww <- stat$ww
         if(!is.null(stat$n)) ww <- stat$n
         yy <- (b2 + (n-2)*seb2)*ww
         yy <- mean(yy)
         if(!is.null(stat$wy)) wy[rownames(stat),1] <- stat$wy
         if(is.null(stat$wy)) wy[rownames(stat),1] <- stat$b*stat$n
         if(any(is.na(wy))) stop("Missing values in wy")
       }
       if( !is.data.frame(stat) && is.list(stat)) {
         nt <- ncol(stat$b)
         trait_names <- colnames(stat$b)
         rsidsLD <- unlist(Glist$rsidsLD)
         b <- wy <- matrix(0,nrow=length(rsidsLD),ncol=nt)
         rownames(b) <- rownames(wy) <- rsidsLD
         rws <- rownames(stat$b)%in%rsidsLD 
         b2 <- (stat$b[rws,])^2
         seb2 <- (stat$seb[rws,])^2
         dfe <- stat$dfe[rws,]
         n <- as.integer(colMeans(stat$dfe)+2)
         ww <- stat$ww[rws,]
         yy <- (b2 + (n-2)*seb2)*ww
         yy <- colMeans(yy)
         wy[rownames(stat$wy[rws,]),] <- stat$wy[rws,]
         if(any(is.na(wy))) stop("Missing values in wy")
       }
       
       
       bm <- dm <- fit <- res <- vector(length=Glist$nchr,mode="list")
       names(bm) <- names(dm) <- names(fit) <- names(res) <- 1:Glist$nchr
       for (chr in chromosomes){
         if(is.null(LD[[chr]])) {
           print(paste("Extract sparse LD matrix for chromosome:",chr))
           LD[[chr]] <- getSparseLD(Glist = Glist, chr = chr, onebased=FALSE)
         } 
         rsidsLD <- names(LD[[chr]]$values)
         clsLD <- match(rsidsLD,Glist$rsids[[chr]])
         LD[[chr]]$values <- lapply(LD[[chr]]$values,function(x){x*n})
         
         # b <- ma[[chr]]$b[rsidsLD,trait]
         # seb <- ma[[chr]]$seb[rsidsLD,trait]
         # ww <- ma[[chr]]$ww[rsidsLD,trait]
         # wy <- ma[[chr]]$wy[rsidsLD,trait]
         # dfe <- ma[[chr]]$dfe[rsidsLD,trait]
         # b2 <- b*b
         # seb2 <- seb*seb
         # yy <- ww*b2 + dfe*seb2*ww;
         # yy <- mean(yy)
         # n <- as.integer(mean(dfe)+2)
         # b <- wy <- rep(0,length(rsidsLD))
         # names(b) <- names(wy) <- rsidsLD
         
         #LD$indices <- lapply(LD$indices,function(x){x-1})
         bmchr <- NULL
         for (trait in 1:nt) {
           print( paste("Fit",methods[method+1] ,"on chromosome:",chr))
           fit[[chr]] <- sbayes_sparse(yy=yy[trait], 
                                       wy=wy[rsidsLD,trait],
                                       b=b[rsidsLD,trait], 
                                       LDvalues=LD[[chr]]$values, 
                                       LDindices=LD[[chr]]$indices, 
                                       method=method, 
                                       nit=nit, 
                                       n=n[trait], 
                                       pi=pi,
                                       nue=nue, 
                                       nub=nub, 
                                       h2=h2, 
                                       lambda=lambda, 
                                       varb=varb, 
                                       vare=vare, 
                                       updateB=updateB, 
                                       updateE=updateE, 
                                       updatePi=updatePi)
           bmchr <- cbind(bmchr, fit[[chr]]$bm)
         }
         colnames(bmchr) <- trait_names
         #res[[chr]] <- data.frame(chr=rep(chr,length(rsidsLD)),rsids=rsidsLD,alleles=Glist$a1[[chr]][clsLD], af=Glist$af[[chr]][clsLD],bmchr)
         res[[chr]] <- data.frame(rsids=rsidsLD,chr=rep(chr,length(rsidsLD)),
                                  pos=Glist$pos[[chr]][clsLD], a1=Glist$a1[[chr]][clsLD],
                                  a2=Glist$a2[[chr]][clsLD], af=Glist$af[[chr]][clsLD],bm=bmchr)
         
         # fit[[chr]] <- sbayes_sparse(yy=yy, 
         #                             wy=wy[rsidsLD],
         #                             b=b[rsidsLD], 
         #                             LDvalues=LD[[chr]]$values, 
         #                             LDindices=LD[[chr]]$indices, 
         #                             method=method, 
         #                             nit=nit, 
         #                             n=n, 
         #                             pi=pi,
         #                             nue=nue, 
         #                             nub=nub, 
         #                             h2=h2, 
         #                             lambda=lambda, 
         #                             varb=varb, 
         #                             vare=vare, 
         #                             updateB=updateB, 
         #                             updateE=updateE, 
         #                             updatePi=updatePi)
         # res[[chr]] <- data.frame(chr=rep(chr,length(rsidsLD)),rsids=rsidsLD,alleles=Glist$a1[[chr]][clsLD], af=Glist$af[[chr]][clsLD],bm=fit[[chr]]$bm)
         rownames(res[[chr]]) <- rsidsLD
         LD[[chr]]$values <- NULL
         LD[[chr]]$indices <- NULL
       }
       res <- do.call(rbind, res)
       rownames(res) <- res$rsids
       fit$stat <- res
       fit$covs <- covs
     }
     
     fit$method <- methods[method+1]
          
     return(fit)
}


bayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                   varg=NULL, varb=NULL, vare=NULL, 
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
  if(method>=4 && is.null(varb)) varb <- (vare*h2)/(m*pi)
  if(is.null(lambda)) lambda <- rep(vare/varb,m)
  if(is.null(varg)) varg <- vare*h2
  
  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (varg/m)
  if(is.null(sse_prior)) sse_prior <- nue*vare
  
  if(algorithm=="default") {
    fit <- .Call("_qgg_bayes",
                 y=y, 
                 W=split(W, rep(1:ncol(W), each = nrow(W))), 
                 b=b,
                 lambda = lambda,
                 pi = pi,
                 varg = varg,
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
    ids <- rownames(W)
    names(fit[[1]]) <- names(fit[[2]]) <- names(fit[[10]]) <- colnames(W)
    fit[[7]] <- crossprod(t(W),fit[[10]])[,1]
    names(fit[[7]]) <- names(fit[[8]]) <- ids
    names(fit) <- c("bm","dm","mu","varb","vare","pi","g","e","param","b")

  } 
  
  return(fit)
  
  
}

sbayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                  varg=NULL, varb=NULL, vare=NULL, 
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
  }    
  m <- ncol(LD)
  
  if(is.null(pi)) pi <- 0.001
  
  if(is.null(h2)) h2 <- 0.5
  
  if(is.null(vare)) vare <- 1
  if(method<4 && is.null(varb)) varb <- (vare*h2)/m
  if(method==4 && is.null(varb)) varb <- (vare*h2)/(m*pi)
  if(is.null(lambda)) lambda <- rep(vare/varb,m)
  if(is.null(varg)) varg <- vare*h2
  
  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (varg/m)
  if(is.null(sse_prior)) sse_prior <- nue*vare
  
  if(is.null(b)) b <- rep(0,m)
  
  
  fit <- .Call("_qgg_sbayes",
               wy=wy, 
               LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
               b = b,
               lambda = lambda,
               yy = yy,
               pi = pi,
               varg = varg,
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
  names(fit) <- c("bm","dm","mu","varb","vare","pi","g","e","param","b")

  return(fit)
  
}


sbayes_dense <- function(yy=NULL, wy=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                        varg=NULL, varb=NULL, vare=NULL, 
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
  if(is.null(varg)) varg <- vare*h2
  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (varg/m)
  if(is.null(sse_prior)) sse_prior <- nue*vare
  if(is.null(b)) b <- rep(0,m)
  
  
  fit <- .Call("_qgg_sbayes",
               wy=wy, 
               LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
               b = b,
               lambda = lambda,
               yy = yy,
               pi = pi,
               varg = varg,
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
  names(fit) <- c("bm","dm","mu","varb","vare","pi","g","e","param","b")
  
  return(fit)
  
}

sbayes_sparse <- function(yy=NULL, wy=NULL, b=NULL, badj=NULL, seb=NULL, 
                          LDvalues=NULL,LDindices=NULL, n=NULL,
                          varg=NULL, varb=NULL, vare=NULL, 
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
  if(is.null(varg)) varg <- vare*h2
  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (varg/m)
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
               varg = varg,
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
  names(fit) <- c("bm","dm","mu","varb","vare","pi","g","e","param","b")
  
  return(fit)
  
}

mtbayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                  varg=NULL, varb=NULL, vare=NULL, 
                  ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                  h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                  nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL) {
  
  if(is.list(y)) nt <- length(y)
  if(!is.list(y)) stop("This is not a multiple trait analysis")
  
  if(method==0) {
    # BLUP and we do not estimate parameters
    updateB=FALSE;
    updateE=FALSE;
  }
  
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
  
  if(is.null(h2)) h2 <- 0.5
  if(is.null(vare)) {
    vare <- diag(sapply(y,var))
  }
  if(method<4 && is.null(varb)) varb <- diag(sapply(y,var)/(m))*h2
  if(method==4 && is.null(varb)) varb <- diag(sapply(y,var)/(m*pi[length(models)]))*h2
  
  if(is.null(varg)) varg <- diag(diag(vare))*h2

  if(is.null(ssb_prior)) ssb_prior <-  (nub-2.0)/nub * (varg/m)
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
  trait_names <- names(y)
  ids <- rownames(W)
  if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
  colnames(fit[[6]]) <- rownames(fit[[6]]) <- trait_names
  colnames(fit[[7]]) <- rownames(fit[[7]]) <- trait_names
  fit[[11]] <- matrix(unlist(fit[[11]]), ncol = nt, byrow = TRUE)
  fit[[12]] <- matrix(unlist(fit[[12]]), ncol = nt, byrow = TRUE)
  colnames(fit[[11]]) <- rownames(fit[[11]]) <- trait_names
  colnames(fit[[12]]) <- rownames(fit[[12]]) <- trait_names
  # add colnames/rownames e, g and gm
  # add colnames/rownames rg and covg
  fit[[13]] <- fit[[13]][[1]]
  fit[[14]] <- fit[[14]][[1]]
  fit[[15]] <- fit[[15]][[1]]
  fit[[16]] <- crossprod(t(W),matrix(unlist(fit[[1]]), ncol=nt))
  fit[[17]] <- cov2cor(fit[[6]])
  fit[[18]] <- cov2cor(fit[[7]])
  fit[[19]] <- cov(fit[[16]])
  fit[[20]] <- cov2cor(fit[[19]])
  colnames(fit[[19]]) <- rownames(fit[[19]]) <- trait_names
  colnames(fit[[20]]) <- rownames(fit[[20]]) <- trait_names
  for(i in 1:nt){
    names(fit[[1]][[i]]) <- colnames(W)
    names(fit[[2]][[i]]) <- colnames(W)
    names(fit[[10]][[i]]) <- colnames(W)
    names(fit[[8]][[i]]) <- ids
    names(fit[[9]][[i]]) <- names(y[[i]])
  }
  names(fit[[1]]) <- trait_names
  names(fit[[2]]) <- trait_names
  names(fit[[4]]) <- trait_names
  names(fit[[5]]) <- trait_names
  names(fit[[8]]) <- trait_names
  names(fit[[9]]) <- trait_names
  rownames(fit[[16]]) <- ids
  colnames(fit[[16]]) <- trait_names
  names(fit[[13]]) <- sapply(models,paste,collapse="_")
  names(fit[[14]]) <- sapply(models,paste,collapse="_")
  names(fit) <- c("bm","dm","mu","Bm","Em","covb","cove","g","e","b","varb","vare","pi","pim","order","gm","rb","re","covg","rg")
  
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



#'
#' @export
#'

checkStat <- function(Glist=NULL, stat=NULL, filename=NULL, excludeMAF=0.01, excludeMAFDIFF=0.05, excludeINFO=0.8, 
                      excludeCGAT=TRUE, excludeINDEL=TRUE, excludeDUPS=TRUE, excludeMHC=FALSE) {
  # effect, effect_se, effect_allele, alternative_allele, effect_allele_freq, nobs
  cpra <- paste(unlist(Glist$chr),unlist(Glist$position),unlist(Glist$a1),unlist(Glist$a2), sep="_")
  df <- data.frame(rsids=unlist(Glist$rsids),cpra,
                   chr=unlist(Glist$chr), position=unlist(Glist$position), 
                   a1=unlist(Glist$a1), a2=unlist(Glist$a2),
                   af=unlist(Glist$af))
  rsidsDUPS <- df$rsids[duplicated(df$rsids)]
  df <- df[!df$rsids%in%rsidsDUPS,]
  rsidsLD <- unlist(Glist$rsidsLD)
  df <- df[df$rsids%in%rsidsLD,]
  rownames(df) <- df$rsids
  
  inGlist <- stat$rsids%in%df$rsids
  message(paste("Number of markers in stat also found in bedfiles:", sum(inGlist)))
  
  stat <- stat[inGlist,]
  df <- df[rownames(stat),]
  
  aligned <- stat$effect_allele==df$a1
  message(paste("Number of effect alleles aligned with first allele in bimfiles:", sum(aligned)))
  message(paste("Number of effect alleles not aligned with first allele in bimfiles:", sum(!aligned)))

  if(!is.null(filename)) png(file=filename)
  
  layout(matrix(1:6,ncol=2,byrow=TRUE))
  
  lm(stat$effect_allele_freq[aligned]~ df$af[aligned])
  plot(stat$effect_allele_freq[aligned],df$af[aligned], ylab="AF in Glist (allele matching)",xlab="AF in stat (allele matching)")
  
  lm(stat$effect_allele_freq[!aligned]~ df$af[!aligned])
  plot(stat$effect_allele_freq[!aligned],df$af[!aligned], ylab="AF in Glist (allele not matching)",xlab="AF in stat (allele not matching)")
  
  stat[!aligned,"effect_allele_freq"] <- 1 - stat[!aligned,"effect_allele_freq"]
  effect <- stat[!aligned,"b"]
  effect_allele <- stat[!aligned,"effect_allele"]
  alternative_allele <- stat[!aligned,"alternative_allele"]
  stat[!aligned,"effect_allele"] <- alternative_allele 
  stat[!aligned,"alternative_allele"] <- effect_allele 
  stat[!aligned,"b"] <- -effect 
  
  lm(stat$effect_allele_freq~ df$af)
  plot(stat$effect_allele_freq,df$af, ylab="AF in Glist",xlab="AF in stat (after allele flipped)")
  
  isDUPS <- duplicated(stat$rsids)
  a1 <- df$a1
  a2 <- df$a2
  isAT <- a1=="A" & a2=="T"
  isTA <- a1=="T" & a2=="A"
  isCG <- a1=="C" & a2=="G"
  isGC <- a1=="G" & a2=="C"
  isCGAT <- isAT | isTA | isCG | isGC
  CGTA <- c("C","G","T","A")
  isINDEL <- !((a1%in%CGTA) & (a2%in%CGTA))
  
  largeMAFDIFF <- abs(df$af-stat$effect_allele_freq) > excludeMAFDIFF
  maf <- stat$effect_allele_freq
  maf[maf>0.5] <- 1-maf[maf>0.5]
  lowMAF <- maf < excludeMAF

  message(paste("Number of markers excluded by low MAF:", sum(lowMAF)))
  message(paste("Number of markers excluded by large difference between MAF difference:", sum(largeMAFDIFF)))

  rsidsQC <- lowMAF | largeMAFDIFF

  if(!is.null(stat$info)) {
    lowINFO <- stat$info < excludeINFO
    rsidsQC <- rsidsQC | lowINFO
    message(paste("Number of markers excluded by low INFO score:", sum(lowINFO)))
  }
  
  if(excludeCGAT) {
    rsidsQC <- rsidsQC | isCGAT
    message(paste("Number of markers excluded by ambiguity (CG or AT):", sum(isCGAT)))
  }
  if(excludeDUPS) {
    rsidsQC <- rsidsQC | isDUPS
    message(paste("Number of markers excluded by duplicated rsids", sum(isDUPS)))
  }
  if(excludeINDEL) {
    rsidsQC <- rsidsQC | isINDEL
    message(paste("Number of markers excluded by being INDEL:", sum(isINDEL)))
  }
  
  rsidsQC <- !rsidsQC
  plot(stat$effect_allele_freq[rsidsQC],df$af[rsidsQC], ylab="AF in Glist",xlab="AF in stat (after qc check)")
  stat <- stat[rsidsQC,]
  maf <- stat$effect_allele_freq
  maf[maf>0.5] <- 1-maf[maf>0.5]
  seb <- stat$seb
  plot(y=seb,x=maf, ylab="SEB",xlab="MAF")
  if(!is.null(filename)) dev.off()
  stat$af <- stat$effect_allele_freq  
  stat$alleles <- stat$effect_allele  
  if(is.null(stat$n)) stat$n <- neff(seb=stat$seb,af=stat$af)
  return(stat)
}




# adjStat <- function(Glist=NULL,stat=NULL,filename=NULL, chr=NULL){
#   chromosomes <- chr
#   if(is.null(chromosomes)) chromosomes <- 1:Glist$nchr
#   badj <- NULL
#   for ( chr in chromosomes) {
#     
#     LD <- getSparseLD(Glist = Glist, chr = chr)
#     rsidsLD <- Glist$rsidsLD[[chr]]
#     
#     zobs <- zpred <- rep(0,length(rsidsLD))
#     names(zobs) <- names(zpred) <- rsidsLD
#     rsidsSTAT <- rownames(stat)[rownames(stat)%in%rsidsLD]
#     zobs[rsidsSTAT] <- stat[rsidsSTAT,"b"]
#     #zobs[rsidsSTAT] <- stat[rsidsSTAT,"b"]/stat[rsidsSTAT,"seb"]
#     for (i in 1:length(LD$indices)){
#       #zsum <- sum(zobs[LD$indices[[i]]]*LD$values[[i]])-zobs[i]
#       #nsum <- sum(abs(LD$values[[i]])) - 1
#       zsum <- sum(zobs[LD$indices[[i]]]*LD$values[[i]])
#       nsum <- sum(abs(LD$values[[i]]))
#       #zpred[i] <- zsum/nsum
#       if(!zobs[i]==0.0) zpred[i] <- zsum/nsum
#     }
#     # quantile normalisation (https://academic.oup.com/bioinformatics/article/19/2/185/372664)
#     zobs_rank <- rank(zobs, ties.method = "min")
#     zpred_rank <- rank(zpred, ties.method = "min")
#     zobs_sort <- sort(zobs)
#     zpred_sort <- sort(zpred)
#     zmean <- (zobs_sort+zpred_sort)/2
#     zobs_adj <- zmean[zobs_rank]
#     zpred_adj <- zmean[zpred_rank]
#     
#     if(!is.null(filename[chr])) {
#       png(filename[chr])
#       layout(matrix(1:4,ncol=2,byrow=TRUE))
#       plot(y=zobs,x=zpred, main=paste("Chr",chr), ylab="Observed Z",xlab="Predicted Z")
#       plot(y=zobs_adj,x=zpred_adj, main=paste("Chr",chr), ylab="Observed Z (normalized)",xlab="Predicted Z (normalized)")
#       plot(y=zobs,x=zobs_adj, main=paste("Chr",chr), ylab="Observed Z", xlab="Observed Z (normalized)")
#       plot(y=zobs,x=zpred_adj, main=paste("Chr",chr), ylab="Observed Z", xlab="Predicted Z (normalized)")
#       dev.off()
#     }
#     print(paste("Finished chr:",chr))
#     badj <- c(badj,zpred_adj)
#   }
#   badj <- badj[names(badj)%in%rownames(stat)]
#   stat <- cbind(stat[names(badj),],badj=badj)
#   return(stat)
# }


#'
#' @export
#'

adjStat <- function(stat = NULL, Glist = NULL, chr=NULL, statistics = "b", 
                  r2 = 0.9, ldSets = NULL, threshold = 1, header=NULL,
                  method = "pruning") {
  if(is.data.frame(stat)) {
    p <- stat$p
    if(is.null(stat$p)) pstat <- pnorm(abs(stat$b/stat$seb),lower.tail=FALSE)
    names(p) <- rownames(stat)
    p <- adjLD(Glist=Glist, stat=p, r2=r2, threshold=threshold)
    p[p>0] <- 1
    if(is.null(header)) header <- c("rsids","chr","pos","allele","a1","a2","af")
    
    if(statistics=="b") {
      b <- stat[rownames(p),"b"]
      badj <- p*stat[rownames(p),"b"]
      colnames(badj) <- paste0("b_",threshold)
      if(any(colnames(stat)%in%header)) statadj <- data.frame(stat[rownames(badj),colnames(stat)%in%header],b,badj)
      if(!any(colnames(stat)%in%header)) statadj <- as.matrix(data.frame(b,badj))
      return(statadj)
    }
    if(statistics=="z") {
      z <- stat[rownames(p),"b"]/stat[rownames(p),"seb"]
      zadj <- p*stat[rownames(p),"b"]/stat[rownames(p),"seb"]
      colnames(zadj) <- paste0("z_",threshold)
      if(any(colnames(stat)%in%header)) statadj <- data.frame(stat[rownames(zadj),colnames(stat)%in%header],z,zadj)
      if(!any(colnames(stat)%in%header)) statadj <- as.matrix(data.frame(z,zadj))
      return(statadj)
    }
  }
}

