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

gbayes <- function(y=NULL, X=NULL, W=NULL, stat=NULL, covs=NULL, trait=NULL, fit=NULL, Glist=NULL, 
                   chr=NULL, rsids=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                   vg=NULL, vb=NULL, ve=NULL, ssg_prior=NULL, ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=TRUE,
                   h2=NULL, pi=0.001, updateB=TRUE, updateE=TRUE, updatePi=TRUE, models=NULL,
                   nug=NULL, nub=4, nue=4, verbose=FALSE,msize=100,
                   GRMlist=NULL, ve_prior=NULL, vg_prior=NULL,tol=0.001,
                   nit=100, nburn=0, nit_local=NULL,nit_global=NULL,
                   method="mixed", algorithm="default") {
  
  methods <- c("blup","mixed","bayesA","blasso","bayesC","ssvs")
  method <- match(method, methods) - 1
  if( !sum(method%in%c(0:5))== 1 ) stop("Method specified not valid") 
  
  nt <- 1
  if(!is.null(y)) {
    if(is.list(y)) nt <- length(y)
    if(is.matrix(y)) nt <- ncol(y)
  }
  if(!is.null(stat)) {
    if(!is.data.frame(stat) && is.list(stat)) nt <- ncol(stat$b)
  }
  
  # define type of analysis
  if(!is.null(GRMlist)) analysis <- "mtmc-mixed"
  
  if(nt==1 && !is.null(y) && !is.null(W) && algorithm=="default") 
    analysis <- "st-blr-individual-level-default"
  
  if(nt==1 && !is.null(y) && !is.null(W) && algorithm=="sbayes") 
    analysis <- "st-blr-individual-level-sbayes"

  if(nt==1 && !is.null(y) && algorithm=="sparse") 
    analysis <- "st-blr-individual-level-sparse-ld"
  
  if( nt==1 && !is.null(y) &&  algorithm=="dense") 
    analysis <- "st-blr-individual-level-dense-ld"
  
  if( nt==1 && is.null(y) && !is.null(stat) && !is.null(Glist)) 
    analysis <- "st-blr-sumstat-sparse-ld"
  
  if(nt>1 && !is.null(y) && !is.null(W))
    analysis <- "mt-blr-individual-level"
  
  if( nt>1 && !is.null(stat) && !is.null(Glist) && algorithm=="default") 
    analysis <- "mt-blr-sumstat-sparse-ld"

  if( nt>1 && !is.null(stat) && !is.null(Glist) && algorithm=="serial") 
    analysis <- "st-blr-sumstat-sparse-ld"
  
  message(paste("Type of analysis chosen:",analysis))  

  # Single and multiple trait BLR based on GRMs    
  if(!is.null(GRMlist)) {
    fit <- bmm(y=y, X=X, W=W, GRMlist=GRMlist,
               vg=vg, ve=ve, nug=nug, nue=nue,
               vg_prior=vg_prior, ve_prior=ve_prior,
               updateG=updateB, updateE=updateE,
               nit=nit, nburn=nburn, tol=tol, verbose=verbose) 
  } 
  
  
  # Single trait BLR using y and W   
  if(nt==1 && !is.null(y) && !is.null(W) && algorithm=="default") {
    
    fit <- bayes(y=y, X=X, W=W, b=b, 
                 badj=badj, seb=seb, LD=LD, n=n,
                 vg=vg, vb=vb, ve=ve, 
                 ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
                 h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
                 nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm)  
  }
  
  
  # Single trait BLR using y and W and sbayes method 
  if(nt==1 && !is.null(y) && !is.null(W) && algorithm=="sbayes") {
    
    fit <- sbayes(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
                  vg=vg, vb=vb, ve=ve, 
                  ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
                  h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
                  nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm)  
  }
  
  # Multiple trait BLR using y and W
  if(nt>1 && !is.null(y) && !is.null(W)) {
    fit <- mtbayes(y=y, X=X, W=W, b=b, badj=badj, seb=seb, LD=LD, n=n,
                   vg=vg, vb=vb, ve=ve, 
                   ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
                   h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
                   nub=nub, nue=nue, nit=nit, method=method, algorithm=algorithm) 
  }
  
  
  
  # Single trait BLR using y and sparse LD provided Glist
  if( nt==1 && !is.null(y) && algorithm=="sparse") {
    
    if(is.null(Glist)) stop("Please provide Glist")
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
    
    if(is.null(nit_local)) nit_local <- nit
    if(is.null(nit_global)) nit_global <- 1
    
    for (it in 1:nit_global) {
      for (chr in chromosomes){
        if(verbose) print(paste("Extract sparse LD matrix for chromosome:",chr))
        LD <- getSparseLD(Glist = Glist, chr = chr, onebased=FALSE)
        LD$values <- lapply(LD$values,function(x){x*n})
        rsidsLD <- names(LD$values)
        clsLD <- match(rsidsLD,Glist$rsids[[chr]])
        wy <- covs[[chr]][rsidsLD,"wy"]
        b <- rep(0,length(wy))
        if(it>1) {
          b <- fit[[chr]]$b
          if(updateB) vb <- fit[[chr]]$param[1]
          if(updateE) ve <- fit[[chr]]$param[2]
          if(updatePi) pi <- fit[[chr]]$param[3]
        }
        if(verbose) print( paste("Fit",methods[method+1] ,"on chromosome:",chr))
        fit[[chr]] <- sbayes_sparse(yy=yy, 
                                    wy=wy,
                                    b=b, 
                                    LDvalues=LD$values, 
                                    LDindices=LD$indices, 
                                    method=method, 
                                    nit=nit_local, 
                                    n=n, 
                                    pi=pi,
                                    nue=nue, 
                                    nub=nub, 
                                    h2=h2, 
                                    lambda=lambda, 
                                    vb=vb, 
                                    ve=ve, 
                                    updateB=updateB, 
                                    updateE=updateE, 
                                    updatePi=updatePi)
        stat[[chr]] <- data.frame(rsids=rsidsLD,chr=rep(chr,length(rsidsLD)),
                                  pos=Glist$pos[[chr]][clsLD], a1=Glist$a1[[chr]][clsLD],
                                  a2=Glist$a2[[chr]][clsLD], af=Glist$af[[chr]][clsLD],bm=fit[[chr]]$bm)
        rownames(stat[[chr]]) <- rsidsLD
      }
    }
    stat <- do.call(rbind, stat)
    rownames(stat) <- stat$rsids
    fit$stat <- stat
    fit$covs <- covs
  }
  
  
  # Single trait BLR using y and dense LD
  if( nt==1 && !is.null(y) &&  algorithm=="dense") {
    
    overlap <- 0
    
    if(is.null(Glist)) stop("Please provide Glist")
    fit <- NULL
    if(is.matrix(y)) ids <- rownames(y)
    if(is.vector(y)) ids <- names(y)
    rws <- match(ids,Glist$ids)
    if(any(is.na(rws))) stop("some elements in names(y) does not match elements in Glist$ids ")       
    n <- length(y)
    
    if(is.null(chr)) chromosomes <- 1:Glist$nchr
    if(!is.null(chr)) chromosomes <- chr
    
    rsids <- unlist(Glist$rsidsLD)
    cls <- lapply(Glist$rsids,function(x) { 
      splitWithOverlap(na.omit(match(rsids,x)),msize,0)})
    vblist <- lapply(sapply(cls,length),function(x) 
    {vector(length=x, mode="numeric")})
    velist <- lapply(sapply(cls,length),function(x) 
    {vector(length=x, mode="numeric")})
    pilist <- lapply(sapply(cls,length),function(x) 
    {vector(length=x, mode="numeric")})
    b <- lapply(Glist$mchr,function(x){rep(0,x)})
    bm <- lapply(Glist$mchr,function(x){rep(0,x)})
    dm <- lapply(Glist$mchr,function(x){rep(0,x)})
    
    if(is.null(nit_local)) nit_local <- nit
    if(is.null(nit_global)) nit_global <- 1
    
    for (it in 1:nit_global) {
      e <- y-mean(y)
      yy <- sum(e**2)
      for (chr in 1:length(Glist$nchr)) {
        for (i in 1:length(cls[[chr]])) {
          wy <- computeWy(y=e,Glist=Glist,chr=chr,cls=cls[[chr]][[i]])
          WW <- computeWW(Glist=Glist, chr=chr, cls=cls[[chr]][[i]], rws=rws)
          if(it>1) {
            if(updateB) vb <- vblist[[chr]][i]
            if(updateE) ve <- velist[[chr]][i]
            if(updatePi) pi <- pilist[[chr]][i]
          }
          fitS <- computeB(wy=wy, yy=yy, WW=WW, n=n,
                           b=b[[chr]][cls[[chr]][[i]]],
                           ve=ve, vb=vb, pi=pi,
                           nub=nub, nue=nue,
                           updateB=updateB, updateE=updateE, updatePi=updatePi,
                           nit=nit, nburn=nburn, method=method) 
          b[[chr]][cls[[chr]][[i]]] <- fitS$b
          bm[[chr]][cls[[chr]][[i]]] <- fitS$bm
          dm[[chr]][cls[[chr]][[i]]] <- fitS$dm
          vblist[[chr]][i] <- fitS$param[1]
          velist[[chr]][i] <- fitS$param[2]
          pilist[[chr]][i] <- fitS$param[3]
          grs <- computeGRS(Glist = Glist, chr = chr, 
                            cls = cls[[chr]][[i]], 
                            b=bm[[chr]][cls[[chr]][[i]]])  
          e <- e - grs[rws,]
        }
      }
    }   
    bm <- unlist(bm)
    dm <- unlist(dm)
    names(bm) <- names(dm) <- unlist(Glist$rsids)
    stat <- data.frame(rsids=rsids,
                       chr=unlist(Glist$chr)[rsids],
                       pos=unlist(Glist$pos)[rsids], 
                       a1=unlist(Glist$a1)[rsids],
                       a2=unlist(Glist$a2)[rsids], 
                       af=unlist(Glist$af)[rsids],
                       bm=bm[rsids])
    fit$stat <- stat
    
  }
  
  # Single trait BLR using summary statistics and sparse LD provided in Glist
#  if( nt==1 && is.null(y) && !is.null(stat) && !is.null(Glist)) {
  if(analysis=="st-blr-sumstat-sparse-ld") {
    # single trait summary statistics
    if(is.data.frame(stat)) {
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
    
    # multiple trait summary statistics
    if( !is.data.frame(stat) && is.list(stat)) {
      nt <- ncol(stat$b)
      trait_names <- colnames(stat$b)
      if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
      rsidsLD <- unlist(Glist$rsidsLD)
      b <- wy <- matrix(0,nrow=length(rsidsLD),ncol=nt)
      rownames(b) <- rownames(wy) <- rsidsLD
      colnames(b) <- colnames(wy) <- trait_names
      rws <- rownames(stat$b)%in%rsidsLD 
      b2 <- (stat$b[rws,])^2
      seb2 <- (stat$seb[rws,])^2
      if(!is.null(stat$dfe)) n <- as.integer(colMeans(stat$dfe)+2)
      if(!is.null(stat$n)) n <- as.integer(colMeans(stat$n))
      ww <- stat$ww[rws,]
      for (i in 1:nt) {
        seb2[,i] <- (n[i]-2)*seb2[,i]
      }
      yy <- (b2 + seb2)*ww
      yy <- colMeans(yy)
      wy[rownames(stat$wy[rws,]),] <- stat$wy[rws,]
      if(any(is.na(wy))) stop("Missing values in wy")
    }
    
    if(is.null(chr)) chromosomes <- 1:Glist$nchr
    if(!is.null(chr)) chromosomes <- chr
    if(is.null(LD)) LD <- vector(length=Glist$nchr,mode="list")
    
    bm <- dm <- fit <- res <- vector(length=Glist$nchr,mode="list")
    names(bm) <- names(dm) <- names(fit) <- names(res) <- 1:Glist$nchr
    for (chr in chromosomes){
      if(is.null(LD[[chr]])) {
        if(verbose) print(paste("Extract sparse LD matrix for chromosome:",chr))
        LD[[chr]] <- getSparseLD(Glist = Glist, chr = chr, onebased=FALSE)
      } 
      rsidsLD <- names(LD[[chr]]$values)
      clsLD <- match(rsidsLD,Glist$rsids[[chr]])
      bmchr <- NULL
      for (trait in 1:nt) {
        if(verbose) print( paste("Fit",methods[method+1], "on chromosome:",chr,"for trait",trait))
        LDvalues <- lapply(LD[[chr]]$values,function(x){x*n[trait]})
        fit[[chr]] <- sbayes_sparse(yy=yy[trait], 
                                    wy=wy[rsidsLD,trait],
                                    b=b[rsidsLD,trait], 
                                    LDvalues=LDvalues, 
                                    LDindices=LD[[chr]]$indices, 
                                    method=method, 
                                    nit=nit, 
                                    n=n[trait], 
                                    pi=pi,
                                    nue=nue, 
                                    nub=nub, 
                                    h2=h2, 
                                    lambda=lambda, 
                                    vb=vb, 
                                    ve=ve, 
                                    updateB=updateB, 
                                    updateE=updateE, 
                                    updatePi=updatePi)
        bmchr <- cbind(bmchr, fit[[chr]]$bm)
      }
      colnames(bmchr) <- trait_names
      res[[chr]] <- data.frame(rsids=rsidsLD,chr=rep(chr,length(rsidsLD)),
                               pos=Glist$pos[[chr]][clsLD], a1=Glist$a1[[chr]][clsLD],
                               a2=Glist$a2[[chr]][clsLD], af=Glist$af[[chr]][clsLD],bm=bmchr)
      rownames(res[[chr]]) <- rsidsLD
      LD[[chr]]$values <- NULL
      LD[[chr]]$indices <- NULL
    }
    res <- do.call(rbind, res)
    rownames(res) <- res$rsids
    fit$stat <- res
    fit$method <- methods[method+1]
  }
  
  # Multi trait BLR using summary statistics and sparse LD provided in Glist
#  if( nt>1 && is.null(y) && !is.null(stat) && !is.null(Glist)) {

  if(analysis=="mt-blr-sumstat-sparse-ld") {
  
    
    # multiple trait summary statistics
    trait_names <- colnames(stat$b)
    if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
    rsidsLD <- unlist(Glist$rsidsLD)
    b <- wy <- matrix(0,nrow=length(rsidsLD),ncol=nt)
    rownames(b) <- rownames(wy) <- rsidsLD
    colnames(b) <- colnames(wy) <- trait_names
    rws <- rownames(stat$b)%in%rsidsLD
    b2 <- (stat$b[rws,])^2
    seb2 <- (stat$seb[rws,])^2
    if(!is.null(stat$dfe)) n <- as.integer(colMeans(stat$dfe)+2)
    if(!is.null(stat$n)) n <- as.integer(colMeans(stat$n))
    ww <- stat$ww[rws,]
    for (i in 1:nt) {
      seb2[,i] <- (n[i]-2)*seb2[,i]
    }
    yy <- (b2 + seb2)*ww
    yy <- colMeans(yy)
    wy[rownames(stat$wy[rws,]),] <- stat$wy[rws,]
    if(any(is.na(wy))) stop("Missing values in wy")
    
    if(is.null(chr)) chromosomes <- 1:Glist$nchr
    if(!is.null(chr)) chromosomes <- chr
    if(is.null(LD)) LD <- vector(length=Glist$nchr,mode="list")
    
    bm <- dm <- fit <- res <- vector(length=Glist$nchr,mode="list")
    names(bm) <- names(dm) <- names(fit) <- names(res) <- 1:Glist$nchr
    for (chr in chromosomes){
      if(is.null(LD[[chr]])) {
        if(verbose) print(paste("Extract sparse LD matrix for chromosome:",chr))
        LD[[chr]] <- getSparseLD(Glist = Glist, chr = chr, onebased=FALSE)
      }
      rsidsLD <- names(LD[[chr]]$values)
      clsLD <- match(rsidsLD,Glist$rsids[[chr]])
      bmchr <- NULL
      fit[[chr]] <- mt_sbayes_sparse(yy=yy,
                                     wy=wy[rsidsLD,],
                                     b=b[rsidsLD,],
                                     LDvalues=LD[[chr]]$values,
                                     LDindices=LD[[chr]]$indices,
                                     n=n,
                                     nit=nit,
                                     pi=pi,
                                     nue=nue,
                                     nub=nub,
                                     h2=h2,
                                     vb=vb,
                                     ve=ve,
                                     ssb_prior=ssb_prior,
                                     sse_prior=sse_prior,
                                     updateB=updateB,
                                     updateE=updateE,
                                     updatePi=updatePi,
                                     method=method)
      res[[chr]] <- data.frame(rsids=rsidsLD,chr=rep(chr,length(rsidsLD)),
                               pos=Glist$pos[[chr]][clsLD], a1=Glist$a1[[chr]][clsLD],
                               a2=Glist$a2[[chr]][clsLD], af=Glist$af[[chr]][clsLD],
                               bm=fit[[chr]]$bm)
      rownames(res[[chr]]) <- rsidsLD
      LD[[chr]]$values <- NULL
      LD[[chr]]$indices <- NULL
    }
    res <- do.call(rbind, res)
    rownames(res) <- res$rsids
    fit$stat <- res
    fit$method <- methods[method+1]
    
  }
  return(fit)
  
}

# Individual level data using W 
bayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                   vg=NULL, vb=NULL, ve=NULL, 
                   ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                   h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                   nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL) {
  ids <- NULL
  if(is.matrix(y)) ids <- rownames(y)
  if(is.vector(y)) ids <- names(y)
  if(scaleY) y <- as.vector(scale(y)) 
  n <- nrow(W)
  m <- ncol(W)
  
  if(is.null(ids)) warning("No names/rownames provided for y")
  if(is.null(rownames(W))) warning("No names/rownames provided for W")
  if(!is.null(ids) & !is.null(rownames(W))) {
    if(any(is.na(match(ids,rownames(W))))) stop("Names/rownames for y does match rownames for W")
  }
  if(is.null(pi)) pi <- 0.001
  if(is.null(h2)) h2 <- 0.5
  e=y-mean(y)
  if(is.null(ve)) ve <- var(e)
  if(method<4 && is.null(vb)) vb <- (ve*h2)/m
  if(method>=4 && is.null(vb)) vb <- (ve*h2)/(m*pi)
  if(is.null(lambda)) lambda <- rep(ve/vb,m)
  if(is.null(vg)) vg <- ve*h2
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m*pi)
  if(is.null(sse_prior)) sse_prior <- nue*ve
  if(is.null(b)) b <- rep(0,m)
  
  if(algorithm=="default") {
    fit <- .Call("_qgg_bayes",
                 y=y, 
                 W=split(W, rep(1:ncol(W), each = nrow(W))), 
                 b=b,
                 lambda = lambda,
                 pi = pi,
                 vg = vg,
                 vb = vb,
                 ve = ve,
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
    names(fit) <- c("bm","dm","coef","vb","ve","pi","g","e","param","b")
  } 
  return(fit)
}

# Individual level data using W 
sbayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                  vg=NULL, vb=NULL, ve=NULL, 
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
  
  if(is.null(ve)) ve <- 1
  if(method<4 && is.null(vb)) vb <- (ve*h2)/m
  if(method>=4 && is.null(vb)) vb <- (ve*h2)/(m*pi)
  if(is.null(lambda)) lambda <- rep(ve/vb,m)
  if(is.null(vg)) vg <- ve*h2
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m*pi)
  if(is.null(sse_prior)) sse_prior <- nue*ve
  if(is.null(b)) b <- rep(0,m)

  fit <- .Call("_qgg_sbayes",
               wy=wy, 
               LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
               b = b,
               lambda = lambda,
               yy = yy,
               pi = pi,
               vg = vg,
               vb = vb,
               ve = ve,
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
  names(fit) <- c("bm","dm","coef","vb","ve","pi","g","e","param","b")

  return(fit)
  
}

sbayes_sparse <- function(yy=NULL, wy=NULL, b=NULL, badj=NULL, seb=NULL, 
                          LDvalues=NULL,LDindices=NULL, n=NULL,
                          vg=NULL, vb=NULL, ve=NULL, 
                          ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                          h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                          nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL) {
  m <- length(LDvalues)
  if(is.null(pi)) pi <- 0.001
  if(is.null(h2)) h2 <- 0.5
  if(is.null(ve)) ve <- 1
  if(method<4 && is.null(vb)) vb <- (ve*h2)/m
  if(method>=4 && is.null(vb)) vb <- (ve*h2)/(m*pi)
  if(is.null(lambda)) lambda <- rep(ve/vb,m)
  if(is.null(vg)) vg <- ve*h2
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m*pi)
  if(is.null(sse_prior)) sse_prior <- nue*ve
  if(is.null(b)) b <- rep(0,m)

  fit <- .Call("_qgg_sbayes_spa",
               wy=wy, 
               LDvalues=LDvalues, 
               LDindices=LDindices, 
               b = b,
               lambda = lambda,
               yy = yy,
               pi = pi,
               vg = vg,
               vb = vb,
               ve = ve,
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
  names(fit) <- c("bm","dm","coef","vb","ve","pi","g","e","param","b")
  return(fit)
}


mt_sbayes_sparse <- function(yy=NULL, wy=NULL, b=NULL, 
                             LDvalues=NULL,LDindices=NULL, n=NULL,
                             vg=NULL, vb=NULL, ve=NULL, 
                             ssb_prior=NULL, sse_prior=NULL, 
                             h2=NULL, pi=NULL, updateB=NULL, 
                             updateE=NULL, updatePi=NULL, models=NULL,
                             nub=NULL, nue=NULL, nit=NULL, method=NULL) {
  
  nt <- length(yy)
  m <- length(LDvalues)
  
  if(method==0) {
    # BLUP and we do not estimate parameters
    updateB=FALSE;
    updateE=FALSE;
  }
  
  if(method==2) stop("Multiple trait not yet implemented for Bayes A") 
  if(method==3) stop("Multiple trait not yet implemented for Bayesian Lasso")
  
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
  if(is.null(ve)) ve <- diag(yy/n,nt)
  if(method<4 && is.null(vb)) vb <- diag(h2/m,nt)
  if(method>=4 && is.null(vb)) vb <- diag(h2/(m*pi[length(models)]),nt)
  if(is.null(vg)) vg <- diag(diag(ve))*h2
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m*0.001)
  if(is.null(sse_prior)) sse_prior <- nue*diag(diag(ve))
  
  fit <- .Call("_qgg_mtsbayes",
               wy=split(wy, rep(1:ncol(wy), each = nrow(wy))),
               yy=yy,
               b = split(b, rep(1:ncol(b), each = nrow(b))),
               LDvalues=LDvalues, 
               LDindices=LDindices, 
               B = vb,
               E = ve,
               ssb_prior=split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
               sse_prior=split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
               models=models,
               pi=pi,
               nub=nub,
               nue=nue,
               updateB = updateB,
               updateE = updateE,
               updatePi = updatePi,
               n=n,
               nit=nit,
               method=as.integer(method))
  
  fit[[6]] <- matrix(unlist(fit[[6]]), ncol = nt, byrow = TRUE)
  fit[[7]] <- matrix(unlist(fit[[7]]), ncol = nt, byrow = TRUE)
  trait_names <- names(yy)
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
  fit[[16]] <- fit[[6]]
  fit[[17]] <- fit[[7]]
  if(sum(diag(fit[[16]]))>0) fit[[16]] <- cov2cor(fit[[16]])
  if(sum(diag(fit[[17]]))>0)  fit[[17]] <- cov2cor(fit[[17]])
  for(i in 1:nt){
    names(fit[[1]][[i]]) <- names(LDvalues)
    names(fit[[2]][[i]]) <- names(LDvalues)
    names(fit[[10]][[i]]) <- names(LDvalues)
  }
  names(fit[[1]]) <- trait_names
  names(fit[[2]]) <- trait_names
  names(fit[[3]]) <- trait_names
  names(fit[[4]]) <- trait_names
  names(fit[[5]]) <- trait_names
  names(fit[[8]]) <- trait_names
  names(fit[[9]]) <- trait_names
  names(fit[[13]]) <- sapply(models,paste,collapse="_")
  names(fit[[14]]) <- sapply(models,paste,collapse="_")
  
  names(fit) <- c("bm","dm","coef","vbs","ves","covb","cove",
                  "wy","r","b","vb","ve","pi","pim","order",
                  "rb","re")
  fit$bm <- as.matrix(as.data.frame(fit$bm))
  fit$dm <- as.matrix(as.data.frame(fit$dm))
  fit$b <- as.matrix(as.data.frame(fit$b))
  return(fit)
}


mtbayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                  vg=NULL, vb=NULL, ve=NULL, 
                  ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                  h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                  nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL) {
  
  if(is.list(y)) nt <- length(y)
  if(!is.matrix(y)) stop("This is not a multiple trait analysis")
  
  if(method==0) {
    # BLUP and we do not estimate parameters
    updateB=FALSE;
    updateE=FALSE;
  }
  
  if(method==2) stop("Multiple trait not yet implemented for Bayes A") 
  if(method==3) stop("Multiple trait not yet implemented for Bayesian Lasso")
  
  n <- nrow(W)
  m <- ncol(W)
  
  if(!is.matrix(y)) stop("y should be a matrix")
  y <- split(y, rep(1:ncol(y), each = nrow(y)))
  
  if(scaleY) y <- lapply(y,function(x){as.vector(scale(x))})
  if(!scaleY) y <- lapply(y,function(x){x-mean(x) })
  
  if(is.null(b)) b <- lapply(1:nt,function(x){rep(0,m)}) 
  if(is.matrix(b)) b <- split(b, rep(1:ncol(b), each = nrow(b)))
  
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
  if(is.null(ve)) ve <- diag(sapply(y,var))
  if(method<4 && is.null(vb)) vb <- diag(sapply(y,var)/(m))*h2
  if(method>=4 && is.null(vb)) vb <- diag(sapply(y,var)/(m*pi[length(models)]))*h2
  if(is.null(vg)) vg <- diag(diag(ve))*h2
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m*0.001)
  if(is.null(sse_prior)) sse_prior <- nue*diag(diag(ve))

  fit <- .Call("_qgg_mtbayes",
               y=y, 
               W=split(W, rep(1:ncol(W), each = nrow(W))), 
               b=b,
               B = vb,
               E = ve,
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
  names(fit) <- c("bm","dm","mus","vbs","ves","covb","cove",
                  "g","e","b","vb","ve","pi","pim","order",
                  "gm","rb","re","covg","rg")
  fit$bm <- as.matrix(as.data.frame(fit$bm))
  fit$dm <- as.matrix(as.data.frame(fit$dm))
  return(fit)
  
}


computeB <- function(wy=NULL, yy=NULL, b=NULL, WW=NULL, n=NULL,
                     vb=NULL, vg=NULL, ve=NULL, lambda=NULL, 
                     ssb_prior=NULL, sse_prior=NULL, 
                     nub=NULL, nue=NULL, 
                     h2=NULL, pi=NULL, 
                     updateB=NULL, updateE=NULL, updatePi=NULL,
                     nit=NULL, nburn=NULL, method=NULL) {
  
  m <- ncol(WW)
  
  if(is.null(pi)) pi <- 0.001
  if(is.null(h2)) h2 <- 0.5
  
  if(is.null(ve)) ve <- 1
  if(method<4 && is.null(vb)) vb <- (ve*h2)/m
  if(method>=4 && is.null(vb)) vb <- (ve*h2)/(m*pi)
  if(is.null(lambda)) lambda <- rep(ve/vb,m)
  if(is.null(vg)) vg <- ve*h2
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m*pi)
  if(is.null(sse_prior)) sse_prior <- nue*ve
  if(is.null(b)) b <- rep(0,m)
  
  fit <- .Call("_qgg_sbayes",
               wy=wy, 
               LD=split(WW, rep(1:ncol(WW), each = nrow(WW))), 
               b = b,
               lambda = lambda,
               yy = yy,
               pi = pi,
               vg = vg,
               vb = vb,
               ve = ve,
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
  
  names(fit[[1]]) <- rownames(WW)
  names(fit) <- c("bm","dm","mus","vbs","ves","pis","wy","r","param","b")
  
  return(fit)
  
}

computeWy <- function(y=NULL, Glist = NULL, chr = NULL, cls = NULL) {
  wy <- cvs(y=y,Glist=Glist,chr=chr,cls=cls)$wy
  return(wy)
}

computeWW <- function(Glist = NULL, chr = NULL, cls = NULL, rws=NULL, scale=TRUE) { 
  W <- getG(Glist=Glist, chr=chr, cls=cls, scale=scale)
  WW <- crossprod(W[rws,])
  return(WW)
}

computeGRS <- function(Glist = NULL, chr = NULL, cls = NULL, b=NULL, scale=TRUE) { 
  af <- Glist$af[[chr]][cls]  
  grs <- .Call("_qgg_mtgrsbed", Glist$bedfiles[chr], 
               Glist$n, cls=cls, af=af, scale=scale, Slist=list(b))
  return(as.matrix(as.data.frame(grs)))
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

qcstat <- function(Glist=NULL, stat=NULL, filename=NULL, 
                   excludeMAF=0.01, excludeMAFDIFF=0.05, excludeINFO=0.8, 
                   excludeCGAT=TRUE, excludeINDEL=TRUE, excludeDUPS=TRUE, excludeMHC=FALSE,
                   excludeMISS=0.05, excludeHWE=1e-12) {


  # stat is a data.frame
  if(!is.data.frame(stat)) stop("stat should be  a data frame")
  if(!is.null(stat$marker)) rownames(stat) <- stat$marker
  if(!is.null(stat$rsids)) rownames(stat) <- stat$rsids
  
  # internal summary statistic column format
  # data.frame(rsids, chr, pos, a1, a2, af, b, seb, stat, p, n)     (single trait)
  # list(marker=(rsids, chr, pos, a1, a2, af), b, seb, stat, p, n)  (multiple trait)

  fm_internal <- c("rsids","chr","pos","a1","a2","af","b","seb")
  fm_external <- c("marker","chromosome", "position", "effect_allele", "non_effect_allele", 
  "effect_allele_freq","effect", "effect_se")
  
  format <- "unknown"
  if(all(fm_internal%in%colnames(stat))) format <- "internal"
  if(all(fm_external%in%colnames(stat))) format <- "external"
  if(format=="unknown") {
    message("Column headings for stat object not found")
    message("Column headings for stat object should be:")
    print(fm_external)
    message("or:")
    print(fm_internal)
    stop("please revised your stat object according to these ")
  }

  # external summary statistic column format
  # optimal format:
  # marker, chromosome, position, effect_allele, non_effect_allele, 
  # effect_allele_freq, effect, effect_se, statistic, p, n
  # (which will allow best quality control)
  #
  # minimal option 1:
  # marker, effect_allele, effect, effect_se, n   (limited quality control)
  #
  # minimal option 2:
  # marker, effect_allele, sign, p, n             (limited quality control)
  
  marker <- data.frame(rsids=unlist(Glist$rsids),unlist(Glist$cpra),
                   chr=unlist(Glist$chr), pos=unlist(Glist$position), 
                   a1=unlist(Glist$a1), a2=unlist(Glist$a2),
                   af=unlist(Glist$af))
  rownames(marker) <- marker$rsids
  
  message("Filtering markers based on information in Glist:")
  message("")
  

  #message("Filtering markers based on qc information in Glist:")
  #message("")
  rsids <-  gfilter(Glist = Glist,
                    excludeMAF=excludeMAF, 
                    excludeMISS=excludeMISS, 
                    excludeCGAT=excludeCGAT, 
                    excludeINDEL=excludeINDEL, 
                    excludeDUPS=excludeDUPS, 
                    excludeHWE=excludeHWE, 
                    excludeMHC=excludeMHC)
  marker <- marker[marker$rsids%in%rsids,]
  message("")
  
  if(!is.null(Glist$rsidsLD)) {
    rsids <- unlist(Glist$rsidsLD)
    message(paste("Number of markers in sparse LD matrices:", sum(marker$rsids%in%rsids)))
    message("")
    marker <- marker[marker$rsids%in%rsids,]
  }
  
  message("Filtering markers based on information in stat:")
  message("")
  
  if(!is.null(stat$rsids)) marker_in_stat <- marker$rsids%in%stat$rsids
  if(!is.null(stat$marker)) marker_in_stat <- marker$rsids%in%stat$marker
  message(paste("Number of markers in stat also found in bimfiles:", sum(marker_in_stat)))
  message("")
  if(sum(marker_in_stat)==0) stop("No marker ids found in bimfiles")
  
  # align marker and stat object
  marker <- marker[marker_in_stat,]
  stat <- stat[marker$rsids,]
  
  if(!is.null(stat$effect_allele)) aligned <- stat$effect_allele==marker$a1
  if(!is.null(stat$a1)) aligned <- stat$a1==marker$a1
  message(paste("Number of effect alleles aligned with first allele in bimfiles:", sum(aligned)))
  message(paste("Number of effect alleles not aligned with first allele in bimfiles:", sum(!aligned)))
  message("")
  

  if(format=="external") {
    #original
    effect <- stat[,"effect"]
    effect_allele <- stat[,"effect_allele"]
    non_effect_allele <- stat[,"non_effect_allele"]
    effect_allele_freq <- stat[,"effect_allele_freq"]
    # aligned
    stat[!aligned,"effect"] <- -effect[!aligned]
    stat[!aligned,"effect_allele"] <- non_effect_allele[!aligned]
    stat[!aligned,"non_effect_allele"] <- effect_allele[!aligned] 
    stat[!aligned,"effect_allele_freq"] <- 1-effect_allele_freq[!aligned]
    #plot(x=stat$effect_allele_freq, y=marker$af, 
    #     ylab="Allele frequency in Glist (after allele matching)", 
    #     xlab="Allele frequency in stat (after allele matching)")
    excludeMAFDIFF <- abs(marker$af-stat$effect_allele_freq) > excludeMAFDIFF
    message(paste("Number of markers excluded by large difference between MAF difference:", sum(excludeMAFDIFF)))
    message("")
    stat <- stat[!excludeMAFDIFF,]
    marker <- marker[!excludeMAFDIFF,]
    if(is.null(stat$n)) stat$n <- neff(seb=stat$effect_se,af=stat$effect_allele_freq)
    colnames(stat)[1:8] <- fm_internal
    
  }  

  if(format=="internal") {
    #original
    effect <- stat[,"b"]
    effect_allele <- stat[,"a1"]
    non_effect_allele <- stat[,"a2"]
    effect_allele_freq <- stat[,"af"]
    # aligned
    stat[!aligned,"b"] <- -effect[!aligned]
    stat[!aligned,"a1"] <- non_effect_allele[!aligned]
    stat[!aligned,"a2"] <- effect_allele[!aligned] 
    stat[!aligned,"af"] <- 1-effect_allele_freq[!aligned]
    #plot(x=stat$af, y=marker$af, 
    #     ylab="Allele frequency in Glist (after allele matching)", 
    #     xlab="Allele frequency in stat (after allele matching)")
    excludeMAFDIFF <- abs(marker$af-stat$af) > excludeMAFDIFF
    message(paste("Number of markers excluded by large difference between MAF difference:", sum(excludeMAFDIFF)))
    message("")
    stat <- stat[!excludeMAFDIFF,]
    marker <- marker[!excludeMAFDIFF,]
    if(is.null(stat$n)) stat$n <- neff(seb=stat$effect_se,af=stat$effect_allele_freq)
  }  
  
  #if(!is.null(filename)) png(file=filename)
  
  if(!is.null(stat$info)) {
    lowINFO <- stat$info < excludeINFO
    message(paste("Number of markers excluded by low INFO score:", sum(lowINFO)))
    message("")
    stat <- stat[!lowINFO,]
  }  
  return(stat)
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
  
  try(lm(stat$effect_allele_freq[aligned]~ df$af[aligned]))
  plot(stat$effect_allele_freq[aligned],df$af[aligned], ylab="AF in Glist (allele matching)",xlab="AF in stat (allele matching)")
  
  try(lm(stat$effect_allele_freq[!aligned]~ df$af[!aligned]))
  plot(stat$effect_allele_freq[!aligned],df$af[!aligned], ylab="AF in Glist (allele not matching)",xlab="AF in stat (allele not matching)")
  
  stat[!aligned,"effect_allele_freq"] <- 1 - stat[!aligned,"effect_allele_freq"]
  effect <- stat[!aligned,"b"]
  effect_allele <- stat[!aligned,"effect_allele"]
  alternative_allele <- stat[!aligned,"alternative_allele"]
  stat[!aligned,"effect_allele"] <- alternative_allele 
  stat[!aligned,"alternative_allele"] <- effect_allele 
  stat[!aligned,"b"] <- -effect 
  
  try(lm(stat$effect_allele_freq~ df$af))
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




################################################################################
# Note on how to obtain a pre-specified prior mode
################################################################################
#
# Inverse Wishart
#
# mode = S/(df+nt+1)	(mode of prior variance estimates)
# => S = mode*(df+nt+1)	(prior S to given prior mode)
# nt = number of traits
# df = prior degrees of freedom
# S = prior scale parameter 
#
# Inverse Chi-square
#
# mode = (S/(df+2))/df		(mode of prior variance estimates)
# => S = (mode*(df+2))/df	(prior S to given prior mode)
# df = prior degrees of freedom
# S = prior scale parameter 
################################################################################

#'
#' @export
#'
#'
#'

    
bmm <- function(y=NULL, X=NULL, W=NULL, GRMlist=NULL,
                vg=NULL, ve=NULL, nug=NULL, nue=NULL,
                vg_prior=NULL, ve_prior=NULL,
                updateG=TRUE, updateE=TRUE,
                nit=500, nburn=0, tol=0.001, verbose=FALSE) {
  
  n <- nrow(y)                            # number of observation
  nt <- ncol(y)                           # number of traits
  tnames <- colnames(y)
  if(is.null(tnames)) tnames <- paste0("T",1:nt)
  e <- matrix(0,nrow=n,ncol=nt)
  mu <- matrix(0,nrow=n,ncol=nt)
  psum <- matrix(0,nrow=n,ncol=nt)
  
  nset <- length(GRMlist)                   # number of sets
  setnames <- names(GRMlist)
  if(is.null(setnames)) setnames <- paste0("Set",1:nset) 
  g <- gm <- NULL
  for ( i in 1:nt) {                        # genetic values (n*nset)
    g[[i]] <- matrix(0,nrow=n,ncol=nset)         
    gm[[i]] <- matrix(0,nrow=n,ncol=nset)
    rownames(g[[i]]) <- rownames(y) 
    colnames(g[[i]]) <- setnames
    rownames(gm[[i]]) <- rownames(y) 
    colnames(gm[[i]]) <- setnames
  }
  names(gm) <- names(g) <- tnames
  
  vgm <- lapply(1:nset,function(x){matrix(0,nt,nt)})
  names(vgm) <- setnames

  vem <- matrix(0,nt,nt)
  
  if(is.null(vg_prior)) vg_prior <- lapply(1:nset,function(x){diag(0.01,nt)})
  if(is.null(ve_prior)) ve_prior <- diag(0.01,nt)
  if(is.null(nug)) nug <- rep(4,nset)
  if(is.null(nue)) nue <- 4
  
  if(!is.null(ve)) ve <- as.matrix(ve)
  if(is.null(ve)) ve <- ve_prior
  if(is.null(vg)) vg <- vg_prior
  
  ves <- matrix(0,nrow=nit,ncol=(nt*(nt+1))/2)
  
  idsG <- rownames(GRMlist[[1]])
  idsT <- rownames(y)
  idsV <- idsG[!idsG%in%idsT]
  train <- match(idsT,idsG)

  # eigen value decomposition of the GRM5 matrices
  U <- D <- vector(length=nset,mode="list")
  vgs <- vector(length=nset,mode="list")
  for ( i in 1:nset ){
    eg <- eigen(GRMlist[[i]][train,train])                    
    ev <- eg$values
    U[[i]] <- eg$vectors[,ev>tol]            # keep eigen vector if ev>tol
    D[[i]] <- eg$values[ev>tol]              # keep eigen value if ev>tol
    vgs[[i]] <- matrix(NA,nrow=nit,ncol=(nt*(nt+1))/2)
  }
  
  for ( i in 1:nit) {                     
    for ( j in 1:nset ) {                   
      rhs <- NULL
      for (t in 1:nt) {
        yadj <- y[,t]-mu[,t]-rowSums(as.matrix(g[[t]][,-j]))
        rhst <- crossprod( U[[j]],yadj ) 
        rhs <- cbind(rhs,rhst)
      }
      Vi <- solve(vg[[j]])
      a <- matrix(0,nrow=nrow(rhs),ncol=nt)
      for (k in 1:nrow(rhs)) {
        iC <- solve( diag(1,nt) + (ve%*%Vi)/D[[j]][k] )      
        ahat <- iC%*%rhs[k,]
        a[k,] <- MASS::mvrnorm(n=1,mu=ahat,Sigma=iC%*%ve) 
      }
      for (t in 1:nt) { 
        g[[t]][,j] <- U[[j]]%*%a[,t] 
      }
      if(i>nburn){
        for (t in 1:nt) { 
          gm[[t]][,j] <- gm[[t]][,j] + g[[t]][,j]  
        } 
      }
      
      # Sample variance components
      df <- nrow(a) + nug[j]
      
      # inverse chisquare
      if (nt==1) {
        scg <- sum((1/D[[j]])*a**2) + (vg_prior[j]*(nug[j]+2))/nug[j]	# => S = (mode*(df+2))/df         
        #scg <- sum((1/D[[j]])*a**2) + (vg[j]*(nug[j]+2))/nug[j]	# => S = (mode*(df+2))/df         
        vg[j] <- scg/rchisq(n=1, df=df, ncp=0)    
        vgs[[j]][i,] <- vg[j]
      }
      
      # inverse wishart
      if (nt>1) {
        S <- t(a*(1/D[[j]]))%*%a + vg_prior[[j]]*(nug[j]+nt+1)		# => S = mode*(df+nt+1)
        #S <- t(a*(1/D[[j]]))%*%a + vg[[j]]*(nug[j]+nt+1)		# => S = mode*(df+nt+1)
        vg[[j]] <- MCMCpack::riwish(df,S)
        vgs[[j]][i,] <- vg[[j]][as.vector(lower.tri(vg[[j]], diag=TRUE))]
      }
    }
    
    # Sample mu
    if(is.null(X)) {
      for (t in 1:nt) {
        yadj <- y[,t]-rowSums(as.matrix(g[[t]])) 
        rhs <- sum(yadj)
        lhs <- (n+ve[t,t]/100000)
        mu[,t] <- rnorm(1,mean=rhs/lhs,sd=1/sqrt(lhs))
      }  
    }
    if(!is.null(X)) {
      for (t in 1:nt) {
        yadj <- y[,t]-rowSums(as.matrix(g[[t]])) 
        mu[,t] <- X%*%solve(t(X)%*%X)%*%t(X)%*%yadj
      }  
    }
    
    for (t in 1:nt) {
      e[,t] <- y[,t]-mu[,t]-rowSums(g[[t]])
      p <-  (1/sqrt(2*pi))*exp(-0.5*(e[,t]**2)) 
      if(i>nburn) psum[,t] <- psum[,t] + 1/p
    }
    Se <- t(e)%*%e + ve_prior*(mean(nue)+nt+1)
    dfe <- nrow(y) + mean(nue)+nt+1
    ve <- MCMCpack::riwish(dfe,Se)
    ves[i,] <- ve[as.vector(lower.tri(ve, diag=TRUE))]
    
    if(nit>nburn) {
      for ( j in 1:nset ) {
        vgm[[j]] <- vgm[[j]] + vg[[j]]  
      }                   
      vem <- vem + ve  
    }
    if(verbose) message(paste("Iteration:",i))
    
  }
  
  for ( j in 1:nset ) {                   
    vgm[[j]] <- vgm[[j]]/(nit-nburn)
    colnames(vgm[[j]]) <- rownames(vgm[[j]]) <- tnames
    for (t in 1:nt) { 
      if(nit>nburn) gm[[t]][,j] <- gm[[t]][,j]/(nit-nburn) 
    } 
  }
  vem <- vem/(nit-nburn)  
  colnames(vem) <- rownames(vem) <- tnames

  # summary of model fit
  logCPO <- NULL
  for (t in 1:nt) {
    logCPO[t] <- sum(log((nit-nburn)*(1/psum[,t])))
  }
  fit <- list(vgs=vgs,ves=ves,vgm=vgm,vem=vem,logCPO=logCPO, gm=gm, g=g,e=e, nit=nit, nburn=nburn)
  vgm <- Reduce(`+`, fit$vgm)
  vem <- fit$vem
  vpm <- vgm + vem
  fit$vpm <- vpm
  fit$h2 <- diag(vgm/vpm)
  fit$h2set <- lapply(fit$vgm,function(x) {diag(x/vpm)})
  fit$rp <- cov2cor(vpm)
  fit$re <- cov2cor(vem)
  fit$rg <- cov2cor(vgm)
  fit$rgset <- lapply(fit$vgm,function(x) {cov2cor(x)})
  
  gset <- NULL
  for ( i in 1:nset) {                        
    gset[[i]] <- matrix(0,nrow=length(idsG),ncol=nt)         
    rownames(gset[[i]]) <- idsG 
    colnames(gset[[i]]) <- tnames
  }
  names(gset) <- setnames
  
  for ( set in 1:nset ){
    for (t1 in 1:nt) {
      gset[[set]][idsT,t1] <- gm[[t1]][,set]
      for (t2 in 1:nt) {
        gset[[set]][idsV,t1] <- gset[[set]][idsV,t1] + (GRMlist[[set]][idsV,idsT]*vgm[t1,t2])%*%gm[[t2]][,set]
      }
    }
  }
  fit$gset <- gset
  fit$g <- Reduce(`+`, fit$gset)
  
  return(fit) # return posterior samples of sigma 
}

