####################################################################################################################
#    Module 6: Bayesian models
####################################################################################################################
#'
#' Bayesian linear regression models
#'
#' @description
#' 
#' Bayesian linear regression (BLR) models:
#' 
#'   - unified mapping of genetic variants, estimation of genetic parameters 
#' (e.g. heritability) and prediction of disease risk) 
#' 
#' - handles different genetic architectures (few large, many small effects)
#' 
#' - scale to large data (e.g. sparse LD)
#' 
#' 
#' In the Bayesian multiple regression model the posterior density of the 
#' model parameters depend on the likelihood of the data given 
#' the parameters and a prior probability for the model parameters
#' 
#' The prior density of marker effects defines whether the model will 
#' induce variable selection and shrinkage or shrinkage only. 
#' Also, the choice of prior will define the extent and type of shrinkage induced.
#' Ideally the choice of prior for the marker effect should reflect the genetic 
#' architecture of the trait, and will vary (perhaps a lot) across traits.
#' 
#' 
#' The following prior distributions are provided:
#' 
#' Bayes N: Assigning a Gaussian prior to marker effects implies that the posterior means are the 
#' BLUP estimates (same as Ridge Regression).
#' 
#' Bayes L: Assigning a double-exponential or Laplace prior is the density used in 
#' the Bayesian LASSO
#' 
#' Bayes A: similar to ridge regression but t-distribution prior (rather than Gaussian) 
#' for the marker effects ; variance comes from an inverse-chi-square distribution instead of being fixed. Estimation 
#' via Gibbs sampling. 
#' 
#' Bayes C: uses a “rounded spike” (low-variance Gaussian) at origin many small 
#' effects can contribute to polygenic component, reduces the dimensionality of 
#' the model (makes Gibbs sampling feasible). 
#' 
#' Bayes R: Hierarchical Bayesian mixture model with 4 Gaussian components, with 
#' variances scaled by 0, 0.0001 , 0.001 , and 0.01 . 
#'
#'
#' @param y is a vector or matrix of phenotypes
#' @param X is a matrix of covariates
#' @param W is a matrix of centered and scaled genotypes
#' @param nit is the number of iterations
#' @param nburn is the number of burnin iterations
#' @param nthin is the thinning parameter
#' @param nit_global is the number of global iterations
#' @param nit_local is the number of local iterations
#' @param pi is the proportion of markers in each marker variance class (e.g. pi=c(0.999,0.001),used if method="ssvs")
#' @param h2 is the trait heritability
#' @param method specifies the methods used (method="bayesN","bayesA","bayesL","bayesC","bayesR")
#' @param algorithm specifies the algorithm
#' @param tol is tolerance, i.e. convergence criteria used in gbayes
#' @param Glist list of information about genotype matrix stored on disk
#' @param stat dataframe with marker summary statistics
#' @param covs is a list of summary statistics (output from internal cvs function)
#' @param fit is a list of results from gbayes
#' @param trait is an integer used for selection traits in covs object
#' @param chr is the chromosome for which to fit BLR models
#' @param rsids is a character vector of rsids
#' @param b is a vector or matrix of marginal marker effects 
#' @param seb is a vector or matrix of standard error of marginal effects
#' @param mask is a vector or matrix of TRUE/FALSE specifying if marker should be ignored 
#' @param bm is a vector or matrix of adjusted marker effects for the BLR model
#' @param LD is a list with sparse LD matrices
#' @param n is a scalar or vector of number of observations for each trait
#' @param vb is a scalar or matrix of marker (co)variances
#' @param vg is a scalar or matrix of genetic (co)variances
#' @param ve is a scalar or matrix of residual (co)variances
#' @param ssb_prior is a scalar or matrix of prior marker (co)variances
#' @param ssg_prior is a scalar or matrix of prior genetic (co)variances
#' @param sse_prior is a scalar or matrix of prior residual (co)variances
#' @param vg_prior is a scalar or matrix of prior genetic (co)variances
#' @param ve_prior is a scalar or matrix of prior residual (co)variances
#' @param nub is a scalar or vector of prior degrees of freedom for marker (co)variances
#' @param nug is a scalar or vector of prior degrees of freedom for prior genetic (co)variances
#' @param nue is a scalar or vector of prior degrees of freedom for prior residual (co)variances
#' @param updateB is a logical for updating marker (co)variances
#' @param updateG is a logical for updating genetic (co)variances
#' @param updateE is a logical for updating residual (co)variances
#' @param updatePi is a logical for updating pi
#' @param adjustE is a logical for adjusting residual variance
#' @param models is a list structure with models evaluated in bayesC
#' @param formatLD is a character specifying LD format (formatLD="dense" is default)
#' @param verbose is a logical; if TRUE it prints more details during iteration
#' @param scaleY is a logical; if TRUE y is centered and scaled 
#' @param msize number of markers used in compuation of sparseld
#' @param lambda is a vector or matrix of lambda values 
#' @param GRMlist is a list providing information about GRM matrix stored in binary files on disk

#' @return Returns a list structure including
#' \item{b}{vector or matrix (mxt) of posterior means for marker effects}
#' \item{d}{vector or matrix (mxt) of posterior means for marker inclusion probabilities}
#' \item{vb}{scalar or vector (t) of posterior means for marker variances}
#' \item{vg}{scalar or vector (t) of posterior means for genomic variances}
#' \item{ve}{scalar or vector (t) of posterior means for residual variances}
#' \item{rb}{matrix (txt) of posterior means for marker correlations}
#' \item{rg}{matrix (txt) of posterior means for genomic correlations}
#' \item{re}{matrix (txt) of posterior means for residual correlations}
#' \item{pi}{vector (1xnmodels) of posterior probabilities for models}
#' \item{h2}{vector (1xt) of posterior means for model probability}
#' \item{param}{ a list current parameters (same information as item listed above) used for restart of the analysis}
#' \item{stat}{matrix (mxt) of marker information and effects used for genomic risk scoring}


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
#' fitM <- gbayes(y=y, W=W, method="bayesN")
#' fitA <- gbayes(y=y, W=W, method="bayesA")
#' fitL <- gbayes(y=y, W=W, method="bayesL")
#' fitC <- gbayes(y=y, W=W, method="bayesC")
#'


#'
#' @export
#'

gbayes <- function(y=NULL, X=NULL, W=NULL, stat=NULL, covs=NULL, trait=NULL, fit=NULL, Glist=NULL, 
                   chr=NULL, rsids=NULL, b=NULL, bm=NULL, seb=NULL, LD=NULL, n=NULL,formatLD="dense",
                   vg=NULL, vb=NULL, ve=NULL, ssg_prior=NULL, ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=TRUE,
                   h2=NULL, pi=0.001, updateB=TRUE, updateG=TRUE, updateE=TRUE, updatePi=TRUE, adjustE=TRUE, models=NULL,
                   nug=4, nub=4, nue=4, verbose=FALSE,msize=100, mask=NULL,
                   GRMlist=NULL, ve_prior=NULL, vg_prior=NULL,tol=0.001,
                   nit=100, nburn=0, nthin=1, nit_local=NULL,nit_global=NULL,
                   method="mixed", algorithm="mcmc") {

    
  # mask
  mask.rsids <- NULL
  if(!is.null(mask)) mask.rsids <- unique(as.vector(apply(mask,2,function(x){as.vector(rownames(mask)[x])})))

  
  # Check methods
  methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
  method <- match(method, methods) - 1
  if( !sum(method%in%c(0:5))== 1 ) stop("Method specified not valid") 

  algorithms <- c("mcmc","em-mcmc")
  algorithm <- match(algorithm, algorithms)
  if(is.na(algorithm)) stop("algorithm argument specified not valid")
  
  if(method==0) {
    # BLUP and we do not estimate parameters
    updateB=FALSE;
    updateE=FALSE;
  }
  
  # Determine number of traits
  nt <- 1
  if(!is.null(y)) {
    if(is.list(y)) nt <- length(y)
    if(is.matrix(y)) nt <- ncol(y)
  }
  if(!is.null(stat)) {
    if(!is.data.frame(stat) && is.list(stat)) nt <- ncol(stat$b)
    if( any(sapply(stat[,-c(1:5)],function(x){any(!is.finite(x))}))) stop("Some elements in stat not finite")
    if( any(sapply(stat[,-c(1:5)],function(x){any(is.na(x))}))) stop("Some elements in stat are NA's")
  }
  
  # Define type of analysis
  if(!is.null(GRMlist)) analysis <- "mtmc-mixed"
  
  if(nt==1 && !is.null(y) && !is.null(W) && formatLD=="dense") 
    analysis <- "st-blr-individual-level-default"
  
  if(nt==1 && !is.null(y) && !is.null(W) && formatLD=="sparse") 
    analysis <- "st-blr-individual-level-sbayes"
  
  if(nt==1 && !is.null(y) && formatLD=="sparse") 
    analysis <- "st-blr-individual-level-sparse-ld"
  
  if( nt==1 && !is.null(y) &&  formatLD=="dense") 
    analysis <- "st-blr-individual-level-dense-ld"
  
  if( nt==1 && is.null(y) && !is.null(stat) && !is.null(Glist)) 
    analysis <- "st-blr-sumstat-sparse-ld"
  
  if(nt>1 && !is.null(y) && !is.null(W))
    analysis <- "mt-blr-individual-level-dense-ld"

  if(nt>1 && !is.null(y) && formatLD=="sparse") 
    analysis <- "mt-blr-individual-level-sparse-ld"
  
  if( nt>1 && !is.null(stat) && !is.null(Glist)) 
    analysis <- "mt-blr-sumstat-sparse-ld"
  
  message(paste("Type of analysis performed:",analysis))  
  
  
  ##############################################################################
  # Different work flows
  ##############################################################################
  
  # Single and multiple trait BLR based on GRMs    
  if(!is.null(GRMlist)) {
    fit <- bmm(y=y, X=X, W=W, GRMlist=GRMlist,
               vg=vg, ve=ve, nug=nug, nue=nue,
               vg_prior=vg_prior, ve_prior=ve_prior,
               updateG=updateB, updateE=updateE,
               nit=nit, nburn=nburn, tol=tol, verbose=verbose) 
  } 
  
  
  # Single trait BLR using y and W   
  if(nt==1 && !is.null(y) && !is.null(W) && formatLD=="dense") {
    
    fit <- bayes(y=y, X=X, W=W, b=b, scaled=TRUE,
                 h2=h2, pi=pi, lambda=lambda, 
                 vg=vg, vb=vb, ve=ve, 
                 nub=nub, nug=nug, nue=nue, 
                 ssb_prior=ssb_prior, ssg_prior=ssg_prior, sse_prior=sse_prior, 
                 updateB=updateB, updateG=updateG, updateE=updateE, updatePi=updatePi, 
                 nit=nit, nburn=nburn, nthin=nthin, method=method)  
  }

  # Multiple trait BLR using y and W
  if(nt>1 && !is.null(y) && !is.null(W)) {
    fit <- mtbayes(y=y, X=X, W=W, b=b, bm=bm, seb=seb, LD=LD, n=n,
                   vg=vg, vb=vb, ve=ve, 
                   ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
                   h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
                   nub=nub, nue=nue, nit=nit, method=method, formatLD=formatLD, algorithm=algorithm, verbose=verbose) 
  }
  
  


  # Single trait BLR using summary statistics and sparse LD provided in Glist
  if(analysis=="st-blr-sumstat-sparse-ld") {
    
    # single trait summary statistics
    if(is.data.frame(stat)) {
      
      nt <- 1
      rsidsLD <- unlist(Glist$rsidsLD)
      m <- length(rsidsLD)
      b <- wy <- ww <- matrix(0,nrow=length(rsidsLD),ncol=nt)
      mask <- matrix(TRUE,nrow=length(rsidsLD),ncol=nt)
      rownames(b) <- rownames(wy) <- rownames(ww) <- rownames(mask) <- rsidsLD
      trait_names <- "bm"     
      
      stat <- stat[rownames(stat)%in%rsidsLD,]
      if(is.null(stat$ww)) stat$ww <- 1/(stat$seb^2 + stat$b^2/stat$n)
      if(is.null(stat$wy)) stat$wy <- stat$b*stat$ww
      if(!is.null(stat$n)) n <- as.integer(median(stat$n))
      ww[rownames(stat),1] <-  stat$ww
      wy[rownames(stat),1] <- stat$wy
      mask[rownames(stat),1] <- FALSE
      if(!is.null(mask.rsids)) mask[rownames(mask)%in%mask.rsids,1] <- TRUE
      if(any(is.na(wy))) stop("Missing values in wy")
      if(any(is.na(ww))) stop("Missing values in ww")
      
      b2 <- stat$b^2
      seb2 <- stat$seb^2
      #yy <- (b2 + (n-2)*seb2)*ww[rownames(stat),]
      yy <- (b2 + (n-2)*seb2)*stat$ww
      yy <- median(yy)
    }
    
    # multiple trait summary statistics
    if( !is.data.frame(stat) && is.list(stat)) {
      
      nt <- ncol(stat$b)
      trait_names <- colnames(stat$b)
      if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
      rsidsLD <- unlist(Glist$rsidsLD)
      m <- length(rsidsLD)
      b <- wy <- ww <- matrix(0,nrow=length(rsidsLD),ncol=nt)
      mask <- matrix(TRUE,nrow=length(rsidsLD),ncol=nt)
      rownames(b) <- rownames(wy) <- rownames(ww) <- rownames(mask) <- rsidsLD
      colnames(b) <- colnames(wy) <- colnames(ww) <- colnames(mask) <- trait_names
      
      rws <- rownames(stat$b)%in%rsidsLD
      if(is.null(stat$ww)) stat$ww <- 1/(stat$seb^2 + stat$b^2/stat$n)
      if(is.null(stat$wy)) stat$wy <- stat$b*stat$ww
      if(!is.null(stat$n)) n <- as.integer(apply(stat$n[rws,],2,median))
      ww[rownames(stat$ww[rws,]),] <- stat$ww[rws,]
      wy[rownames(stat$wy[rws,]),] <- stat$wy[rws,]
      mask[rownames(stat$wy[rws,]),] <- FALSE
      if(!is.null(mask.rsids)) mask[rownames(mask)%in%mask.rsids,] <- TRUE
      
      if(any(is.na(wy))) stop("Missing values in wy")
      if(any(is.na(ww))) stop("Missing values in ww")
      
      b2 <- (stat$b[rws,])^2
      seb2 <- (stat$seb[rws,])^2
      for (i in 1:nt) {
        seb2[,i] <- (n[i]-2)*seb2[,i]
      }
      yy <- (b2 + seb2)*stat$ww[rws,]
      yy <- apply(yy,2,median)
      
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
      bmchr <- dmchr <- NULL
      for (trait in 1:nt) {
        if(verbose) print( paste("Fit",methods[method+1], "on chromosome:",chr,"for trait",trait))
        #LDvalues <- LD[[chr]]$values
        fit[[chr]] <- sbayes_sparse(yy=yy[trait], 
                                    wy=wy[rsidsLD,trait],
                                    ww=ww[rsidsLD,trait],
                                    b=b[rsidsLD,trait], 
                                    LDvalues=LD[[chr]]$values, 
                                    LDindices=LD[[chr]]$indices,
                                    mask=mask[rsidsLD,trait],
                                    method=method, 
                                    nit=nit, 
                                    nburn=nburn,
                                    nthin=nthin,
                                    n=n[trait],
                                    m=m,
                                    pi=pi,
                                    nue=nue, 
                                    nub=nub, 
                                    h2=h2[trait], 
                                    lambda=lambda, 
                                    vb=vb, 
                                    ve=ve,
                                    ssb_prior=ssb_prior, 
                                    sse_prior=sse_prior,
                                    updateB=updateB, 
                                    updateE=updateE, 
                                    updatePi=updatePi,
                                    updateG=updateG,
                                    adjustE=adjustE,
                                    algorithm=algorithm)
        bmchr <- cbind(bmchr, fit[[chr]]$bm)
        dmchr <- cbind(dmchr, fit[[chr]]$dm)
      }
      colnames(bmchr) <- trait_names
      res[[chr]] <- data.frame(rsids=rsidsLD,chr=rep(chr,length(rsidsLD)),
                               pos=Glist$pos[[chr]][clsLD], ea=Glist$a1[[chr]][clsLD],
                               nea=Glist$a2[[chr]][clsLD], eaf=Glist$af[[chr]][clsLD],
                               bm=bmchr, dm=dmchr, stringsAsFactors = FALSE)
      rownames(res[[chr]]) <- rsidsLD
      LD[[chr]]$values <- NULL
      LD[[chr]]$indices <- NULL
    }
    res <- do.call(rbind, res)
    rownames(res) <- res$rsids
    fit$stat <- res
    fit$stat$vm <- 2*(1-fit$stat$eaf)*fit$stat$eaf*fit$stat$bm^2
    fit$method <- methods[method+1]
    
    fit$ves <- lapply(fit[chromosomes],function(x){x$ves})
    fit$vgs <- lapply(fit[chromosomes],function(x){x$vgs})
    fit$vbs <- lapply(fit[chromosomes],function(x){x$vbs})
    fit$pis <- lapply(fit[chromosomes],function(x){x$pis})
    fit$pim <- lapply(fit[chromosomes],function(x){x$pim})
    fit$param <- lapply(fit[chromosomes],function(x){x$param})
    
    fit$mask <- mask
    zve <- sapply(fit$ves[chromosomes],function(x){coda::geweke.diag(x[nburn:length(x)])$z})
    zvg <- sapply(fit$vgs[chromosomes],function(x){coda::geweke.diag (x[nburn:length(x)])$z})
    zvb <- sapply(fit$vbs[chromosomes],function(x){coda::geweke.diag(x[nburn:length(x)])$z})
    zpi <- sapply(fit$pis[chromosomes],function(x){coda::geweke.diag(x[nburn:length(x)])$z})
    ve <- sapply(fit$ves[chromosomes],function(x){mean(x[nburn:length(x)])})
    vg <- sapply(fit$vgs[chromosomes],function(x){mean(x[nburn:length(x)])})
    vb <- sapply(fit$vbs[chromosomes],function(x){mean(x[nburn:length(x)])})
    pi <- sapply(fit$pim[chromosomes],function(x){1-x[1]})
    fit$conv <- data.frame(zve=zve,zvg=zvg, zvb=zvb, zpi=zpi)  
    fit$post <- data.frame(ve=ve,vg=vg, vb=vb,pi=pi)  
    fit$ve <- mean(ve)
    fit$vg <- sum(vg)

  }
  

  # Multi trait BLR using summary statistics and sparse LD provided in Glist
  #  if( nt>1 && is.null(y) && !is.null(stat) && !is.null(Glist)) {
  
  if(analysis=="mt-blr-sumstat-sparse-ld") {
    
    # multiple trait summary statistics
    
    nt <- ncol(stat$b)
    trait_names <- colnames(stat$b)
    if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
    rsidsLD <- unlist(Glist$rsidsLD)
    b <- wy <- ww <- matrix(0,nrow=length(rsidsLD),ncol=nt)
    rownames(b) <- rownames(wy) <- rownames(ww) <- rsidsLD
    colnames(b) <- colnames(wy) <- colnames(ww) <- trait_names
    
    rws <- rownames(stat$b)%in%rsidsLD
    #if(is.null(stat$ww)) stat$ww <- 1/(stat$seb^2 + (stat$b^2)/stat$n)
    #if(is.null(stat$wy)) stat$wy <- stat$b*stat$n
    if(!is.null(stat$n)) n <- as.integer(colMeans(stat$n[rws,]))
    if(!is.null(stat$ww)) ww[rownames(stat$ww[rws,]),] <- stat$ww[rws,]
    if(is.null(stat$ww)) ww[rownames(stat$n[rws,]),] <- stat$n[rws,]
    if(!is.null(stat$wy)) wy[rownames(stat$wy[rws,]),] <- stat$wy[rws,]
    if(is.null(stat$wy)) wy[rownames(stat$b[rws,]),] <- stat$b[rws,]*stat$n[rws,]
    if(any(is.na(wy))) stop("Missing values in wy")
    if(any(is.na(ww))) stop("Missing values in ww")
    
    b2 <- (stat$b[rws,])^2
    seb2 <- (stat$seb[rws,])^2
    for (i in 1:nt) {
      seb2[,i] <- (n[i]-2)*seb2[,i]
    }
    yy <- (b2 + seb2)*stat$ww[rws,]
    yy <- colMeans(yy)
    
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
                                     ww=ww[rsidsLD,],
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
                                     models=models,
                                     method=method,
                                     verbose=verbose)
      res[[chr]] <- data.frame(rsids=rsidsLD,chr=rep(chr,length(rsidsLD)),
                               pos=Glist$pos[[chr]][clsLD], ea=Glist$a1[[chr]][clsLD],
                               nea=Glist$a2[[chr]][clsLD], eaf=Glist$af[[chr]][clsLD],
                               bm=fit[[chr]]$bm, dm=fit[[chr]]$dm,  stringsAsFactors = FALSE)
      rownames(res[[chr]]) <- rsidsLD
      LD[[chr]]$values <- NULL
      LD[[chr]]$indices <- NULL
    }
    res <- do.call(rbind, res)
    rownames(res) <- res$rsids
    fit$stat <- res
    fit$method <- methods[method+1]
    fit$stat$vm <- 2*(1-fit$stat$eaf)*fit$stat$eaf*fit$stat$bm^2
    
  }
  return(fit)
  
}

##############################################################################
# Core functions used in work flows
##############################################################################

# Single trait BLR based on individual level data 
bayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, scaled=TRUE,
                  h2=0.5, pi=0.001, lambda=NULL, 
                  vb=NULL, vg=NULL, ve=NULL, 
                  nub=4, nug=4, nue=4, 
                  ssb_prior=NULL, ssg_prior=NULL, sse_prior=NULL, 
                  updateB=NULL, updateG=NULL, updateE=NULL, updatePi=NULL, 
                  nit=500, nburn=100, nthin=1, method=NULL) {
  ids <- NULL
  if(is.matrix(y)) ids <- rownames(y)
  if(is.vector(y)) ids <- names(y)
  if(scaled) y <- as.vector(scale(y)) 
  n <- nrow(W)
  m <- ncol(W)
  
  if(!is.null(ids) & !is.null(rownames(W))) {
    if(any(is.na(match(ids,rownames(W))))) stop("Names/rownames for y does match rownames for W")
  }
  vy <- var(y)
  if(is.null(pi)) pi <- 0.001
  if(is.null(h2)) h2 <- 0.5
  if(is.null(vg)) vg <- vy*h2
  if(is.null(ve)) ve <- vy*(1-h2)
  if(method<4 && is.null(vb)) vb <- vg/m
  if(method>=4 && is.null(vb)) vb <- vg/(m*pi)
  if(is.null(lambda)) lambda <- rep(ve/vb,m)
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/(m*pi))
  if(is.null(sse_prior)) sse_prior <- ((nue-2.0)/nue)*ve
  if(is.null(ssg_prior)) ssg_prior=((nug-2.0)/nug)*vg;
  if(is.null(b)) b <- rep(0,m)

  pi <- c(1-pi,pi)
  gamma <- c(0,1.0)
  if(method==5) pi <- c(0.95,0.02,0.02,0.01)
  if(method==5) gamma <- c(0,0.01,0.1,1.0)
  
  seed <- sample.int(.Machine$integer.max, 1)
  
  fit <- .Call("_qgg_bayes",
               y=y, 
               W=as.list(as.data.frame(W)), 
               b=b,
               lambda = lambda,
               pi = pi,
               gamma = gamma,
               vb = vb,
               vg = vg,
               ve = ve,
               ssb_prior=ssb_prior,
               ssg_prior=ssg_prior,
               sse_prior=sse_prior,
               nub=nub,
               nug=nug,
               nue=nue,
               updateB = updateB,
               updateG = updateG,
               updateE = updateE,
               updatePi = updatePi,
               nit=nit,
               nburn=nburn,
               nthin=nthin,
               method=as.integer(method),
               seed=seed) 
  ids <- rownames(W)
  names(fit[[1]]) <- names(fit[[2]]) <- names(fit[[11]]) <- colnames(W)
  names(fit[[9]]) <- ids
  names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pis","pim","g","b","d","param")
  
  return(fit)
}



# Single trait BLR using summary statistics and sparse LD provided in Glist 
sbayes_sparse <- function(yy=NULL, wy=NULL, ww=NULL, 
                          LDvalues=NULL,LDindices=NULL, b=NULL,  
                          n=NULL, m=NULL, 
                          h2=NULL, pi=NULL, lambda=NULL, mask=NULL,
                          vb=NULL, vg=NULL, ve=NULL, 
                          nub=4, nug=4, nue=4, 
                          ssb_prior=NULL, ssg_prior=NULL, sse_prior=NULL, 
                          updateB=NULL, updateG=NULL, updateE=NULL, updatePi=NULL, adjustE=NULL, 
                          nit=NULL, nburn=NULL, nthin=1, 
                          method=NULL, algorithm=NULL) {

  if(is.null(m)) m <- length(LDvalues)
  vy <- yy/(n-1)
  if(is.null(pi)) pi <- 0.001
  if(is.null(h2)) h2 <- 0.5
  if(is.null(ve)) ve <- vy*(1-h2)
  if(is.null(vg)) vg <- vy*h2
  if(method<4 && is.null(vb)) vb <- vg/m
  if(method>=4 && is.null(vb)) vb <- vg/(m*pi)
  if(is.null(lambda)) lambda <- rep(ve/vb,m)
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/(m*pi))
  if(is.null(ssg_prior)) ssg_prior <- ((nug-2.0)/nug)*vg
  if(is.null(sse_prior)) sse_prior <- ((nue-2.0)/nue)*ve
  if(is.null(b)) b <- rep(0,m)
  
  pi <- c(1-pi,pi)
  gamma <- c(0,1.0)
  if(method==5) pi <- c(0.95,0.02,0.02,0.01)
  if(method==5) gamma <- c(0,0.01,0.1,1.0)

  seed <- sample.int(.Machine$integer.max, 1)
  
  fit <- .Call("_qgg_sbayes",
               yy=yy,
               wy=wy, 
               ww=ww, 
               LDvalues=LDvalues, 
               LDindices=LDindices, 
               b = b,
               lambda = lambda,
               mask=mask,
               pi = pi,
               gamma = gamma,
               vb = vb,
               vg = vg,
               ve = ve,
               ssb_prior=ssb_prior,
               ssg_prior=ssg_prior,
               sse_prior=sse_prior,
               nub=nub,
               nug=nug,
               nue=nue,
               updateB = updateB,
               updateG = updateG,
               updateE = updateE,
               updatePi = updatePi,
               adjustE = adjustE,
               n=n,
               nit=nit,
               nburn=nburn,
               nthin=nthin,
               method=as.integer(method),
               algo=as.integer(algorithm),
               seed=seed)
  
  names(fit[[1]]) <- names(LDvalues)
  names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pis","pim","r","b","param")
  return(fit)
}



# Single trait BLR using summary statistics and sparse LD provided in Glist
blr <- function(yy=NULL, Xy=NULL, XX=NULL, n=NULL,
                     mask=NULL, lambda=NULL,
                     vg=NULL, vb=NULL, ve=NULL, h2=NULL, pi=NULL,
                     ssb_prior=NULL, ssg_prior=NULL, sse_prior=NULL, nub=4, nug=4, nue=4,
                     updateB=TRUE, updateG=TRUE, updateE=TRUE, updatePi=TRUE, 
                     adjustE=TRUE, models=NULL,
                     nit=500, nburn=100, nthin=1, method="bayesC", algorithm="mcmc", verbose=FALSE) {
  
  # Check methods and parameter settings
  methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
  method <- match(method, methods) - 1
  if( !sum(method%in%c(0:5))== 1 ) stop("method argument specified not valid")
  algorithms <- c("mcmc","em-mcmc")
  algorithm <- match(algorithm, algorithms)
  if(is.na(algorithm)) stop("algorithm argument specified not valid")
  
  if( is.null(n) ) stop("Please provide n")
  if( is.null(yy) ) stop("Please provide yy")
  if( is.null(Xy) ) stop("Please provide Xy matrix")
  if( is.null(XX) ) stop("Please provide XX matrix")
  
  xx <- diag(XX)  
  m <- length(xx)
  
  XX <- cov2cor(XX)  
  XXvalues <- as.list(as.data.frame(XX)) 
  XXindices <- lapply(1:ncol(XX),function(x) { (1:ncol(XX))-1 } )
  
  b <- rep(0, m)
  mask <- rep(FALSE, m)

  # Prepare starting parameters
  vy <- yy/(n-1)
  if(is.null(pi)) pi <- 0.001
  if(is.null(h2)) h2 <- 0.5
  if(is.null(ve)) ve <- vy*(1-h2)
  if(is.null(vg)) vg <- vy*h2
  if(method<4 && is.null(vb)) vb <- vg/m
  if(method>=4 && is.null(vb)) vb <- vg/(m*pi)
  if(is.null(lambda)) lambda <- rep(ve/vb,m)
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/(m*pi))
  if(is.null(ssg_prior)) ssg_prior <- ((nug-2.0)/nug)*vg
  if(is.null(sse_prior)) sse_prior <- ((nue-2.0)/nue)*ve
  if(is.null(b)) b <- rep(0,m)
  
  pi <- c(1-pi,pi)
  gamma <- c(0,1.0)
  if(method==5) pi <- c(0.95,0.02,0.02,0.01)
  if(method==5) gamma <- c(0,0.01,0.1,1.0)
  
  seed <- sample.int(.Machine$integer.max, 1)
  
  fit <- .Call("_qgg_sbayes",
               yy=yy,
               wy=Xy, 
               ww=xx, 
               LDvalues=XXvalues, 
               LDindices=XXindices, 
               b = b,
               lambda = lambda,
               mask=mask,
               pi = pi,
               gamma = gamma,
               vb = vb,
               vg = vg,
               ve = ve,
               ssb_prior=ssb_prior,
               ssg_prior=ssg_prior,
               sse_prior=sse_prior,
               nub=nub,
               nug=nug,
               nue=nue,
               updateB = updateB,
               updateG = updateG,
               updateE = updateE,
               updatePi = updatePi,
               adjustE = adjustE,
               n=n,
               nit=nit,
               nburn=nburn,
               nthin=nthin,
               method=as.integer(method),
               algo=as.integer(algorithm),
               seed=seed)
  names(fit[[1]]) <- names(XXvalues)
  names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pis","pim","r","b","param")
  return(fit)
}


################################################################################
# Single trait fine-mapping BLR using summary statistics and sparse LD provided in Glist 
# gmap full version

#' Finemapping using Bayesian Linear Regression Models
#'
#' This function implements Bayesian linear regression models to provide unified mapping of 
#' genetic variants, estimate genetic parameters (e.g. heritability), and predict disease risk.
#' It is designed to handle various genetic architectures and scale efficiently with large datasets.
#'
#' @description
#' In the Bayesian multiple regression model, the posterior density of the model parameters depends
#' on the likelihood of the data given the parameters and a prior probability for the model parameters.
#' The choice of the prior for marker effects can influence the type and extent of shrinkage induced in the model.
#'
#' @param Glist A list containing information on genotypic data, including SNPs, chromosomes, positions, and optionally, LD matrices.
#' @param stat A data frame or list of summary statistics including effect sizes, standard errors, sample sizes, etc.
#' @param sets Optional list specifying sets of SNPs for mapping.
#' @param models Optional list of predefined models for Bayesian regression.
#' @param rsids Vector of SNP identifiers.
#' @param ids Vector of sample identifiers.
#' @param vb Initial value for the marker effect variance (default: NULL).
#' @param vg Initial value for the genetic variance (default: NULL).
#' @param ve Initial value for the residual variance (default: NULL).
#' @param vy Initial value for the phenotypic variance (default: NULL).
#' @param pi Vector of initial values for pi parameters in the model (default of pi=c(0.999,0.001) for bayesC and pi=c(0.994,0.003,0.002,0.001).
#' @param gamma Vector of initial values for gamma parameters in the model (default of gamma=c(0,1) for bayesC and gamma=c(0,0.01,0.1,1).
#' @param h2 Heritability estimate (default: 0.5).
#' @param mc Number of potentiel genome-wide causal markers for the trait analysed - only used for specification of ssb_prior (default: 5000).
#' @param lambda Vector of initial values for penalty parameters in the model.
#' @param mask Logical matrix indicating SNPs to exclude from analysis.
#' @param nub,nug,nue Degrees of freedom parameters for the priors of marker, genetic, and residual variances, respectively.
#' @param ssb_prior,ssg_prior,sse_prior Priors for the marker, genetic, and residual variances.
#' @param vb_prior,vg_prior,ve_prior Additional priors for marker, genetic, and residual variances (default: NULL).
#' @param updateB,updateG,updateE,updatePi Logical values specifying whether to update marker effects, genetic variance, residual variance, and inclusion probabilities, respectively.
#' @param updateMH Logical values specifying whether to update marker effects using a Metropolis-Hasting algorithm.
#' @param formatLD Format of LD matrix ("dense" by default).
#' @param checkLD Logical, whether to check the LD matrix for inconsistencies (default: FALSE).
#' @param shrinkLD,shrinkCor Logical, whether to apply shrinkage to the LD or correlation matrices (default: FALSE).
#' @param pruneLD Logical, whether to prune LD matrix (default: FALSE).
#' @param checkConvergence Logical, whether to check for convergence of the Gibbs sampler (default: FALSE).
#' @param critVe,critVg,critVb,critPi,critB,critB1,critB2 Convergence criteria for residual, genetic, and marker variances, inclusion probabilities, and marker effects.
#' @param eigen_threshold Threshold for eigen value decomposition (default: eigen_threshold=0.995, other examples: eigen_threshold=c(0.995,0.99) )
#' @param cs_threshold,cs_r2 PIP and r2 thresholds credible set construction (default: cs_threshold=0.9, cs_r2=0.5)
#' @param mh_p,mh_r2 Probability and r2 thresholds used in MH step (default: mh_p=0.05, mh_r2=0.95)
#' @param verbose Logical, whether to print detailed output for debugging (default: FALSE).
#' @param eigen_threshold Threshold for eigenvalues in eigen decomposition (default: 0.995).
#' @param nit Number of iterations in the MCMC sampler (default: 5000).
#' @param nburn Number of burn-in iterations (default: 500).
#' @param nthin Thinning interval for MCMC (default: 5).
#' @param nrun Number of parallel runs MCMC (default: 1).
#' @param method The regression method to use, options include "blup", "bayesN", "bayesA", "bayesL", "bayesC", "bayesR".
#' @param algorithm Algorithm for MCMC sampling, options include "mcmc", "em-mcmc", "mcmc-eigen".
#' @param output Level of output, options include "summary", "full".
#' @param seed Random seed for reproducibility (default: 10).
#'
#' @return Returns a list structure including the following components:
#'
#' @author Peter Sørensen
#'
#' @export
#'
gmap <- function(Glist=NULL, stat=NULL, sets=NULL, models=NULL,
                 rsids=NULL, ids=NULL, mask=NULL, lambda=NULL,  
                 vb=NULL, vg=NULL, ve=NULL, vy=NULL,
                 pi=NULL, gamma=NULL, mc=5000, h2=0.5,   
                 nub=4, nug=4, nue=4, 
                 ssb_prior=NULL, ssg_prior=NULL, sse_prior=NULL,
                 vb_prior=NULL, vg_prior=NULL, ve_prior=NULL,
                 updateB=TRUE, updateG=TRUE, updateE=TRUE, updatePi=TRUE, updateMH=FALSE,
                 formatLD="dense", checkLD=FALSE, shrinkLD=FALSE, shrinkCor=FALSE, pruneLD=FALSE, 
                 checkConvergence=FALSE, critVe=3, critVg=3, critVb=3, critPi=3, 
                 critB=3, critB1=0.5, critB2=3, 
                 verbose=FALSE, eigen_threshold=0.995, 
                 cs_threshold=0.9, cs_r2=0.5, cs_method="CS2",
                 mh_p=0.05, mh_r2=0.95,
                 nit=1000, nburn=100, nthin=1, nrun=1, output="summary",
                 method="bayesR", algorithm="mcmc", seed=10) {
  
  
  # Check methods and parameter settings
  methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
  method <- match(method, methods) - 1
  if( !sum(method%in%c(0:5))== 1 ) stop("method argument specified not valid")
  algorithms <- c("mcmc","em-mcmc", "mcmc-eigen")
  algorithm <- match(algorithm, algorithms)
  if(is.na(algorithm)) stop("algorithm argument specified not valid")
  
  # Check input data
  if(!is.data.frame(stat)) stop("Stat is not a data frame")
  if( any(sapply(stat[,-c(1:5)],function(x){any(!is.finite(x))}))) stop("Some elements in stat not finite")
  if( any(sapply(stat[,-c(1:5)],function(x){any(is.na(x))}))) stop("Some elements in stat NA")
  if(is.null(sets)) stop("Please provide sets")
  
  m <- sum(stat$rsids%in%unlist(Glist$rsids))
  
  yy <- median((stat$b^2 + (stat$n-2)*stat$seb^2)*stat$n)
  n <- median(stat$n)
  if(is.null(stat[["ww"]])) stat$ww <- (yy/n)/(stat$seb^2 + stat$b^2/stat$n)
  if(is.null(stat[["wy"]])) stat$wy <- stat$b*stat$ww
  
  
  # Prepare starting parameters
  
  if(method==4 && is.null(pi)) pi <- c(1-0.001,0.001)
  if(method==5 && is.null(pi)) pi <- c(0.990,0.006,0.003,0.001)
  if(method==4 && is.null(gamma)) gamma <- c(0,1.0)
  if(method==5 && is.null(gamma)) gamma <- c(0,0.01,0.1,1.0)
  
  #vy <- median(2*stat$eaf*(1-stat$eaf)*(stat$n*stat$seb^2 + stat$b^2))
  vy <- yy/(n-1)
  if(is.null(ve)) ve <- vy*(1-h2)
  if(is.null(vg)) vg <- vy*h2
  if(method>=4 && is.null(vb)) vb <- vg/(mc*sum(pi*gamma))
  if(method>=4 && is.null(ssb_prior))  ssb_prior <- ((nub-2.0)/nub)*(vg/(mc*sum(pi*gamma)))
  ssg_prior <-  ((nug-2.0)/nug)*vg
  sse_prior <- ((nue-2.0)/nue)*ve
  
  
  sets <- mapSets(sets=sets, rsids=stat$rsids, index=FALSE)
  if(any(sapply(sets,function(x){any(is.na(x))}))) stop("NAs in sets detected - please remove these")
  
  chr <- as.numeric(unlist(Glist$chr))
  chrSets <- sapply(mapSets(sets = sets, Glist = Glist, index = TRUE), function(x) unique(chr[x]))
  #if (length(Glist$bedfiles) == 1) chrSets <- setNames(rep(1, length(chrSets)), names(chrSets))
  if (length(Glist$bedfiles) == 1) chrSets <- setNames(rep(1, length(sets)), names(sets))
  lsets <- sapply(chrSets,length)
  sets <- sets[lsets==1]
  if(any(lsets>1)) stop(paste("Following marker sets mapped to multiple chromosome:",paste(which(lsets>1),collapse=",")))
  if(any(lsets==0)) stop(paste("Following marker sets not mapped to any chromosome:",paste(which(lsets==0),collapse=",")))
  
  # Prepare output
  bm <- dm <- vector(mode="list",length=length(sets))
  ves <- vgs <- vbs <- pis <- conv <- vector(mode="list",length=length(sets))
  bs <- ds <- prob <- vector(mode="list",length=length(sets))
  pim <- vector(mode="list",length=length(sets))
  logcpo <- rep(0,length(sets))
  fdr <- csets <- vector(mode="list",length=length(sets))
  names(bm) <- names(dm) <- names(pim) <- names(sets)     
  names(ves) <- names(vgs) <- names(pis) <- names(vbs) <- names(conv) <- names(sets)     
  names(bs)  <- names(ds) <- names(prob) <- names(sets)
  names(logcpo) <- names(fdr) <- names(csets) <- names(sets)
  attempts <- rep(1, length=length(sets))
  
  if(is.null(ids)) ids <- Glist$idsLD
  if(is.null(ids)) ids <- Glist$ids
  
  # Loop over sets
  for (i in 1:length(sets)) {
    
    chr <- chrSets[[i]]
    rsids <- sets[[i]]
    rws <- match(rsids,stat$rsids)
    message(paste("Processing region:",i))
    
    pos <- getPos(Glist=Glist, chr=chr, rsids=rsids)
    message(paste("Region size in Mb:",round((max(pos)-min(pos))/1000000,2)))
    if(!is.null(Glist$map)) map <- getMap(Glist=Glist, chr=chr, rsids=rsids)
    if(!is.null(Glist$map)) message(paste("Region size in cM:",round(max(map)-min(map),2)))
    
    W <- getG(Glist=Glist, chr=chr, rsids=rsids, ids=ids, scale=TRUE)
    B <- crossprod(scale(W))/(nrow(W)-1)
    
    if(shrinkLD) B <- corpcor::cor.shrink(W)
    
    if(algorithm==3) eig <- eigen(B, symmetric=TRUE)
    
    for (j in 1:length(eigen_threshold)) {
      
      if(algorithm==1) {
        LD <- NULL
        LDvalues <- as.list(as.data.frame(B))
        LDindices <- lapply(1:ncol(B),function(x) { (1:nrow(B))-1 } )
        names(LDvalues) <- names(LDindices) <- rsids
      }
      if(algorithm==3) {
        ww <- stat[rws,"n"]
        scaleb <- sqrt(1/(stat[rws, "n"]*stat[rws, "seb"]+stat[rws, "b"]^2))
        keep <- cumsum(eig$values)/sum(eig$values) < eigen_threshold[j]
        z <- t(eig$vectors[,keep]) %*% (stat[rws, "b"]*scaleb)
        wy <- z / sqrt(eig$values[keep])
        Q <- diag(sqrt(eig$values[keep]))%*%t(eig$vectors[,keep])
        colnames(Q) <- colnames(B)
        LD <- NULL
        LDvalues <- as.list(as.data.frame(Q))
        LDindices <- lapply(1:ncol(Q),function(x) { (1:nrow(Q))-1 } )
        names(LDvalues) <- names(LDindices) <- rsids
      }
      
      n <- mean(stat[rws,"n"])
      m <- length(rsids)
      
      b <- rep(0, m)
      mask <- rep(FALSE, m)
      lambda <- rep(ve/vb,m)

      for (run in 1:nrun) {
        
        if (run > 1) {
          seed[run] <- sample.int(.Machine$integer.max, 1)
        }
        if(algorithm==1) {
          fit <- .Call("_qgg_sbayes_reg",
                       yy=yy,
                       wy=stat$wy[rws],
                       ww=stat$ww[rws],
                       LDvalues=LDvalues,
                       LDindices=LDindices,
                       b = b,
                       lambda = lambda,
                       mask=mask,
                       pi = pi,
                       gamma = gamma,
                       vb = vb,
                       vg = vg,
                       ve = ve,
                       ssb_prior=ssb_prior,
                       ssg_prior=ssg_prior,
                       sse_prior=sse_prior,
                       nub=nub,
                       nug=nug,
                       nue=nue,
                       mh_p=mh_p,
                       mh_r2=mh_r2,
                       updateB = updateB,
                       updateG = updateG,
                       updateE = updateE,
                       updatePi = updatePi,
                       updateMH = updateMH,
                       n=n,
                       nit=nit,
                       nburn=nburn,
                       nthin=nthin,
                       method=as.integer(method),
                       algo=as.integer(algorithm),
                       seed=as.integer(seed[run]))
        }   
        if(algorithm==3) {
          fit <- .Call("_qgg_sbayes_reg_eigen",
                       wy=wy,
                       ww=ww,
                       LDvalues=LDvalues,
                       LDindices=LDindices,
                       b = b,
                       lambda = lambda,
                       mask=mask,
                       pi = pi,
                       gamma = gamma,
                       vb = vb,
                       vg = vg,
                       ve = ve,
                       ssb_prior=ssb_prior,
                       ssg_prior=ssg_prior,
                       sse_prior=sse_prior,
                       nub=nub,
                       nug=nug,
                       nue=nue,
                       mh_p=mh_p,
                       mh_r2=mh_r2,
                       updateB = updateB,
                       updateG = updateG,
                       updateE = updateE,
                       updatePi = updatePi,
                       updateMH = updateMH,
                       n=n,
                       nit=nit,
                       nburn=nburn,
                       nthin=nthin,
                       method=as.integer(method),
                       algo=as.integer(algorithm),
                       seed=as.integer(seed[run]))
        }   
        names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pis","pim","r","b","param","bs","ds","prob")
        if(nrun>1) {
          if (run==1) {
            rbm <- rep(0,length(fit$bm)) 
            rdm <- rep(0,length(fit$dm)) 
            rvbs <- rep(0,length(fit$vbs)) 
            rvgs <- rep(0,length(fit$vgs)) 
            rves <- rep(0,length(fit$ves)) 
            rpis <- rep(0,length(fit$pis)) 
            rpim <- rep(0,length(fit$pim)) 
            rbs <- rep(0,length(fit$bs)) 
            rds <- rep(0,length(fit$ds)) 
            rprob <- rep(0,length(fit$prob)) 
          }
          
          
          rbm <- rbm + fit$bm/nrun
          rdm <- rdm + fit$dm/nrun
          rvbs <- rvbs + fit$vbs/nrun
          rvgs <- rvgs + fit$vgs/nrun
          rves <- rves + fit$ves/nrun
          rpis <- rpis + fit$pis/nrun
          rpim <- rpim + fit$pim/nrun
          rbs <- rbs + fit$bs/nrun
          rds <- rds + fit$ds/nrun
          rprob <- rprob + fit$prob/nrun
          if(run==nrun) {
            fit$bm <- rbm
            fit$dm <- rdm
            fit$vbs <- rvbs
            fit$vgs <- rvgs
            fit$ves <- rves
            fit$pis <- rpis
            fit$pim <- rpim
            fit$bs <- rbs
            fit$ds <- rds
            fit$prob <- rprob
          }
        }
      }
      
      if(algorithm==3) fit$bm <- fit$bm/scaleb
      names(fit$bm) <- names(fit$dm) <- names(fit$b) <- names(LDvalues)
      fit$bs <- matrix(fit$bs,nrow=length(fit$bm))
      fit$ds <- matrix(fit$ds,nrow=length(fit$bm))
      fit$prob <- matrix(fit$prob,nrow=length(fit$bm))
      rownames(fit$bs) <- rownames(fit$ds) <- rownames(fit$prob) <- names(LDvalues)
      colnames(fit$bs) <- colnames(fit$ds) <- colnames(fit$prob) <- 1:(nit+nburn)
      
      if(algorithm==3) {
        for (k in 1:nrow(fit$bs)) {
          fit$bs[k,] <- fit$bs[k,]/scaleb[k]
        }
      }
      
      # Check convergence            
      critve <- critvg <- critvb <- critpi <- critb <- FALSE
      if(!updateE) critve <- TRUE
      if(!updateG) critvg <- TRUE
      if(!updateB) critvb <- TRUE
      if(!updatePi) critpi <- TRUE
      zve <- coda::geweke.diag(fit$ves[nburn:(nburn+nit)])$z
      zvg <- coda::geweke.diag(fit$vgs[nburn:(nburn+nit)])$z
      zvb <- coda::geweke.diag(fit$vbs[nburn:(nburn+nit)])$z
      zpi <- coda::geweke.diag(fit$pis[nburn:(nburn+nit)])$z
      zb <- coda::geweke.diag(apply(fit$bs[,nburn:(nburn+nit)],2,var))$z
      if(!is.na(zve)) critve <- abs(zve)<critVe
      if(!is.na(zvg)) critvg <- abs(zvg)<critVg
      if(!is.na(zvb)) critvb <- abs(zvb)<critVb
      if(!is.na(zpi)) critpi <- abs(zpi)<critPi
      if(!is.na(zb)) critb <- abs(zb)<critB
      
      critb1 <- critb2 <- FALSE
      
      brws <- fit$dm>0
      
      # Check divergence        
      if(sum(brws)>1 && all(c(critve, critvg, critvb, critpi, critb))) {
        
        tstat <- fit$bs[brws,]/stat[rws, "b"][brws]
        pdiv <- apply(tstat, 1, function(x) {
          sum(x[nburn:length(x)] > -critB1 & x[nburn:length(x)] <= 1+critB1)
        })
        pdiv <- pdiv/length(nburn:ncol(tstat))
        pdiv <- pdiv[is.finite(pdiv)]
        if(any(pdiv<0.95) && verbose ) plot(pdiv)
        critb1 <- !any(pdiv<0.95)    # FALSE if any pdiv is less than 0.95
        if(!critb1) message(paste("Convergence not reached for critB1 "))
        critb <- critb1
      }
      
      # Check mismatch
      if(checkLD && sum(brws)>1 && !all(c(critve, critvg, critvb, critpi, critb))) {
        # Identify mismatch between LD and summary statistics 
        bout <- checkb(B=B[brws,brws],
                       b=stat[rws, "b"][brws],
                       seb=stat[rws, "seb"][brws],
                       critB=critB2, verbose=verbose)
        critb2 <- !any(bout$outliers)    # FALSE if there are any outliers
        if(critb2) message(paste("Convergence not reached for critB2 "))
        critb <- critb2
      }
      
      converged <- critve & critvg & critvb & critpi & critb
      if(!checkConvergence) converged <- TRUE
      # Make plots to monitor convergence
      if(verbose) {
        layout(matrix(1:4,ncol=2))
        pipsets <- splitWithOverlap(1:length(rsids),100,99)
        pip <- fit$dm
        plot(pip, ylim=c(0,max(pip)), ylab="PIP",xlab="Position", frame.plot=FALSE)
        plot(-log10(stat[rws,"p"]), ylab="-log10(P)",xlab="Position", frame.plot=FALSE)
        hist(fit$ves, main="Ve", xlab="")
        plot(y=fit$bm, x=stat[rws,"b"], ylab="Adjusted",xlab="Marginal", frame.plot=FALSE)
        abline(h=0,v=0, lwd=2, col=2, lty=2)
      }
      attempts[i] <- j      
      if(!converged) {
        message(paste("Convergence not reached using eigen_threshold:",eigen_threshold[j]))
        criteria_names <- c("Variance of errors (critve)", 
                            "Genetic variance (critvg)", 
                            "Marker variance (critvb)", 
                            "Inclusion probability (critpi)", 
                            "Posterior mean (critb)")
        criteria_status <- c(critve, critvg, critvb, critpi, critb)
        message("Convergence criteria:")
        for (k in seq_along(criteria_names)) {
          message(sprintf("  %s: %s", criteria_names[k], ifelse(criteria_status[k], "Met", "Not Met")))
        }
        if (checkConvergence & algorithm == 1) {
          message("Serious convergence issues detected. Please check your summary statistics and LD reference data.")
          message("If using in-sample LD and summary statistics, consider increasing the default values for the following parameters: critVe=3, critVg=3, critVb=3, critPi=3, critB=3.")
          message('If using external LD and summary statistics, consider using algorithm="mcmc-eigen" as a more robust method.')
          message("Program terminated.")
          return(NULL)  
        }
      }
      # Exit outer loop if convergence is reached
      if (converged) {
        if(verbose & algorithm==1) message("Convergence reached")
        if(verbose & algorithm==3) message(paste("Convergence reached using eigen_threshold:",eigen_threshold[j]))
        break
      }
    }
    if(verbose) message("Computing FDR")

    cutoffs <- seq(0.01, 0.99, by = 0.01)  # Generate 1:99 as fractions
    cutoff_indices <- lapply(cutoffs, function(cutoff) fit$dm > cutoff)
    bfdrs <- sapply(cutoff_indices, function(rws) {
      if (any(rws)) {
        fdrs <- rowMeans(1 - fit$prob[rws, , drop = FALSE], na.rm = TRUE)
        c(mean = mean(fdrs, na.rm = TRUE), quantile(fdrs, c(0.025, 0.975), na.rm = TRUE))
      } else {
        c(NA, NA, NA)
      }
    })
    bfdrs <- t(bfdrs)
    rownames(bfdrs) <- round(cutoffs, 2)
    

    # Save results
    if(verbose) message("Saving results")
    bm[[i]] <- fit$bm
    dm[[i]] <- fit$dm
    pim[[i]] <- fit$pim
    ves[[i]] <- fit$ves
    vbs[[i]] <- fit$vbs
    vgs[[i]] <- fit$vgs
    pis[[i]] <- fit$pis
    conv[[i]] <- c(zve,zvg,zvb,zpi,zb) 
    if(output=="full") {
      bs[[i]] <- fit$bs
      ds[[i]] <- fit$ds
      prob[[i]] <- fit$prob
    }
    fdr[[i]] <- bfdrs
    logcpo[i] <- fit$param[4]
    if(verbose) message("Compute credible sets")
    if(sum(fit$dm)>cs_threshold) csets[[i]] <- crs(prob=fit$dm, B=B, 
                                                   threshold=cs_threshold, 
                                                   r2=cs_r2, 
                                                   method=cs_method)
    names(bm[[i]]) <- names(dm[[i]]) <- rsids
  }
  fit <- NULL
  fit$bm <- bm
  fit$dm <- dm
  fit$pim <- pim
  fit$ves <- ves
  fit$vbs <- vbs
  fit$vgs <- vgs
  fit$pis <- pis
  if(output=="full") {
    fit$bs <- bs
    fit$ds <- ds
    fit$prob <- prob
  }
  fit$fdr <- fdr
  fit$logcpo <- logcpo
  fit$cs <- csets
  
  pip <- sapply(fit$dm,sum)
  minb <- sapply(fit$bm,min)
  maxb <- sapply(fit$bm,max)
  m <- sapply(fit$bm,length)
  
  bm <- unlist(unname(fit$bm))
  dm <- unlist(unname(fit$dm))
  marker <- data.frame(rsids=unlist(Glist$rsids),
                       chr=unlist(Glist$chr), pos=unlist(Glist$pos), 
                       ea=unlist(Glist$a1), nea=unlist(Glist$a2),
                       eaf=unlist(Glist$af),stringsAsFactors = FALSE)
  marker <- marker[marker$rsids%in%names(bm),]
  fit$stat <- data.frame(marker,bm=bm[marker$rsids],
                         dm=dm[marker$rsids], stringsAsFactors = FALSE)
  fit$stat$vm <- 2*(1-fit$stat$eaf)*fit$stat$eaf*fit$stat$bm^2
  fit$method <- methods[method+1]
  fit$mask <- mask
  ve <- sapply(fit$ves,function(x){mean(x[nburn:length(x)])})
  vg <- sapply(fit$vgs,function(x){mean(x[nburn:length(x)])})
  vb <- sapply(fit$vbs,function(x){mean(x[nburn:length(x)])})
  pi <- sapply(fit$pis,function(x){mean(x[nburn:length(x)])})
  ve_ci <- t(sapply(fit$ves,function(x){quantile(x[nburn:length(x)], c(0.025,0.975))}))
  vg_ci <- t(sapply(fit$vgs,function(x){quantile(x[nburn:length(x)], c(0.025,0.975))}))
  vb_ci <- t(sapply(fit$vbs,function(x){quantile(x[nburn:length(x)], c(0.025,0.975))}))
  pi_ci <- t(sapply(fit$pis,function(x){quantile(x[nburn:length(x)], c(0.025,0.975))}))
  
  fit$ci <- list(ve=cbind(mean=ve,ve_ci),
                 vg=cbind(mean=vg,vg_ci), 
                 vb=cbind(mean=vb,vb_ci), 
                 pi=cbind(mean=pi,pi_ci))  
  
  if(!is.null(Glist$map)) map <- unlist(Glist$map)
  pos <- unlist(Glist$pos)
  sets <- lapply(fit$bm,names)
  setsindex <- mapSets(sets=sets, rsids=unlist(Glist$rsids))
  if(!is.null(Glist$map)) cm <- sapply(setsindex, function(x){ max(map[x])-min(map[x]) })
  mb <- sapply(setsindex, function(x){ (max(pos[x])-min(pos[x]))/1000000 })
  minmb <- sapply(setsindex, function(x){ min(pos[x]) })
  maxmb <- sapply(setsindex, function(x){ max(pos[x]) })
  
  chr <- unlist(Glist$chr)
  chr <- sapply(setsindex,function(x){as.numeric(unique(chr[x]))[1]})
  
  b <- stat[fit$stat$rsids,"b"]
  
  conv <- t(as.data.frame(conv))
  colnames(conv) <- c("zve","zvg","zvb","zpi","zb")
  fit$conv <- data.frame(conv,ntrials=attempts, cutoff=eigen_threshold[attempts])
  if(is.null(Glist$map)) fit$post <- data.frame(ve=ve,vg=vg, vb=vb, pi=pi, pip=pip, minb=minb, maxb=maxb, m=m, mb=mb, chr=chr, minmb=minmb, maxmb=maxmb)  
  if(!is.null(Glist$map)) fit$post <- data.frame(ve=ve,vg=vg, vb=vb, pi=pi, pip=pip, minb=minb, maxb=maxb, m=m, mb=mb, cm=cm, chr=chr, minmb=minmb, maxmb=maxmb)  
  rownames(fit$conv) <- rownames(fit$post) <- names(sets) 
  
  fit$ve <- mean(ve)
  fit$vg <- sum(vg)
  fit$b <- b
  fit$seed <- seed
  return(fit)
}

cpo <- function(yobs=NULL, ypred=NULL, nit=NULL, nburn=nburn) {
  psum <- rep(0,nrow(ypred))
  for (i in 1:ncol(ypred)) {
    x <- (yobs-ypred[,i])[,1]
    p <- (1/sqrt(2*base::pi))*exp(-0.5*(x^2))
    psum <- psum + 1/p
  }
  sum(log((nit-nburn)*(1/psum)))
}


crs <- function(prob = NULL, B = NULL, threshold = 0.9, r2 = 0.5, keep = FALSE, cutoff=0.001, method="CS2") {
  
  # Input validation
  if (is.null(prob) || is.null(B)) stop("Both 'prob' and 'B' must be provided.")
  if (!is.vector(prob) || !is.matrix(B)) stop("'prob' must be a vector and 'B' must be a matrix.")
  if (is.null(names(prob)) || is.null(rownames(B)) || is.null(colnames(B))) {
    stop("'prob' and 'B' must have names for proper indexing.")
  }
  
  
  # Step 1: Sort PIPs in descending order
  prob[prob<cutoff] <- 0 
  dsorted <- sort(prob, decreasing = TRUE)
  
  if(method=="CS3") {
    bprob <- B^2%*%prob
    o <- order(bprob, decreasing=TRUE)
    dsorted <- prob[o]
  }
  
  credible_sets <- list()  # Initialize as an empty list
  
  # Identify LD sets and map to sorted PIPs
  sets <- sapply(names(dsorted), function(x) {
    colnames(B)[B[x, ]^2 > r2]  # Identify SNPs in LD
  }, simplify = FALSE)  
  
  sets <- mapSets(sets = sets, rsids = names(dsorted), index = FALSE)  # Ensure mapSets is properly defined
  
  # Step 2: Identify credible sets
  for (j in seq_along(sets)) {
    if (length(sets[[j]]) == 0) next  # Skip empty sets
    cumulative_pip <- sum(dsorted[names(dsorted)%in%sets[[j]]], na.rm = TRUE)  # Avoid NA issues
    
    if (cumulative_pip >= threshold) {
      dset <- dsorted[names(dsorted)%in%sets[[j]]]
      dset <- dset[dset > 0]  # Keep only non-zero PIPs
      
      if (length(dset) > 0 && any(cumsum(dset) >= threshold)) {
        crset <- names(dset)[1:which(cumsum(dset) >= threshold)[1]]
      } else {
        crset <- names(dset)  # Default to full set if threshold isn't met
      }
      
      credible_sets[[length(credible_sets) + 1]] <- crset
      names(credible_sets)[length(credible_sets)] <- paste0("Set", length(credible_sets))
      
      # Set prob values to zero for markers in the credible set
      dsorted[sets[[j]]] <- 0
    }
  }
  
  # Return results
  if (length(credible_sets) == 0) return(list(NULL))
  return(credible_sets)
}


# crs <- function(prob = NULL, B = NULL, threshold = 0.8, r2 = 0.5, keep = FALSE) {
#   # Input validation
#   if (is.null(prob) || is.null(B)) stop("Both 'prob' and 'B' must be provided.")
#   if (!is.vector(prob) || !is.matrix(B)) stop("'prob' must be a vector and 'B' must be a matrix.")
#   if (is.null(names(prob)) || is.null(rownames(B)) || is.null(colnames(B))) {
#     stop("'prob' and 'B' must have names for proper indexing.")
#   }
#   
#   # Step 1: Sort PIPs in descending order
#   dsorted <- sort(prob, decreasing = TRUE)
#   credible_sets <- NULL
#   
#   # Step 2: Identify credible sets of size 1
#   high_pip_markers <- names(dsorted)[dsorted > threshold]
#   
#   if (length(high_pip_markers) > 0) {
#     # Add markers with PIP > threshold as credible sets of size 1
#     credible_sets <- as.list(high_pip_markers)
#     names(credible_sets) <- paste0("Set", seq_along(credible_sets))
#     
#     # Remove these markers from further analysis
#     dsorted <- dsorted[!(names(dsorted) %in% high_pip_markers)]
#   }
#   
#   # Step 3: Identify credible sets of size > 1
#   remaining_markers <- names(dsorted)
#   k <- 1
#   while (sum(dsorted)>= threshold) {
#     #lead_marker <- remaining_markers[1]
#     lead_marker <- remaining_markers[1:min(length(remaining_markers), k)]
#     
#     # Step 3.1: Identify LD friends of the lead marker
#     #ld_friends <- colnames(B)[B[lead_marker, ]^2 > r2]
#     #ld_friends <- intersect(ld_friends, remaining_markers[-1])  # Only include unprocessed markers
#     ld_friends <- sapply(lead_marker, function(x) {
#       colnames(B)[B[x, ]^2 > r2]
#     })
#     ld_friends <- unique(unlist(ld_friends))
#     # Only include unprocessed markers
#     ld_friends <- intersect(ld_friends, remaining_markers[-c(1:k)])  
#     
#     #if (length(ld_friends) > 0) {
#     #  cumulative_pip <- sum(dsorted[c(lead_marker, ld_friends)])
#     #} else {
#     #  cumulative_pip <- dsorted[lead_marker]
#     #}
#     # Compute cumulative PIP for the lead marker and LD friends
#     relevant_markers <- intersect(names(dsorted), c(lead_marker, ld_friends))
#     cumulative_pip <- sum(dsorted[relevant_markers])
#     
#     # Step 3.2: Check if cumulative PIP exceeds threshold
#     if (cumulative_pip >= threshold) {
#       #dset <- dsorted[names(dsorted) %in% c(lead_marker, ld_friends)]
#       #crset <- names(dset)[1:which(cumsum(dset) >= threshold)[1]]
#       dset <- dsorted[names(dsorted) %in% relevant_markers]
#       # Ensure a valid credible set index before indexing
#       if (any(cumsum(dset) >= threshold)) {
#         crset <- names(dset)[1:which(cumsum(dset) >= threshold)[1]]
#       } else {
#         crset <- names(dset)  # Default to full set if no threshold is met
#       }
#       credible_sets[[length(credible_sets) + 1]] <- crset
#       names(credible_sets)[length(credible_sets)] <- paste0("Set", length(credible_sets))
# 
#       # Remove credible set markers from further analysis
#       if (keep) {
#         remaining_markers <- setdiff(remaining_markers, crset)
#       } else {
#         remaining_markers <- setdiff(remaining_markers, c(lead_marker, ld_friends))
#       }
#       k <- 1  # Reset k after finding a set
#     } else {
#       # If cumulative PIP does not exceed threshold, move to the next marker
#       #remaining_markers <- remaining_markers[-1]
#       k <- k + 1
#     }
#     
#     # Update dsorted dynamically based on remaining markers
#     dsorted <- dsorted[remaining_markers]
#   }
#   
#   # Return results
#   if(is.null(credible_sets)) return(list(NULL))
#   return(credible_sets)
# }

bfdr <- function(prob = NULL, probs = NULL, cutoff = 0.1, ql = 0.025, qu = 0.975, verbose = TRUE) {
  # prob is a vector mx1 of posterior means of inclusion probability for m variables
  # probs is a mxnit matrix of inclusion probabilities for m variables over nit iterations
  # cutoff is a chosen cutoff (between 0 and 1) that defines the discovery set
  # ql and qu is the lower and upper quantiles for fdr
  # Validate inputs
  if (is.null(prob) || is.null(probs)) {
    stop("Both 'prob' and 'probs' must be provided.")
  }
  if (!is.numeric(prob) || !is.matrix(probs)) {
    stop("'prob' must be a numeric vector and 'probs' must be a matrix.")
  }
  if (length(prob) != nrow(probs)) {
    stop("The length of 'prob' must match the number of rows in 'probs'.")
  }
  if (cutoff <= 0 || cutoff >= 1) {
    stop("'cutoff' must be between 0 and 1.")
  }
  
  # Identify discovery set
  rws <- prob > cutoff
  if (sum(rws) == 0) {
    warning("No variables selected with the given cutoff.")
    return(NULL)
  }
  
  # Compute FDRs
  fdrs <- apply(1 - probs[rws, , drop = FALSE], 2, mean)
  
  # Plot results if requested
  if (verbose) {
    layout(matrix(1:2, ncol = 2))
    matplot(t(probs[rws, , drop = FALSE]), frame.plot = FALSE, 
            ylab = "Posterior Probability", xlab = "Iteration", type = "l")
    hist(fdrs, breaks = 30, xlab = "McMC-Bayes FDR", 
         main = NULL, freq = TRUE, col = "lightblue")
  }
  
  # Return summary statistics
  list(fdr=c(mean = mean(fdrs, na.rm = TRUE), quantile(fdrs, c(ql, qu), na.rm = TRUE)),
       set=names(prob)[rws])
}

lcpo <- function(yobs=NULL, ypred=NULL, nit=NULL, nburn=NULL) {
  psum <- rep(0,nrow(ypred))
  for (i in 1:ncol(ypred)) {
    x <- (yobs-ypred[,i])[,1]
    p <- (1/sqrt(2*base::pi))*exp(-0.5*(x^2))
    psum <- psum + 1/p
  }
  sum(log((nit-nburn)*(1/psum)))
}

check_divergence <- function(bm, bs, ci_level = 0.95) {
  # Ensure bm is a vector and bs is a matrix
  if (!is.vector(bm) || !is.matrix(bs)) {
    stop("bm must be a vector and bs must be a matrix")
  }
  
  # Check dimensions
  if (length(bm) != nrow(bs)) {
    stop("Length of bm must match the number of rows in bs")
  }
  
  # Initialize results
  n <- length(bm)
  ci_divergence <- logical(n)
  condition1 <- logical(n)
  condition2 <- logical(n)
  
  # Loop through each row to compute confidence intervals and conditions
  for (i in 1:n) {
    # Calculate confidence interval for row i of bs
    alpha <- (1 - ci_level) / 2
    ci_lower <- quantile(bs[i, ], alpha)
    ci_upper <- quantile(bs[i, ], 1 - alpha)
    
    # Check if bm[i] falls within the confidence interval
    ci_divergence[i] <- !(ci_lower <= bm[i] & bm[i] <= ci_upper)
    
    # Conditions based on CI
    condition1[i] <- bm[i] < 0 & ci_upper > 0    # bm[i] < 0 and upper bound > 0
    condition2[i] <- bm[i] > 0 & ci_lower < 0    # bm[i] > 0 and lower bound < 0
  }
  
  # Combine all divergence checks
  overall_divergence <- condition1 | condition2
  
  # Return results as a list
  results <- list(
    ci_divergence = ci_divergence,
    condition1 = condition1,
    condition2 = condition2,
    overall_divergence = overall_divergence
  )
  
  return(results)
}



#' Get Credible Sets
#'
#' This function calculates credible sets for a list of marker sets using probabilistic measures.
#'
#' @param Glist A list containing genomic data (e.g., `chr` and `bedfiles`).
#' @param fit A model object containing SNP statistics (`stat` with `rsids` and `dm` columns).
#' @param sets A list of marker sets (SNP identifiers).
#' @return A list of credible sets for each marker set.
#' 
#' @export
#'
getCredibleSets <- function(Glist=NULL, fit=NULL, sets=NULL) {
  # Validate inputs
  if (is.null(Glist) || is.null(fit) || is.null(sets)) {
    stop("All inputs (Glist, fit, sets) must be provided.")
  }
  
  if (!all(c("rsids", "dm") %in% colnames(fit$stat))) {
    stop("fit$stat must contain 'rsids' and 'dm' columns.")
  }
  
  sets <- mapSets(sets = sets, rsids = fit$stat$rsids, index = FALSE)
  rownames(fit$stat) <- fit$stat$rsids
  if (any(sapply(sets, function(x) { any(is.na(x)) }))) {
    stop("NAs detected in sets - please remove these.")
  }
  
  # Chromosome mapping
  chr <- as.character(unlist(Glist$chr))
  chrSets <- sapply(mapSets(sets = sets, Glist = Glist, index = TRUE), 
                    function(x) unique(chr[x]))
  if (length(Glist$bedfiles) == 1) {
    chrSets <- setNames(rep(1, length(sets)), names(sets))
  }
  
  # Validate chromosome mapping
  lsets <- sapply(chrSets, length)
  if (any(lsets > 1)) {
    stop(paste("The following marker sets map to multiple chromosomes:", 
               paste(which(lsets > 1), collapse = ", ")))
  }
  if (any(lsets == 0)) {
    stop(paste("The following marker sets do not map to any chromosome:", 
               paste(which(lsets == 0), collapse = ", ")))
  }
  
  if (length(sets) == 0) {
    warning("No valid marker sets provided. Returning an empty list.")
    return(list())
  }
  
  # Compute credible sets
  cs <- lapply(seq_along(sets), function(i) {
    rsids <- sets[[i]]
    prob <- fit$stat[rsids, "dm"]
    names(prob) <- rsids
    W <- getG(Glist = Glist, chr = chrSets[i], rsids = rsids, scale = TRUE)
    B <- crossprod(scale(W)) / (nrow(W) - 1)
    qgg:::crs(prob = prob, B = B, threshold = 0.8)
  })
  names(cs) <- names(sets)
  return(cs)
}


# Multiple trait BLR using summary statistics and sparse LD provided in Glist
mtblr <- function(yy=NULL, Xy=NULL, XX=NULL, n=NULL,
                  b=NULL,h2=NULL, pi=0.001, models=NULL, pimodels=NULL,
                  vg=NULL, vb=NULL, ve=NULL,
                  ssb_prior=NULL, sse_prior=NULL,
                  updateB=TRUE, updateE=TRUE, updatePi=TRUE,
                  nub=4, nue=4, nit=1000, nburn=500, nthin=1, method="bayesC", verbose=NULL) {
  
  ww <- lapply(XX,diag)
  wy <- lapply(Xy, as.vector)
  m <- mean(sapply(wy,length))
  nt <- length(wy)
  XXvalues <- lapply(XX, function(x){as.list(as.data.frame(x))})
  XXindices <- lapply(1:m,function(x) { (1:m)-1 } )

    
  #if(!method=="bayesC") stop("Only method==bayesC is allowed")
  #method <- 4
  if(method=="bayesC") method <- 4
  if(method=="bayesR") method <- 5
  
  if(is.null(b)) b <- lapply(1:nt,function(x){rep(0,m)})
  if(is.matrix(b)) b <- split(b, rep(1:ncol(b), each = nrow(b)))
  
  if(is.null(models)) {
    models <- rep(list(0:1), nt)
    models <- t(do.call(expand.grid, models))
    models <- split(models, rep(1:ncol(models), each = nrow(models)))
  }
  if(is.character(models)) {
    if(models=="restrictive") {
      models <- list(rep(0,nt),rep(1,nt))
    }
  }
  if(is.null(pimodels)) pimodels <- c(1-pi,rep(pi/(length(models)-1),length(models)-1)) 
  
  vy <- diag(yy/(n-1),nt)
  if(is.null(h2)) h2 <- 0.5
  if(is.null(vg)) vg <- diag(diag(vy)*h2)
  if(is.null(ve)) ve <- diag(diag(vy)*(1-h2))
  if(method<4 && is.null(vb)) vb <- diag((diag(vy)*h2)/m)
  if(method>=4 && is.null(vb)) vb <- diag((diag(vy)*h2)/(m*pi))
  #if(method>=4 && is.null(vb)) vb <- diag((diag(vy)*h2)/(m*pi[length(models)]))
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  diag(((nub-2.0)/nub)*(diag(vg)/m))
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  diag(((nub-2.0)/nub)*(diag(vg)/(m*pi)))
  #if(method>=4 && is.null(ssb_prior))  ssb_prior <-  diag(((nub-2.0)/nub)*(diag(vg)/(m*pi[length(models)])))
  if(is.null(sse_prior)) sse_prior <- diag(((nue-2.0)/nue)*diag(ve))
  
  seed <- sample.int(.Machine$integer.max, 1)
  
  fit <- .Call("_qgg_mtblr",
               wy=wy,
               ww=ww,
               yy=yy,
               b = b,
               XXvalues=XXvalues,
               XXindices=XXindices,
               B = vb,
               E = ve,
               ssb_prior=split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
               sse_prior=split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
               models=models,
               pi=pimodels,
               nub=nub,
               nue=nue,
               updateB = updateB,
               updateE = updateE,
               updatePi = updatePi,
               n=n,
               nit=nit,
               nburn=nburn,
               nthin=nthin,
               seed=seed,
               method=as.integer(method))

  # fit[[6]] <- matrix(unlist(fit[[6]]), ncol = nt, byrow = TRUE)
  # fit[[7]] <- matrix(unlist(fit[[7]]), ncol = nt, byrow = TRUE)
  # trait_names <- names(yy)
  # if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
  # colnames(fit[[6]]) <- rownames(fit[[6]]) <- trait_names
  # colnames(fit[[7]]) <- rownames(fit[[7]]) <- trait_names
  # fit[[11]] <- matrix(unlist(fit[[11]]), ncol = nt, byrow = TRUE)
  # fit[[12]] <- matrix(unlist(fit[[12]]), ncol = nt, byrow = TRUE)
  # colnames(fit[[11]]) <- rownames(fit[[11]]) <- trait_names
  # colnames(fit[[12]]) <- rownames(fit[[12]]) <- trait_names
  # # add colnames/rownames e, g and gm
  # # add colnames/rownames rg and covg
  # fit[[13]] <- fit[[13]][[1]]
  # fit[[14]] <- fit[[14]][[1]]
  # fit[[15]] <- fit[[15]][[1]]
  # fit[[16]] <- fit[[6]]
  # fit[[17]] <- fit[[7]]
  # if(sum(diag(fit[[16]]))>0) fit[[16]] <- cov2cor(fit[[16]])
  # if(sum(diag(fit[[17]]))>0)  fit[[17]] <- cov2cor(fit[[17]])
  # for(i in 1:nt){
  #   names(fit[[1]][[i]]) <- names(XXvalues)
  #   names(fit[[2]][[i]]) <- names(XXvalues)
  #   names(fit[[10]][[i]]) <- names(XXvalues)
  # }
  # names(fit[[1]]) <- trait_names
  # names(fit[[2]]) <- trait_names
  # names(fit[[3]]) <- trait_names
  # names(fit[[4]]) <- trait_names
  # names(fit[[5]]) <- trait_names
  # names(fit[[8]]) <- trait_names
  # names(fit[[9]]) <- trait_names
  # names(fit[[13]]) <- sapply(models,paste,collapse="_")
  # names(fit[[14]]) <- sapply(models,paste,collapse="_")
  # 
  # names(fit) <- c("bm","dm","coef","vbs","ves","covb","cove",
  #                 "wy","r","b","covg","ve","pi","pim","order",
  #                 "rb","re")
  # fit$bm <- as.matrix(as.data.frame(fit$bm))
  # fit$dm <- as.matrix(as.data.frame(fit$dm))
  # fit$b <- as.matrix(as.data.frame(fit$b))
  names(fit) <- c("bm","dm","wy","r","b","d","o",
                  "vbs","vgs","ves",
                  "covb","covg","cove",
                  "vb","vg","ve",
                  "pi","pim","pitrait","pimarker")
  
  trait_names <- names(yy)
  if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
  variable_names <- names(XXvalues[[1]])
  if(is.null(variable_names)) variable_names <- paste0("V",1:length(XXvalues[[1]]))
  
  for(i in 1:7){
    fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
    rownames(fit[[i]]) <- variable_names
    colnames(fit[[i]]) <- trait_names 
  }
  for(i in 8:10){
    fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
    rownames(fit[[i]]) <- paste0("Iter",1:nrow(fit[[i]]))
    colnames(fit[[i]]) <- trait_names 
  }
  
  for(i in 11:16){
    fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
    colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
  }
  fit[[17]] <- fit[[17]][[1]]
  fit[[18]] <- fit[[18]][[1]]
  names(fit[[17]]) <- sapply(models,paste,collapse="_")
  names(fit[[18]]) <- sapply(models,paste,collapse="_")
  if(method==4) {
    fit <- fit[1:18]
  }
  if(method==5) {
    fit[[19]] <- as.matrix(as.data.frame(fit[[19]]))
    rownames(fit[[19]]) <- c("0","0.01","0.1","1.0")
    colnames(fit[[19]]) <- trait_names
    fit[[20]] <- fit[[20]][[1]]
    names(fit[[20]]) <- c("0","1")
  }
  if(sum(diag(fit$covb))>0) fit$rb <- cov2cor(fit$covb)
  if(sum(diag(fit$covg))>0) fit$rg <- cov2cor(fit$covg)
  if(sum(diag(fit$cove))>0) fit$re <- cov2cor(fit$cove)
  
  return(fit)
}


# Multiple trait BLR using summary statistics and sparse LD provided in Glist
mtsblr <- function(stat=NULL, LD=NULL, n=NULL, vy=NULL, scaled=TRUE,
                   b=NULL,h2=0.1, pi=0.001, models=NULL,pimodels=NULL,
                   vg=NULL, vb=NULL, ve=NULL,
                   ssb_prior=NULL, sse_prior=NULL,
                   updateB=TRUE, updateE=TRUE, updatePi=TRUE,
                   nub=4, nue=4, nit=1000, nburn=500, nthin=1, method="bayesC", verbose=NULL) {
  # Check method
  methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
  method <- match(method, methods) - 1
  if( !sum(method%in%c(0:5))== 1 ) stop("Method specified not valid")
  
  if(method=="bayesC") method <- 4
  if(method=="bayesR") method <- 5
  
    
  # Prepare summary statistics input
  if( is.null(stat) ) stop("Please provide summary statistics")
  
  if( is.null(stat$n) ) stop("Please provide summary statistics that include n")
  
  nt <- ncol(stat$b)
  m <- nrow(stat$b)
  trait_names <- colnames(stat$b)
  if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
  b <- matrix(0,nrow=m,ncol=nt)
  rownames(b) <- rownames(stat$b) 
  colnames(b) <- trait_names
  
  # If genotypes are centered and y standardized to unit variance
  if(scaled) ww <- 1/(stat$seb^2 + stat$b^2/stat$n)
  
  # If genotypes are centered and y not standardized
  if(!scaled ) {
    if(is.null(vy)) stop("Please provide phenotypic variance for each trait")
    ww <- t(t(1/(stat$seb^2 + stat$b^2/stat$n))*vy)
  }
  # If genotypes are centered and scaled
  #ww <- stat$n
  
  wy <- stat$b*ww
  
  yy <- (stat$b^2 + (stat$n-2)*stat$seb^2)*ww
  yy <- apply(yy,2,median)
  n <- as.integer(apply(stat$n,2,median))
  
  if(is.null(models)) {
    models <- rep(list(0:1), nt)
    models <- t(do.call(expand.grid, models))
    models <- split(models, rep(1:ncol(models), each = nrow(models)))
  }
  if(is.character(models)) {
    if(models=="restrictive") {
      models <- list(rep(0,nt),rep(1,nt))
    }
  }
  if(is.null(pimodels)) pimodels <- c(1-pi,rep(pi/(length(models)-1),length(models)-1)) 
  
  seed <- sample.int(.Machine$integer.max, 1)
  
  fit <- mt_sbayes_sparse(yy=yy,
                          ww=ww,
                          wy=wy,
                          b=b,
                          LDvalues=LD$values,
                          LDindices=LD$indices,
                          n=n,
                          nit=nit,
                          nburn=nburn,
                          nthin=nthin,
                          seed=seed,
                          pi=pi,
                          pimodels=pimodels,
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
                          models=models,
                          method=method,
                          verbose=verbose)
  
  return(fit)
}

checkb <- function(B = NULL, b = NULL, seb = NULL, critB = 3.0, 
                   shrink = 0.0001, verbose=TRUE) {
  # Validate inputs
  if (is.null(B) || is.null(b) || is.null(seb)) {
    stop("B, b, and seb must be provided.")
  }
  if (ncol(B) != length(b)) {
    stop("Dimensions of B and b do not match.")
  }
  
  m <- length(b)
  bpred <- numeric(m)  # Preallocate for efficiency
  
  # Precompute diagonal shrinkage for speed
  shrink_diag <- diag(shrink, m - 1)
  
  for (i in seq_len(m)) {
    # Avoid repeatedly recalculating B[-i, -i]
    B_sub <- B[-i, -i] + shrink_diag
    Bi <- chol2inv(chol(B_sub))  # Efficient inversion
    bpred[i] <- sum(B[i, -i] %*% Bi %*% b[-i])  # Predicted value
  }
  
  # Calculate outliers
  abs_diff <- abs((bpred - b) / seb)
  outliers <- abs_diff > critB
  
  # Plot outliers
  if(verbose) {
    # Set up a 1x2 layout for two side-by-side plots
    layout(matrix(1:2, ncol = 2))  # Create two panels side by side
    
    # Scatter plot of predicted vs observed values
    plot(
      bpred, b,
      xlab = "Predicted Beta Values",  # Add meaningful x-axis label
      ylab = "Observed Beta Values",      # Add meaningful y-axis label
      main = "Predicted vs Observed",    # Add a descriptive title
      pch = 19,                          # Use solid points for better visibility
      col = "blue",                      # Use a distinct color
      frame.plot = FALSE                 # Remove the box around the plot
    )
    abline(a = 0, b = 1, col = "red", lty = 2)  # Add y = x line for reference
    
    # Plot of scaled residuals
    plot(
      (bpred - b) / seb,
      xlab = "Index",                     # Label x-axis as "Index"
      ylab = "Z statistics",          # Add meaningful y-axis label
      main = "Outlier Statistics",          # Add a descriptive title
      pch = 19,                           # Use solid points
      col = ifelse(abs((bpred - b) / seb) > critB, "red", "black"),  # Highlight outliers
      frame.plot = FALSE                  # Remove the box around the plot
    )
    abline(h = c(-3, 3), col = "red", lty = 2)  # Add horizontal lines at ±3
  }
  # Return results
  return(list(
    bpred = bpred,          # Predicted values
    outliers = outliers,    # Logical vector indicating outliers
    abs_diff = abs_diff     # Scaled differences for diagnostic purposes
  ))
}

# adjustB <- function(b=NULL, LD = NULL, msize=NULL, overlap=NULL, shrink=0.001, threshold=1e-8) {
#   m <- length(b)
#   badj <- rep(0,m)
#   sets <- splitWithOverlap(1:m,msize,overlap)
#   for( i in 1:length(sets) ) {
#     mset <- length(sets[[i]])
#     bset <- b[sets[[i]]]
#     B <- LD[sets[[i]],sets[[i]]]
#     for (j in 1:mset) {
#       Bi <- chol2inv(chol(B[-j,-j]+diag(shrink,mset-1)))
#       badj[sets[[i]][j]] <- sum(B[j,-j]*Bi%*%bset[-j])
#     }
#     plot(y=badj[sets[[i]]],x=b[sets[[i]]],ylab="Predicted", xlab="Observed",  frame.plot=FALSE)
#     abline(0,1, lwd=2, col=2, lty=2)
#   }
#   return(b=badj)
# }

#' Adjust B-values
#'
#' This function adjusts the B-values based on the LD structure and other parameters.
#' The adjustment is done in subsets, and a plot of observed vs. predicted values is produced for each subset.
#'
#' @param b A numeric vector containing the B-values to be adjusted. If NULL (default), no adjustments are made.
#' @param LD A matrix representing the linkage disequilibrium (LD) structure.
#' @param msize An integer specifying the size of the subsets.
#' @param overlap An integer specifying the overlap size between consecutive subsets.
#' @param shrink A numeric value used for shrinkage. Default is 0.001.
#' @param threshold A numeric value specifying the threshold. Default is 1e-8.
#'
#' @return A list containing the adjusted B-values.
#' 
#' @keywords internal
#' @export

adjustB <- function(b=NULL, seb=NULL, LD = NULL, 
                       msize=100, overlap=50, nrep=20, shrink=0.001, threshold=1) {
  sets <- lapply(1:nrep,function(x){sample(1:length(b),msize)})
  #sets <- splitWithOverlap(1:length(b),msize,overlap)
  nsets <- length(sets)
  if(length(sets[[nsets]])<10) {
    sets[[nsets-1]] <- c(sets[[nsets-1]],sets[[nsets]])
    sets <- sets[1:(nsets-1)]
    nsets <- length(sets)
  }
  for( i in 1:length(sets) ) {
    badj <- b
    mset <- length(sets[[i]])
    bset <- b[sets[[i]]]
    B <- LD[sets[[i]],sets[[i]]]
    for (j in 1:mset) {
      Bi <- chol2inv(chol(B[-j,-j]+diag(shrink,mset-1)))
      bj <- sum(B[j,-j]*Bi%*%bset[-j])
      bset[j]  <- 0.5*bj+0.5*b[sets[[i]][j]]
      #bset[j]  <- bj
      b[sets[[i]][j]] <- bset[j] 
    }
    print(sum((b-badj)^2))
    plot(y=badj[sets[[i]]],x=b[sets[[i]]],ylab="Predicted", xlab="Observed",  frame.plot=FALSE)
    abline(0,1, lwd=2, col=2, lty=2)
  }
  return(b=b)
}

adjLDregion <- function(LD=NULL, p=NULL, r2=0.5, thold=1) {
  rsids <- colnames(LD)
  o <- order(p, decreasing = FALSE)
  m <- length(rsids)
  indx1 <- rep(T, m)
  indx2 <- rep(F, m)
  message(paste("Pruning using r2 threshold:",r2))
  message(paste("Pruning using p-value threshold:",thold))
  ldSets <- apply(LD,1,function(x){colnames(LD)[x>r2]})
  ldSets <- mapSets(sets = ldSets, rsids = rsids)
  for (j in o) {
    if (p[j] <= thold) {
      if (indx1[j]) {
        rws <- ldSets[[j]]
        indx1[rws] <- F
        indx2[j] <- T
      }
    }
  }
  return(rsids[!indx2])
}


# Single trait BLR using summary statistics and sparse LD provided in Glist
sblr <- function(stat=NULL, b=NULL, seb=NULL, n=NULL, vy=1,
                 LD=NULL, LDvalues=NULL,LDindices=NULL,
                 mask=NULL, lambda=NULL,
                 vg=NULL, vb=NULL, ve=NULL, h2=NULL, pi=NULL,
                 ssb_prior=NULL, ssg_prior=NULL, sse_prior=NULL, nub=4,  nug=4, nue=4,
                 updateB=TRUE, updateE=TRUE, updatePi=TRUE, updateG=TRUE,
                 adjustE=TRUE, models=NULL,
                 nit=500, nburn=100, nthin=1, method="bayesC", algorithm=1, verbose=FALSE) {
  
  # Check method
  methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
  method <- match(method, methods) - 1
  if( !sum(method%in%c(0:5))== 1 ) stop("Method specified not valid")
  
  # Prepare summary statistics input
  if( is.null(stat) ) stop("Please provide summary statistics")
  m <- nrow(stat)
  if(is.null(mask)) mask <- rep(FALSE, m)
  if(is.null(stat$n)) stat$n <- stat$dfe
  if( is.null(stat$n) ) stop("Please provide summary statistics that include n")
  
  n <- as.integer(median(stat$n))
  
  #ww <- 1/(stat$seb^2 + stat$b^2/stat$n)
  # use this if y is not scaled
  ww <- vy/(stat$seb^2 + stat$b^2/stat$n)
  wy <- stat$b*ww
  if(!is.null(stat$ww)) ww <- stat$ww
  if(!is.null(stat$wy)) wy <- stat$wy
  
  b2 <- stat$b^2
  seb2 <- stat$seb^2
  yy <- (b2 + (stat$n-2)*seb2)*ww
  yy <- median(yy)
  
  # Prepare sparse LD matrix
  if( is.null(LD) ) stop("Please provide LD matrix")
  #LDvalues <- split(LD, rep(1:ncol(LD), each = nrow(LD)))
  #LDvalues=as.list(as.data.frame(LD)) 
  #LDindices <- lapply(1:ncol(LD),function(x) { (1:ncol(LD))-1 } )
  
  # Prepare starting parameters
  vy <- yy/(n-1)
  if(is.null(pi)) pi <- 0.001
  if(is.null(h2)) h2 <- 0.5
  if(is.null(ve)) ve <- vy*(1-h2)
  if(is.null(vg)) vg <- vy*h2
  if(method<4 && is.null(vb)) vb <- vg/m
  if(method>=4 && is.null(vb)) vb <- vg/(m*pi)
  if(is.null(lambda)) lambda <- rep(ve/vb,m)
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/(m*pi))
  if(is.null(sse_prior)) sse_prior <- ((nue-2.0)/nue)*ve
  if(is.null(ssg_prior)) ssg_prior <- ((nug-2.0)/nug)*vg
  if(is.null(b)) b <- rep(0,m)
  
  pi <- c(1-pi,pi)
  gamma <- c(0,1.0)
  if(method==5) pi <- c(0.95,0.02,0.02,0.01)
  if(method==5) gamma <- c(0,0.01,0.1,1.0)
  
  seed <- sample.int(.Machine$integer.max, 1)
  
  fit <- .Call("_qgg_sbayes",
               yy=yy,
               wy=wy,
               ww=ww,
               LDvalues=LD$values,
               LDindices=LD$indices,
               b = b,
               lambda = lambda,
               mask=mask,
               pi = pi,
               gamma = gamma,
               vb = vb,
               vg = vg,
               ve = ve,
               ssb_prior=ssb_prior,
               ssg_prior=ssg_prior,
               sse_prior=sse_prior,
               nub=nub,
               nug=nug,
               nue=nue,
               updateB = updateB,
               updateG = updateG,
               updateE = updateE,
               updatePi = updatePi,
               adjustE = adjustE,
               n=n, 
               nit=nit,
               nburn=nburn,
               nthin=nthin,
               method=as.integer(method),
               algo=as.integer(algorithm),
               seed=seed)
  names(fit[[1]]) <- names(LDvalues)
  names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pis","pim","r","b","param")
  return(fit)
}


# Multiple trait BLR using summary statistics and sparse LD provided in Glist 
mt_sbayes_sparse <- function(yy=NULL, ww=NULL, wy=NULL, b=NULL, 
                             LDvalues=NULL,LDindices=NULL, n=NULL,
                             vg=NULL, vb=NULL, ve=NULL, 
                             ssb_prior=NULL, sse_prior=NULL, 
                             h2=NULL, pi=NULL, pimodels=NULL, updateB=NULL, 
                             updateE=NULL, updatePi=NULL, models=NULL,
                             nub=NULL, nue=NULL, nit=1000, nburn=500, nthin=1, 
                             seed=NULL, method=NULL, verbose=NULL) {
  
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
  if(is.matrix(b)) b <- split(b, rep(1:ncol(b), each = nrow(b)))
  
  if(is.null(models)) {
    models <- rep(list(0:1), nt)
    models <- t(do.call(expand.grid, models))
    models <- split(models, rep(1:ncol(models), each = nrow(models)))
  }
  if(is.character(models)) {
    if(models=="restrictive") {
      models <- list(rep(0,nt),rep(1,nt))
    }
  }
  if(is.null(pimodels)) pimodels <- c(1-pi,rep(pi/(length(models)-1),length(models)-1)) 
  
  vy <- diag(yy/(n-1),nt)
  if(is.null(h2)) h2 <- 0.5
  if(is.null(vg)) vg <- diag(diag(vy)*h2)
  if(is.null(ve)) ve <- diag(diag(vy)*(1-h2))
  if(method<4 && is.null(vb)) vb <- diag((diag(vy)*h2)/m)
  #if(method>=4 && is.null(vb)) vb <- diag((diag(vy)*h2)/(m*pi[length(models)]))
  if(method>=4 && is.null(vb)) vb <- diag((diag(vy)*h2)/(m*pi))
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  diag(((nub-2.0)/nub)*(diag(vg)/m))
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  diag(((nub-2.0)/nub)*(diag(vg)/(m*pi)))
  #if(method>=4 && is.null(ssb_prior))  ssb_prior <-  diag(((nub-2.0)/nub)*(diag(vg)/(m*pi[length(models)])))
  if(is.null(sse_prior)) sse_prior <- diag(((nue-2.0)/nue)*diag(ve))

  fit <- .Call("_qgg_mtsbayes",
               wy=split(wy, rep(1:ncol(wy), each = nrow(wy))),
               ww=split(ww, rep(1:ncol(ww), each = nrow(ww))),
               yy=yy,
               b = b,
               LDvalues=LDvalues, 
               LDindices=LDindices, 
               B = vb,
               E = ve,
               ssb_prior=split(ssb_prior, rep(1:ncol(ssb_prior), each = nrow(ssb_prior))),
               sse_prior=split(sse_prior, rep(1:ncol(sse_prior), each = nrow(sse_prior))),
               models=models,
               pi=pimodels,
               nub=nub,
               nue=nue,
               updateB = updateB,
               updateE = updateE,
               updatePi = updatePi,
               n=n,
               nit=nit,
               nburn=nburn,
               nthin=nthin,
               seed=seed,
               method=as.integer(method))
  
  # fit[[6]] <- matrix(unlist(fit[[6]]), ncol = nt, byrow = TRUE)
  # fit[[7]] <- matrix(unlist(fit[[7]]), ncol = nt, byrow = TRUE)
  # trait_names <- names(yy)
  # if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
  # colnames(fit[[6]]) <- rownames(fit[[6]]) <- trait_names
  # colnames(fit[[7]]) <- rownames(fit[[7]]) <- trait_names
  # fit[[11]] <- matrix(unlist(fit[[11]]), ncol = nt, byrow = TRUE)
  # fit[[12]] <- matrix(unlist(fit[[12]]), ncol = nt, byrow = TRUE)
  # colnames(fit[[11]]) <- rownames(fit[[11]]) <- trait_names
  # colnames(fit[[12]]) <- rownames(fit[[12]]) <- trait_names
  # # add colnames/rownames e, g and gm
  # # add colnames/rownames rg and covg
  # fit[[13]] <- fit[[13]][[1]]
  # fit[[14]] <- fit[[14]][[1]]
  # fit[[15]] <- fit[[15]][[1]]
  # fit[[16]] <- fit[[6]]
  # fit[[17]] <- fit[[7]]
  # fit[[18]] <- fit[[11]]
  # if(sum(diag(fit[[16]]))>0) fit[[16]] <- cov2cor(fit[[16]])
  # if(sum(diag(fit[[17]]))>0)  fit[[17]] <- cov2cor(fit[[17]])
  # if(sum(diag(fit[[18]]))>0)  fit[[18]] <- cov2cor(fit[[18]])
  # for(i in 1:nt){
  #   names(fit[[1]][[i]]) <- names(LDvalues)
  #   names(fit[[2]][[i]]) <- names(LDvalues)
  #   names(fit[[10]][[i]]) <- names(LDvalues)
  # }
  # names(fit[[1]]) <- trait_names
  # names(fit[[2]]) <- trait_names
  # names(fit[[3]]) <- trait_names
  # names(fit[[4]]) <- trait_names
  # names(fit[[5]]) <- trait_names
  # names(fit[[8]]) <- trait_names
  # names(fit[[9]]) <- trait_names
  # names(fit[[13]]) <- sapply(models,paste,collapse="_")
  # names(fit[[14]]) <- sapply(models,paste,collapse="_")
  # 
  # names(fit) <- c("bm","dm","coef","vbs","ves","covb","cove",
  #                 "wy","r","b","covg","ve","pi","pim","order",
  #                 "rb","re","rg")
  # fit$bm <- as.matrix(as.data.frame(fit$bm))
  # fit$dm <- as.matrix(as.data.frame(fit$dm))
  # fit$b <- as.matrix(as.data.frame(fit$b))
  # names(fit) <- c("bm","dm","wy","r","b","d","o",
  #                 "vbs","vgs","ves",
  #                 "covb","covg","cove",
  #                 "vb","vg","ve",
  #                 "pi","pim")
  # 
  # trait_names <- names(yy)
  # if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
  # variable_names <- names(LDvalues)
  # if(is.null(variable_names)) variable_names <- paste0("V",1:length(LDvalues))
  # 
  # for(i in 1:7){
  #   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  #   rownames(fit[[i]]) <- variable_names
  #   colnames(fit[[i]]) <- trait_names 
  # }
  # for(i in 8:10){
  #   fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
  #   rownames(fit[[i]]) <- paste0("Iter",1:nrow(fit[[i]]))
  #   colnames(fit[[i]]) <- trait_names 
  # }
  # 
  # for(i in 11:16){
  #   fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
  #   colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
  # }
  # fit[[17]] <- fit[[17]][[1]]
  # fit[[18]] <- fit[[18]][[1]]
  # names(fit[[17]]) <- sapply(models,paste,collapse="_")
  # names(fit[[18]]) <- sapply(models,paste,collapse="_")
  # if(sum(diag(fit$covb))>0) fit$rb <- cov2cor(fit$covb)
  # if(sum(diag(fit$covg))>0) fit$rg <- cov2cor(fit$covg)
  # if(sum(diag(fit$cove))>0) fit$re <- cov2cor(fit$cove)
  names(fit) <- c("bm","dm","wy","r","b","d","o",
                  "vbs","vgs","ves",
                  "covb","covg","cove",
                  "vb","vg","ve",
                  "pi","pim","pitrait","pimarker")
  
  trait_names <- names(yy)
  if(is.null(trait_names)) trait_names <- paste0("T",1:nt)
  variable_names <- names(LDvalues)
  if(is.null(variable_names)) variable_names <- paste0("V",1:length(LDvalues))
  
  for(i in 1:7){
    fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
    rownames(fit[[i]]) <- variable_names
    colnames(fit[[i]]) <- trait_names 
  }
  for(i in 8:10){
    fit[[i]] <- as.matrix(as.data.frame(fit[[i]]))
    rownames(fit[[i]]) <- paste0("Iter",1:nrow(fit[[i]]))
    colnames(fit[[i]]) <- trait_names 
  }
  
  for(i in 11:16){
    fit[[i]] <- matrix(unlist(fit[[i]]), ncol = nt, byrow = TRUE)
    colnames(fit[[i]]) <- rownames(fit[[i]]) <- trait_names
  }
  fit[[17]] <- fit[[17]][[1]]
  fit[[18]] <- fit[[18]][[1]]
  names(fit[[17]]) <- sapply(models,paste,collapse="_")
  names(fit[[18]]) <- sapply(models,paste,collapse="_")
  if(method==4) {
    fit <- fit[1:18]
  }
  if(method==5) {
    fit[[19]] <- as.matrix(as.data.frame(fit[[19]]))
    rownames(fit[[19]]) <- c("0","0.01","0.1","1.0")
    colnames(fit[[19]]) <- trait_names
    fit[[20]] <- fit[[20]][[1]]
    names(fit[[20]]) <- c("0","1")
  }
  if(sum(diag(fit$covb))>0) fit$rb <- cov2cor(fit$covb)
  if(sum(diag(fit$covg))>0) fit$rg <- cov2cor(fit$covg)
  if(sum(diag(fit$cove))>0) fit$re <- cov2cor(fit$cove)
  
  return(fit)
}


# Multiple trait BLR based on individual level data based on fast algorithm  
mtbayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, bm=NULL, seb=NULL, LD=NULL, n=NULL,
                    vg=NULL, vb=NULL, ve=NULL, formatLD=NULL,
                    ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
                    h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
                    nub=NULL, nue=NULL, nit=NULL, method=NULL, algorithm=NULL, verbose=NULL) {
  
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
  
  vy <- diag(sapply(y,var))
  if(is.null(h2)) h2 <- 0.5
  if(is.null(ve)) ve <- diag(diag(vy)*(1-h2))
  if(method<4 && is.null(vb)) vb <- diag((diag(vy)*h2)/m)
  if(method>=4 && is.null(vb)) vb <- diag((diag(vy)*h2)/(m*pi[length(models)])) 
  if(is.null(vg)) vg <- diag(diag(vy)*h2)
  if(method<4 && is.null(ssb_prior))  ssb_prior <-  diag(((nub-2.0)/nub)*(diag(vg)/m))
  if(method>=4 && is.null(ssb_prior))  ssb_prior <-  diag(((nub-2.0)/nub)*(diag(vg)/(m*pi[length(models)])))
  if(is.null(sse_prior)) sse_prior <- diag(((nue-2.0)/nue)*diag(ve))

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
  names(fit[[3]]) <- trait_names
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
  fit$mus <- as.data.frame(fit$mus)
  fit$vbs <- as.data.frame(fit$vbs)
  fit$ves <- as.data.frame(fit$ves)
  fit$g <- as.data.frame(fit$g)
  fit$e <- as.data.frame(fit$e)
  fit$b <- as.data.frame(fit$b)
  fit$gm <- as.data.frame(fit$gm)
  return(fit)
  
}



cvarspm <- function( spm ) {
  stopifnot( methods::is( spm, "dgCMatrix" ) )
  ans <- sapply( seq.int(spm@Dim[2]), function(j) {
    if( spm@p[j+1] == spm@p[j] ) { return(0) } # all entries are 0: var is 0
    mean <- base::sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) +
      mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) } ) / ( spm@Dim[1] - 1 )
  names(ans) <- spm@Dimnames[[2]]
  ans
}

cmeanspm <- function( spm ) {
  stopifnot( methods::is( spm, "dgCMatrix" ) )
  ans <- sapply( seq.int(spm@Dim[2]), function(j) {
    if( spm@p[j+1] == spm@p[j] ) { return(0) } # all entries are 0: var is 0
    base::sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]} ) 
  names(ans) <- spm@Dimnames[[2]]
  ans
}

#' @importFrom Matrix sparseMatrix
computeStat <- function(X=NULL, y=NULL, scale=FALSE) {

  if (!inherits(X, "sparseMatrix")) {
    stop("The provided matrix is not sparse.")
  }  
  if(is.vector(y)) y <- as.matrix(y)
  if(is.data.frame(y)) y <- as.matrix(y)
  #selected <- intersect(rownames(y),rownames(X))
  
  #if(length(selected)) stop("Number of matching rows in y and X is less than 10")
  #y <- y[selected,]
  #X <- X[selected,]
  #if(scale) y <- scale(y)
  #mu <- Matrix::colMeans(X)
  mu <- cmeanspm(X)
  #sigma <- rep(0,ncol(X))
  #for(i in 1:ncol(X)) { sigma[i] <- var(X[,i]) }
  sigma <- cvarspm(X)
  sigma <- sqrt(sigma)
  
  n <- nrow(X)
  mu_crossprod <- n * crossprod(t(mu))
  XX <- Matrix::crossprod(X)
  XX <- as.matrix((XX - mu_crossprod) / outer(sigma, sigma))
  if(ncol(y)==1) {
    Xy <- Matrix::crossprod(X,y)
    Xy <- as.vector(Xy)
    Xy <- (Xy - mu*sum(y))/sigma
    yy <- sum(y^2)
    n <- length(y)
  }
  if(ncol(y)>1) {
    Xy <- Matrix::crossprod(X,y)
    Xy <- as.matrix(Xy)
    for(i in 1:ncol(y)) {
      Xy[,i] <- (Xy[,i] - mu*sum(y[,i]))/sigma
    }
    Xy <- as.list(as.data.frame(Xy)) 
    yy <- colSums(y^2)
    n <- rep(nrow(y),ncol(y))
  }
  list(XX=XX, Xy=Xy, yy=yy, n=n)
}

#' @importFrom Matrix sparseMatrix
designMatrix <- function(sets=NULL, values=NULL, rowids=NULL, format="sparse") {
  if(format=="sparse") {
    # Compute design matrix for marker sets in sparse format
    is <- qgg::mapSets(sets=sets, rsids=rowids, index=TRUE)
    js <- rep(1:length(is),times=sapply(is,length))
    x <- rep(1,length(js))
    if(!is.null(values)) {
      if(any(!names(values)==rowids)) stop("Mismatch between names(values) and rowids")
      x <- values[unlist(is)]
    }
    W <- sparseMatrix(unlist(is),as.integer(js),x=x)
    indx <- 1:max(sapply(is,max))
    rowids <- rowids[indx]
    colnames(W) <- names(is)
    rownames(W) <- rowids
  }
  if(format=="dense") {
    sets <- mapSets(sets=sets,rsids=rowids, index=TRUE)
    W <- matrix(0,nrow=length(rowids), ncol=length(sets))
    for(i in 1:length(sets)) {
      W[sets[[i]],i] <- 1
    }
    colnames(W) <- names(sets)
    rownames(W) <- rowids
  }
  return(W)
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
#' @param causal indices for "causal" markers 
#' @param chr chromosome to plot 
#' @param what character fro what to plot (e.g. "trace", "bpm", "ve", "vg", "vb") 
#' @keywords internal

#'
#' @export
#'

plotBayes <- function(fit=NULL, causal=NULL, what="bm", chr=1) {
  # if(!is.list(fit[[1]])) {
  #   layout(matrix(1:6,nrow=3,ncol=2))
  #   plot(fit[[1]],ylab="Posterior", xlab="Marker", main="Marker effect", frame.plot=FALSE)  
  #   if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch=4, cex=2, lwd=3 )
  #   #plot(fit[[2]],ylab="Posterior mean", xlab="Marker", pch="✈",,main="Marker indicator", frame.plot=FALSE)  
  #   plot(fit[[2]],ylab="Posterior mean", xlab="Marker", main="Marker indicator", ylim=c(0,1), frame.plot=FALSE)  
  #   if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch=4, cex=2, lwd=3 )
  #   #if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch="✈", cex=2, lwd=3 )
  #   hist(fit[[6]],xlab="Posterior", main="Pi")  
  #   hist(fit[[4]],xlab="Posterior", main="Marker variance")  
  #   hist(fit[[5]],xlab="Posterior", main="Residual variance")  
  #   plot(fit[[6]],xlab="Sample", ylab="Pi", frame.plot=FALSE)  
  # } 
  # if(is.list(fit[[1]])) {
  #   layout(matrix(1:4,2,2))
  #   matplot(as.data.frame(fit[[1]]),ylab="Marker effect", frame.plot=FALSE)  
  #   if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="green", pch=4, cex=2, lwd=3 )
  #   matplot(as.data.frame(fit[[2]]),ylab="Marker indicator", frame.plot=FALSE)  
  #   if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="green", pch=4, cex=2, lwd=3 )
  #   matplot(as.data.frame(fit[[4]]),ylab="Marker variance", frame.plot=FALSE)  
  #   matplot(as.data.frame(fit[[5]]),ylab="Residual variance", frame.plot=FALSE)  
  # } 
  if(length(chr)==1) {
    if(what=="bm") plot(fit[[chr]]$bm,ylab="Posterior Mean Effect", xlab="Marker", main="Adjusted Marker Effect", frame.plot=FALSE)  
    
    if(what=="dm") {
      if(fit$method=="bayesC") {
        plot(fit[[chr]]$dm, ylab="PIP", xlab="Marker", main="Posterior Inclusion Probability", ylim=c(0,1), frame.plot=FALSE)
      }
      if(fit$method=="bayesR") {
        plot(fit[[chr]]$dm, 
             xlab="Marker", 
             ylab="Marker Class",
             frame.plot=FALSE, ylim=c(0,3), main="Posterior Mean Marker Class")
        abline(h=1, lty=2, col=2, lwd=2)
        abline(h=2, lty=2, col=2, lwd=2)
        abline(h=3, lty=2, col=2, lwd=2)
      } 

    }
    if(what=="ve") hist(fit[[chr]]$ves,xlab="Posterior", main="Residual Variance")
    if(what=="vb") hist(fit[[chr]]$vbs,xlab="Posterior", main="Marker Variance")
    if(what=="vg") hist(fit[[chr]]$vgs,xlab="Posterior", main="Genetic Variance")
    fit[[chr]]$h2 <- fit[[chr]]$vgs/(fit[[chr]]$vgs+fit[[chr]]$ves)
    if(what=="h2") hist(fit[[chr]]$h2,xlab="Posterior", main="Heritability")  
  }
  if(length(chr)==1 & what=="trace") {
    layout(matrix(1:4,ncol=2))
    plot(fit[[chr]]$ves, ylab="Ve", xlab="Iteration", main="Residual Variance")
    plot(fit[[chr]]$vbs, ylab="Vb", xlab="Iteration", main="Marker Variance")
    plot(fit[[chr]]$vgs, ylab="Vg", xlab="Iteration", main="Genetic Variance")
    fit[[chr]]$h2 <- fit[[chr]]$vgs/(fit[[chr]]$vgs+fit[[chr]]$ves)
    plot(fit[[chr]]$h2, ylab="h2", xlab="Iteration", main="Heritability")  
  }
  
  if(length(chr)>1) {
    if(what=="Ve") {
      x <- sapply(fit[chr], function(x) {mean(x$ves)})
      sd <- sapply(fit[chr], function(x) {sd(x$ves)})
      main <- "Residual Variance"
    }
    if(what=="Vg") {
      x <- sapply(fit[chr], function(x) {mean(x$vgs)})
      sd <- sapply(fit[chr], function(x) {sd(x$vgs)})
      main <- "Genetic Variance"
    }
    if(what=="Vb") {
      x <- sapply(fit[chr], function(x) {mean(x$vbs)})
      sd <- sapply(fit[chr], function(x) {sd(x$vbs)})
      main <- "Marker Variance"
    }
    if(what=="h2") {
      x <- sapply(fit[chr], function(x) {mean(x$vgs/(x$vgs+x$ves))})
      sd <- sapply(fit[chr], function(x) {sd(x$vgs/(x$vgs+x$ves))})
      main <- "Heritability"
    }
    names(x) <- names(sd) <- paste("Chr",chr)
    plotForest(x=x,sd=sd, reorder=FALSE, xlab=what, main=main) 
  }
  
}

plotCvs <- function(fit=NULL, causal=NULL) {
  
  layout(matrix(1:length(fit$p),ncol=1))
  for (i in 1:length(fit$p)) {
    plot(-log10(fit$p[[i]]),ylab="mlogP", xlab="Marker", main="Marker effect", frame.plot=FALSE)  
    if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch=4, cex=2, lwd=3 )
  }
}


#' Split Vector with Overlapping Segments
#'
#' Splits a vector into segments of a specified length with a specified overlap.
#' The function returns a list where each element contains a segment of the vector.
#'
#' @param vec A numeric or character vector to be split.
#' @param seg.length An integer specifying the length of each segment.
#' @param overlap An integer specifying the number of overlapping elements between consecutive segments.
#'
#' @return A list where each element is a segment of the input vector. The segments can overlap based on the specified overlap.
#'
#' @keywords internal
#' @export
splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

#' Adjust Linkage Disequilibrium (LD) Using Map Information
#'
#' Adjusts the linkage disequilibrium (LD) values based on map information, effective population size,
#' and sample size used for map construction.
#'
#' @param LD A matrix representing the linkage disequilibrium (LD) structure.
#' @param map A numeric vector containing the map information.
#' @param neff An integer representing the effective population size. Default is 11600.
#' @param nmap An integer representing the sample size used for map construction. Default is 186.
#' @param threshold A numeric value specifying the threshold for setting LD to zero. Currently unused in the function. Default is 0.001.
#'
#' @return A matrix where each value is the adjusted LD based on the input map, effective population size, and map construction sample size.
#'
#' @keywords internal
#' @export
adjustMapLD <- function(LD = NULL, map=NULL, neff=11600, nmap=186, threshold=0.001) {
  # neff: effective population size
  # nmap: sample size used for map contruction
  # threshold: used for setting LD to zero
  rho <- sapply(map,function(x){map-x})
  rho <- 4*neff*abs(rho)
  shrink <- exp(-rho/(2*nmap))
  return(LD*shrink)
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

bmm <- function(y=NULL, X=NULL, W=NULL, GRMlist=NULL,
                vg=NULL, ve=NULL, nug=NULL, nue=NULL,
                vg_prior=NULL, ve_prior=NULL,
                updateG=TRUE, updateE=TRUE,
                nit=500, nburn=0, tol=0.001, verbose=FALSE) {
  if(is.vector(y)) y <- as.matrix(y)
  n <- nrow(y)                            # number of observation
  nt <- ncol(y)                           # number of traits
  tnames <- colnames(y)
  if(is.null(tnames)) tnames <- paste0("T",1:nt)
  e <- matrix(0,nrow=n,ncol=nt)
  mu <- matrix(0,nrow=n,ncol=nt)
  mus <- matrix(0,nrow=n,ncol=nt)
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

  if(any(!idsT%in%idsG)) stop("Some ids in y not in GRM")
  
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
        scg <- sum((1/D[[j]])*a**2) + (vg_prior[[j]]*(nug[j]+2))/nug[j]	# => S = (mode*(df+2))/df         
        #scg <- sum((1/D[[j]])*a**2) + (vg[j]*(nug[j]+2))/nug[j]	# => S = (mode*(df+2))/df         
        vg[j] <- scg/rchisq(n=1, df=df, ncp=0)    
        vgs[[j]][i,] <- vg[[j]]
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
        if(i>nburn) mus[,t] <- mus[,t] + mu[,t]
      }  
    }
    if(!is.null(X)) {
      for (t in 1:nt) {
        yadj <- y[,t]-rowSums(as.matrix(g[[t]])) 
        mu[,t] <- X%*%solve(t(X)%*%X)%*%t(X)%*%yadj
        if(i>nburn) mus[,t] <- mus[,t] + mu[,t]
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
    
    if(i>nburn) {
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
      if(i>nburn) gm[[t]][,j] <- gm[[t]][,j]/(nit-nburn) 
    } 
  }
  mus <- mus/(nit-nburn)
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
  
  # nt = (sqrt(1 + 8 * ncol(ves)) - 1) / 2
  # mat <- matrix(0,nt,nt)
  # upperindex <- upper.tri(mat, diag=TRUE)
  # lowerindex <- lower.tri(mat, diag=TRUE)
  # fit$res <- apply(fit$ves,1,function(x) {
  #   mat[upperindex] <- x
  #   mat[lowerindex] <- x
  #   rgmat <- cov2cor(mat)
  #   rgmat[upper.tri(mat, diag=FALSE)]
  # })
  # fit$rgs <- apply(fit$vgs,1,function(x) {
  #   mat[upperindex] <- x
  #   mat[lowerindex] <- x
  #   rgmat <- cov2cor(mat)
  #   rgmat[upper.tri(mat, diag=FALSE)]
  # })
  
  
  gset <- NULL
  for ( i in 1:nset) {                        
    gset[[i]] <- matrix(0,nrow=length(idsG),ncol=nt)         
    rownames(gset[[i]]) <- idsG 
    colnames(gset[[i]]) <- tnames
  }
  names(gset) <- setnames
  
  for ( set in 1:nset ){
    for (t1 in 1:nt) {
      isDups <- duplicated(idsT)
      ghat <- gm[[t1]][!isDups,set]
      names(ghat) <- idsT[!isDups]
      gset[[set]][names(ghat),t1] <- ghat
      #gset[[set]][idsT,t1] <- gm[[t1]][,set]
      if(length(idsV)>0) {
        for (t2 in 1:nt) {
          gset[[set]][idsV,t1] <- gset[[set]][idsV,t1] + (GRMlist[[set]][idsV,idsT]*vgm[t1,t2])%*%gm[[t2]][,set]
        }
      }
    }
  }
  fit$gset <- gset
  fit$g <- Reduce(`+`, fit$gset)
  fit$mu <- colMeans(mus) 
  return(fit) # return posterior samples of sigma 
}


# X <- matrix(rnorm(100000),nrow=1000)
# set <- sample(1:ncol(X),5)
# g <- rowSums(X[,set])
# e <- rnorm(nrow(X),mean=0,sd=1)
# y <- g + e
# y <- y - mean(y)
# XX <- crossprod(X)
# Xy <- crossprod(X,y)
# 
# tol <- 0.0001
# eg <- eigen(XX)                    
# ev <- eg$values
# U <- eg$vectors[,ev>tol]            
# D <- eg$values[ev>tol]              
# 
# nit <- 500
# vbs <- ves <- rep(0,nit)
# b <- r <- rm <- rep(0,ncol(XX))
# nt <- 1
# vb <- vb_prior <- 0.0001
# vb <- 25
# nub <- 4
# ve <- 10000
# for ( i in 1:nit ) {                   
#   radj <- Xy-r
#   rhs <- crossprod(U,radj) 
#   for (k in 1:nrow(rhs)) {
#     iC <- solve( diag(1,nt) + (ve/vb)/D[k])      
#     bhat <- iC%*%rhs[k]
#     b[k] <- MASS::mvrnorm(n=1,mu=bhat,Sigma=iC%*%ve) 
#   }
#   r <- U%*%b 
#   rm <- rm + r/nit  
#   
#   # Sample variance components
#   df <- length(b) + nub
#   
#   # inverse chisquare
#   scb<- sum((1/D)*b**2) + (vb_prior*(nub+2))/nub	# => S = (mode*(df+2))/df         
#   vb <- scb/rchisq(n=1, df=df, ncp=0)    
#   vbs[i] <- vb
#   ve <- var(r)
#   ves[i] <- ve
# }
# layout(matrix(1:3,ncol=3))
# plot(rm)
# plot(vbs)
# plot(ves)
# gmap1 <- function(Glist=NULL, stat=NULL, sets=NULL, models=NULL,
#                   rsids=NULL, ids=NULL, mask=NULL, lambda=NULL,  
#                   vb=NULL, vg=NULL, ve=NULL, pi=0.001, h2=0.5, 
#                   nub=4, nug=4, nue=4, 
#                   ssb_prior=NULL, ssg_prior=NULL, sse_prior=NULL,
#                   vb_prior=NULL, vg_prior=NULL, ve_prior=NULL,
#                   updateB=TRUE, updateG=TRUE, updateE=TRUE, updatePi=TRUE,
#                   formatLD="dense", checkLD=FALSE, shrinkLD=FALSE, shrinkCor=FALSE, pruneLD=FALSE, 
#                   checkConvergence=FALSE, critVe=3, critVg=3, critVb=3, critPi=3, 
#                   critB=3, critB1=0.5, critB2=3, 
#                   verbose=FALSE, eigen_threshold=0.995, cs_threshold=0.9, cs_r2=0.5,
#                   nit=1000, nburn=100, nthin=1, output="summary",
#                   method="bayesR", algorithm="mcmc-eigen", seed=10) {
#   
#   
#   # Check methods and parameter settings
#   methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
#   method <- match(method, methods) - 1
#   if( !sum(method%in%c(0:5))== 1 ) stop("method argument specified not valid")
#   algorithms <- c("mcmc","em-mcmc", "mcmc-eigen")
#   algorithm <- match(algorithm, algorithms)
#   if(is.na(algorithm)) stop("algorithm argument specified not valid")
#   
#   # check this again
#   if(is.data.frame(stat)) {
#     if( any(sapply(stat[,-c(1:5)],function(x){any(!is.finite(x))}))) stop("Some elements in stat not finite")
#     if( any(sapply(stat[,-c(1:5)],function(x){any(is.na(x))}))) stop("Some elements in stat NA")
#     nt <- 1
#     m <- sum(stat$rsids%in%unlist(Glist$rsids))
#     if(!is.null(Glist$rsidsLD)) m <- sum(stat$rsids%in%unlist(Glist$rsidsLD))
#   }
#   
#   # Prepare summary statistics
#   if(nt==1) {
#     yy <- median((stat$b^2 + (stat$n-2)*stat$seb^2)*stat$n)
#     n <- median(stat$n)
#   }
#   if(is.null(stat[["ww"]])) stat$ww <- (yy/n)/(stat$seb^2 + stat$b^2/stat$n)
#   if(is.null(stat[["wy"]])) stat$wy <- stat$b*stat$ww
#   if(nt>1) {
#     yy <- (stat$b^2 + (stat$n-2)*stat$seb^2)*stat$ww
#     yy <- apply(yy,2,median)
#     n <- apply(stat$n,2,median)
#   }
#   
#   # Prepare input
#   b <- matrix(0, nrow=length(stat$rsids), ncol=nt)
#   if(is.null(mask)) mask <- matrix(FALSE, nrow=length(rsids), ncol=nt)
#   
#   vy <- yy/(n-1)
#   if(is.null(ve)) ve <- vy*(1-h2)
#   if(is.null(vg)) vg <- vy*h2
#   mc <- min(c(5000,m))
#   if(method>=4 && is.null(vb)) vb <- vg/(mc*pi)
#   if(method>=4 && is.null(ssb_prior))  ssb_prior <- ((nub-2.0)/nub)*(vg/(mc*pi))
#   
#   if(!is.null(sets) && algorithm==3)  { 
#     
#     sets <- mapSets(sets=sets, rsids=stat$rsids, index=FALSE)
#     if(any(sapply(sets,function(x){any(is.na(x))}))) stop("NAs in sets detected - please remove these")
#     
#     chr <- as.numeric(unlist(Glist$chr))
#     chrSets <- sapply(mapSets(sets = sets, Glist = Glist, index = TRUE), function(x) unique(chr[x]))
#     if (length(Glist$bedfiles) == 1) chrSets <- setNames(rep(1, length(chrSets)), names(chrSets))
#     lsets <- sapply(chrSets,length)
#     sets <- sets[lsets==1]
#     if(any(lsets>1)) stop(paste("Following marker sets mapped to multiple chromosome:",paste(which(lsets>1),collapse=",")))
#     if(any(lsets==0)) stop(paste("Following marker sets not mapped to any chromosome:",paste(which(lsets==0),collapse=",")))
#     
#     # Prepare output
#     bm <- dm <- vector(mode="list",length=length(sets))
#     ves <- vgs <- vbs <- pis <- conv <- vector(mode="list",length=length(sets))
#     bs <- ds <- prob <- vector(mode="list",length=length(sets))
#     pim <- vector(mode="list",length=length(sets))
#     logcpo <- rep(0,length(sets))
#     fdr <- csets <- vector(mode="list",length=length(sets))
#     names(bm) <- names(dm) <- names(pim) <- names(sets)     
#     names(ves) <- names(vgs) <- names(pis) <- names(vbs) <- names(conv) <- names(sets)     
#     names(bs)  <- names(ds) <- names(prob) <- names(sets)
#     names(logcpo) <- names(fdr) <- names(csets) <- names(sets)
#     attempts <- rep(1, length=length(sets))
#     
#     
#     if(is.null(ids)) ids <- Glist$idsLD
#     if(is.null(ids)) ids <- Glist$ids
#     
#     # Compute phenotypic 
#     vy <- median(2*stat$eaf*(1-stat$eaf)*(stat$n*stat$seb^2 + stat$b^2))
#     
#     # BLR model for each set
#     for (i in 1:length(sets)) {
#       
#       chr <- chrSets[[i]]
#       rsids <- sets[[i]]
#       rws <- match(rsids,stat$rsids)
#       message(paste("Processing region:",i))
#       
#       pos <- getPos(Glist=Glist, chr=chr, rsids=rsids)
#       message(paste("Region size in Mb:",round((max(pos)-min(pos))/1000000,2)))
#       if(!is.null(Glist$map)) map <- getMap(Glist=Glist, chr=chr, rsids=rsids)
#       if(!is.null(Glist$map)) message(paste("Region size in cM:",round(max(map)-min(map),2)))
#       
#       # Prepare input
#       W <- getG(Glist=Glist, chr=chr, rsids=rsids, ids=ids, scale=TRUE)
#       B <- crossprod(scale(W))/(nrow(W)-1)
#       
#       if(shrinkLD) B <- corpcor::cor.shrink(W)
#       
#       eig <- eigen(B, symmetric=TRUE)
#       
#       for (j in 1:length(eigen_threshold)) {
#         
#         keep <- cumsum(eig$values)/sum(eig$values) < eigen_threshold[j]
#         
#         z <- t(eig$vectors[,keep]) %*% stat[rws, "b"]
#         
#         scaleb <- sqrt(1/(stat[rws, "n"]*stat[rws, "seb"]+stat[rws, "b"]^2))
#         z <- t(eig$vectors[,keep]) %*% (stat[rws, "b"]*scaleb)
#         
#         # Scale each element by the inverse square root of the corresponding eigenvalue
#         w <- z / sqrt(eig$values[keep])
#         
#         Q <- diag(sqrt(eig$values[keep]))%*%t(eig$vectors[,keep])
#         
#         colnames(Q) <- colnames(B)
#         
#         
#         LD <- NULL
#         LDvalues <- as.list(as.data.frame(Q))
#         LDindices <- lapply(1:ncol(Q),function(x) { (1:nrow(Q))-1 } )
#         rsids <- colnames(Q)
#         names(LDvalues) <- rsids
#         names(LDindices) <- rsids
#         
#         n <- mean(stat[rws,"n"])
#         xx <- stat[rws,"n"]
#         m <- ncol(Q)
#         
#         b <- rep(0, m)
#         
#         pi <- c(0.992,0.005,0.003,0.001)
#         gamma <- c(0,0.01,0.1,1)
#         
#         ve <- vy*(1-h2)
#         vg <- vy*h2
#         vb <- vg/(m*sum(pi*gamma))
#         
#         if(is.null(ssb_prior)) ssb_prior <-  ((nub-2.0)/nub)*(vg/(m*sum(pi*gamma)))
#         ssg_prior <-  ((nug-2.0)/nug)*vg
#         sse_prior <- ((nue-2.0)/nue)*ve
#         
#         
#         lambda <- rep(ve/vb,m)
#         mask <- rep(FALSE, m)
#         
#         fit <- .Call("_qgg_sbayes_reg_eigen",
#                      wy=w,
#                      ww=xx,
#                      LDvalues=LDvalues,
#                      LDindices=LDindices,
#                      b = b,
#                      lambda = lambda,
#                      mask=mask,
#                      pi = pi,
#                      gamma = gamma,
#                      vb = vb,
#                      vg = vg,
#                      ve = ve,
#                      ssb_prior=ssb_prior,
#                      ssg_prior=ssg_prior,
#                      sse_prior=sse_prior,
#                      nub=nub,
#                      nug=nug,
#                      nue=nue,
#                      updateB = updateB,
#                      updateE = updateE,
#                      updatePi = updatePi,
#                      updateG = updateG,
#                      n=n,
#                      nit=nit,
#                      nburn=nburn,
#                      nthin=nthin,
#                      method=as.integer(method),
#                      algo=as.integer(algorithm),
#                      seed=seed)
#         names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pis","pim","r","b","param","bs","ds","prob")
#         fit$bm <- fit$bm/scaleb
#         names(fit$bm) <- names(fit$dm) <- names(fit$b) <- names(LDvalues)
#         fit$bs <- matrix(fit$bs,nrow=length(fit$bm))
#         fit$ds <- matrix(fit$ds,nrow=length(fit$bm))
#         fit$prob <- matrix(fit$prob,nrow=length(fit$bm))
#         rownames(fit$bs) <- rownames(fit$ds) <- rownames(fit$prob) <- names(LDvalues)
#         colnames(fit$bs) <- colnames(fit$ds) <- colnames(fit$prob) <- 1:(nit+nburn)
#         # Re-scale betas
#         for (k in 1:nrow(fit$bs)) {
#           fit$bs[k,] <- fit$bs[k,]/scaleb[k]
#         }
#         
#         # Check convergence            
#         critve <- critvg <- critvb <- critpi <- critb <- FALSE
#         if(!updateE) critve <- TRUE
#         if(!updateG) critvg <- TRUE
#         if(!updateB) critvb <- TRUE
#         if(!updatePi) critpi <- TRUE
#         zve <- coda::geweke.diag(fit$ves[nburn:(nburn+nit)])$z
#         zvg <- coda::geweke.diag(fit$vgs[nburn:(nburn+nit)])$z
#         zvb <- coda::geweke.diag(fit$vbs[nburn:(nburn+nit)])$z
#         zpi <- coda::geweke.diag(fit$pis[nburn:(nburn+nit)])$z
#         zb <- coda::geweke.diag(apply(fit$bs[,nburn:(nburn+nit)],2,var))$z
#         if(!is.na(zve)) critve <- abs(zve)<critVe
#         if(!is.na(zvg)) critvg <- abs(zvg)<critVg
#         if(!is.na(zvb)) critvb <- abs(zvb)<critVb
#         if(!is.na(zpi)) critpi <- abs(zpi)<critPi
#         if(!is.na(zb)) critb <- abs(zb)<critB
#         
#         critb1 <- critb2 <- FALSE
#         
#         brws <- fit$dm>0
#         
#         # Check divergence        
#         if(sum(brws)>1 && all(c(critve, critvg, critvb, critpi, critb))) {
#           
#           tstat <- fit$bs[brws,]/stat[rws, "b"][brws]
#           pdiv <- apply(tstat, 1, function(x) {
#             sum(x[nburn:length(x)] > -critB1 & x[nburn:length(x)] <= 1+critB1)
#           })
#           pdiv <- pdiv/length(nburn:ncol(tstat))
#           pdiv <- pdiv[is.finite(pdiv)]
#           if(any(pdiv<0.95)) plot(pdiv)
#           critb1 <- !any(pdiv<0.95)    # FALSE if any pdiv is less than 0.95
#           if(!critb1) message(paste("Convergence not reached for critB1 "))
#           critb <- critb1
#         }
#         
#         # Check mismatch
#         if(checkLD && sum(brws)>1 && !all(c(critve, critvg, critvb, critpi, critb))) {
#           # Identify mismatch between LD and summary statistics 
#           bout <- checkb(B=B[brws,brws],
#                          b=stat[rws, "b"][brws],
#                          seb=stat[rws, "seb"][brws],
#                          critB=critB2, verbose=verbose)
#           critb2 <- !any(bout$outliers)    # FALSE if there are any outliers
#           if(critb2) message(paste("Convergence not reached for critB2 "))
#           critb <- critb2
#         }
#         
#         converged <- critve & critvg & critvb & critpi & critb
#         
#         # Make plots to monitor convergence
#         if(verbose) {
#           layout(matrix(1:4,ncol=2))
#           pipsets <- splitWithOverlap(1:length(rsids),100,99)
#           pip <- fit$dm
#           plot(pip, ylim=c(0,max(pip)), ylab="PIP",xlab="Position", frame.plot=FALSE)
#           plot(-log10(stat[rws,"p"]), ylab="-log10(P)",xlab="Position", frame.plot=FALSE)
#           hist(fit$ves, main="Ve", xlab="")
#           plot(y=fit$bm, x=stat[rws,"b"], ylab="Adjusted",xlab="Marginal", frame.plot=FALSE)
#           abline(h=0,v=0, lwd=2, col=2, lty=2)
#         }
#         attempts[i] <- j      
#         if(!converged) {
#           message(paste("Convergence not reached using eigen_threshold:",eigen_threshold[j]))
#           criteria_names <- c("Variance of errors (critve)", 
#                               "Genetic variance (critvg)", 
#                               "Marker variance (critvb)", 
#                               "Inclusion probability (critpi)", 
#                               "Posterior mean (critb)")
#           criteria_status <- c(critve, critvg, critvb, critpi, critb)
#           
#           message("Convergence criteria:")
#           for (k in seq_along(criteria_names)) {
#             message(sprintf("  %s: %s", criteria_names[k], ifelse(criteria_status[k], "Met", "Not Met")))
#           }
#         }
#         # Exit outer loop if convergence is reached
#         if (converged) {
#           if(verbose) message(paste("Convergence reached using eigen_threshold:",eigen_threshold[j]))
#           break
#         }
#         
#       }
#       
#       cutoffs <- seq(0.01, 0.99, by = 0.01)  # Generate 1:99 as fractions
#       cutoff_indices <- lapply(cutoffs, function(cutoff) fit$dm > cutoff)
#       bfdrs <- sapply(cutoff_indices, function(rws) {
#         if (any(rws)) {
#           fdrs <- rowMeans(1 - fit$prob[rws, , drop = FALSE], na.rm = TRUE)
#           c(mean = mean(fdrs, na.rm = TRUE), quantile(fdrs, c(0.025, 0.975), na.rm = TRUE))
#         } else {
#           c(NA, NA, NA)
#         }
#       })
#       bfdrs <- t(bfdrs)
#       rownames(bfdrs) <- round(cutoffs, 2)
#       
#       # Save results
#       bm[[i]] <- fit$bm
#       dm[[i]] <- fit$dm
#       pim[[i]] <- fit$pim
#       ves[[i]] <- fit$ves
#       vbs[[i]] <- fit$vbs
#       vgs[[i]] <- fit$vgs
#       pis[[i]] <- fit$pis
#       conv[[i]] <- c(zve,zvg,zvb,zpi,zb) 
#       if(output=="full") {
#         bs[[i]] <- fit$bs
#         ds[[i]] <- fit$ds
#         prob[[i]] <- fit$prob
#       }
#       fdr[[i]] <- bfdrs
#       logcpo[i] <- fit$param[4]
#       if(sum(fit$dm)>cs_threshold) csets[[i]] <- crs(prob=fit$dm, B=B, threshold=cs_threshold, r2=cs_r2)
#       names(bm[[i]]) <- names(dm[[i]]) <- rsids
#     }
#     fit <- NULL
#     fit$bm <- bm
#     fit$dm <- dm
#     fit$pim <- pim
#     fit$ves <- ves
#     fit$vbs <- vbs
#     fit$vgs <- vgs
#     fit$pis <- pis
#     if(output=="full") {
#       fit$bs <- bs
#       fit$ds <- ds
#       fit$prob <- prob
#     }
#     fit$fdr <- fdr
#     fit$logcpo <- logcpo
#     fit$cs <- csets
#   }  
#   
#   pip <- sapply(fit$dm,sum)
#   minb <- sapply(fit$bm,min)
#   maxb <- sapply(fit$bm,max)
#   m <- sapply(fit$bm,length)
#   
#   bm <- unlist(unname(fit$bm))
#   dm <- unlist(unname(fit$dm))
#   marker <- data.frame(rsids=unlist(Glist$rsids),
#                        chr=unlist(Glist$chr), pos=unlist(Glist$pos), 
#                        ea=unlist(Glist$a1), nea=unlist(Glist$a2),
#                        eaf=unlist(Glist$af),stringsAsFactors = FALSE)
#   marker <- marker[marker$rsids%in%names(bm),]
#   fit$stat <- data.frame(marker,bm=bm[marker$rsids],
#                          dm=dm[marker$rsids], stringsAsFactors = FALSE)
#   fit$stat$vm <- 2*(1-fit$stat$eaf)*fit$stat$eaf*fit$stat$bm^2
#   fit$method <- methods[method+1]
#   fit$mask <- mask
#   ve <- sapply(fit$ves,function(x){mean(x[nburn:length(x)])})
#   vg <- sapply(fit$vgs,function(x){mean(x[nburn:length(x)])})
#   vb <- sapply(fit$vbs,function(x){mean(x[nburn:length(x)])})
#   pi <- sapply(fit$pis,function(x){mean(x[nburn:length(x)])})
#   ve_ci <- t(sapply(fit$ves,function(x){quantile(x[nburn:length(x)], c(0.025,0.975))}))
#   vg_ci <- t(sapply(fit$vgs,function(x){quantile(x[nburn:length(x)], c(0.025,0.975))}))
#   vb_ci <- t(sapply(fit$vbs,function(x){quantile(x[nburn:length(x)], c(0.025,0.975))}))
#   pi_ci <- t(sapply(fit$pis,function(x){quantile(x[nburn:length(x)], c(0.025,0.975))}))
#   
#   fit$ci <- list(ve=cbind(mean=ve,ve_ci),
#                  vg=cbind(mean=vg,vg_ci), 
#                  vb=cbind(mean=vb,vb_ci), 
#                  pi=cbind(mean=pi,pi_ci))  
#   
#   if(!is.null(Glist$map)) map <- unlist(Glist$map)
#   pos <- unlist(Glist$pos)
#   sets <- lapply(fit$bm,names)
#   setsindex <- mapSets(sets=sets, rsids=unlist(Glist$rsids))
#   if(!is.null(Glist$map)) cm <- sapply(setsindex, function(x){ max(map[x])-min(map[x]) })
#   mb <- sapply(setsindex, function(x){ (max(pos[x])-min(pos[x]))/1000000 })
#   minmb <- sapply(setsindex, function(x){ min(pos[x]) })
#   maxmb <- sapply(setsindex, function(x){ max(pos[x]) })
#   
#   chr <- unlist(Glist$chr)
#   chr <- sapply(setsindex,function(x){as.numeric(unique(chr[x]))})
#   
#   b <- stat[fit$stat$rsids,"b"]
#   
#   conv <- t(as.data.frame(conv))
#   colnames(conv) <- c("zve","zvg","zvb","zpi","zb")
#   fit$conv <- data.frame(conv,ntrials=attempts, cutoff=eigen_threshold[attempts])
#   if(is.null(Glist$map)) fit$post <- data.frame(ve=ve,vg=vg, vb=vb, pi=pi, pip=pip, minb=minb, maxb=maxb, m=m, mb=mb, chr=chr, minmb=minmb, maxmb=maxmb)  
#   if(!is.null(Glist$map)) fit$post <- data.frame(ve=ve,vg=vg, vb=vb, pi=pi, pip=pip, minb=minb, maxb=maxb, m=m, mb=mb, cm=cm, chr=chr, minmb=minmb, maxmb=maxmb)  
#   rownames(fit$conv) <- rownames(fit$post) <- names(sets) 
#   
#   fit$ve <- mean(ve)
#   fit$vg <- sum(vg)
#   fit$b <- b
#   return(fit)
# }
# 
# 
# 
# gmap0 <- function(y=NULL, X=NULL, W=NULL, stat=NULL, trait=NULL, sets=NULL, fit=NULL, Glist=NULL,
#                   chr=NULL, rsids=NULL, ids=NULL, b=NULL, bm=NULL, seb=NULL, mask=NULL, LD=NULL, n=NULL,
#                   vg=NULL, vb=NULL, ve=NULL, ssg_prior=NULL, ssb_prior=NULL, sse_prior=NULL,
#                   lambda=NULL, scaleY=TRUE, shrinkLD=FALSE, shrinkCor=FALSE, formatLD="dense", pruneLD=TRUE, 
#                   r2=0.05, checkLD=TRUE,
#                   h2=NULL, pi=0.001, updateB=TRUE, updateG=TRUE, updateE=TRUE, updatePi=TRUE,
#                   adjustE=TRUE, models=NULL,
#                   checkConvergence=FALSE, critVe=3, critVg=5, critVb=5, critPi=3, ntrial=1,
#                   nug=4, nub=4, nue=4, verbose=FALSE,msize=100,threshold=NULL,
#                   ve_prior=NULL, vg_prior=NULL,tol=0.001,
#                   nit=100, nburn=50, nit_local=NULL,nit_global=NULL,
#                   method="bayesC", algorithm="mcmc") {
#   
#   
#   # Check methods and parameter settings
#   methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
#   method <- match(method, methods) - 1
#   if( !sum(method%in%c(0:5))== 1 ) stop("method argument specified not valid")
#   algorithms <- c("mcmc","em-mcmc")
#   algorithm <- match(algorithm, algorithms)
#   if(is.na(algorithm)) stop("algorithm argument specified not valid")
#   if(shrinkLD) {
#     if(is.null(Glist$map)) {
#       warning("No map information in Glist - LD matrix shrinkage turned off")
#       shrinkLD <- FALSE
#     }
#   } 
#   
#   # check this again
#   if(is.data.frame(stat)) {
#     if( any(sapply(stat[,-c(1:5)],function(x){any(!is.finite(x))}))) stop("Some elements in stat not finite")
#     if( any(sapply(stat[,-c(1:5)],function(x){any(is.na(x))}))) stop("Some elements in stat NA")
#     nt <- 1
#     rsids <- stat$rsids
#     m <- sum(rsids%in%unlist(Glist$rsids))
#     if(!is.null(Glist$rsidsLD)) m <- sum(rsids%in%unlist(Glist$rsidsLD))
#     stat$b <- as.matrix(stat$b)
#     stat$seb <- as.matrix(stat$seb)
#     stat$n <- as.matrix(stat$n)
#     stat$p <- as.matrix(stat$p)
#     rownames(stat$b) <- rownames(stat$seb) <- rsids
#     rownames(stat$n) <- rownames(stat$p) <- rsids
#     if(!is.null(stat[["ww"]])) {
#       stat$ww <- as.matrix(stat$ww)
#       rownames(stat$ww) <- rsids
#     }
#     if(!is.null(stat[["wy"]])) {
#       stat$wy <- as.matrix(stat$wy)
#       rownames(stat$wy) <- rsids
#     }
#   }
#   if(!is.data.frame(stat) && is.list(stat)) {
#     nt <- ncol(stat$b)
#     rsids <- rownames(stat$b)
#     m <- sum(rsids%in%unlist(Glist$rsids))
#     if(!is.null(Glist$rsidsLD)) m <- sum(rsids%in%unlist(Glist$rsidsLD))
#   }
#   
#   
#   # Prepare summary statistics
#   if(nt==1) {
#     yy <- median((stat$b^2 + (stat$n-2)*stat$seb^2)*stat$n)
#     n <- median(stat$n)
#   }
#   if(is.null(stat[["ww"]])) stat$ww <- (yy/n)/(stat$seb^2 + stat$b^2/stat$n)
#   #if(is.null(stat[["ww"]])) stat$ww <- 1/(stat$seb^2 + stat$b^2/stat$n)
#   if(is.null(stat[["wy"]])) stat$wy <- stat$b*stat$ww
#   # if(nt==1) {
#   #   yy <- median((stat$b^2 + (stat$n-2)*stat$seb^2)*stat$ww)
#   #   n <- median(stat$n)
#   # }
#   if(nt>1) {
#     yy <- (stat$b^2 + (stat$n-2)*stat$seb^2)*stat$ww
#     yy <- apply(yy,2,median)
#     n <- apply(stat$n,2,median)
#   }
#   
#   # Prepare input
#   b <- matrix(0, nrow=length(rsids), ncol=nt)
#   if(is.null(mask)) mask <- matrix(FALSE, nrow=length(rsids), ncol=nt)
#   rownames(b) <- rownames(mask) <- rsids
#   
#   vy <- yy/(n-1)
#   if(is.null(pi)) pi <- 0.001
#   if(is.null(h2)) h2 <- 0.5
#   if(is.null(ve)) ve <- vy*(1-h2)
#   if(is.null(vg)) vg <- vy*h2
#   mc <- min(c(5000,m))
#   if(method>=4 && is.null(vb)) vb <- vg/(mc*pi)
#   if(method>=4 && is.null(ssb_prior))  ssb_prior <- ((nub-2.0)/nub)*(vg/(mc*pi))
#   
#   if(is.null(trait)) trait <- 1
#   message(paste("Processing trait:",trait))
#   
#   
#   if(!is.null(sets))  { 
#     
#     sets <- mapSets(sets=sets, rsids=stat$rsids, index=FALSE)
#     if(any(sapply(sets,function(x){any(is.na(x))}))) stop("NAs in sets detected - please remove these")
#     
#     chr <- unlist(Glist$chr)
#     chrSets <- mapSets(sets=sets, Glist=Glist, index=TRUE)
#     chrSets <- sapply(chrSets,function(x){as.numeric(unique(chr[x]))})
#     lsets <- sapply(chrSets,length)
#     sets <- sets[lsets==1]
#     if(any(lsets>1)) {
#       warning(paste("Marker sets mapped to multiple chromosome:",paste(which(lsets>1),collapse=",")))
#     }
#     if(any(lsets==0)) {
#       warning(paste("Marker sets mapped to multiple chromosome:",paste(which(lsets==0),collapse=",")))
#     }
#     
#     # Prepare output
#     bm <- dm <- vector(mode="list",length=length(sets))
#     ves <- vgs <- vbs <- pis <- bs <- ds <- vector(mode="list",length=length(sets))
#     pim <- vector(mode="list",length=length(sets))
#     names(bm) <- names(dm) <- names(ves) <- names(vgs) <- names(pis) <- names(bs)  <- names(ds) <- names(sets)     
#     names(pim) <- names(bs)  <- names(ds) <- names(sets)     
#     attempts <- rep(0, length=length(sets))
#     
#     
#     if(is.null(ids)) ids <- Glist$idsLD
#     if(is.null(ids)) ids <- Glist$ids
#     
#     #message(paste("Processing chromosome:",chr))
#     if(formatLD=="sparse") {
#       sparseLD <- getSparseLD(Glist=Glist,chr=chr)
#     }
#     # BLR model for each set
#     for (i in 1:length(sets)) {
#       
#       chr <- chrSets[[i]]
#       rsids <- sets[[i]]
#       message(paste("Processing region:",i,"on chromosome:",chr))
#       
#       pos <- getPos(Glist=Glist, chr=chr, rsids=rsids)
#       message(paste("Region size in Mb:",round((max(pos)-min(pos))/1000000,2)))
#       if(!is.null(Glist$map)) map <- getMap(Glist=Glist, chr=chr, rsids=rsids)
#       if(!is.null(Glist$map)) message(paste("Region size in cM:",round(max(map)-min(map),2)))
#       
#       if(formatLD=="dense") {
#         W <- getG(Glist=Glist, chr=chr, rsids=rsids, ids=ids, scale=TRUE)
#         B <- crossprod(scale(W))/(length(ids)-1)
#         if(shrinkCor) B <- corpcor::cor.shrink(W)
#         if(shrinkLD) B <- adjustMapLD(LD = B, map=map)
#         LD <- NULL
#         #LD$values <- split(B, rep(1:ncol(B), each = nrow(B)))
#         LD$values <- as.list(as.data.frame(B))
#         LD$indices <- lapply(1:ncol(B),function(x) { (1:ncol(B))-1 } )
#         rsids <- colnames(B)
#         names(LD$values) <- rsids
#         names(LD$indices) <- rsids
#         msize_set <- length(rsids)
#       }
#       
#       
#       if(formatLD=="sparse") {
#         B <- regionLD(sparseLD = sparseLD, onebased=FALSE, rsids=rsids, format="dense")
#         if(shrinkCor) B <- corpcor::cor.shrink(W)
#         if(shrinkLD) B <- adjustMapLD(LD = B, map=map)
#         LD <- NULL
#         #LD$values <- split(B, rep(1:ncol(B), each = nrow(B)))
#         LD$values <- as.list(as.data.frame(B))
#         LD$indices <- lapply(1:ncol(B),function(x) { (1:ncol(B))-1 } )
#         rsids <- colnames(B)
#         names(LD$values) <- rsids
#         names(LD$indices) <- rsids
#         msize_set <- length(rsids)
#       }
#       
#       #ntrial <- 5
#       converged <- FALSE
#       
#       updateB_reg <- updateB
#       updatePi_reg <- updatePi
#       pi_reg <- pi
#       r2_reg <- r2
#       
#       for (trial in 1:ntrial) {
#         
#         if (!converged) {
#           
#           if(pruneLD) {
#             message("Adjust summary statistics using pruning")
#             pruned <- adjLDregion(LD=B, p=stat$p[rsids,trait], r2=r2_reg, thold=1) 
#             mask[pruned,trait] <- TRUE
#           }
#           
#           attempts[i] <- trial
#           
#           fit <- sbayes_region(yy=yy[trait],
#                                wy=stat$wy[rsids,trait],
#                                ww=stat$ww[rsids,trait],
#                                b=b[rsids,trait],
#                                mask=mask[rsids,trait],
#                                LDvalues=LD$values,
#                                LDindices=LD$indices,
#                                method=method,
#                                algorithm=algorithm,
#                                nit=nit,
#                                nburn=nburn,
#                                n=n[trait],
#                                m=msize_set,
#                                pi=pi_reg,
#                                nue=nue,
#                                nub=nub,
#                                ssb_prior=ssb_prior,
#                                updateB=updateB_reg,
#                                updateE=updateE,
#                                updatePi=updatePi_reg,
#                                updateG=updateG,
#                                adjustE=adjustE)
#           
#           # Check convergence            
#           critve <- critvg <- critvb <- critpi <- FALSE
#           if(!updateE) critve <- TRUE
#           if(!updateG) critvg <- TRUE
#           if(!updateB) critvb <- TRUE
#           if(!updatePi) critpi <- TRUE
#           zve <- coda::geweke.diag(fit$ves[nburn:length(fit$ves)])$z
#           zvg <- coda::geweke.diag(fit$vgs[nburn:length(fit$vgs)])$z
#           zvb <- coda::geweke.diag(fit$vbs[nburn:length(fit$vbs)])$z
#           zpi <- coda::geweke.diag(fit$pis[nburn:length(fit$pis)])$z
#           if(!is.na(zve)) critve <- abs(zve)<critVe
#           if(!is.na(zvg)) critvg <- abs(zvg)<critVg
#           if(!is.na(zvb)) critvb <- abs(zvb)<critVb
#           if(!is.na(zpi)) critpi <- abs(zpi)<critPi
#           
#           critb1 <- fit$dm>0.01 & fit$bm>0 & fit$bm>stat$b[rsids,trait]
#           critb2 <- fit$dm>0.01 & fit$bm<0 & fit$bm<stat$b[rsids,trait]
#           critb <- !any(critb1 | critb2)
#           converged <- critve & critvg & critvb & critpi & critb
#           
#           
#           if (!converged & checkConvergence) {
#             message("")
#             message(paste("Region not converged in attempt:",trial))
#             if(!critve) message(paste("Zve:",zve))
#             if(!critvg) message(paste("Zvg:",zvg))
#             if(!critvb) message(paste("Zvb:",zvb))
#             if(!critpi) message(paste("Zpi:",zpi))
#             W <- getG(Glist=Glist, chr=chr, rsids=rsids, ids=ids, scale=TRUE)
#             B <- crossprod(scale(W))/(length(ids)-1)
#             if(shrinkCor) B <- corpcor::cor.shrink(W)
#             if(shrinkLD) B <- adjustMapLD(LD = B, map=map)
#             if(checkLD) { 
#               message("Adjust summary statistics using imputation")
#               badj <- adjustB(b=stat$b[rsids,trait], LD = B, 
#                               msize=200, overlap=50, shrink=0.001, threshold=1e-8) 
#               # badj <- adjustB(b=stat$b[rsids,trait], LD = B, 
#               #                       msize=500, overlap=100, shrink=0.001, threshold=1e-8) 
#               # z <- (badj-stat$b[rsids,trait])/stat$seb[rsids,trait]
#               # outliers <- names(z[abs(z)>1.96])
#               #mask[outliers,trait] <- TRUE
#               #stat$b[outliers,trait] <- badj[abs(z)>1.96]
#               #stat$ww[outliers,trait] <- 1/(stat$seb[outliers,trait]^2 + stat$b[outliers,trait]^2/stat$n[outliers,trait])
#               #stat$wy[outliers,trait] <- stat$b[outliers,trait]*stat$ww[outliers,trait]
#               stat$b[rsids,trait] <- badj
#               stat$ww[rsids,trait] <- 1/(stat$seb[rsids,trait]^2 + stat$b[rsids,trait]^2/stat$n[rsids,trait])
#               stat$wy[rsids,trait] <- stat$b[rsids,trait]*stat$ww[rsids,trait]
#             }
#             if(pruneLD) {
#               #if(pruneLD) {
#               message("Adjust summary statistics using pruning")
#               pruned <- adjLDregion(LD=B, p=stat$p[rsids,trait], r2=r2, thold=1) 
#               mask[pruned,trait] <- TRUE
#             }
#             # if(trial==3) {
#             #   message("Set updateB and updatePi to FALSE")
#             #   updateB_reg <- FALSE 
#             #   updatePi_reg <- FALSE 
#             # }
#             if(trial>0) {
#               message("Decrease r2 by a factor 10")
#               #updateB_reg <- FALSE 
#               updatePi_reg <- FALSE
#               r2_reg <- r2_reg*0.1
#               #pi_reg <- pi_reg*0.1
#             }
#           }
#         }
#       }
#       
#       # Make plots to monitor convergence
#       if(verbose) {
#         layout(matrix(1:4,ncol=2))
#         pipsets <- splitWithOverlap(1:length(rsids),100,99)
#         pip <- fit$dm
#         plot(pip, ylim=c(0,max(pip)), ylab="PIP",xlab="Position", frame.plot=FALSE)
#         plot(-log10(stat$p[rsids,trait]), ylab="-log10(P)",xlab="Position", frame.plot=FALSE)
#         hist(fit$ves, main="Ve", xlab="")
#         plot(y=fit$bm, x=stat$b[rsids,trait], ylab="Adjusted",xlab="Marginal", frame.plot=FALSE)
#         abline(h=0,v=0, lwd=2, col=2, lty=2)
#       }
#       
#       # Save results
#       bm[[i]] <- fit$bm
#       dm[[i]] <- fit$dm
#       pim[[i]] <- fit$pim
#       ves[[i]] <- fit$ves
#       vbs[[i]] <- fit$vbs
#       vgs[[i]] <- fit$vgs
#       pis[[i]] <- fit$pis
#       selected <- NULL
#       if(!is.null(threshold)) selected <- fit$dm>=threshold
#       if(any(selected)) {
#         bs[[i]] <- matrix(fit$bs,nrow=length(rsids))
#         ds[[i]] <- matrix(fit$ds,nrow=length(rsids))
#         rownames(bs[[i]]) <- rownames(ds[[i]]) <- rsids
#         colnames(bs[[i]]) <- colnames(ds[[i]]) <- 1:(nit+nburn)
#         bs[[i]] <- bs[[i]][selected,]
#         ds[[i]] <- ds[[i]][selected,]
#       }
#       names(bm[[i]]) <- names(dm[[i]]) <- rsids
#     }
#     fit <- NULL
#     fit$bm <- bm
#     fit$dm <- dm
#     fit$pim <- pim
#     fit$ves <- ves
#     fit$vbs <- vbs
#     fit$vgs <- vgs
#     fit$pis <- pis
#     fit$attempts <- attempts
#     if(!is.null(threshold)) fit$bs <- bs
#     if(!is.null(threshold)) fit$ds <- ds
#     
#   }  
#   
#   if(is.null(sets))  { 
#     
#     
#     # Prepare output
#     bm <- dm <- vector(mode="list",length=22)
#     ves <- vgs <- vbs <- pis <- bs <- ds <- vector(mode="list",length=22)
#     pim <- attempts <- vector(mode="list",length=22)
#     
#     chromosomes <- 1:22
#     if(!is.null(chr)) chromosomes <- chr 
#     
#     if(is.null(ids)) ids <- Glist$idsLD
#     if(is.null(ids)) ids <- Glist$ids
#     
#     for (chr in chromosomes) {
#       
#       message(paste("Processing chromosome:",chr))
#       rsidsLD <- Glist$rsidsLD[[chr]]
#       rsidsLD <- rsidsLD[rsidsLD%in%rownames(b)]
#       sets <- split(rsidsLD, ceiling(seq_along(rsidsLD) / msize))
#       #if(is.null(ssb_prior)) {
#       #  h2 <- 0.5
#       #  pi <- 0.001
#       #  vy <- 1
#       #  vg <- h2*vy
#       #  nub <- 4
#       #  ww <- 1/(stat$seb^2 + stat$b/stat$n)
#       #  mx <- sum(ww/mean(stat$n))
#       #  ssb_prior <- vy*h2*(nub+2)/mx/pi
#       #}
#       
#       if(formatLD=="sparse") {
#         sparseLD <- getSparseLD(Glist=Glist,chr=chr)
#       }
#       
#       # BLR model for each set
#       for (i in 1:length(sets)) {
#         
#         message(paste("Processing region:",i))
#         rsids <- sets[[i]]
#         pos <- getPos(Glist=Glist, chr=chr, rsids=rsids)
#         message(paste("Region size in Mb:",round((max(pos)-min(pos))/1000000,2)))
#         if(!is.null(Glist$map)) map <- getMap(Glist=Glist, chr=chr, rsids=rsids)
#         if(!is.null(Glist$map)) message(paste("Region size in cM:",round(max(map)-min(map),2)))
#         
#         if(formatLD=="dense") {
#           W <- getG(Glist=Glist, chr=chr, rsids=rsids, ids=ids, scale=TRUE)
#           B <- crossprod(scale(W))/(length(ids)-1)
#           if(shrinkCor) B <- corpcor::cor.shrink(W)
#           if(shrinkLD) B <- adjustMapLD(LD = B, map=map)
#           LD <- NULL
#           #LD$values <- split(B, rep(1:ncol(B), each = nrow(B)))
#           LD$values <- as.list(as.data.frame(B))
#           LD$indices <- lapply(1:ncol(B),function(x) { (1:ncol(B))-1 } )
#           rsids <- colnames(B)
#           names(LD$values) <- rsids
#           names(LD$indices) <- rsids
#           msize_set <- length(rsids)
#         }
#         
#         
#         
#         if(formatLD=="sparse") {
#           B <- regionLD(sparseLD = sparseLD, onebased=FALSE, rsids=rsids, format="dense")
#           if(shrinkLD) B <- adjustMapLD(LD = B, map=map)
#           LD <- NULL
#           #LD$values <- split(B, rep(1:ncol(B), each = nrow(B)))
#           LD$values <- as.list(as.data.frame(B))
#           LD$indices <- lapply(1:ncol(B),function(x) { (1:ncol(B))-1 } )
#           rsids <- colnames(B)
#           names(LD$values) <- rsids
#           names(LD$indices) <- rsids
#           msize_set <- length(rsids)
#           #LD <- regionLD(sparseLD = sparseLD, onebased=FALSE, rsids=rsids, format="sparse")
#           #rsids <- LD$rsids
#           #msize <- length(rsids)
#         }
#         
#         #ntrial <- 5
#         converged <- FALSE
#         
#         updateB_reg <- updateB
#         updatePi_reg <- updatePi
#         pi_reg <- pi
#         
#         for (trial in 1:ntrial) {
#           
#           if (!converged) {
#             
#             attempts[[chr]][[i]] <- trial
#             
#             
#             fit <- sbayes_region(yy=yy[trait],
#                                  wy=stat$wy[rsids,trait],
#                                  ww=stat$ww[rsids,trait],
#                                  b=b[rsids,trait],
#                                  mask=mask[rsids,trait],
#                                  LDvalues=LD$values,
#                                  LDindices=LD$indices,
#                                  method=method,
#                                  algorithm=algorithm,
#                                  nit=nit,
#                                  nburn=nburn,
#                                  n=n[trait],
#                                  m=msize_set,
#                                  pi=pi,
#                                  nue=nue,
#                                  nub=nub,
#                                  ssb_prior=ssb_prior,
#                                  updateB=updateB_reg,
#                                  updateE=updateE,
#                                  updatePi=updatePi_reg,
#                                  updateG=updateG,
#                                  adjustE=adjustE)
#             #   }
#             
#             # Check convergence            
#             critve <- critvg <- critvb <- critpi <- FALSE
#             if(!updateE) critve <- TRUE
#             if(!updateG) critvg <- TRUE
#             if(!updateB) critvb <- TRUE
#             if(!updatePi) critpi <- TRUE
#             zve <- coda::geweke.diag(fit$ves[nburn:length(fit$ves)])$z
#             zvg <- coda::geweke.diag(fit$vgs[nburn:length(fit$vgs)])$z
#             zvb <- coda::geweke.diag(fit$vbs[nburn:length(fit$vbs)])$z
#             zpi <- coda::geweke.diag(fit$pis[nburn:length(fit$pis)])$z
#             if(!is.na(zve)) critve <- abs(zve)<critVe
#             if(!is.na(zvg)) critvg <- abs(zvg)<critVg
#             if(!is.na(zvb)) critvb <- abs(zvb)<critVb
#             if(!is.na(zpi)) critpi <- abs(zpi)<critPi
#             
#             critb1 <- fit$dm>0.01 & fit$bm>0 & fit$bm>stat$b[rsids,trait]
#             critb2 <- fit$dm>0.01 & fit$bm<0 & fit$bm<stat$b[rsids,trait]
#             critb <- !any(critb1 | critb2)
#             converged <- critve & critvg & critvb & critpi & critb
#             
#             if (!converged & checkConvergence) {
#               message("")
#               message(paste("Region not converged in attempt:",trial))
#               if(!critve) message(paste("Zve:",zve))
#               if(!critvg) message(paste("Zvg:",zvg))
#               if(!critvb) message(paste("Zvb:",zvb))
#               if(!critpi) message(paste("Zpi:",zpi))
#               W <- getG(Glist=Glist, chr=chr, rsids=rsids, ids=ids, scale=TRUE)
#               B <- crossprod(scale(W))/(length(ids)-1)
#               if(shrinkCor) B <- corpcor::cor.shrink(W)
#               if(shrinkLD) B <- adjustMapLD(LD = B, map=map)
#               if(checkLD) { 
#                 # message("Adjust summary statistics using imputation")
#                 # badj <- adjustB(b=stat$b[rsids,trait], LD = B, 
#                 #                       msize=500, overlap=100, shrink=0.001, threshold=1e-8) 
#                 # z <- (badj-stat$b[rsids,trait])/stat$seb[rsids,trait]
#                 # outliers <- names(z[abs(z)>1.96])
#                 #mask[outliers,trait] <- TRUE
#                 message("Adjust summary statistics using imputation")
#                 badj <- adjustB(b=stat$b[rsids,trait], LD = B, 
#                                 msize=200, overlap=50, shrink=0.001, threshold=1e-8) 
#                 #z <- (badj-stat$b[rsids,trait])/stat$seb[rsids,trait]
#                 #outliers <- names(z[abs(z)>1.96])
#                 #mask[outliers,trait] <- TRUE
#                 #stat$b[outliers,trait] <- badj[abs(z)>1.96]
#                 #stat$ww[outliers,trait] <- 1/(stat$seb[outliers,trait]^2 + stat$b[outliers,trait]^2/stat$n[outliers,trait])
#                 #stat$wy[outliers,trait] <- stat$b[outliers,trait]*stat$ww[outliers,trait]
#                 stat$b[rsids,trait] <- badj
#                 stat$ww[rsids,trait] <- 1/(stat$seb[rsids,trait]^2 + stat$b[rsids,trait]^2/stat$n[rsids,trait])
#                 stat$wy[rsids,trait] <- stat$b[rsids,trait]*stat$ww[rsids,trait]
#               }
#               #if(pruneLD) {
#               if(pruneLD) {
#                 message("Adjust summary statistics using pruning")
#                 pruned <- adjLDregion(LD=B, p=stat$p[rsids,trait], r2=r2, thold=1) 
#                 mask[pruned,trait] <- TRUE
#               }
#               # if(trial==3) {
#               #   message("Set updateB and updatePi to FALSE")
#               #   updateB_reg <- FALSE 
#               #   updatePi_reg <- FALSE 
#               # }
#               if(trial>1) {
#                 message("Decrease Pi by a factor 10")
#                 #updateB_reg <- FALSE 
#                 updatePi_reg <- FALSE
#                 pi_reg <- pi_reg*0.1
#               }
#             }
#           }
#         }
#         
#         # Make plots to monitor convergence
#         if(verbose) {
#           layout(matrix(1:4,ncol=2))
#           pipsets <- splitWithOverlap(1:length(rsids),100,99)
#           #pip <- sapply(pipsets,function(x){sum(fit$dm[x])})
#           pip <- fit$dm
#           plot(pip, ylim=c(0,max(pip)), ylab="PIP",xlab="Position", frame.plot=FALSE)
#           plot(-log10(stat$p[rsids,trait]), ylab="-log10(P)",xlab="Position", frame.plot=FALSE)
#           hist(fit$ves, main="Ve", xlab="")
#           plot(y=fit$bm, x=stat$b[rsids,trait], ylab="Adjusted",xlab="Marginal", frame.plot=FALSE)
#           abline(h=0,v=0, lwd=2, col=2, lty=2)
#         }
#         
#         # Save results
#         bm[[chr]][[i]] <- fit$bm
#         dm[[chr]][[i]] <- fit$dm
#         pim[[chr]][[i]] <- fit$pim
#         ves[[chr]][[i]] <- fit$ves
#         vbs[[chr]][[i]] <- fit$vbs
#         vgs[[chr]][[i]] <- fit$vgs
#         pis[[chr]][[i]] <- fit$pis
#         #bs[[chr]][[i]] <- matrix(fit$bs,nrow=length(rsids))
#         #ds[[chr]][[i]] <- matrix(fit$ds,nrow=length(rsids))
#         #rownames(bs[[chr]][[i]]) <- rownames(ds[[chr]][[i]]) <- rsids
#         #colnames(bs[[chr]][[i]]) <- colnames(ds[[chr]][[i]]) <- 1:nit
#         selected <- NULL
#         if(!is.null(threshold)) selected <- fit$dm>=threshold
#         if(any(selected)) {
#           bs[[chr]][[i]] <- matrix(fit$bs,nrow=length(rsids))
#           ds[[chr]][[i]] <- matrix(fit$ds,nrow=length(rsids))
#           rownames(bs[[chr]][[i]]) <- rownames(ds[[chr]][[i]]) <- rsids
#           colnames(bs[[chr]][[i]]) <- colnames(ds[[chr]][[i]]) <- 1:(nit+nburn)
#           bs[[chr]][[i]] <- bs[[chr]][[i]][selected,]
#           ds[[chr]][[i]] <- ds[[chr]][[i]][selected,]
#         }
#         names(bm[[chr]][[i]]) <- names(dm[[chr]][[i]]) <- rsids
#       }
#     }
#     
#     fit <- NULL
#     fit$bm <- unlist(bm, recursive=FALSE)
#     fit$dm <- unlist(dm, recursive=FALSE)
#     fit$pim <- unlist(pim, recursive=FALSE)
#     fit$ves <- unlist(ves, recursive=FALSE)
#     fit$vbs <- unlist(vbs, recursive=FALSE)
#     fit$vgs <- unlist(vgs, recursive=FALSE)
#     fit$pis <- unlist(pis, recursive=FALSE)
#     fit$attempts <- unlist(attempts, recursive=TRUE)
#     if(!is.null(threshold)) fit$bs <- unlist(bs, recursive=FALSE)
#     if(!is.null(threshold)) fit$ds <- unlist(ds, recursive=FALSE)
#   }
#   
#   
#   pip <- sapply(fit$dm,sum)
#   minb <- sapply(fit$bm,min)
#   maxb <- sapply(fit$bm,max)
#   m <- sapply(fit$bm,length)
#   
#   bm <- unlist(unname(fit$bm))
#   dm <- unlist(unname(fit$dm))
#   #selected <- dm>0
#   #bm <- bm[selected]
#   #dm <- dm[selected]
#   marker <- data.frame(rsids=unlist(Glist$rsids),
#                        chr=unlist(Glist$chr), pos=unlist(Glist$pos), 
#                        ea=unlist(Glist$a1), nea=unlist(Glist$a2),
#                        eaf=unlist(Glist$af),stringsAsFactors = FALSE)
#   marker <- marker[marker$rsids%in%names(bm),]
#   fit$stat <- data.frame(marker,bm=bm[marker$rsids],
#                          dm=dm[marker$rsids], stringsAsFactors = FALSE)
#   fit$stat$vm <- 2*(1-fit$stat$eaf)*fit$stat$eaf*fit$stat$bm^2
#   fit$method <- methods[method+1]
#   fit$mask <- mask
#   zve <- sapply(fit$ves,function(x){coda::geweke.diag(x[nburn:length(x)])$z})
#   zvg <- sapply(fit$vgs,function(x){coda::geweke.diag(x[nburn:length(x)])$z})
#   zvb <- sapply(fit$vbs,function(x){coda::geweke.diag(x[nburn:length(x)])$z})
#   zpi <- sapply(fit$pis,function(x){coda::geweke.diag(x[nburn:length(x)])$z})
#   ve <- sapply(fit$ves,function(x){mean(x[nburn:length(x)])})
#   vg <- sapply(fit$vgs,function(x){mean(x[nburn:length(x)])})
#   vb <- sapply(fit$vbs,function(x){mean(x[nburn:length(x)])})
#   pi <- sapply(fit$pim,function(x){1-x[1]})
#   
#   if(!is.null(Glist$map)) map <- unlist(Glist$map)
#   pos <- unlist(Glist$pos)
#   sets <- lapply(fit$bm,names)
#   setsindex <- mapSets(sets=sets, rsids=unlist(Glist$rsids))
#   if(!is.null(Glist$map)) cm <- sapply(setsindex, function(x){ max(map[x])-min(map[x]) })
#   mb <- sapply(setsindex, function(x){ (max(pos[x])-min(pos[x]))/1000000 })
#   minmb <- sapply(setsindex, function(x){ min(pos[x]) })
#   maxmb <- sapply(setsindex, function(x){ max(pos[x]) })
#   
#   chr <- unlist(Glist$chr)
#   chr <- sapply(setsindex,function(x){as.numeric(unique(chr[x]))})
#   
#   b <- stat[fit$stat$rsids,"b"]
#   
#   #fit$region <-  NULL
#   fit$conv <- data.frame(zve=zve,zvg=zvg, zvb=zvb, zpi=zpi)  
#   if(is.null(Glist$map)) fit$post <- data.frame(ve=ve,vg=vg, vb=vb, pi=pi, pip=pip, minb=minb, maxb=maxb, m=m, mb=mb, chr=chr, minmb=minmb, maxmb=maxmb)  
#   if(!is.null(Glist$map)) fit$post <- data.frame(ve=ve,vg=vg, vb=vb, pi=pi, pip=pip, minb=minb, maxb=maxb, m=m, mb=mb, cm=cm, chr=chr, minmb=minmb, maxmb=maxmb)  
#   rownames(fit$conv) <- rownames(fit$post) <- names(sets) 
#   fit$ve <- mean(ve)
#   fit$vg <- sum(vg)
#   fit$b <- b
#   return(fit)
# }
# 
# 
# 
# # Single trait fine-mapping BLR using summary statistics and sparse LD provided in Glist 
# sbayes_region <- function(yy=NULL, wy=NULL, ww=NULL, b=NULL, bm=NULL, mask=NULL, seb=NULL, 
#                           LDvalues=NULL,LDindices=NULL, n=NULL, m=NULL,
#                           vg=NULL, vb=NULL, ve=NULL, 
#                           ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
#                           h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, 
#                           updateG=NULL, adjustE=NULL, models=NULL,
#                           nub=NULL, nue=NULL, nit=NULL, nburn=NULL, method=NULL, algorithm=NULL, verbose=NULL) {
#   
#   if(is.null(m)) m <- length(LDvalues)
#   vy <- yy/(n-1)
#   if(is.null(pi)) pi <- 0.001
#   if(is.null(h2)) h2 <- 0.5
#   if(is.null(ve)) ve <- vy*(1-h2)
#   if(is.null(vg)) vg <- vy*h2
#   if(method<4 && is.null(vb)) vb <- vg/m
#   if(method>=4 && is.null(vb)) vb <- vg/(m*pi)
#   if(is.null(lambda)) lambda <- rep(ve/vb,m)
#   if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
#   if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/(m*pi))
#   if(is.null(sse_prior)) sse_prior <- ((nue-2.0)/nue)*ve
#   if(is.null(b)) b <- rep(0,m)
#   
#   pi <- c(1-pi,pi)
#   gamma <- c(0,1.0)
#   if(method==5) pi <- c(0.95,0.02,0.02,0.01)
#   if(method==5) gamma <- c(0,0.01,0.1,1.0)
#   if(is.null(algorithm)) algorithm <- 0
#   
#   fit <- .Call("_qgg_sbayes_reg",
#                wy=wy, 
#                ww=ww, 
#                LDvalues=LDvalues, 
#                LDindices=LDindices, 
#                b = b,
#                lambda = lambda,
#                mask = mask,
#                yy = yy,
#                pi = pi,
#                gamma = gamma,
#                vg = vg,
#                vb = vb,
#                ve = ve,
#                ssb_prior=ssb_prior,
#                sse_prior=sse_prior,
#                nub=nub,
#                nue=nue,
#                updateB = updateB,
#                updateE = updateE,
#                updatePi = updatePi,
#                updateG = updateG,
#                adjustE = adjustE,
#                n=n,
#                nit=nit,
#                nburn=nburn,
#                method=as.integer(method),
#                algo=as.integer(algorithm))
#   names(fit[[1]]) <- names(LDvalues)
#   names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pis","pim","r","b","param","bs","ds")
#   return(fit)
# }
# Single trait BLR using y and dense LD
#if( nt==1 && !is.null(y) &&  algorithm=="dense") {
# if( nt==1 && !is.null(y) &&  formatLD=="dense") {
#   
#   overlap <- 0
#   
#   if(is.null(Glist)) stop("Please provide Glist")
#   fit <- NULL
#   if(is.matrix(y)) ids <- rownames(y)
#   if(is.vector(y)) ids <- names(y)
#   rws <- match(ids,Glist$ids)
#   if(any(is.na(rws))) stop("some elements in names(y) does not match elements in Glist$ids ")       
#   n <- length(y)
#   
#   if(is.null(chr)) chromosomes <- 1:Glist$nchr
#   if(!is.null(chr)) chromosomes <- chr
#   
#   rsids <- unlist(Glist$rsidsLD)
#   cls <- lapply(Glist$rsids,function(x) { 
#     splitWithOverlap(na.omit(match(rsids,x)),msize,0)})
#   vblist <- lapply(sapply(cls,length),function(x) 
#   {vector(length=x, mode="numeric")})
#   velist <- lapply(sapply(cls,length),function(x) 
#   {vector(length=x, mode="numeric")})
#   pilist <- lapply(sapply(cls,length),function(x) 
#   {vector(length=x, mode="numeric")})
#   b <- lapply(Glist$mchr,function(x){rep(0,x)})
#   bm <- lapply(Glist$mchr,function(x){rep(0,x)})
#   dm <- lapply(Glist$mchr,function(x){rep(0,x)})
#   
#   if(is.null(nit_local)) nit_local <- nit
#   if(is.null(nit_global)) nit_global <- 1
#   
#   for (it in 1:nit_global) {
#     e <- y-mean(y)
#     yy <- sum(e**2)
#     for (chr in 1:length(Glist$nchr)) {
#       for (i in 1:length(cls[[chr]])) {
#         wy <- computeWy(y=e,Glist=Glist,chr=chr,cls=cls[[chr]][[i]])
#         WW <- computeWW(Glist=Glist, chr=chr, cls=cls[[chr]][[i]], rws=rws)
#         if(it>1) {
#           if(updateB) vb <- vblist[[chr]][i]
#           if(updateE) ve <- velist[[chr]][i]
#           if(updatePi) pi <- pilist[[chr]][i]
#         }
#         fitS <- computeB(wy=wy, yy=yy, WW=WW, n=n,
#                          b=b[[chr]][cls[[chr]][[i]]],
#                          ve=ve, vb=vb, pi=pi,
#                          nub=nub, nue=nue,
#                          updateB=updateB, updateE=updateE, updatePi=updatePi,
#                          nit=nit, nburn=nburn, method=method) 
#         b[[chr]][cls[[chr]][[i]]] <- fitS$b
#         bm[[chr]][cls[[chr]][[i]]] <- fitS$bm
#         dm[[chr]][cls[[chr]][[i]]] <- fitS$dm
#         vblist[[chr]][i] <- fitS$param[1]
#         velist[[chr]][i] <- fitS$param[2]
#         pilist[[chr]][i] <- fitS$param[3]
#         grs <- computeGRS(Glist = Glist, chr = chr, 
#                           cls = cls[[chr]][[i]], 
#                           b=bm[[chr]][cls[[chr]][[i]]])  
#         e <- e - grs[rws,]
#       }
#     }
#   }   
#   bm <- unlist(bm)
#   dm <- unlist(dm)
#   names(bm) <- names(dm) <- unlist(Glist$rsids)
#   rsids2rws <- match(rsids,unlist(Glist$rsids))
#   stat <- data.frame(rsids=rsids,
#                      chr=unlist(Glist$chr)[rsids2rws],
#                      pos=unlist(Glist$pos)[rsids2rws], 
#                      ea=unlist(Glist$a1)[rsids2rws],
#                      nea=unlist(Glist$a2)[rsids2rws], 
#                      eaf=unlist(Glist$af)[rsids2rws],
#                      bm=bm[rsids],
#                      dm=dm[rsids], stringsAsFactors = FALSE)
#   fit$stat <- stat
#   fit$stat$vm <- 2*(1-fit$stat$eaf)*fit$stat$eaf*fit$stat$bm^2
#   fit$method <- methods[method+1]
#   
# }

# # Single trait BLR using y and W and sbayes method 
# if(nt==1 && !is.null(y) && !is.null(W) && formatLD=="sparse") {
#   
#   fit <- sbayes_wy(y=y, X=X, W=W, b=b, bm=bm, seb=seb, LD=LD, n=n,
#                 vg=vg, vb=vb, ve=ve, 
#                 ssb_prior=ssb_prior, sse_prior=sse_prior, lambda=lambda, scaleY=scaleY,
#                 h2=h2, pi=pi, updateB=updateB, updateE=updateE, updatePi=updatePi, models=models,
#                 nub=nub, nue=nue, nit=nit, method=method, formatLD=formatLD, algorithm=algorithm)  
# }

# # Single trait BLR using y and sparse LD provided Glist
# if( nt==1 && !is.null(y) && formatLD=="sparse") {
#   
#   if(is.null(Glist)) stop("Please provide Glist")
#   fit <- NULL
#   if(is.matrix(y)) ids <- rownames(y)
#   if(is.vector(y)) ids <- names(y)
#   rws <- match(ids,Glist$ids)
#   if(any(is.na(rws))) stop("some elements in names(y) does not match elements in Glist$ids ")       
#   n <- length(y)
#   
#   if(is.null(chr)) chromosomes <- 1:Glist$nchr
#   if(!is.null(chr)) chromosomes <- chr
#   
#   bm <- dm <- fit <- stat <- vector(length=Glist$nchr,mode="list")
#   names(bm) <- names(dm) <- names(fit) <- names(stat) <- 1:Glist$nchr
#   
#   yy <- sum((y-mean(y))**2)
#   
#   if(is.null(covs)) {
#     covs <- vector(length=Glist$nchr,mode="list")
#     names(covs) <- 1:Glist$nchr
#     for (chr in chromosomes){
#       print(paste("Computing summary statistics for chromosome:",chr))
#       covs[[chr]] <- cvs(y=y,Glist=Glist,chr=chr)
#     }
#   } 
#   
#   if(is.null(nit_local)) nit_local <- nit
#   if(is.null(nit_global)) nit_global <- 1
#   
#   for (it in 1:nit_global) {
#     for (chr in chromosomes){
#       if(verbose) print(paste("Extract sparse LD matrix for chromosome:",chr))
#       LD <- getSparseLD(Glist = Glist, chr = chr, onebased=FALSE)
#       LD$values <- lapply(LD$values,function(x){x*n})
#       rsidsLD <- names(LD$values)
#       clsLD <- match(rsidsLD,Glist$rsids[[chr]])
#       wy <- covs[[chr]][rsidsLD,"wy"]
#       b <- rep(0,length(wy))
#       if(it>1) {
#         b <- fit[[chr]]$b
#         if(updateB) vb <- fit[[chr]]$param[1]
#         if(updateE) ve <- fit[[chr]]$param[2]
#         if(updatePi) pi <- fit[[chr]]$param[3]
#       }
#       stop("Need to add ww and mask to sbayes_sparse - use cvs function to get it")
#       if(verbose) print( paste("Fit",methods[method+1] ,"on chromosome:",chr))
#       fit[[chr]] <- sbayes_sparse(yy=yy, 
#                                   wy=wy,
#                                   b=b, 
#                                   LDvalues=LD$values, 
#                                   LDindices=LD$indices, 
#                                   method=method, 
#                                   nit=nit_local, 
#                                   nburn=nburn, 
#                                   n=n, 
#                                   pi=pi,
#                                   nue=nue, 
#                                   nub=nub, 
#                                   h2=h2, 
#                                   lambda=lambda, 
#                                   vb=vb, 
#                                   ve=ve, 
#                                   updateB=updateB, 
#                                   updateE=updateE, 
#                                   updatePi=updatePi,
#                                   algorithm=algorithm)
#       stat[[chr]] <- data.frame(rsids=rsidsLD,chr=rep(chr,length(rsidsLD)),
#                                 pos=Glist$pos[[chr]][clsLD], ea=Glist$a1[[chr]][clsLD],
#                                 nea=Glist$a2[[chr]][clsLD], eaf=Glist$af[[chr]][clsLD],
#                                 bm=fit[[chr]]$bm,dm=fit[[chr]]$dm,stringsAsFactors = FALSE)
#       rownames(stat[[chr]]) <- rsidsLD
#     }
#   }
#   stat <- do.call(rbind, stat)
#   rownames(stat) <- stat$rsids
#   fit$stat <- stat
#   fit$stat$vm <- 2*(1-fit$stat$eaf)*fit$stat$eaf*fit$stat$bm^2
#   fit$method <- methods[method+1]
#   fit$covs <- covs
# }

# # Single trait BLR based on individual level data based on fast algorithm  
# sbayes_wy <- function(y=NULL, X=NULL, W=NULL, b=NULL, bm=NULL, seb=NULL, LD=NULL, n=NULL,
#                    vg=NULL, vb=NULL, ve=NULL, 
#                    ssb_prior=NULL, sse_prior=NULL, lambda=NULL, scaleY=NULL,
#                    h2=NULL, pi=NULL, updateB=NULL, updateE=NULL, updatePi=NULL, models=NULL,
#                    nub=NULL, nue=NULL, nit=NULL, method=NULL, formatLD=NULL, algorithm=NULL) {
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
#   }    
#   m <- ncol(LD)
#   
#   if(is.null(pi)) pi <- 0.001
#   if(is.null(h2)) h2 <- 0.5
#   
#   if(is.null(ve)) ve <- 1
#   if(method<4 && is.null(vb)) vb <- (ve*h2)/m
#   if(method>=4 && is.null(vb)) vb <- (ve*h2)/(m*pi)
#   if(is.null(lambda)) lambda <- rep(ve/vb,m)
#   if(is.null(vg)) vg <- ve*h2
#   if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
#   if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m*pi)
#   if(is.null(sse_prior)) sse_prior <- ((nue-2.0)/nue)*ve
#   if(is.null(b)) b <- rep(0,m)
#   
#   fit <- .Call("_qgg_sbayes",
#                wy=wy, 
#                #LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
#                LD=as.list(as.data.frame(LD)), 
#                b = b,
#                lambda = lambda,
#                yy = yy,
#                pi = pi,
#                vg = vg,
#                vb = vb,
#                ve = ve,
#                ssb_prior=ssb_prior,
#                sse_prior=sse_prior,
#                nub=nub,
#                nue=nue,
#                updateB = updateB,
#                updateE = updateE,
#                updatePi = updatePi,
#                n=n,
#                nit=nit,
#                method=as.integer(method))
#   names(fit[[1]]) <- rownames(LD)
#   if(!is.null(W)) fit[[7]] <- crossprod(t(W),fit[[10]])[,1]
#   names(fit[[7]]) <- ids
#   stop("check fit names")
#   names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pi","r","param","b")
#   
#   return(fit)
#   
# }

# # Single trait BLR using summary statistics and sparse LD provided in Glist
# sbayes <- function(stat=NULL, b=NULL, seb=NULL, n=NULL,
#                    LD=NULL, LDvalues=NULL,LDindices=NULL,
#                    mask=NULL, lambda=NULL,
#                    vg=NULL, vb=NULL, ve=NULL, h2=NULL, pi=NULL,
#                    ssb_prior=NULL, sse_prior=NULL, nub=4, nue=4,
#                    updateB=TRUE, updateE=TRUE, updatePi=TRUE, updateG=TRUE,
#                    adjustE=TRUE, models=NULL,
#                    nit=500, nburn=100, nthin=1, method="bayesC", algorithm=1, verbose=FALSE) {
#   
#   # Check method
#   methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
#   method <- match(method, methods) - 1
#   if( !sum(method%in%c(0:5))== 1 ) stop("Method specified not valid")
#   
#   # Prepare summary statistics input
#   if( is.null(stat) ) stop("Please provide summary statistics")
#   m <- nrow(stat)
#   if(is.null(mask)) mask <- rep(FALSE, m)
#   if(is.null(stat$n)) stat$n <- stat$dfe
#   if( is.null(stat$n) ) stop("Please provide summary statistics that include n")
#   
#   n <- as.integer(median(stat$n))
#   ww <- 1/(stat$seb^2 + stat$b^2/stat$n)
#   wy <- stat$b*ww
#   if(!is.null(stat$ww)) ww <- stat$ww
#   if(!is.null(stat$wy)) wy <- stat$wy
#   
#   b2 <- stat$b^2
#   seb2 <- stat$seb^2
#   yy <- (b2 + (stat$n-2)*seb2)*ww
#   yy <- median(yy)
#   
#   # Prepare sparse LD matrix
#   if( is.null(LD) ) stop("Please provide LD matrix")
#   #LDvalues <- split(LD, rep(1:ncol(LD), each = nrow(LD)))
#   LDvalues=as.list(as.data.frame(LD)) 
#   LDindices <- lapply(1:ncol(LD),function(x) { (1:ncol(LD))-1 } )
#   
#   # Prepare starting parameters
#   vy <- yy/(n-1)
#   if(is.null(pi)) pi <- 0.001
#   if(is.null(h2)) h2 <- 0.5
#   if(is.null(ve)) ve <- vy*(1-h2)
#   if(is.null(vg)) vg <- vy*h2
#   if(method<4 && is.null(vb)) vb <- vg/m
#   if(method>=4 && is.null(vb)) vb <- vg/(m*pi)
#   if(is.null(lambda)) lambda <- rep(ve/vb,m)
#   if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
#   if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/(m*pi))
#   if(is.null(sse_prior)) sse_prior <- ((nue-2.0)/nue)*ve
#   if(is.null(b)) b <- rep(0,m)
#   
#   pi <- c(1-pi,pi)
#   gamma <- c(0,1.0)
#   if(method==5) pi <- c(0.95,0.02,0.02,0.01)
#   if(method==5) gamma <- c(0,0.01,0.1,1.0)
#   
#   seed <- sample.int(.Machine$integer.max, 1)
#   
#   fit <- .Call("_qgg_sbayes_spa",
#                wy=wy,
#                ww=ww,
#                LDvalues=LDvalues,
#                LDindices=LDindices,
#                b = b,
#                lambda = lambda,
#                mask=mask,
#                yy = yy,
#                pi = pi,
#                gamma = gamma,
#                vg = vg,
#                vb = vb,
#                ve = ve,
#                ssb_prior=ssb_prior,
#                sse_prior=sse_prior,
#                nub=nub,
#                nue=nue,
#                updateB = updateB,
#                updateE = updateE,
#                updatePi = updatePi,
#                updateG = updateG,
#                adjustE = adjustE,
#                n=n,
#                nit=nit,
#                nburn=nburn,
#                nthin=nthin,
#                method=as.integer(method),
#                algo=as.integer(algorithm),
#                seed=seed)
#   names(fit[[1]]) <- names(LDvalues)
#   names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pis","pim","r","b","param")
#   return(fit)
# }
# 
# # Single trait BLR using summary statistics and sparse LD provided in Glist
# sbayesXy <- function(yy=NULL, Xy=NULL, XX=NULL, n=NULL,
#                      mask=NULL, lambda=NULL,
#                      vg=NULL, vb=NULL, ve=NULL, h2=NULL, pi=NULL,
#                      ssb_prior=NULL, sse_prior=NULL, nub=4, nue=4,
#                      updateB=TRUE, updateE=TRUE, updatePi=TRUE, updateG=TRUE,
#                      adjustE=TRUE, models=NULL,
#                      nit=500, nburn=100, nthin=1, method="bayesC", algorithm="mcmc", verbose=FALSE) {
#   
#   # Check methods and parameter settings
#   methods <- c("blup","bayesN","bayesA","bayesL","bayesC","bayesR")
#   method <- match(method, methods) - 1
#   if( !sum(method%in%c(0:5))== 1 ) stop("method argument specified not valid")
#   algorithms <- c("mcmc","em-mcmc")
#   algorithm <- match(algorithm, algorithms)
#   if(is.na(algorithm)) stop("algorithm argument specified not valid")
#   
#   if( is.null(n) ) stop("Please provide n")
#   if( is.null(yy) ) stop("Please provide yy")
#   if( is.null(Xy) ) stop("Please provide Xy matrix")
#   if( is.null(XX) ) stop("Please provide XX matrix")
#   
#   xx <- diag(XX)  
#   m <- length(xx)
#   
#   XX <- cov2cor(XX)  
#   #XXvalues <- split(XX, rep(1:ncol(XX), each = nrow(XX)))
#   XXvalues <- as.list(as.data.frame(XX)) 
#   
#   XXindices <- lapply(1:ncol(XX),function(x) { (1:ncol(XX))-1 } )
#   
#   b <- rep(0, m)
#   mask <- rep(FALSE, m)
#   
#   
#   
#   # Prepare starting parameters
#   vy <- yy/(n-1)
#   if(is.null(pi)) pi <- 0.001
#   if(is.null(h2)) h2 <- 0.5
#   if(is.null(ve)) ve <- vy*(1-h2)
#   if(is.null(vg)) vg <- vy*h2
#   if(method<4 && is.null(vb)) vb <- vg/m
#   if(method>=4 && is.null(vb)) vb <- vg/(m*pi)
#   if(is.null(lambda)) lambda <- rep(ve/vb,m)
#   if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
#   if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/(m*pi))
#   if(is.null(sse_prior)) sse_prior <- ((nue-2.0)/nue)*ve
#   if(is.null(b)) b <- rep(0,m)
#   
#   pi <- c(1-pi,pi)
#   gamma <- c(0,1.0)
#   if(method==5) pi <- c(0.95,0.02,0.02,0.01)
#   if(method==5) gamma <- c(0,0.01,0.1,1.0)
#   
#   seed <- sample.int(.Machine$integer.max, 1)
#   
#   fit <- .Call("_qgg_sbayes_spa",
#                wy=Xy,
#                ww=xx,
#                LDvalues=XXvalues,
#                LDindices=XXindices,
#                b = b,
#                lambda = lambda,
#                mask=mask,
#                yy = yy,
#                pi = pi,
#                gamma = gamma,
#                vg = vg,
#                vb = vb,
#                ve = ve,
#                ssb_prior=ssb_prior,
#                sse_prior=sse_prior,
#                nub=nub,
#                nue=nue,
#                updateB = updateB,
#                updateE = updateE,
#                updatePi = updatePi,
#                updateG = updateG,
#                adjustE = adjustE,
#                n=n,
#                nit=nit,
#                nburn=nburn,
#                nthin=nthin,
#                method=as.integer(method),
#                algo=as.integer(algorithm),
#                seed=seed)
#   names(fit[[1]]) <- names(XXvalues)
#   names(fit) <- c("bm","dm","coef","vbs","vgs","ves","pis","pim","r","b","param")
#   return(fit)
# }

# computeB <- function(wy=NULL, yy=NULL, b=NULL, WW=NULL, n=NULL,
#                      vb=NULL, vg=NULL, ve=NULL, lambda=NULL, 
#                      ssb_prior=NULL, sse_prior=NULL, 
#                      nub=NULL, nue=NULL, 
#                      h2=NULL, pi=NULL, 
#                      updateB=NULL, updateE=NULL, updatePi=NULL,
#                      nit=NULL, nburn=NULL, method=NULL) {
#   
#   m <- ncol(WW)
#   
#   if(is.null(pi)) pi <- 0.001
#   if(is.null(h2)) h2 <- 0.5
#   
#   if(is.null(ve)) ve <- 1
#   if(method<4 && is.null(vb)) vb <- (ve*h2)/m
#   if(method>=4 && is.null(vb)) vb <- (ve*h2)/(m*pi)
#   if(is.null(lambda)) lambda <- rep(ve/vb,m)
#   if(is.null(vg)) vg <- ve*h2
#   if(method<4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m)
#   if(method>=4 && is.null(ssb_prior))  ssb_prior <-  ((nub-2.0)/nub)*(vg/m*pi)
#   if(is.null(sse_prior)) sse_prior <- nue*ve
#   if(is.null(b)) b <- rep(0,m)
#   
#   fit <- .Call("_qgg_sbayes",
#                wy=wy, 
#                LD=split(WW, rep(1:ncol(WW), each = nrow(WW))), 
#                b = b,
#                lambda = lambda,
#                yy = yy,
#                pi = pi,
#                vg = vg,
#                vb = vb,
#                ve = ve,
#                ssb_prior=ssb_prior,
#                sse_prior=sse_prior,
#                nub=nub,
#                nue=nue,
#                updateB = updateB,
#                updateE = updateE,
#                updatePi = updatePi,
#                n=n,
#                nit=nit,
#                method=as.integer(method))
#   
#   names(fit[[1]]) <- rownames(WW)
#   stop("Check names of fit again")
#   names(fit) <- c("bm","dm","mus","vbs","ves","pis","wy","r","param","b")
#   
#   return(fit)
#   
# }
# 
# computeWy <- function(y=NULL, Glist = NULL, chr = NULL, cls = NULL) {
#   wy <- cvs(y=y,Glist=Glist,chr=chr,cls=cls)$wy
#   return(wy)
# }
# 
# computeWW <- function(Glist = NULL, chr = NULL, cls = NULL, rws=NULL, scale=TRUE) { 
#   W <- getG(Glist=Glist, chr=chr, cls=cls, scale=scale)
#   WW <- crossprod(W[rws,])
#   return(WW)
# }
