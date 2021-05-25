####################################################################################################################
#    Module 7: LDSC
####################################################################################################################
#'
#' LDSC regression analysis
#'
#'
#' @description
#' The ldsc function is used for LDSC analysis
#'
#' @param y vector or matrix of phenotypes
#' @param X design matrix of fixed effects
#' @param W matrix of centered and scaled genotypes
#' @param Glist list of information about genotype matrix stored on disk
#' @param rsids vector of marker rsids used in the analysis
#' @param ids vector of individuals used in the analysis
#' @param lambda overall shrinkage factor
#' @param weights vector of single marker weights used in BLUP
#' @param method used in solver (currently only methods="gsru": gauss-seidel with resiudal update)
#' @param maxit maximum number of iterations used in the Gauss-Seidel procedure
#' @param tol tolerance, i.e. the maximum allowed difference between two consecutive iterations of the solver to declare convergence
#' @param sets	list containing marker sets rsids
#' @param scale logical if TRUE the genotypes in Glist will be scaled to mean zero and variance one
#' @param ncores number of cores used in the analysis


#' @author Peter Soerensen

#' @examples
#'
#' # Simulate data
#' W <- matrix(rnorm(1000000), ncol = 1000)
#' 	colnames(W) <- as.character(1:ncol(W))
#' 	rownames(W) <- as.character(1:nrow(W))
#' m <- ncol(W)
#' causal <- sample(1:ncol(W),50)
#' y <- rowSums(W[,causal]) + rnorm(nrow(W),sd=sqrt(50))
#'
#' X <- model.matrix(y~1)
#'
#' Sg <- 50
#' Se <- 50
#' h2 <- Sg/(Sg+Se)
#' lambda <- Se/(Sg/m)
#' lambda <- m*(1-h2)/h2
#'
#' # BLUP of single marker effects and total genomic effects based on Gauss-Seidel procedure
#' fit <- gsolve( y=y, X=X, W=W, lambda=lambda)
#'



#' @export
ldscore <- function(Glist=NULL, chr=NULL, onebased=TRUE, nbytes=4) {
  
  chromosomes <- chr
  
  if(is.null(chr)) chromosomes <- 1:22  
  
  lscore2 <- vector(length=length(chromosomes),mode="list")
  
  for (chr in chromosomes) {
    message(paste("Processing chromosome:",chr))
    m <- Glist$mchr[chr]
    msize <- Glist$msize
    rsids <- Glist$rsids[[chr]]
    
    # LD indexes
    k1 <- rep(1, m)
    k1[1:msize] <- msize - 1:msize + 2
    k2 <- rep((2 * msize + 1), m)
    k2[(m - msize + 1):m] <- msize + m - ((m - msize + 1):m) + 1
    
    ldchr <- rep(0,m)
    names(ldchr) <- rsids
    
    fnLD <- Glist$fnLD[[chr]]
    bfLD <- file(fnLD, "rb")
    
    for (j in 1:m) {
      rwsLD <- k1[j]:k2[j]
      ld <- readBin(bfLD, "numeric", n = (2*msize+1), size = nbytes, endian = "little")
      ldchr[j] <- sum(ld[rwsLD]**2)
    }
    close(bfLD)
    lscore2[[chr]] <- ldchr
  }
  
  return(unlist(lscore2))
}


#' @export
ldsc <- function(Glist=NULL, lscore=NULL, z=NULL, b=NULL, seb=NULL, p=NULL, n=NULL, sets=NULL, intercept=TRUE, what="h2") {
  if(!is.null(Glist) & is.null(lscore) ) lscore <- unlist(Glist$lscore)
  
  if(!is.null(z)) nt <- ncol(z)
  
  if(!is.null(b)) {
    nt <- ncol(b)
    for (t in 1:nt) {
      z <- cbind(z,(b[,t]/seb[,t]))
    }
  }  

  lscore <- lscore[rownames(z)]
  if(is.null(n)) {
    n <- NULL
    if(is.null(seb)) stop("Please provide n or alternatively seb")
    for ( t in 1:nt) {
      n <- c(n,neff(seb[,t],af[,t]))
    }
  }
  maxZ2 <- max(0.001 * max(n), 80)
  h2 <- NULL
  for ( t in 1:nt) {
    #z2 <- (b[,t]/seb[,t])**2
    z2 <- z[,t]**2
    z2 <- z2[!is.na(z2)] 
    z2 <- z2[z2<800]
    if(intercept) X <- cbind(1,n[t]*lscore[names(z2)]/length(z2))
    if(!intercept) X <- matrix(n[t]*lscore[names(z2)]/length(z2),ncol=1)
    y <- z2
    XtX <- crossprod(X)
    Xy <- crossprod(X,y)
    h2_ldsc <- solve(XtX, Xy)
    if(intercept){
      if(h2_ldsc[2]<0) h2_ldsc[2] <- NA
    }
    if(!intercept){
      if(h2_ldsc[1]<0) h2_ldsc[1] <- NA
    }
    h2 <- c(h2,h2_ldsc)
  }
  if(intercept) {
    h2 <- matrix(h2,ncol=2,byrow=TRUE)
    rownames(h2) <- colnames(z)
    colnames(h2) <- c("intercept","h2")
  }
  if(!intercept) {
    h2 <- matrix(h2,ncol=1,byrow=TRUE)
    rownames(h2) <- colnames(z)
    colnames(h2) <- "h2"
  }
  result <- h2
  
  #---------------------------------#
  # Block Jackknife to estimate h2 SE
  if(SE==TRUE){
    if(intercept==TRUE){
      
    }
  }
      
  
  
  
  if(what=="correlation") {
    rg <- matrix(0,nt,nt)
    rownames(rg) <- colnames(rg) <- colnames(b)
    for (t1 in 1:nt) {
      for (t2 in t1:nt) {
        #Z1 <- (b[,t1]/seb[,t1])
        #Z2 <- (b[,t2]/seb[,t2])
        #Z2_1 <- (b[,t1]/seb[,t1])**2
        #Z2_2 <- (b[,t2]/seb[,t2])**2
        Z1 <- z[,t1]
        Z2 <- z[,t2]
        Z2_1 <- Z1**2
        Z2_2 <- Z2**2
        rws1 <- Z2_1<maxZ2
        rws2 <- Z2_2<maxZ2
        rws <- rws1 & rws2
        m <- sum(rws)
        X <- cbind(1,sqrt(n[t1])*sqrt(n[t2])*lscore[rws]/m)
        y <- Z1[rws]*Z2[rws]
        XtX <- crossprod(X)
        Xy <- crossprod(X,y)
        if(intercept & !any(is.na(h2[c(t1,t2),2]))) rg[t1,t2] <- solve(XtX, Xy)[2]/(sqrt(h2[t1,2])*sqrt(h2[t2,2]))
        if(!intercept & !any(is.na(h2[c(t1,t2)]))) rg[t1,t2] <- solve(XtX, Xy)[2]/(sqrt(h2[t1])*sqrt(h2[t2]))
      }
      if(intercept) rg[t1,t1] <- h2[t1,2]
      if(!intercept) rg[t1,t1] <- h2[t1]
    }
    rownames(rg) <- colnames(rg) <- colnames(z)
    result <- rg
  }
  return(result)
}

neff <- function(seb=NULL,af=NULL,Vy=1) {
  seb2 <- seb**2
  vaf <- 2*af*(1-af)
  neff <- round(median(Vy/(vaf*seb2)))
  return(neff)
}
