####################################################################################################################
#    Module 2: Marker set tests
####################################################################################################################
#' 
#' Multi marker association analysis based on marker set tests
#'
#' @description
#' Genetic marker set tests based on summing the single genetic marker test statistics.
#' The single marker test statistics can be obtained from GBLUP and GFBLUP model fits or from standard GWAS. 

#' @details 
#' The sum test is powerful if the genomic feature harbors many genetic markers having small to moderate effects. 
#' The distribution of this test statistic under the null hypothesis (associated markers are picked at random from the total 
#' number of tested genetic markers) is difficult to describe in terms of exact or approximate 
#' distributions, and an empirical distribution is required.
#' Genetic marker set tests based on the covariance statistics for a set of genetic markers.
#' The covariance test statistic is derived from a GBLUP (or GFBLUP) model fit. It is a measure of covariance between the total genomic effect for all markers 
#' and the genomic effect for the genetic markers in the genomic feature. It also relates to the explained sums of
#' squares for the genetic markers. 
#' The distribution of these test statistics under the null hypothesis is difficult to describe in terms of exact or approximate 
#' distributions, and an empirical distribution is required.
#' Genetic marker set tests based on the hyperG test statistics for a set of genetic markers.
#' The hyperG marker set test tests a predefined set of markers (i.e. those within a particular genomic feature)
#' for an association with the trait phenotype.
#' Under the null hypothesis (associated markers are picked at random from the total number of tested 
#' genetic markers) it is assumed that the observed count statistic is a realization from a hypergeometric 
#' distribution.
#' This hypothesis can be formulated and tested in a number of ways. Here we consider a test statistic based 
#' on counting the number of genetic markers in the feature that are associated to trait phenotype. 
#' A test based on the count test statistic is likely to have high power to detect association if the genomic 
#' feature harbours genetic markers with large effects. 
#'                        
#' @param stat vector of single marker statistics (e.g. marker effects, t-stat, p-values)
#' @param sets list of marker sets - names corresponds to rownames in stat
#' @param nperm number of permutations
#' @param ncores number of cores
#' @param W matrix of centered and scaled genotypes (used if method = cvat or score)
#' @param fit is the fit object obtained from a linear mixed model fit using the greml function
#' @param g vector (or list) of genetic effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param s vector (or list) of single marker effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param method including sum, cvat, hyperG, score
#' @param threshold used if method = hyperG
#' @return Returns a dataframe including 
#' \item{setT}{marker set test statistics} 
#' \item{nset}{number of markers in the set}
#' \item{p}{p-value for marker set}
#' @author Peter Soerensen
#' @references Rohde, P. D., Demontis, D., Cuyabano, B. C. D., Boerglum, A. D., & Soerensen, P. (2016). Covariance Association Test (CVAT) Identifies Genetic Markers Associated with Schizophrenia in Functionally Associated Biological Processes. Genetics, 203(4), 1901-1913.
#' @references Rohde, P. D., Edwards S. M., Sarup P., Soerensen, P. (August, 2014). Gene-based Association Approach Identify Genes Across Stress Traits in Fruit Flies. Poster presented at the 10th World Congress of Genetics Applied to Livestock Production (WCGALP), Vancouver, Canada.
#' @examples
#'
#'  
#'  # Simulate data
#'  W <- matrix(rnorm(20000000), ncol = 10000)
#'  colnames(W) <- as.character(1:ncol(W))
#'  rownames(W) <- as.character(1:nrow(W))
#'  y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#'  
#'  # Create model
#'  data <- data.frame(y = y, mu = 1)
#'  fm <- y ~ 0 + mu
#'  X <- model.matrix(fm, data = data)
#'  
#'  # Compute GRM
#'  GB <- computeGRM(W = W)
#'  
#'  # REML analyses and single marker association test
#'  fit <- greml(y = y, X = X, GRM = list(GB), verbose = TRUE)
#'  ma <- mlma(fit = fit, W = W)
#'  ma <- as.matrix(ma)
#'  
#'  f <- factor(rep(1:100,each=100), levels=1:100) 
#'  sets <- split(as.character(1:10000),f=f) 
#'
#' # Set test based on sums 
#' res <- mma(stat = ma[,3]**2, sets = sets, method = "sum", nperm = 10000)
#' 
#' 
#' # Set test based on hyperG 
#' res <- mma(stat = ma[,4], sets = sets, method = "hyperG", threshold = 0.05)
#' 


#' @export

mma <- function(stat=NULL,sets=NULL,ncores=1, nperm=1000, method="sum") {
     
     m <- length(stat)
     #rws <- 1:m 
     #names(rws) <- names(stat)
     #sets <- lapply(sets, function(x) {rws[x]}) 
     if (is.matrix(stat)) sets <- mapSets(sets=sets,rsids=rownames(stat),index=TRUE)
     if (is.vector(stat)) sets <- mapSets(sets=sets,rsids=names(stat),index=TRUE)
     
     nsets <- length(sets)
     msets <- sapply(sets, length)
     
     if (is.matrix(stat)) { 
          p <- apply(stat,2, function(x) { gsets(stat=x,sets=sets, ncores=ncores, np=nperm) } )
          setstat <- apply(stat,2,function(x) { sapply(sets, function(y) {sum(x[y])}) })
          rownames(setstat) <- rownames(p) <- names(msets) <- names(sets)  
          res <- list(m=msets,stat=setstat,p=p)
     }
     if (is.vector(stat)) { 
          setstat <- sapply(sets, function(x) {sum(stat[x])})
          p <- gsets(stat=stat,sets=sets, ncores=ncores, np=nperm, method=method) 
          res <- cbind(m=msets,stat=setstat,p=p)
          rownames(res) <- names(sets)  
     }     
     res
     
}

#' @export

setTest <- function(stat = NULL, W = NULL, sets = NULL, nperm = NULL, method = "sum", threshold = 0.05) {
     
  if (method == "sum") setT <- sumTest(stat = stat, sets = sets, nperm = nperm) 
  if (method == "cvat") setT <- cvat(s = stat, W = W, sets = sets, nperm = nperm) 
  if (method == "hyperG") setT <- hgTest(p = stat, sets = sets, threshold = threshold) 
  if (method == "score") setT <- scoreTest(e = e, W = W, sets = sets, nperm = nperm)
     
  return(setT)

}

sumTest <- function(stat = NULL, sets = NULL, nperm = NULL, method = "sum") {
     
  if (method == "mean") setT <- sapply(sets, function(x) {mean(stat[x])})
  if (method == "sum") setT <- sapply(sets, function(x) {sum(stat[x])})
  if (method == "max") setT <- sapply(sets, function(x) {max(stat[x])})
  if (!is.null(nperm)) {
    p <- rep(0, length(sets)) 
    n <- length(stat)
    nset <- sapply(sets, length)
    rws <- 1:n
    names(rws) <- names(stat)
    sets <- lapply(sets, function(x) {rws[x]}) 
    for (i in 1:nperm) {
      rws <- sample(1:n, 1)
      o <- c(rws:n, 1:(rws - 1))
      pstat <- stat[o]
      if (method == "mean") setTP <- sapply(sets, function(x) {mean(pstat[x])})
      if (method == "sum") setTP <- sapply(sets, function(x) {sum(pstat[x])})
      if (method == "max") setTP <- sapply(sets, function(x) {max(pstat[x])})
      p <- p + as.numeric(setT > setTP) 
    }  
    p <- 1 - p / nperm
    setT <- data.frame(setT, nset, p)
  }
     
  return(setT)

}

msetTest <- function(stat = NULL, sets = NULL, nperm = NULL, method = "sum") {
     
     setT <- apply(stat, 2, function(x) {setTest(stat = x, sets = sets, nperm = nperm, method = method)})
     names(setT) <- colnames(stat)
     setT  

} 

gsett <- function(stat = NULL, W = NULL, sets = NULL, nperm = NULL, method = "sum", threshold = 0.05) {
     
     m <- length(stat)
     rws <- 1:m
     names(rws) <- names(stat)
     sets <- lapply(sets, function(x) {rws[x]}) 
     
     setT <- sapply(sets, function(x) {sum(stat[x])})
     
     if (!is.null(nperm)) {
          p <- rep(0, length(sets)) 
          n <- length(stat)
          nset <- sapply(sets, length)
          sets <- lapply(sets, function(x) {rws[x]}) 
          for (i in 1:nperm) {
               rws <- sample(1:n, 1)
               o <- c(rws:n, 1:(rws - 1))
               ostat <- stat[o]
               setTP <- sapply(sets, function(x) {sum(ostat[x])})
               p <- p + as.numeric(setT > setTP) 
          }  
          p <- 1 - p / nperm
          setT <- data.frame(setT, nset, p)
     }
     return(setT)
}

#' @export

cvat <- function(fit = NULL, s = NULL, g = NULL, W = NULL, sets = NULL, nperm = 100) {
     if (!is.null(fit)) {s <- crossprod(W/ncol(W),fit$Py)*fit$theta[1]}
     Ws <- t(t(W) * as.vector(s))
     if (is.null(g)) g <- W %*% s   
     cvs <- colSums(as.vector(g) * Ws)
     #setT <- setTest(stat = cvs, sets = sets, nperm = nperm, method = "sum")$p
     #names(setT) <- names(sets)
     setT <- setTest(stat = cvs, sets = sets, nperm = nperm, method = "sum")
     if (!is.null(names(sets))) rownames(setT) <- names(sets)
     return(setT)

}

scoreTest <- function(e = NULL, W = NULL, sets = NULL, nperm = 100) {
     
     we2 <- as.vector((t(W) %*% e)**2)   
     names(we2) <- colnames(W)       
     setT <- setTest(stat = we2, sets = sets, nperm = nperm, method = "sum")$p
     return(setT)

}


#' @export

hgTest <- function(p = NULL, sets = NULL, threshold = 0.05) {
     
     N <- length(p)
     Na <- sum(p < threshold)
     Nna <- N - Na
     Nf <- sapply(sets, length)
     Naf <- sapply(sets, function(x) {sum(p[x] < threshold)})
     Nnaf <- Nf - Naf
     Nanf <- Na - Naf
     Nnanf <- Nna - Nnaf
     phyperg <- 1 - phyper(Naf - 1, Nf, N - Nf, Na)
     phyperg

}



#' @export

gsets <- function(stat=NULL,sets=NULL,ncores=1, np=1000, method="sum") {
     
     m <- length(stat)
     nsets <- length(sets)
     msets <- sapply(sets, length)
     setstat <- sapply(sets, function(x) {sum(stat[x])})
   

     
     res <- .Fortran("psets", 
                     m = as.integer(m),
                     stat = as.double(stat),
                     nsets = as.integer(nsets),
                     setstat = as.double(setstat),
                     msets = as.integer(msets),
                     p = as.integer(rep(0,nsets)),
                     np = as.integer(np),
                     ncores = as.integer(ncores),
                     PACKAGE = 'qgg'
     )
     
     p <- res$p/np
}


#' @export

mapSets <- function( sets=NULL, rsids=NULL, Glist=NULL, index=TRUE ) { 
     if(!is.null(Glist)) rsids <- unlist(Glist$rsids)
     nsets <- sapply(sets,length)
     rs <- rep(names(sets),times=nsets)
     rsSets <- unlist(sets,use.names=FALSE)
     rsSets <- match(rsSets,rsids)
     inW <- !is.na(rsSets)
     rsSets <- rsSets[inW]
     if(!index) rsSets <- rsids[rsSets]
     rs <-  rs[inW]
     rs <- factor(rs, levels=unique(rs))  
     rsSets <- split(rsSets,f=rs)
     return(rsSets)
}

