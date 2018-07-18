####################################################################################################################
#    Module 2: SetTests
####################################################################################################################
#' 
#' Genetic marker set tests based on sum statistics
#'
#' @description
#' Set test based on summing the single genetic marker test statistics.
#' The sum test is powerful if the genomic feature harbors many genetic markers having small to moderate effects. 
#'                       
#' @details
#' The singler marker test statistics can be obtained from GBLUP and GFBLUP model fits or from standard GWAS. 
#' The distribution of this test statistic under the null hypothesis (associated markers are picked at random from the total 
#' number of tested genetic markers) is difficult to describe in terms of exact or approximate 
#' distributions, and an empirical distribution is required.
#'                        
#' @param stat vector of single marker statistics (e.g. marker effects, t-stat, p-values)
#' @param sets list of marker sets - names corresponds to rownames in stat
#' @param nperm number of permutations
#' @param W matrix of centered and scaled genotypes (used if method = cvat or score)
#' @param method including sum, cvat, hyperG, score
#' @param threshold used if method = hyperG
#' @return Returns a dataframe including 
#' \item{setT}{marker set test statistics} 
#' \item{nset}{number of markers in the set}
#' \item{p}{p-value for marker set}
#' @author Peter Sørensen
#’ @references Rohde, P. D., Demontis, D., Cuyabano, B. C. D., Børglum, A. D., & Sørensen, P. (2016). Covariance Association Test (CVAT) Identifies Genetic Markers Associated with Schizophrenia in Functionally Associated Biological Processes. Genetics, 203(4), 1901-1913.
#’ @references Rohde, P. D., Edwards S. M., Sarup P., Sørensen, P. (August, 2014). Gene-based Association Approach Identify Genes Across Stress Traits in Fruit Flies. Poster presented at the 10th World Congress of Genetics Applied to Livestock Production (WCGALP), Vancouver, Canada.
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(2000000), ncol = 10000)
#'   colnames(W) <- as.character(1:ncol(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#' 
#' # REML analyses 
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ mu
#' fit <- gfm(fm = fm, W = W, sets = list(colnames(W)), data = data)
#' fit$df <- 10
#' fit$p <- pt(fit$s / sqrt(fit$vs), df = fit$df, lower.tail = FALSE) 
#' 
#' sets <- list(A = as.character(1:100), B = as.character(101:1000), C = as.character(1001:5000), D = as.character(5001:10000))
#'
#' # Set test based on sums 
#' res <- setTest(stat = fit$s**2, sets = sets, method = "sum", nperm = 100)
#' 
#' # Set test based on cvat 
#' res <- setTest(stat = fit$s, W = W, sets = sets, method = "cvat", nperm = 100)
#' 
#' # Set test based on hyperG 
#' res <- setTest(stat = fit$p, sets = sets, method = "hyperG", threshold = 0.05)
#' 
#' @export
#'

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

####################################################################################################################
#' 
#' Genetic marker set tests based on the covariance test
#'
#' @description
#' Genetic marker set tests based on the covariance statistics for a set of genetic markers.
#' 
#' @details
#' The covariance test statistic is derived from a GBLUP (or GFBLUP) model fit. It is a measure of covariance between the total genomic effect for all markers 
#' and the genomic effect for the genetic markers in the genomic feature. It also relates to the explained sums of
#' squares for the genetic markers. 
#' The distribution of this test statistic under the null hypothesis is difficult to describe in terms of exact or approximate 
#' distributions, and an empirical distribution is required.
#' 
#' @param fit is the fit object obtained from a linear mixed model fit using the greml function
#' @param g vector (or list) of genetic effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param s vector (or list) of single marker effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param W matrix of centered and scaled genotypes (n x m)
#' @param sets list of marker sets corresponding to columns in W
#' @param nperm number of permutations
#' @return Returns a dataframe including 
#' \item{setT}{covariance test statistics} 
#' \item{nset}{number of markers in the set}
#' \item{p}{p-value}
#' @author Peter Sørensen
#’ @references Rohde, P. D., Demontis, D., Cuyabano, B. C. D., Børglum, A. D., & Sørensen, P. (2016). Covariance Association Test (CVAT) Identifies Genetic Markers Associated with Schizophrenia in Functionally Associated Biological Processes. Genetics, 203(4), 1901-1913.
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(20000000), ncol = 10000)
#' 	colnames(W) <- as.character(1:ncol(W))
#' 	rownames(W) <- as.character(1:nrow(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#'
#' # Create model
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ 0 + mu
#' X <- model.matrix(fm, data = data)
#'
#' # Create framework for lists
#' setsGB <- list(A = colnames(W)) # gblup model
#' setsGF <- list(C1 = colnames(W)[1:1000], C2 = colnames(W)[1001:2000], C3 = colnames(W)[2000:10000]) # gfblup model
#' setsGT <- list(C1 = colnames(W)[1:10], C2 = colnames(W)[1001:1010], C3 = colnames(W)[1:10000]) # true model
#'
#' # Compute G
#' G <- computeG(W = W)
#' GB <- lapply(setsGB, function(x) {computeG(W = W[, x])})
#' GF <- lapply(setsGF, function(x) {computeG(W = W[, x])})
#' GT <- lapply(setsGT, function(x) {computeG(W = W[, x])})
#'
#' # REML analyses and multi marker association (set) test
#' fitGB <- greml(y = y, X = X, G = GB, verbose = TRUE)
#'
#' # Use fit object as input
#' cvat(fit = fitGB, W = W, sets = setsGF, nperm = 1000)
#' cvat(fit = fitGB, W = W, sets = setsGT, nperm = 1000)
#'
#' # Use single coefficients as input 
#' s <- crossprod(W / ncol(W), fitGB$Py) * fitGB$theta[1]
#' cvat(s = s, W = W, sets = setsGF, nperm = 1000)
#' cvat(s = s, W = W, sets = setsGT, nperm = 1000)
#' 
#' @export
#'

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



####################################################################################################################
#' 
#' Genetic marker set tests based on the hyperG test
#' 
#' @description
#' Genetic marker set tests based on the hyperG test statistics for a set of genetic markers.
#'
#' @details
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
#' @param p vector of single marker p-values (e.g. based on a t-stat)
#' @param sets list of marker sets - names corresponds to rownames in stat
#' @param threshold single marker p-value cut-off
#' @return Returns vector of p values with length equal to the number of sets 
#' @author Peter Sørensen
#’ @references Rohde, P. D., Edwards S. M., Sarup P., Sørensen, P. (August, 2014). Gene-based Association Approach Identify Genes Across Stress Traits in Fruit Flies. Poster presented at the 10th World Congress of Genetics Applied to Livestock Production (WCGALP), Vancouver, Canada.
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(2000000), ncol = 10000)
#'   colnames(W) <- as.character(1:ncol(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#' 
#' # REML analyses 
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ mu
#' fit <- gfm(fm = fm, W = W, sets = list(colnames(W)), data = data)
#' 
#' # hyperG set test 
#' sets <- list(A = as.character(1:100), B = as.character(101:1000), C = as.character(1001:5000), D = as.character(5001:10000))
#' p <- pt(fit$s / fit$vs)
#' res <- hgTest(p = p, sets = sets, threshold = 0.05)
#' 
#' @export
#'

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

gsets <- function(stat=NULL,sets=NULL,ncores=1, np=1000) {
     
     m <- length(stat)
     #rws <- 1:m 
     #names(rws) <- names(stat)
     #sets <- lapply(sets, function(x) {rws[x]}) 

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
     
     res$p/np
     
}
