####################################################################################################################
#    Module 8: GSIM
####################################################################################################################
#'
#' Genomic simulation
#' 
#'
#' @description
#' The gsolve function is used for solving of linear mixed model equations. The algorithm used to solve the equation
#' system is based on a Gauss-Seidel (GS) method (matrix-free with residual updates) that handles large data sets.
#'
#' The linear mixed model fitted can account for multiple traits, multiple genetic factors (fixed or random genetic
#' marker effects), adjust for complex family relationships or population stratification, and adjust for other
#' non-genetic factors including lifestyle characteristics. Different genetic architectures (infinitesimal,
#' few large and many small effects) is accounted for by modeling genetic markers in different sets as fixed or
#' random effects and by specifying individual genetic marker weights.

#'
#' @param y vector or matrix of phenotypes
#' @param X design matrix of fixed effects
#' @param W matrix of centered and scaled genotypes
#' @param Glist list of information about genotype matrix stored on disk


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




#'
#' @export
#'

gsim <- function(nt=1,W=NULL,n=1000,m=1000) {
  if(is.null(W)) {
    W <- matrix(runif(n),ncol=1)
    for (i in 2:m) {
      W <- cbind(W,scale(W[,i-1]) + runif(n))  
    }
  }
  n <- nrow(W)
  m <- ncol(W)
  if(is.null(colnames(W))) colnames(W) <- paste0("m",1:m)
  if(is.null(rownames(W))) rownames(W) <- paste0("id",1:n)
  
  y <- e <- vector(length=nt,mode="list")
  names(y) <- paste0("D",1:nt)
  set0 <- sample(1:ncol(W),2)
  set1 <- b1 <- g1 <- vector(length=nt,mode="list")
  g <- NULL
  for (i in 1:nt){
    b0 <- sample(c(0.25,-0.25,0.5,-0.5),2)
    g0 <- W[,set0]%*%b0
    set1[[i]] <- sample(1:ncol(W),2)
    b1[[i]] <- sample(c(0.25,-0.25,0.5,-0.5),length(set1[[i]]))
    g1[[i]] <- W[,set1[[i]]]%*%b1[[i]]
    e[[i]] <- rnorm(nrow(W),mean=0,sd=1)
    y[[i]] <- as.vector(g0+g1[[i]]+e[[i]])
    names(y[[i]]) <- rownames(W)
    g <- cbind(g,g0+g1[[i]])
  }
  colnames(g) <- paste0("D",1:nt) 
  if(nt==1) return( list( y=y[[1]],W=W, e=e[[1]],g=g,b0=b0,b1=b1,set0=set0,set1=set1,causal=c(set0,unlist(set1))))
  if(nt>1) return( list( y=y,W=W, e=e,g=g,b0=b0,b1=b1,set0=set0,set1=set1,causal=c(set0,unlist(set1))))
}

