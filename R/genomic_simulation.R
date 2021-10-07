####################################################################################################################
#    Module 8: GSIM
####################################################################################################################
#'
#' LD score regression
#' @description
#' The ldsc function is used for LDSC analysis
#'
#' @param Glist list of information about genotype matrix stored on disk
#' @param ldscores vector of LD scores (optional as LD scores are stored within Glist)
#' @param z matrix of z statistics for n traits
#' @param b matrix of marker effects for n traits if z matrix not is given
#' @param seb matrix of standard errors of marker effects for n traits if z matrix not is given
#' @param n vector of sample sizes for the traits (element i corresponds to column vector i in z matrix)
#' @param interept logical if TRUE the LD score regression includes intercept
#' @param what either computation of heritability (what="h2") or genetic correlation between traits (what="rg")
#' @param SE.h2 logical if TRUE standard errors and significance for the heritability estimates are computed using a block jackknife approach
#' @param SE.rg logical if TRUE standard errors and significance for the genetic correlations are computed using a block jackknife approach
#' @param blk numeric size of the blocks used in the jackknife estimation of standard error (default = 200)
#'
#' @return Returns a matrix of heritability estimates when what="h2", and if SE.h2=TRUE standard errors (SE) and significance levels (P) are returned. 
#'         If what="rg" an n-by-n matrix of correlations is returned where the diagonal elements being h2 estimates. 
#'         If SE.rg=TRUE a list is returned with n-by-n matrices of genetic correlations, estimated standard errors and significance levels.
#'
#' @author Peter Soerensen
#' @author Palle Duun Rohde
#'
#' @examples
#'
#' #Simulate data
#' W1 <- getG(Glist, chr=1, scale=TRUE)
#' W2 <- getG(Glist, chr=2, scale=TRUE)
#'
#' W <- cbind(W1,W2)
#' causal <- sample(1:ncol(W),5)
#'
#' b1 <- rnorm(length(causal))
#' b2 <- rnorm(length(causal))
#' y1 <- W[, causal]%*%b1 + rnorm(nrow(W))
#' y2 <- W[, causal]%*%b2 + rnorm(nrow(W))
#'
# # Create model
#' data1 <- data.frame(y = y1, mu = 1)
#' data2 <- data.frame(y = y2, mu = 1)
#' X1 <- model.matrix(y ~ 0 + mu, data = data1)
#' X2 <- model.matrix(y ~ 0 + mu, data = data2)
#'
#' # Linear model analyses and single marker association test
#' maLM1 <- lma(y=y1, X=X1,W = W)
#' maLM2 <- lma(y=y2,X=X2,W = W)
#' 
#' # Compute heritability and genetic correlations for trait 1 and 2
#' z1 <- maLM1[,"stat"]
#' z2 <- maLM2[,"stat"]
#'
#'z <- cbind(z1=z1,z2=z2)
#'
#'h2 <- ldsc(Glist, z=z, n=c(500,500), what="h2")
#'rg <- ldsc(Glist, z=z, n=c(500,500), what="rg")
#'
#' @export 

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

