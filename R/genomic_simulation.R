####################################################################################################################
#    Module 8: Genomic simulation
####################################################################################################################
#'
#' Genomic simulation
#' 
#'
#' Simulate Genotype and Phenotype Data
#'
#' This function simulates genotype and phenotype data based on the `Glist`, which is 
#' information about the genotype matrix. 
#'
#' @param Glist A list of information about the genotype matrix. Default is `NULL`.
#' @param chr The chromosome(s) being used in the simulation. Default is 1.
#' @param nt Number of traits. Default is 1.
#' @param W Matrix of centered and scaled genotypes. Default is `NULL`.
#' @param n Number of individuals. Default is 1000.
#' @param m Number of markers. Default is 1000.
#' @param rsids A character vector of rsids. Default is `NULL`.
#'
#' @return 
#' A list containing:
#' \itemize{
#'   \item \code{y}: Phenotypes.
#'   \item \code{W}: Matrix of centered and scaled genotypes.
#'   \item \code{e}: Errors.
#'   \item \code{g}: Genotype effect.
#'   \item \code{b0}, \code{b1}: Coefficients.
#'   \item \code{set0}, \code{set1}: Selected markers.
#'   \item \code{causal}: Causal markers.
#' }
#' 
#' @examples
#' ## Plink bed/bim/fam files
#' bedfiles <- system.file("extdata", paste0("sample_chr",1:2,".bed"), package = "qgg")
#' bimfiles <- system.file("extdata", paste0("sample_chr",1:2,".bim"), package = "qgg")
#' famfiles <- system.file("extdata", paste0("sample_chr",1:2,".fam"), package = "qgg")
#'
#' # Summarize bed/bim/fam files
#' Glist <- gprep(study="Example", bedfiles=bedfiles, bimfiles=bimfiles, famfiles=famfiles)
#'
#' # Simulate phenotype
#' sim <- gsim(Glist=Glist, chr=1, nt=1)
#' head(sim$y)
#' head(sim$e)
#' head(sim$causal)
#'
#' @author Peter Soerensen
#' @export

gsim <- function(Glist=NULL, chr=1, nt=1,W=NULL, n=1000, m=1000, rsids=NULL) {
  
  
  if(!is.null(Glist)) {
    if(is.null(rsids)) rsids <- Glist$rsidsLD[[chr]]
    W <- getG(Glist=Glist, chr=chr, rsids=rsids)
  }
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
  b0 <- sample(c(0.25,-0.25,0.5,-0.5),2)
  for (i in 1:nt){
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
  if(nt==1) return( list(y=y[[1]],W=W,e=e[[1]],g=g,b0=b0,b1=b1,set0=set0,set1=set1,causal=c(set0,unlist(set1))))
  if(nt>1) return( list(y=as.matrix(as.data.frame(y)),W=W,e=as.matrix(as.data.frame(e)),g=g,b0=b0,b1=b1,set0=set0,set1=set1,causal=c(set0,unlist(set1))))
}

