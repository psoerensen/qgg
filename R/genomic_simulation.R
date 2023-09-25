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


#' Simulate Genetic Data Based on Given Parameters
#'
#' This function simulates phenotype data by random sampling of markers available on `Glist`.
#' Default parameters for the simulated phenotype reflect the genetic architecture assumed by BayesC prior(Habier et al., 2011).
#' This function is under active development.
#'
#' @param Glist A list containing genetic data. If NULL, the function will stop with an error.
#' @param h2 Heritability. If NULL, heritability of 0.5 is assumed.
#' @param m Number of causal markers. The values for either `m` or `prp.cau` should be provided at any given time. 
#' If the list of quality controlled markers is not available then list of raw markers is used.  
#' If `m` is NULL and `prp.cau` is also NULL, `prp.cau` will default to 0.001.
#' @param prp.cau Proportion of causal markers. The values for either `m` or `prp.cau` should be provided at any given time. 
#' @param n Number of individuals randomly sampled from `Glist`. If NULL, all the individuals on `Glist` is used.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{y}: Vector of simulated phenotypes.
#'   \item \code{g}: Vector of simulated genetic values.
#'   \item \code{e}: Vector of simulated residual effects.
#'   \item \code{b}: Vector of effect sizes of the simulated causal markers.
#'   \item \code{causal}: Vector of ids for the simulated causal markers.
#'   \item \code{h2}: Estimated heritability of the simulated phenotype.
#' }
#'
#' @examples
#' @author Peter Soerensen 
#' @author Merina Shrestha
#' @keywords internal
#' @export

gsimC <- function (Glist = NULL, h2 = NULL, m = NULL, prp.cau = NULL, n = NULL)
{if(is.null(Glist)){stop("Error:Glist is NULL")}
  y <- g <- e <- NULL
  if(is.null(h2)){hsnp = 0.5} else {hsnp = h2}
  print(paste("heritability :",hsnp,sep=""))
  if(is.null(m) && is.null(prp.cau)) {
    prp.cau = 0.001
  } else if(!is.null(m) && is.null(prp.cau)) {
    snp.cau = m 
  } else if(is.null(m) && !is.null(prp.cau)){
    if(!is.null(Glist$rsidsLD)) { 
      tot.snps <- length(unlist(Glist$rsidsLD))   
    } else if(is.null(Glist$rsidsLD)) {
      tot.snps <- length(unlist(Glist$rsids))  
    }
    snp.cau <- tot.snps * prp.cau  
  }
  print(paste("Number of causal snps :",snp.cau,sep=""))
  pop.var <- hsnp/snp.cau 
  pop.sd <- sqrt(pop.var) 
  b0 <- rnorm(n=snp.cau, mean=0, sd=pop.sd)
  if(!is.null(Glist$rsidsLD)) { 
    rsidsLD <- unlist(Glist$rsidsLD)
    print(paste("QC'd SNPs used")) 
  } else if(is.null(Glist$rsidsLD)){
    rsidsLD <- unlist(Glist$rsids)
    print(paste("Raw SNPs used"))
  }
  rsids.cau <- sample(rsidsLD,snp.cau)
  if(is.null(n)){
    ids.1 <- NULL
    print(paste("Number of samples:",length(Glist$ids)))
  } else {
    ids.1 <- sample(x=Glist$ids, size=n, replace = FALSE) 
    print(paste("Number of samples:",n))
  }
  var.x = 1
  tpb <- txtProgressBar(min=0, max=length(rsids.cau), style=3)
  for (snp.num in 1:length(rsids.cau)){  
    setTxtProgressBar(tpb,snp.num) 
    w.snp <- rsids.cau[snp.num]
    for(chr.ms in 1:length(Glist$bedfiles)){ 
      for (j in 1:length(Glist$rsids[[chr.ms]])){ 
        if(w.snp == Glist$rsids[[chr.ms]][j]){
          cls <- j
          w <- NULL
          w <- getG(Glist = Glist, chr = chr.ms, bedfiles = NULL, bimfiles = NULL,
                    famfiles = NULL, ids = ids.1, rsids = NULL, rws = NULL, cls = cls,
                    impute = TRUE, scale = TRUE)
          b1 <- as.matrix(b0[snp.num])
          g1 <- w %*% b1
          ifelse( var.x == 1, g1.1 <- g1, g1.1 <- g1.1 + g1) # this is the way to add matrix
          var.x <- var.x + 1
        }}}} 
  rm(var.x)
  close(tpb) 
  g <- cbind(g, g1.1)
  res.var <- var(g)*((1/hsnp)-1) 
  res.sd <- sqrt(res.var)
  if(is.null(n)) {
    if(length(Glist$ids) != 0){
      e <- rnorm(length(Glist$ids), mean = 0, sd = res.sd)
      y <- as.vector(g + e)
      names(y) <- Glist$ids 
    } else {stop("Error: no ids in Glist")}
  } else if(!is.null(n)){
    e <- rnorm(n, mean = 0, sd = res.sd) 
    y <- as.vector(g + e)
    names(y) <- ids.1   
  }
  return(list(y = y, g = g, e = e, b = b0, causal = rsids.cau, h2 = round(var(g)/(var(g)+var(e)),2)))
} 


#' Simulate Genetic Data Based on Given Parameters
#'
#' This function simulates phenotype data by random sampling of markers available on `Glist`.
#' Default parameters for the simulated phenotype reflect the genetic architecture assumed by BayesR prior(Erbe et al., 2012).
#' Marker effect are sampled from mixture distributions leading to small, moderate and large effect sizes. 
#' This function is under active development.
#' 
#' @param Glist A list containing genetic data. If NULL, the function will stop with an error.
#' @param h2 Heritability. If NULL, heritability of 0.5 is assumed.
#' @param m Number of causal markers. The values for either `m` or `prp.cau` should be provided at any given time. 
#' If the list of quality controlled markers is not available then list of raw markers is used.  
#' If `m` is NULL and `prp.cau` is also NULL, `prp.cau` will default to 0.001.
#' @param prp.cau Proportion of causal markers. The values for either `m` or `prp.cau` should be provided at any given time. 
#' @param n Number of individuals randomly sampled from `Glist`. If NULL, all the individuals on `Glist` is used.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{y}: Vector of simulated phenotypes.
#'   \item \code{g}: Vector of simulated genetic values.
#'   \item \code{e}: Vector of simulated residual effects.
#'   \item \code{b}: list of vectors for effect sizes of the simulated causal markers in three classes.
#'   \item \code{causal}: list of vectors for ids for the simulated causal markers in three classes.
#'   \item \code{h2}: Estimated heritability of the simulated phenotype.
#' }
#' 
#' @examples
#' @author Peter Soerensen 
#' @author Merina Shrestha
#' @keywords internal
#' @export

gsimR <- function (Glist = NULL, h2 = NULL, m = NULL, prp.cau = NULL, n = NULL)  
{if(is.null(Glist)){stop("Error:Glist is NULL")}
  y <- g <- e <- NULL
  if(is.null(h2)){hsnp = 0.5} else {hsnp = h2}
  print(paste("heritability :",hsnp,sep=""))
  if(is.null(m) && is.null(prp.cau)) {
    prp.cau = 0.001
    m = NULL
  }
  if(!is.null(m) && is.null(prp.cau)) {
    snp.cau = m
    prp.cau = NULL
  }
  if(is.null(m) && !is.null(prp.cau)){
    if(!is.null(Glist$rsidsLD)) { 
      tot.snps <- length(unlist(Glist$rsidsLD))} else if(is.null(Glist$rsidsLD)) {tot.snps <- length(unlist(Glist$rsids))}
    snp.cau <- tot.snps * prp.cau  #### number of causal snps
  }
  print(paste("Number of causal snps :",snp.cau,sep=""))
  
  snp.cau.1 <- round(snp.cau * 0.94, digits=0)
  snp.cau.2 <- round(snp.cau * 0.05, digits=0)
  snp.cau.3 <- round(snp.cau * 0.01, digits=0)
  if(is.null(snp.cau.1)){snp.cau.1 <- 1}
  if(is.null(snp.cau.2)){snp.cau.2 <- 1}
  if(is.null(snp.cau.3)){snp.cau.3 <- 1}
  lapply(1:3, function(i) print(paste("Number of causal snps in class", i, ":", get(paste("snp.cau.", i, sep="")))))
  
  hsnp.1 <- (snp.cau.1 * 0.001)/hsnp
  hsnp.2 <- (snp.cau.2 * 0.01)/hsnp
  hsnp.3 <- (snp.cau.3 * 0.1)/hsnp
  
  if(is.null(Glist$rsidsLD)){
    rsidsLD <- unlist(Glist$rsids)
  } else { rsidsLD <- unlist(Glist$rsidsLD)} 
  rsids.cau.1 <- sample(rsidsLD,snp.cau.1) 
  rsidsLD <- rsidsLD[!(rsidsLD %in% rsids.cau.1)]
  rsids.cau.2 <- sample(rsidsLD,snp.cau.2) 
  rsidsLD <- rsidsLD[!(rsidsLD %in% rsids.cau.2)]
  rsids.cau.3 <- sample(rsidsLD,snp.cau.3) 
  if(is.null(n)){
    ids.1 <- NULL
    print(paste("Number of samples:",length(Glist$ids)))
  } else {
    ids.1 <- sample(x=Glist$ids, size=n, replace = FALSE)
    print(paste("Number of samples:",n))
  }
  
  g.eff <- function(h_snp, snp_cau, rsids_cau, Glist) {
    pop.var <- h_snp/snp_cau
    pop.sd <- sqrt(pop.var)
    b <- rnorm(n=snp_cau, mean=0, sd=pop.sd)
    g_final <- 0
    var.x = 1
    tpb <- txtProgressBar(min=0, max=length(rsids_cau), style=3)
    for (snp.num in 1:length(rsids_cau)){ 
      setTxtProgressBar(tpb,snp.num)
      w.snp <- rsids_cau[snp.num]
      for(chr.ms in 1:length(Glist$bedfiles)){ 
        for (j in 1:length(Glist$rsids[[chr.ms]])){  
          if(w.snp == Glist$rsids[[chr.ms]][j]){
            cls <- j
            w <- getG(Glist = Glist, chr = chr.ms, bedfiles = NULL, bimfiles = NULL, 
                      famfiles = NULL, ids = ids.1, rsids = NULL, rws = NULL, cls = cls, 
                      impute = TRUE, scale = TRUE) 
            b_temp <- as.matrix(b[snp.num])  
            g <- w %*% b_temp 
            ifelse(var.x == 1, g_final <- g, g_final <- g_final + g)
            var.x = var.x + 1
          }}}}
    rm(var.x)
    close(tpb)
    return(list(g_final = g_final, b = b))
  }
  print("Class 1")
  res1 <- g.eff(hsnp.1, snp.cau.1, rsids.cau.1, Glist)
  print("Class 2")
  res2 <- g.eff(hsnp.2, snp.cau.2, rsids.cau.2, Glist)
  print("Class 3")
  res3 <- g.eff(hsnp.3, snp.cau.3, rsids.cau.3, Glist)
  g <- cbind(g, res1$g_final + res2$g_final + res3$g_final) 
  res.var <- var(g) * ((1/hsnp)-1) 
  res.sd <- sqrt(res.var)
  if(is.null(n)) {
    if(length(Glist$ids) != 0){
      e <- rnorm(length(Glist$ids), mean = 0, sd = res.sd)
      y <- as.vector(g + e)
      names(y) <- Glist$ids 
    } else {stop("Error: no ids in Glist")}
  } else if(!is.null(n)){
    e <- rnorm(n, mean = 0, sd = res.sd)
    y <- as.vector(g + e)
    names(y) <- ids.1   
  }
  return(list(y = y, e = e, g = g, b = list(class1=res1$b,class2=res2$b,class3=res3$b), 
              causal=list(class1=rsids.cau.1,class2=rsids.cau.2,class3=rsids.cau.3),
              h2 = round(var(g)/(var(g)+var(e)),2)))
}

