####################################################################################################################
#    Module 7: LDSC
####################################################################################################################
#'
#' LD score regression
#' @description
#' The ldsc function is used for LDSC analysis
#'
#' @param Glist list of information about genotype matrix stored on disk
#' @param stat dataframe with marker summary statistics
#' @param ldscores vector of LD scores (optional as LD scores are stored within Glist)
#' @param z matrix of z statistics for n traits
#' @param b matrix of marker effects for n traits if z matrix not is given
#' @param seb matrix of standard errors of marker effects for n traits if z matrix not is given
#' @param af vector of allele frequencies 
#' @param n vector of sample sizes for the traits (element i corresponds to column vector i in z matrix)
#' @param intercept logical if TRUE the LD score regression includes intercept
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
#'
#' # Plink bed/bim/fam files
#'  #bedfiles <- system.file("extdata", paste0("sample_chr",1:2,".bed"), package = "qgg")
#'  #bimfiles <- system.file("extdata", paste0("sample_chr",1:2,".bim"), package = "qgg")
#'  #famfiles <- system.file("extdata", paste0("sample_chr",1:2,".fam"), package = "qgg")
#'  #
#'  ## Summarize bed/bim/fam files
#'  #Glist <- gprep(study="Example", bedfiles=bedfiles, bimfiles=bimfiles, famfiles=famfiles)
#'
#'  #
#'  ## Filter rsids based on MAF, missingness, HWE
#'  #rsids <-  gfilter(Glist = Glist, excludeMAF=0.05, excludeMISS=0.05, excludeHWE=1e-12) 
#'  #
#'  ## Compute sparse LD (msize=size of LD window)
#'  ##ldfiles <- system.file("extdata", paste0("sample_chr",1:2,".ld"), package = "qgg")
#'  ##Glist <- gprep(Glist, task="sparseld", msize=200, rsids=rsids, ldfiles=ldfiles, overwrite=TRUE)
#'  #
#'  #
#'  ##Simulate data
#'  #W1 <- getG(Glist, chr=1, scale=TRUE)
#'  #W2 <- getG(Glist, chr=2, scale=TRUE)
#'
#'  #W <- cbind(W1,W2)
#'  #causal <- sample(1:ncol(W),5)
#'
#'  #b1 <- rnorm(length(causal))
#'  #b2 <- rnorm(length(causal))
#'  #y1 <- W[, causal]%*%b1 + rnorm(nrow(W))
#'  #y2 <- W[, causal]%*%b2 + rnorm(nrow(W))
#'
# # Create model
#'  #data1 <- data.frame(y = y1, mu = 1)
#'  #data2 <- data.frame(y = y2, mu = 1)
#'  #X1 <- model.matrix(y ~ 0 + mu, data = data1)
#'  #X2 <- model.matrix(y ~ 0 + mu, data = data2)
#'
#'  ## Linear model analyses and single marker association test
#'  #maLM1 <- lma(y=y1, X=X1,W = W)
#'  #maLM2 <- lma(y=y2,X=X2,W = W)
#'  #
#'  ## Compute heritability and genetic correlations for trait 1 and 2
#'  #z1 <- maLM1[,"stat"]
#'  #z2 <- maLM2[,"stat"]
#'
#'  #z <- cbind(z1=z1,z2=z2)
#'
#'  #h2 <- ldsc(Glist, z=z, n=c(500,500), what="h2")
#'  #rg <- ldsc(Glist, z=z, n=c(500,500), what="rg")
#'
#'
#'
#' @export 

ldsc <- function(Glist=NULL, ldscores=NULL, z=NULL, b=NULL, seb=NULL, af=NULL, stat=NULL, tol=1e-8,
                 n=NULL, intercept=TRUE, what="h2", maxZ2=NULL, SE.h2=FALSE, SE.rg=FALSE, blk=200) {
  
  if(!is.null(Glist) & is.null(ldscores) ) ldscores <- unlist(Glist$ldscores)
  ldscores <- unlist(ldscores)
  
  if(any(is.na(ldscores))) stop("Missing values in ldscores")
  if(is.null(names(ldscores))) stop("Missing names in ldscores")
  
  if(!is.null(stat)) {
    if(is.data.frame(stat)) {
      z <- as.matrix(stat$b/stat$seb)
      rownames(z) <- stat$marker
      if(is.null(n)) {
        if("n"%in%colnames(stat)) n <- mean(stat$n,na.rm=T)
        if(!"n"%in%colnames(stat)) n <- neff(seb=stat$seb, af=stat$eaf)
      }
      nt <- 1
    }  
  }
  
  if(!is.null(z)) nt <- ncol(z)
  
  if(!is.null(b)) {
    nt <- ncol(b)
    for (t in 1:nt) {
      z <- cbind(z,(b[,t]/seb[,t]))
    }
    colnames(z) <- colnames(b)
    rownames(z) <- rownames(b)
  }  
  z <- as.matrix(z[rownames(z)%in%names(ldscores),,drop=FALSE])
  ldscores <- ldscores[rownames(z)]
  if(is.null(n)) {
    n <- NULL
    if(is.null(seb)) stop("Please provide n or alternatively seb")
    if(is.null(af)) stop("Please provide af")
    for ( t in 1:nt) {
      n <- c(n,neff(seb[,t],af[,t]))
    }
  }
  if(is.null(maxZ2)) maxZ2 <- max(0.001 * max(n), 80)
  h2 <- NULL
  for ( t in 1:nt) {
    z2 <- z[,t]**2
    z2 <- z2[!is.na(z2)] 
    z2 <- z2[z2<maxZ2]
    if(intercept) X <- cbind(1,n[t]*ldscores[names(z2)]/length(z2))
    if(!intercept) X <- matrix(n[t]*ldscores[names(z2)]/length(z2),ncol=1)
    y <- z2
    XtX <- crossprod(X)
    Xy <- crossprod(X,y)
    h2_ldsc <- solve(XtX, Xy)
    if(intercept){
      if(is.na(h2_ldsc[2])) h2_ldsc[2] <- 0
      if(h2_ldsc[2]<0) h2_ldsc[2] <- 0
      if(h2_ldsc[2]>1) h2_ldsc[2] <- 1
    }
    if(!intercept){
      if(is.na(h2_ldsc[1])) h2_ldsc[1] <- 0
      if(h2_ldsc[1]<0) h2_ldsc[1] <- 0
      if(h2_ldsc[1]>1) h2_ldsc[1] <- 1
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
  if(intercept)  result <- h2[,2]
  isNA <- result<=tol | is.na(result)
  result[isNA] <- tol
  
  
  #---------------------------------#
  # Block Jackknife to estimate h2 SE
  if(SE.h2==TRUE){
    if(intercept==TRUE){
      P <- SE <- matrix(0,ncol=1,nrow=nt)
      for ( t in 1:nt) {
        z2 <- z[,t]**2
        z2 <- z2[!is.na(z2)]
        z2 <- z2[z2<maxZ2]
        X <- cbind(1,n[t]*ldscores[names(z2)]/length(z2))
        y <- z2
        XtX <- crossprod(X)
        Xy <- crossprod(X,y)
        
        idx <- split(1:length(y), cut(1:length(y), blk, labels=F))
        
        h2.jack <- numeric(blk)
        for(i in 1:length(idx)){
          XtX.idx <- crossprod(X[idx[[i]],])
          Xy.idx <- crossprod(X[idx[[i]],],y[idx[[i]]])
          h2.jack[i] <- solve(XtX-XtX.idx, Xy-Xy.idx)[2]
        }
        SE[t,1] <- sqrt(mean((h2.jack - mean(h2.jack))^2)*(blk - 1))
        P[t,1] <- pchisq((h2[t,2]/SE[t,])^2, df = 1, lower.tail = FALSE)
      }
      h2 <- cbind(h2,SE,P)
      colnames(h2)[3:4] <- c("SE","P")
      h2[h2[,2]<0,2] <- h2[h2[,2]<0,3] <-  h2[h2[,2]<0,4] <- NA
      result <- h2
    }
    if(intercept==FALSE){
      P <- SE <- matrix(0,ncol=1,nrow=nt)
      for ( t in 1:nt) {
        z2 <- z[,t]**2
        z2 <- z2[!is.na(z2)]
        z2 <- z2[z2<maxZ2]
        X <- matrix(n[t]*ldscores[names(z2)]/length(z2),ncol=1)
        y <- z2
        XtX <- crossprod(X)
        Xy <- crossprod(X,y)
        
        idx <- split(1:length(y), cut(1:length(y), blk, labels=F))
        
        h2.jack <- numeric(blk)
        for(i in 1:length(idx)){
          XtX.idx <- crossprod(X[idx[[i]],])
          Xy.idx <- crossprod(X[idx[[i]],],y[idx[[i]]])
          h2.jack[i] <- solve(XtX-XtX.idx, Xy-Xy.idx)[1]
        }
        SE[t,1] <- sqrt(mean((h2.jack - mean(h2.jack))^2)*(blk - 1))
        P[t,1] <- pchisq((h2[t,1]/SE[t,])^2, df = 1, lower.tail = FALSE)
      }
      h2 <- cbind(h2,SE,P)
      colnames(h2)[2:3] <- c("SE","P")
      h2[h2[,1]<0,1] <- h2[h2[,1]<0,2] <-  h2[h2[,1]<0,3] <- NA
      result <- h2
    }
  }
  
  if(what=="rg") {
    rg <- matrix(0,nt,nt)
    rownames(rg) <- colnames(rg) <- colnames(z)
    for (t1 in 1:nt) {
      for (t2 in t1:nt) {
        Z1 <- z[,t1]
        Z2 <- z[,t2]
        Z2_1 <- Z1**2
        Z2_2 <- Z2**2
        rws1 <- Z2_1<maxZ2
        rws2 <- Z2_2<maxZ2
        rws <- rws1 & rws2
        m <- sum(rws)
        X <- cbind(1,sqrt(n[t1])*sqrt(n[t2])*ldscores[rws]/m)
        y <- Z1[rws]*Z2[rws]
        XtX <- crossprod(X)
        Xy <- crossprod(X,y)
        if(intercept && !any(is.na(h2[c(t1,t2),2]))) rg[t1,t2] <- solve(XtX, Xy)[2]/(sqrt(h2[t1,2])*sqrt(h2[t2,2]))
        #if(!intercept & !any(is.na(h2[c(t1,t2)]))) rg[t1,t2] <- solve(XtX,Xy)[2]/(sqrt(h2[t1])*sqrt(h2[t2]))
        if(!intercept && !any(is.na(h2[c(t1,t2)]))) rg[t1,t2] <- (Xy[2]/XtX[2,2])/(sqrt(h2[t1])*sqrt(h2[t2]))
      }
      if(intercept) rg[t1,t1] <- h2[t1,2]
      if(!intercept) rg[t1,t1] <- h2[t1]
      #rg[t2,t1] <- rg[t1,t2]
    }
    rownames(rg) <- colnames(rg) <- colnames(z)
    result <- NULL
    result$h2 <- diag(rg)
    for (i in 1:ncol(rg)) {
      for (j in i:ncol(rg)) {
        rg[j,i] <- rg[i,j]
      }
    }
    diag(rg) <- 1 
    result$rg <- rg
    result$rg[result$rg > 1] <- 1
    result$rg[result$rg < -1] <- -1
    
    isNA <- result$h2<=tol | is.na(result$h2)
    result$h2[isNA] <- tol
    result$rg[isNA,] <- 0
    result$rg[,isNA] <- 0
    diag(result$rg) <- 1
    
    
    #---------------------------------#
    # Block Jackknife to estimate rg SE
    if(SE.rg==TRUE){
      if(intercept==TRUE){
        P <- SE <- matrix(NA,ncol=nt,nrow=nt)
        rownames(SE) <- colnames(SE) <- rownames(P) <- colnames(P) <- colnames(z)
        
        for (t1 in 1:nt) {
          for (t2 in t1:nt) {
            if(t1!=t2){
              Z1 <- z[,t1]
              Z2 <- z[,t2]
              Z2_1 <- Z1**2
              Z2_2 <- Z2**2
              rws1 <- Z2_1<maxZ2
              rws2 <- Z2_2<maxZ2
              rws <- rws1 & rws2
              m <- sum(rws)
              
              X1 <- cbind(1,n[t1]*ldscores[rws]/sum(rws))
              y1 <- Z2_1[rws]
              X2 <- cbind(1,n[t2]*ldscores[rws]/sum(rws))
              y2 <- Z2_2[rws]
              
              XtX1 <- crossprod(X1)
              Xy1 <- crossprod(X1,y1)
              XtX2 <- crossprod(X2)
              Xy2 <- crossprod(X2,y2)
              
              X <- cbind(1,sqrt(n[t1])*sqrt(n[t2])*ldscores[rws]/m)
              y <- Z1[rws]*Z2[rws]
              XtX <- crossprod(X)
              Xy <- crossprod(X,y)
              
              idx <- split(1:length(y), cut(1:length(y), blk, labels=F))
              
              rg.jack <- numeric(blk)
              for(i in 1:length(idx)){
                
                XtX.idx1 <- crossprod(X1[idx[[i]],])
                Xy.idx1 <- crossprod(X1[idx[[i]],],y1[idx[[i]]])
                h2.1 <- solve(XtX1-XtX.idx1, Xy1-Xy.idx1)[2]
                h2.1[h2.1<0] <- 0.0001
                
                XtX.idx2 <- crossprod(X2[idx[[i]],])
                Xy.idx2 <- crossprod(X2[idx[[i]],],y2[idx[[i]]])
                h2.2 <- solve(XtX2-XtX.idx2, Xy2-Xy.idx2)[2]
                h2.2[h2.2<0] <- 0.0001
                
                XtX.idx <- crossprod(X[idx[[i]],])
                Xy.idx <- crossprod(X[idx[[i]],],y[idx[[i]]])
                rg.jack[i] <- solve(XtX-XtX.idx, Xy-Xy.idx)[2]/(sqrt(h2.1)*sqrt(h2.2))
              }
              
              SE[t1,t2] <- sqrt(mean((rg.jack - mean(rg.jack))^2)*(blk - 1))
              P[t1,t2] <- pchisq((rg[t1,t2]/SE[t1,t2])^2, df = 1, lower.tail = FALSE)
              
              SE[t2,t1] <- SE[t1, t2]
              P[t2,t1] <- P[t1, t2]
              message("SE for trait combination: ",t1, "-", t2, " completed")
            }
          }
        }
        return(list(rg=result,SE=SE, P=P))                         
      }
      
      if(intercept==FALSE){
        P <- SE <- matrix(NA,ncol=nt,nrow=nt)
        rownames(SE) <- colnames(SE) <- rownames(P) <- colnames(P) <- colnames(b)
        
        for (t1 in 1:nt) {
          for (t2 in t1:nt) {
            if(t1!=t2){
              Z1 <- z[,t1]
              Z2 <- z[,t2]
              Z2_1 <- Z1**2
              Z2_2 <- Z2**2
              rws1 <- Z2_1<maxZ2
              rws2 <- Z2_2<maxZ2
              rws <- rws1 & rws2
              m <- sum(rws)
              
              X1 <- matrix(n[t1]*ldscores[rws]/sum(rws),ncol=1)
              y1 <- Z2_1[rws]
              X2 <- matrix(n[t2]*ldscores[rws]/sum(rws),ncol=1)
              y2 <- Z2_2[rws]
              
              XtX1 <- crossprod(X1)
              Xy1 <- crossprod(X1,y1)
              XtX2 <- crossprod(X2)
              Xy2 <- crossprod(X2,y2)
              
              X <- matrix(1,sqrt(n[t1])*sqrt(n[t2])*ldscores[rws]/m,ncol=1)
              y <- Z1[rws]*Z2[rws]
              XtX <- crossprod(X)
              Xy <- crossprod(X,y)
              
              idx <- split(1:length(y), cut(1:length(y), blk, labels=F))
              
              rg.jack <- numeric(blk)
              for(i in 1:length(idx)){
                
                XtX.idx1 <- crossprod(X1[idx[[i]],])
                Xy.idx1 <- crossprod(X1[idx[[i]],],y1[idx[[i]]])
                h2.1 <- solve(XtX1-XtX.idx1, Xy1-Xy.idx1)[1]
                h2.1[h2.1<0] <- 0.0001
                
                XtX.idx2 <- crossprod(X2[idx[[i]],])
                Xy.idx2 <- crossprod(X2[idx[[i]],],y2[idx[[i]]])
                h2.2 <- solve(XtX2-XtX.idx2, Xy2-Xy.idx2)[1]
                h2.2[h2.2<0] <- 0.0001
                
                XtX.idx <- crossprod(X[idx[[i]],])
                Xy.idx <- crossprod(X[idx[[i]],],y[idx[[i]]])
                rg.jack[i] <- solve(XtX-XtX.idx, Xy-Xy.idx)[1]/(sqrt(h2.1)*sqrt(h2.2))
              }
              
              SE[t1,t2] <- sqrt(mean((rg.jack - mean(rg.jack))^2)*(blk - 1))
              P[t1,t2] <- pchisq((rg[t1,t2]/SE[t1,t2])^2, df = 1, lower.tail = FALSE)
              
              SE[t2,t1] <- SE[t1, t2]
              P[t2,t1] <- P[t1, t2]
              
            }
          }
        }
        return(list(rg=result,SE=SE, P=P))     
      }
    }
  }
  
  
  return(result)
}




neff <- function(seb=NULL,af=NULL,Vy=1) {
  seb2 <- seb**2
  vaf <- 2*af*(1-af)
  neff <- round(median(Vy/(vaf*seb2)))
  return(neff)
}


