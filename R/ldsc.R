####################################################################################################################
#    Module 7: LDSC
####################################################################################################################
#'
#' LD score regression
#' @description
#' The ldsc function is used for LDSC analysis
#'
#' @param Glist list of information about genotype matrix stored on disk
#' @param lscore vector of LD scores (optional)
#' @param z matrix of z statistics for n traits
#' @param b matrix of marker effects for n traits
#' @param seb matrix of standard errors of marker effects for n traits

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
#' @author Palle Duun Rohde

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

ldsc <- function(Glist=NULL, lscore=NULL, z=NULL, b=NULL, seb=NULL, p=NULL, n=NULL, sets=NULL, intercept=TRUE, what="h2", SE.h2=FALSE, SE.rg=FALSE, blk=200) {
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
               z2 <- z[,t]**2
               z2 <- z2[!is.na(z2)] 
               z2 <- z2[z2<maxZ2]
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
          if(SE.h2==TRUE){
               if(intercept==TRUE){
                    P <- SE <- matrix(0,ncol=1,nrow=nt)
                    for ( t in 1:nt) {
                         z2 <- z[,t]**2
                         z2 <- z2[!is.na(z2)]
                         z2 <- z2[z2<maxZ2]
                         X <- cbind(1,n[t]*lscore[names(z2)]/length(z2))
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
                    h2[h2[,2]>1,2] <- h2[h2[,2]>1,3] <-  h2[h2[,2]>1,4] <- NA
                    result <- h2
               }
               if(intercept==FALSE){
                    P <- SE <- matrix(0,ncol=1,nrow=nt)
                    for ( t in 1:nt) {
                         z2 <- z[,t]**2
                         z2 <- z2[!is.na(z2)]
                         z2 <- z2[z2<maxZ2]
                         X <- matrix(n[t]*lscore[names(z2)]/length(z2),ncol=1)
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
                    h2[h2[,1]>1,1] <- h2[h2[,1]>1,2] <-  h2[h2[,1]>1,3] <- NA
                    result <- h2
               }
          }
         

     
     if(what=="correlation") {
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
                    X <- cbind(1,sqrt(n[t1])*sqrt(n[t2])*lscore[rws]/m)
                    y <- Z1[rws]*Z2[rws]
                    XtX <- crossprod(X)
                    Xy <- crossprod(X,y)
                    if(intercept & !any(is.na(h2[c(t1,t2),2]))) rg[t1,t2] <- solve(XtX, Xy)[2]/(sqrt(h2[t1,2])*sqrt(h2[t2,2]))
                    if(!intercept & !any(is.na(h2[c(t1,t2)]))) rg[t1,t2] <- solve(XtX, Xy)[2]/(sqrt(h2[t1])*sqrt(h2[t2]))
               }
               if(intercept) rg[t1,t1] <- h2[t1,2]
               if(!intercept) rg[t1,t1] <- h2[t1]
               rg[t2,t1] <- rg[t1,t2]
          }
          rownames(rg) <- colnames(rg) <- colnames(z)
          result <- rg
          
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
                         
                         X1 <- cbind(1,n[t1]*lscore[rws]/sum(rws))
                         y1 <- Z2_1[rws]
                         X2 <- cbind(1,n[t2]*lscore[rws]/sum(rws))
                         y2 <- Z2_2[rws]
                         
                         XtX1 <- crossprod(X1)
                         Xy1 <- crossprod(X1,y1)
                         XtX2 <- crossprod(X2)
                         Xy2 <- crossprod(X2,y2)
                         
                         X <- cbind(1,sqrt(n[t1])*sqrt(n[t2])*lscore[rws]/m)
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
                              
                              X1 <- matrix(n[t1]*lscore[rws]/sum(rws),ncol=1)
                              y1 <- Z2_1[rws]
                              X2 <- matrix(n[t2]*lscore[rws]/sum(rws),ncol=1)
                              y2 <- Z2_2[rws]
                              
                              XtX1 <- crossprod(X1)
                              Xy1 <- crossprod(X1,y1)
                              XtX2 <- crossprod(X2)
                              Xy2 <- crossprod(X2,y2)
                              
                              X <- matrix(1,sqrt(n[t1])*sqrt(n[t2])*lscore[rws]/m,ncol=1)
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
     return(result)
}


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

neff <- function(seb=NULL,af=NULL,Vy=1) {
     seb2 <- seb**2
     vaf <- 2*af*(1-af)
     neff <- round(median(Vy/(vaf*seb2)))
     return(neff)
}


