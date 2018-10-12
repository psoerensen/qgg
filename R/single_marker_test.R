####################################################################################################################
#    Module 5: LMM marker association test
####################################################################################################################
#'
#' Linear mixed model marker association test
#'
#' @description
#' Marker test using linear mixed model (LMM) to test for an association of single markers with a phenotype.
#'
#' @details
#' Linear mixed model single marker association (LMMA) statistics are based on either exact or approximate methods.
#' Exact methods estimate variance components and effects of single markers jointly.
#' Approximate methods estimate single marker effects conditionally.
#'
#' @param fit list of information about linear model fit (output from greml)
#' @param W matrix of centered and scaled genotypes (n x m)
#' @param m is the total number of markers in W 
#' @param statistic is the single marker test statistic used. Default is the "mastor", alternatives include "gblup", "bolt-lmm" and "grammar-gamma"
#' @return Returns a dataframe including 
#' \item{coef}{single marker coefficients} 
#' \item{se}{standard error of coefficients}
#' \item{stat}{single marker test statistic}
#' \item{p}{p-value}
#' @author Peter Sørensen
#' @references Chen, W. M., & Abecasis, G. R. (2007). Family-based association tests for genomewide association scans. The American Journal of Human Genetics, 81(5), 913-926.
#' @references Loh, P. R., Tucker, G., Bulik-Sullivan, B. K., Vilhjalmsson, B. J., Finucane, H. K., Salem, R. M., ... & Patterson, N. (2015). Efficient Bayesian mixed-model analysis increases association power in large cohorts. Nature genetics, 47(3), 284-290.
#' @references Kang, H. M., Sul, J. H., Zaitlen, N. A., Kong, S. Y., Freimer, N. B., Sabatti, C., & Eskin, E. (2010). Variance component model to account for sample structure in genome-wide association studies. Nature genetics, 42(4), 348-354.
#' @references Lippert, C., Listgarten, J., Liu, Y., Kadie, C. M., Davidson, R. I., & Heckerman, D. (2011). FaST linear mixed models for genome-wide association studies. Nature methods, 8(10), 833-835.
#' @references Listgarten, J., Lippert, C., Kadie, C. M., Davidson, R. I., Eskin, E., & Heckerman, D. (2012). Improved linear mixed models for genome-wide association studies. Nature methods, 9(6), 525-526.
#' @references Listgarten, J., Lippert, C., & Heckerman, D. (2013). FaST-LMM-Select for addressing confounding from spatial structure and rare variants. Nature Genetics, 45(5), 470-471.
#' @references Lippert, C., Quon, G., Kang, E. Y., Kadie, C. M., Listgarten, J., & Heckerman, D. (2013). The benefits of selecting phenotype-specific variants for applications of mixed models in genomics. Scientific reports, 3.
#' @references Zhou, X., & Stephens, M. (2012). Genome-wide efficient mixed-model analysis for association studies. Nature genetics, 44(7), 821-824.
#' @references Svishcheva, G. R., Axenovich, T. I., Belonogova, N. M., van Duijn, C. M., & Aulchenko, Y. S. (2012). Rapid variance components-based method for whole-genome association analysis. Nature genetics, 44(10), 1166-1170.
#' @references Yang, J., Zaitlen, N. A., Goddard, M. E., Visscher, P. M., & Price, A. L. (2014). Advantages and pitfalls in the application of mixed-model association methods. Nature genetics, 46(2), 100-106.
#' @references Bulik-Sullivan, B. K., Loh, P. R., Finucane, H. K., Ripke, S., Yang, J., Patterson, N., ... & Schizophrenia Working Group of the Psychiatric Genomics Consortium. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nature genetics, 47(3), 291-295.
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
#' # REML analyses and single marker association test
#' fitGB <- greml(y = y, X = X, G = GB, verbose = TRUE)
#' maGB <- lmma(fit = fitGB, W = W)
#'
#' @export
#'

lmma <- function( fit=NULL, W=NULL, m=NULL, statistic="mastor") {
  
  if (is.null(m)) m <- ncol(W)
  n <- nrow(W)
  # get linear model fit results
  Sg <- fit$theta[1]
  Py <- fit$Py
  Vy <- fit$Vy
  yVy <- fit$yVy
  trPG <- fit$trPG[1]
  trVG <- fit$trVG[1]
  
  # compute single marker, coefficients, test statistics and p-values
  nW <- colSums(!W==0)
  WPy <- crossprod(W,Py)
  WVy <- crossprod(W,Vy)
  ww <- colSums(W**2)
  
  if (statistic=="mastor") {
    coef <- (WPy/(trPG/m))/m
    se <- (1/(trPG/m))/m
    stat <- WPy**2/sum(Py**2)
    p <- pchisq(stat,df=1,ncp=0, lower.tail=FALSE)
  }
  
  mma <- data.frame(coef=coef,se=se,stat=stat,p=p)
  rownames(mma) <- colnames(W)
  
  return(mma)
  
}

#' @export

plotma <- function(ma=NULL,chr=NULL,rsids=NULL,thresh=5) {
  
  mlogObs <- -log10(ma$p)
  m <- length(mlogObs)
  rwsSNP <- rsids
  if(!is.null(thresh)) { rwsSNP <- mlogObs > thresh }  
  if(is.null(chr)) { chr <- rep(1,m) }  
  
  layout(matrix(1:2,ncol=1))
  o <- order(mlogObs,decreasing=TRUE)
  mlogExp <- -log10( (1:m)/m ) 
  plot(y=mlogObs[o],x=mlogExp, col=2, pch="+",
       frame.plot=FALSE, main="", xlab="Expected -log10(p)", ylab="Observed -log10(p)")
  abline(a=0,b=1)
  
  colSNP <- "red"
  
  is.even <- function(x) x %% 2 == 0
  colCHR <- rep(gray(0.3),m)
  colCHR[is.even(chr)] <- gray(0.9)
  
  plot(y=mlogObs,x=1:m,ylab="Observed -log10(p)",xlab="Position",col=colCHR,   
       pch=".",frame.plot=FALSE)
  points(y=mlogObs[rwsSNP],x=(1:m)[rwsSNP],col=colSNP)  
  abline(h=thresh,col=2,lty=2)
}   


#' @export

lma <- function( y=NULL, X=NULL, W=NULL, Glist=NULL, ids=NULL, rsids=NULL, msize=100, scaled=TRUE) {
     if(is.vector(y)) y <- matrix(y,ncol=1, dimnames= list(names(y),"trait"))
     ids <- rownames(y)
     if(!is.null(W)) {
          if(any(!ids==rownames(W))) stop("Some names of y does not match rownames of W")
          if(!is.null(X)) y <- residuals(lm(y~X))
          res <- smlm(y=y,X=X,W=W)
          #if(is.null(colnames(y))) colnames(y) <- paste("t",1:ncol(y),sep="")
          #lapply(res, function(x){colnames(x) <- colnames(y)})
          #lapply(res, function(x){rownames(x) <- colnames(W)})
     }
     if(!is.null(Glist)) {
          if(any(!ids%in%Glist$ids)) stop("Some names of y does not match names in Glist$ids")
          if(!is.null(X)) y <- as.matrix(residuals(lm(y~X)))
          nt <- ncol(y) 
          m <- Glist$m
          n <- Glist$n
          cls <- 1:m
          if(!is.null(rsids)) cls <- match(rsids,Glist$rsids)
          rws <- 1:n
          if(!is.null(ids)) rws <- match(ids,Glist$ids)
          s <- se <- stat <- p <- matrix(NA,nrow=m,ncol=nt)
          rownames(s) <- rownames(se) <- rownames(stat) <- rownames(p) <- Glist$rsids 
          colnames(s) <- colnames(se) <- colnames(stat) <- colnames(p) <- colnames(y) 
          sets <- split(cls, ceiling(seq_along(cls)/msize))
          nsets  <-  length(sets)
          for (i in 1:nsets) {
               cls <- sets[[i]]
               W <- readbed(Glist=Glist,rws=rws, cls=cls, scaled=scaled)
               #W <- W[rws,]
               res <- smlm(y=y,X=X,W=W)
               s[cls,] <- res[[1]]
               se[cls,] <- res[[2]]
               stat[cls,] <- res[[3]]
               p[cls,] <- res[[4]]
               print(paste("Finished block",i,"out of",nsets,"blocks"))
          }
          cls <- unlist(sets)
          res <- list(s=s[cls,],se=se[cls,],stat=stat[cls,],p=p[cls,])
          if (nt==1) res <- as.data.frame(res)
     }
     return(res)
}


#' @export

smlm <- function( y=NULL, X=NULL, W=NULL) {
     if(is.vector(y)) y <- matrix(y,ncol=1)
     nt <- ncol(y) 
     ones <- matrix(1,nrow=nrow(y),ncol=nt)
     ones[is.na(y)] <- 0
     y[is.na(y)] <- 0
     m <- ncol(W)
     n <- nrow(W)
     Wy <- crossprod(W,y)
     #wwadj <- matrix((crossprod(W,ones)**2)/colSums(ones),nrow=m,ncol=nt,byrow=TRUE)
     W2 <- W**2
     ww <- crossprod(W2,ones) 
     yy <- matrix(colSums((y**2)*ones),nrow=m,ncol=nt,byrow=TRUE)
     sse <- yy-(Wy**2)/ww
     sse[is.na(sse)] <- 0
     coef <- Wy*(1/ww)
     coef[is.na(coef)] <- 0
     dfe <- colSums(ones)-2
     dfe <- matrix(dfe,nrow=m,ncol=nt,byrow=TRUE)
     se <- sqrt(sse/dfe)/sqrt(ww) 
     tt <- coef/se
     ptt <- 2*pt(-abs(tt),df=dfe)
     res <- list(s=coef,se=se,stat=tt,p=ptt)
     return(res)
}



