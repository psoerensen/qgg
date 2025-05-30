####################################################################################################################
#    Module 1: LMM marker association test
####################################################################################################################
#'
#' Single marker association analysis using linear models or linear mixed models
#'
#' @description
#' The function glma performs single marker association analysis between genotype markers and the phenotype
#' either based on linear model analysis (LMA) or mixed linear model analysis (MLMA).
#'
#' The basic MLMA approach involves 1) building a genetic relationship matrix (GRM) that models genome-wide
#' sample structure, 2) estimating the contribution of the GRM to phenotypic variance using a random effects model
#' (with or without additional fixed effects) and 3) computing association statistics that account for this component
#' on phenotypic variance.
#'
#' MLMA methods are the method of choice when conducting association mapping in the presence of sample structure,
#' including geographic population structure, family relatedness and/or cryptic relatedness. MLMA methods prevent
#' false positive associations and increase power. The general recommendation when using MLMA is to exclude candidate
#' markers from the GRM. This can be efficiently implemented via a leave-one-chromosome-out analysis.
#' Further, it is recommend that analyses of randomly ascertained quantitative traits should include all markers
#' (except for the candidate marker and markers in LD with the candidate marker) in the GRM, except as follows.
#' First, the set of markers included in the GRM can be pruned by LD to reduce running time (with association
#' statistics still computed for all markers). Second, genome-wide significant markers of large effect should be
#' conditioned out as fixed effects or as an additional random effect (if a large number of associated markers).
#' Third, when population stratification is less of a concern, it may be useful using the top associated markers
#' selected based on the global maximum from out-of sample predictive accuracy.
#'
#'
#'
#' @param y vector or matrix of phenotypes
#' @param X design matrix for factors modeled as fixed effects
#' @param fit list of information about linear mixed model fit (output from greml)
#' @param Glist list of information about genotype matrix stored on disk
#' @param W matrix of centered and scaled genotypes
#' @param rsids vector of marker rsids used in the analysis
#' @param ids vector of individuals used in the analysis
#' @param statistic single marker test statistic used (currently based on the "mastor" statistics).
#' @param msize number of genotype markers used for batch processing
#' @param scale logical if TRUE the genotypes have been scaled to mean zero and variance one
#' @param verbose is a logical; if TRUE it prints more details during optimization
#' @param chr chromosome for which summary statistics are computed
#'
#' @return Returns a dataframe (if number of traits = 1) else a list including
#' \item{coef}{single marker coefficients}
#' \item{se}{standard error of coefficients}
#' \item{stat}{single marker test statistic}
#' \item{p}{p-value}
#' @author Peter Soerensen
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
#' W <- matrix(rnorm(1000000), ncol = 1000)
#' 	colnames(W) <- as.character(1:ncol(W))
#' 	rownames(W) <- as.character(1:nrow(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 501:510]) + rnorm(nrow(W))
#'
#' # Create model
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ 0 + mu
#' X <- model.matrix(fm, data = data)
#' 
#' # Linear model analyses and single marker association test
#' stat <- glma(y=y,X=X,W = W)
#'
#' head(stat)
#'
#' \donttest{
#' # Compute GRM
#' GRM <- grm(W = W)
#'
#' # Estimate variance components using REML analysis
#' fit <- greml(y = y, X = X, GRM = list(GRM), verbose = TRUE)
#'
#' # Single marker association test
#' stat <- glma(fit = fit, W = W)
#'
#' head(stat)
#'
#' }
#' 
#' 

#' @export
#'
#' 
glma <- function(y = NULL, X = NULL, W = NULL, Glist = NULL, chr=NULL,  fit = NULL, verbose=FALSE,
                 statistic = "mastor", ids = NULL, rsids = NULL, msize = 100, scale = TRUE) {
  if (is.null(fit)) {
    ma <- sma(y = y, X = X, W = W, Glist = Glist, chr=chr, ids = ids, rsids = rsids, msize = msize, scale = scale, verbose=verbose)
    return(ma)
  }
  if (!is.null(fit)) {
    ma <- mlma(y = y, X = X, fit = fit, W = W, statistic = statistic)
    return(ma)
  }
}

mlma <- function(y = NULL, X = NULL, fit = NULL, W = NULL, m = NULL, statistic = "mastor") {
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
  nW <- colSums(!W == 0)
  WPy <- crossprod(W, Py)
  WVy <- crossprod(W, Vy)
  ww <- colSums(W**2)
  
  if (statistic == "mastor") {
    coef <- (WPy / (trPG / m)) / m
    se <- (1 / (trPG / m)) / m
    stat <- WPy**2 / sum(Py**2)
    p <- pchisq(stat, df = 1, ncp = 0, lower.tail = FALSE)
  }
  
  mma <- data.frame(coef = coef, se = se, stat = stat, p = p)
  rownames(mma) <- colnames(W)
  
  return(as.matrix(mma))
}


plotma <- function(ma = NULL, chr = NULL, rsids = NULL, thresh = 5) {
  mlogObs <- -log10(ma$p)
  m <- length(mlogObs)
  rwsSNP <- rsids
  if (!is.null(thresh)) {
    rwsSNP <- mlogObs > thresh
  }
  if (is.null(chr)) {
    chr <- rep(1, m)
  }
  
  layout(matrix(1:2, ncol = 1))
  o <- order(mlogObs, decreasing = TRUE)
  mlogExp <- -log10((1:m) / m)
  plot(
    y = mlogObs[o], x = mlogExp, col = 2, pch = "+",
    frame.plot = FALSE, main = "", xlab = "Expected -log10(p)", ylab = "Observed -log10(p)"
  )
  abline(a = 0, b = 1)
  
  colSNP <- "red"
  
  is.even <- function(x) x %% 2 == 0
  colCHR <- rep(gray(0.3), m)
  colCHR[is.even(chr)] <- gray(0.9)
  
  plot(y = mlogObs, x = 1:m, ylab = "Observed -log10(p)", xlab = "Position", col = colCHR,
       pch = ".", frame.plot = FALSE
  )
  points(y = mlogObs[rwsSNP], x = (1:m)[rwsSNP], col = colSNP)
  abline(h = thresh, col = 2, lty = 2)
}


sma <- function(y = NULL, X = NULL, W = NULL, Glist = NULL, chr=NULL, ids = NULL, rsids = NULL,
                msize = 100, scale = TRUE, verbose=FALSE) {
  if (is.list(y)) {
    if(is.null(names(y))) names(y) <- paste0("T",1:length(y)) 
    y <- as.matrix(as.data.frame(y))
  }
  
  if (is.vector(y)) y <- matrix(y, ncol = 1, dimnames = list(names(y), "trait"))
  ids <- rownames(y)
  
  nt <- ncol(y)
  if (!is.null(W)) {
    if (any(!ids == rownames(W))) stop("Some names of y does not match rownames of W")
    if (!is.null(X)) y <- residuals(lm(y ~ X))
    if (is.null(X)) y <- residuals(lm(y ~ 1))
    ma <- smlm(y = y, X = X, W = W)
    if (nt == 1) ma <- as.data.frame(ma)
  }
  if (!is.null(Glist)) {
    if (any(!ids %in% Glist$ids)) stop("Some names of y does not match names in Glist$ids")
    if (!is.null(X) && !any(is.na(y))) y <- as.matrix(residuals(lm(y ~ X)))
    if (!is.null(X) && any(is.na(y))) stop("Need to fix missing data lm")
    if (is.null(X) && !any(is.na(y))) y <- as.matrix(residuals(lm(y ~ 1)))
    if (is.null(X)) {y <- apply(y,2,function(x) { 
      x[!is.na(x)] <- x[!is.na(x)] - mean(x, na.rm=TRUE)
      x})}
    nt <- ncol(y)
    ma <- vector(length=Glist$nchr,mode="list")
    if(is.null(chr)) chromosomes <- 1:Glist$nchr
    if(!is.null(chr)) chromosomes <- chr
    message("Type of analysis: glma")
    for (chr in  chromosomes) {
      message(paste("Processing chromosome:", chr, "out of", Glist$nchr, "chromosomes"))
      m <- Glist$mchr[chr]
      cls <- 1:m
      n <- Glist$n
      rws <- 1:n
      if (!is.null(ids)) rws <- match(ids, Glist$ids)
      s <- se <- stat <- p <- matrix(NA, nrow = m, ncol = nt)
      dfe <- ww <- wy <- matrix(NA, nrow = m, ncol = nt)
      rownames(s) <- rownames(se) <- rownames(stat) <- rownames(p) <- Glist$rsids[[chr]]
      colnames(s) <- colnames(se) <- colnames(stat) <- colnames(p) <- colnames(y)
      rownames(dfe) <- rownames(ww) <- rownames(wy) <- Glist$rsids[[chr]]
      colnames(dfe) <- colnames(ww) <- colnames(wy) <- colnames(y)
      if (!is.null(rsids)) cls <- match(rsids, Glist$rsids[[chr]])
      rsidsChr <- Glist$rsids[[chr]]
      if (!is.null(rsids)) {
        rsidsChr <- rsidsChr[rsidsChr%in%rsids]
        cls <- match(rsidsChr,Glist$rsids[[chr]])
      }
      sets <- split(cls, ceiling(seq_along(cls) / msize))
      # Check if the last set has length 1
      if (length(sets[[length(sets)]]) == 1 && length(sets) > 1) {
        # Merge the last element with the second-last
        sets[[length(sets) - 1]] <- c(sets[[length(sets) - 1]], sets[[length(sets)]])
        # Remove the last element
        sets <- sets[-length(sets)]
      }
      nsets <- length(sets)
      for (i in 1:nsets) {
        cls <- sets[[i]]
        W <- getG(Glist, chr=chr, scale=scale, cls=cls)
        #W <- getW(Glist = Glist, rws = rws, cls = cls, scale = scale)
        res <- smlm(y = y, X = X, W = W[rws,])
        s[cls, ] <- res[[1]]
        se[cls, ] <- res[[2]]
        stat[cls, ] <- res[[3]]
        p[cls, ] <- res[[4]]
        dfe[cls, ] <- res[[5]]
        ww[cls, ] <- res[[6]]
        wy[cls, ] <- res[[7]]
        if(verbose) message(paste("Finished chromosome segment", i, "out of", nsets))
      }
      cls <- unlist(sets)
      ma[[chr]] <- list(b = s[cls, ], seb = se[cls, ], stat = stat[cls, ], p = p[cls, ],
                        n = dfe[cls, ], ww = ww[cls, ], wy = wy[cls, ])
      if (nt == 1) ma[[chr]] <- as.data.frame(ma[[chr]])
    }
    if (nt == 1) {
      ma <- do.call(rbind, ma)
      if(!is.null(Glist)) {
        rsids <- rownames(ma)
        ma$rsids <- rsids
        rsids2rws <- match(rsids,unlist(Glist$rsids))
        ma$chr <- unlist(Glist$chr)[rsids2rws]
        ma$pos <- unlist(Glist$pos)[rsids2rws]
        #ma$a1 <- unlist(Glist$a1)[rsids2rws]
        #ma$a2 <- unlist(Glist$a2)[rsids2rws]
        #ma$af <- unlist(Glist$af)[rsids2rws]
        ma$ea <- unlist(Glist$a1)[rsids2rws]
        ma$nea <- unlist(Glist$a2)[rsids2rws]
        ma$eaf <- unlist(Glist$af)[rsids2rws]
        ma <- ma[,c(8:13,1:7)]
        ma <- na.omit(ma)
      }
    }
    if(nt>1) {
      b <- do.call(rbind, lapply(ma,function(x){x$b}))
      seb <- do.call(rbind, lapply(ma,function(x){x$seb}))
      stat <- do.call(rbind, lapply(ma,function(x){x$stat}))
      p <- do.call(rbind, lapply(ma,function(x){x$p}))
      n <- do.call(rbind, lapply(ma,function(x){x$n}))
      ww <- do.call(rbind, lapply(ma,function(x){x$ww}))
      wy <- do.call(rbind, lapply(ma,function(x){x$wy}))
      ma <- list(b = b, seb = seb, stat = stat, p = p,
                 n = n, ww = ww, wy = wy)
      if(!is.null(Glist)) {
        rsids <- rownames(ma$b)
        rsids2rws <- match(rsids,unlist(Glist$rsids))
        ma$marker <- data.frame( rsids=rsids, 
                                 chr=unlist(Glist$chr)[rsids2rws],
                                 pos=unlist(Glist$pos)[rsids2rws],
                                 #a1=unlist(Glist$a1)[rsids2rws],
                                 #a2=unlist(Glist$a2)[rsids2rws],
                                 #af=unlist(Glist$af)[rsids2rws])
                                 ea=unlist(Glist$a1)[rsids2rws],
                                 nea=unlist(Glist$a2)[rsids2rws],
                                 eaf=unlist(Glist$af)[rsids2rws])
      }
      
    }
  }
  return(ma)
  }
smlm <- function(y = NULL, X = NULL, W = NULL) {
  if (is.vector(y)) y <- matrix(y, ncol = 1)
  nt <- ncol(y)
  ones <- matrix(1, nrow = nrow(y), ncol = nt)
  ones[is.na(y)] <- 0
  y[is.na(y)] <- 0
  m <- ncol(W)
  n <- nrow(W)
  Wy <- crossprod(W, y)
  W2 <- W**2
  ww <- crossprod(W2, ones)
  yy <- matrix(colSums((y**2) * ones), nrow = m, ncol = nt, byrow = TRUE)
  sse <- yy - (Wy**2) / ww
  sse[is.na(sse)] <- 0
  coef <- Wy * (1 / ww)
  coef[is.na(coef)] <- 0
  dfe <- colSums(ones) - 2
  dfe <- matrix(dfe, nrow = m, ncol = nt, byrow = TRUE)
  se <- sqrt(sse / dfe) / sqrt(ww)
  tt <- coef / se
  ptt <- 2 * pt(-abs(tt), df = dfe)
  res <- list(b = coef, seb = se, stat = tt, p = ptt, dfe = dfe, ww = ww, wy = Wy)
  return(res)
}


cvs <- function(y=NULL, Glist = NULL, chr = NULL, bedfiles = NULL, bimfiles = NULL, famfiles = NULL, ids = NULL, rsids = NULL,
                rws = NULL, cls = NULL, impute = TRUE, scale = TRUE) {
  
  if(is.null(cls)) cls <- 1:Glist$mchr[chr]
  af <- Glist$af[[chr]][cls]
  ids <- NULL
  if(is.matrix(y)) ids <- rownames(y)
  if(is.vector(y)) ids <- names(y)
  if(is.vector(y)) y <- as.matrix(y)
  nt <- ncol(y)
  for (i in 1:ncol(y)) {
    y[,i] <- y[,i] - mean(y[,i])
  }
  if(is.null(ids)) stop("No names/rownames provided for y")
  if(!is.null(ids)) {
    if(any(is.na(match(ids,Glist$ids))))  stop("Names/rownames for y does match rownames for W")
    if(any(is.na(Glist$ids%in%ids)))  stop("Names/rownames for y does match rownames for W")
  }
  rws <- match(ids,Glist$ids)
  if(any(is.na(rws))) stop("Some ids in y does not match individuals in Glist$ids")
  dfe <- nrow(y)-1
  ylist <- lapply(1:nt,function(x) {rep(0,Glist$n)})
  weights <- lapply(1:nt,function(x) {rep(0,Glist$n)})
  for (i in 1:nt) {
    ylist[[i]][rws] <- y[,i]
    weights[[i]][rws] <- 1.0
  }
  if(!file.exists(Glist$bedfiles[chr])) stop(paste("Bedfile:", Glist$bedfiles[chr],"does not exist"))
  covs <- .Call("_qgg_summarybed", 
                Glist$bedfiles[chr], 
                n=Glist$n,
                cls=cls,
                af=af,
                weights=weights,
                y=ylist)
  
  res <- vector(length=nt,mode="list")
  names(res) <- colnames(y)
  for ( i in 1:nt) {
    rsids <- Glist$rsids[[chr]][cls]
    ww <- covs[[1]][[i]]
    wy <- covs[[2]][[i]]
    b <- (covs[[2]][[i]]/covs[[1]][[i]])
    seb <- 1/sqrt(covs[[1]][[i]])
    tstat <- (covs[[2]][[i]]/covs[[1]][[i]])*sqrt(covs[[1]][[i]])
    p <- 2 * pt(-abs(tstat), df = dfe - 2)
    res[[i]] <- data.frame(rsids,ww,wy,b,seb,tstat,p)
    rownames(res[[i]]) <- rsids
  }
  if(nt==1) res <- res[[1]]
  return(res)
}


