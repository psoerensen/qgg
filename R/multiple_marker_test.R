####################################################################################################################
#    Module 2: Marker set tests
####################################################################################################################
#'
#' Gene set enrichment analysis
#'
#' @description
#' The function gsea can perform several different gene set enrichment analyses. The general procedure is to obtain
#' single marker statistics (e.g. summary statistics), from which it is possible to compute and evaluate a test statistic
#' for a set of genetic markers that measures a joint degree of association between the marker set and the phenotype.
#' The marker set is defined by a genomic feature such as genes, biological pathways, gene interactions,
#' gene expression profiles etc.
#'
#' Currently, four types of gene set enrichment analyses can be conducted with gsea; sum-based, count-based,
#' score-based, and our own developed method, the covariance association test (CVAT). For details and comparisons of
#' test statistics consult doi:10.1534/genetics.116.189498.
#'
#' The sum test is based on the sum of all marker summary statistics located within the feature set. The single marker
#' summary statistics can be obtained from linear model analyses (from PLINK or using the qgg glma approximation),
#' or from single or multiple component REML analyses (GBLUP or GFBLUP) from the greml function. The sum test is powerful
#' if the genomic feature harbors many genetic markers that have small to moderate effects.
#'
#' The count-based method is based on counting the number of markers within a genomic feature that show association
#' (or have single marker p-value below a certain threshold) with the phenotype. Under the null hypothesis (that the
#' associated markers are picked at random from the total number of markers, thus, no enrichment of markers in any
#' genomic feature) it is assumed that the observed count statistic is a realization from a hypergeometric distribution.
#'
#' The score-based approach is based on the product between the scaled genotypes in a genomic feature and the residuals
#' from the liner mixed model (obtained from greml).
#'
#' The covariance association test (CVAT) is derived from the fit object from greml (GBLUP or GFBLUP), and measures
#' the covariance between the total genomic effects for all markers and the genomic effects of the markers within the
#' genomic feature.
#'
#' The distribution of the test statistics obtained from the sum-based, score-based and CVAT is unknown, therefore
#' a circular permutation approach is used to obtain an empirical distribution of test statistics.
#'
#' @param stat vector or matrix of single marker statistics (e.g. coefficients, t-statistics, p-values)
#' @param sets list of marker sets - names corresponds to row names in stat
#' @param nperm number of permutations used for obtaining an empirical p-value
#' @param ncores number of cores used in the analysis
#' @param Glist list providing information about genotypes stored on disk
#' @param W matrix of centered and scaled genotypes (used if method = cvat or score)
#' @param fit list object obtained from a linear mixed model fit using the greml function
#' @param g vector (or matrix) of genetic effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param e vector (or matrix) of residual effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param method including sum, cvat, hyperg, score
#' @param threshold used if method='hyperg' (threshold=0.05 is default)

#' @return Returns a dataframe or a list including
#' \item{stat}{marker set test statistics}
#' \item{m}{number of markers in the set}
#' \item{p}{enrichment p-value for marker set}

#' @author Peter Soerensen


#' @examples
#'
#'
#'  # Simulate data
#'  W <- matrix(rnorm(1000000), ncol = 1000)
#'  colnames(W) <- as.character(1:ncol(W))
#'  rownames(W) <- as.character(1:nrow(W))
#'  y <- rowSums(W[, 1:10]) + rowSums(W[, 501:510]) + rnorm(nrow(W))
#'
#'  # Create model
#'  data <- data.frame(y = y, mu = 1)
#'  fm <- y ~ 0 + mu
#'  X <- model.matrix(fm, data = data)
#'
#'  # Single marker association analyses
#'  stat <- glma(y=y,X=X,W=W)
#'
#'  # Create marker sets
#'  f <- factor(rep(1:100,each=10), levels=1:100)
#'  sets <- split(as.character(1:1000),f=f)
#'
#'  # Set test based on sums
#'  b2 <- stat[,"stat"]**2
#'  names(b2) <- rownames(stat)
#'  mma <- gsea(stat = b2, sets = sets, method = "sum", nperm = 100)
#'  head(mma)
#'
#'  # Set test based on hyperG
#'  p <- stat[,"p"]
#'  names(p) <- rownames(stat)
#'  mma <- gsea(stat = p, sets = sets, method = "hyperg", threshold = 0.05)
#'  head(mma)
#'
#' \donttest{
#'  G <- grm(W=W)
#'  fit <- greml(y=y, X=X, GRM=list(G=G), theta=c(10,1))
#'
#'  # Set test based on cvat
#'  mma <- gsea(W=W,fit = fit, sets = sets, nperm = 1000, method="cvat")
#'  head(mma)
#'
#'  # Set test based on score
#'  mma <- gsea(W=W,fit = fit, sets = sets, nperm = 1000, method="score")
#'  head(mma)
#'
#' }


#' @export
gsea <- function(stat = NULL, sets = NULL, Glist = NULL, W = NULL, fit = NULL, g = NULL, e = NULL, threshold = 0.05, method = "sum", nperm = 1000, ncores = 1) {
  if(is.data.frame(stat)) {
    colstat <- !colnames(stat)%in%c("rsids","chr","pos","ea","nea","eaf")
    if(any(!rownames(stat)==stat$rsids)) stop("Row names of stat does not match stat$rsids")
    if(any(colstat)) stat <- as.matrix(stat[,colstat])**2
    if(!any(colstat)) stat <- as.matrix(stat[,colstat])
  }
  if (method == "sum") {
    #m <- length(stat)
    if (is.matrix(stat)) sets <- mapSets(sets = sets, rsids = rownames(stat), index = TRUE)
    if (is.vector(stat)) sets <- mapSets(sets = sets, rsids = names(stat), index = TRUE)
    nsets <- length(sets)
    msets <- sapply(sets, length)
    if (is.matrix(stat)) {
      p <- apply(stat, 2, function(x) {
        gsets(stat = x, sets = sets, ncores = ncores, np = nperm)
      })
      setstat <- apply(stat, 2, function(x) {
        sapply(sets, function(y) {
          sum(x[y])
        })
      })
      rownames(setstat) <- rownames(p) <- names(msets) <- names(sets)
      res <- list(m = msets, stat = setstat, p = p)
    }
    if (is.vector(stat)) {
      setstat <- sapply(sets, function(x) {
        sum(stat[x])
      })
      p <- gsets(stat = stat, sets = sets, ncores = ncores, np = nperm, method = method)
      res <- cbind(m = msets, stat = setstat, p = p)
      rownames(res) <- names(sets)
      res <- as.data.frame(res)
    }
    return(res)
  }
  if (method == "cvat") {
    if (!is.null(W)) res <- cvat(fit = fit, g = g, W = W, sets = sets, nperm = nperm)
    if (!is.null(Glist)) {
      sets <- mapSets(sets = sets, rsids = Glist$rsids, index = TRUE)
      nsets <- length(sets)
      msets <- sapply(sets, length)
      ids <- fit$ids
      Py <- fit$Py
      g <- as.vector(fit$g)
      Sg <- fit$theta[1]
      stat <- gstat(method = "cvat", Glist = Glist, g = g, Sg = Sg, Py = Py, ids = ids)
      setstat <- sapply(sets, function(x) {
        sum(stat[x])
      })
      p <- gsets(stat = stat, sets = sets, ncores = ncores, np = nperm, method = "sum")
      res <- cbind(m = msets, stat = setstat, p = p)
      rownames(res) <- names(sets)
    }
    return(res)
  }
  if (method == "score") {
    if (!is.null(W)) res <- scoretest(e = fit$e, W = W, sets = sets, nperm = nperm)
    if (!is.null(Glist)) {
        sets <- mapSets(sets = sets, rsids = Glist$rsids, index = TRUE)
      nsets <- length(sets)
      msets <- sapply(sets, length)
      ids <- fit$ids
      e <- fit$e
      stat <- gstat(method = "score", Glist = Glist, e = e, ids = ids)
      setstat <- sapply(sets, function(x) {
        sum(stat[x])
      })
      p <- gsets(stat = stat, sets = sets, ncores = ncores, np = nperm, method = "sum")
      res <- cbind(m = msets, stat = setstat, p = p)
      rownames(res) <- names(sets)
    }
    return(res)
  }
  if (method == "hyperg") {
    res <- hgtest(p = stat, sets = sets, threshold = threshold)
    return(as.data.frame(res))
  }
}


gsets <- function(stat = NULL, sets = NULL, ncores = 1, np = 1000, method = "sum") {
  m <- length(stat)
  nsets <- length(sets)
  msets <- sapply(sets, length)
  setstat <- sapply(sets, function(x) {
    sum(stat[x])
  })
  p <- .Call("_qgg_psets", msets = msets,
              setstat = setstat,
              stat = stat,
              np = np)
  p <- p/np 

  return(p)
}


#' Map Sets to RSIDs
#'
#' This function maps sets to rsids. If a `Glist` is provided, `rsids` are extracted from the `Glist`.
#' It returns a list of matched RSIDs for each set.
#'
#' @param sets A list of character vectors where each vector represents a set of items. If the names
#'   of the sets are not provided, they are named as "Set1", "Set2", etc.
#' @param rsids A character vector of RSIDs. If `Glist` is provided, this parameter is ignored.
#' @param Glist A list containing an element `rsids` which is a character vector of RSIDs.
#' @param index A logical. If `TRUE` (default), it returns indices of RSIDs; otherwise, it returns the RSID names.
#' 
#' @return A list where each element represents a set and contains matched RSIDs or their indices.
#' 
#' @keywords internal
#' @export
mapSets <- function(sets = NULL, rsids = NULL, Glist = NULL, index = TRUE) {
  if (!is.null(Glist)) rsids <- unlist(Glist$rsids)
  nsets <- sapply(sets, length)
  if(is.null(names(sets))) names(sets) <- paste0("Set",1:length(sets))
  rs <- rep(names(sets), times = nsets)
  rsSets <- unlist(sets, use.names = FALSE)
  rsSets <- match(rsSets, rsids)
  inW <- !is.na(rsSets)
  rsSets <- rsSets[inW]
  if (!index) rsSets <- rsids[rsSets]
  rs <- rs[inW]
  rs <- factor(rs, levels = unique(rs))
  rsSets <- split(rsSets, f = rs)
  return(rsSets)
}




gstat <- function(method = NULL, Glist = NULL, g = NULL, Sg = NULL, Py = NULL, e = NULL, msize = 100, rsids = NULL,
                  impute = TRUE, scale = TRUE, ids = NULL, ncores = 1) {
  n <- Glist$n
  rws <- match(ids, Glist$ids)
  if (any(is.na(rws))) stop("Some ids in fit object not found in Glist")
  nr <- length(rws)
  nbytes <- ceiling(n / 4)
  cls <- 1:Glist$m
  if (!is.null(rsids)) cls <- match(rsids, Glist$rsids)
  nc <- length(cls)
  cls <- split(cls, ceiling(seq_along(cls) / msize))
  msets <- sapply(cls, length)
  nsets <- length(msets)
  setstat <- NULL
  fnRAW <- Glist$fnRAW
  for (j in 1:nsets) {
    nc <- length(cls[[j]])
    direction <- rep(1, nc)

# Check this again - should be getW using a bedfile
      W <- getW(Glist = Glist, rws = rws, cls = cls, scale = scale)
# End check this again

    if (method == "cvat") {
      s <- crossprod(W / nc, Py) * Sg
      Ws <- t(t(W) * as.vector(s))
      setstat <- c(setstat, colSums(g * Ws))
    }
    if (method == "score") {
      we2 <- as.vector((t(W) %*% e)**2)
      setstat <- c(setstat, we2)
    }
    message(paste("Finished block", j, "out of", nsets, "blocks"))
  }
  return(setstat)
}




settest <- function(stat = NULL, W = NULL, sets = NULL, nperm = NULL, method = "sum", threshold = 0.05) {
  if (method == "sum") setT <- sumtest(stat = stat, sets = sets, nperm = nperm)
  if (method == "hyperG") setT <- hgtest(p = stat, sets = sets, threshold = threshold)
  return(setT)
}



sumtest <- function(stat = NULL, sets = NULL, nperm = NULL, method = "sum") {
  if (method == "mean") {
    setT <- sapply(sets, function(x) {
      mean(stat[x])
    })
  }
  if (method == "sum") {
    setT <- sapply(sets, function(x) {
      sum(stat[x])
    })
  }
  if (method == "max") {
    setT <- sapply(sets, function(x) {
      max(stat[x])
    })
  }
  if (!is.null(nperm)) {
    p <- rep(0, length(sets))
    n <- length(stat)
    nset <- sapply(sets, length)
    rws <- 1:n
    names(rws) <- names(stat)
    sets <- lapply(sets, function(x) {
      rws[x]
    })
    for (i in 1:nperm) {
      rws <- sample(1:n, 1)
      o <- c(rws:n, 1:(rws - 1))
      pstat <- stat[o]
      if (method == "mean") {
        setTP <- sapply(sets, function(x) {
          mean(pstat[x])
        })
      }
      if (method == "sum") {
        setTP <- sapply(sets, function(x) {
          sum(pstat[x])
        })
      }
      if (method == "max") {
        setTP <- sapply(sets, function(x) {
          max(pstat[x])
        })
      }
      p <- p + as.numeric(setT > setTP)
    }
    p <- 1 - p / nperm
    setT <- data.frame(setT, nset, p)
  }

  return(setT)
}



cvat <- function(fit = NULL, s = NULL, g = NULL, W = NULL, sets = NULL, nperm = 100) {
  if (!is.null(fit)) {
    s <- crossprod(W / ncol(W), fit$Py) * fit$theta[1]
  }
  Ws <- t(t(W) * as.vector(s))
  if (is.null(g)) g <- W %*% s
  cvs <- colSums(as.vector(g) * Ws)
  # setT <- setTest(stat = cvs, sets = sets, nperm = nperm, method = "sum")$p
  # names(setT) <- names(sets)
  setT <- settest(stat = cvs, sets = sets, nperm = nperm, method = "sum")
  if (!is.null(names(sets))) rownames(setT) <- names(sets)
  return(setT)
}


scoretest <- function(e = NULL, W = NULL, sets = NULL, nperm = 100) {
  we2 <- as.vector((t(W) %*% e)**2)
  names(we2) <- colnames(W)
  setT <- settest(stat = we2, sets = sets, nperm = nperm, method = "sum")$p
  return(setT)
}


hgtest <- function(p = NULL, sets = NULL, threshold = 0.05) {
  population_size <- length(p)
  sample_size <- sapply(sets, length)
  n_successes_population <- sum(p < threshold)
  n_successes_sample <- sapply(sets, function(x) {
    sum(p[x] < threshold)
  })
  phyperg <- rep(1,length(sets))
  names(phyperg) <- names(sets)
  for (i in 1:length(sets)) {
    phyperg[i] <- 1.0-phyper(n_successes_sample[i]-1, n_successes_population,
                             population_size-n_successes_population,
                             sample_size[i])
  }
  # Calculate enrichment factor
  ef <- (n_successes_sample/sample_size)/
    (n_successes_population/population_size)
  
  # Create data frame for table
  res <- data.frame(ng = sample_size,
                   nag = n_successes_sample,
                   ef=ef,
                   p = phyperg)
  #colnames(res) <- c("Feature", "Number of Genes",
  #                  "Number of Associated Genes",
  #                  "Enrichment Factor",
  #                  "p")
  res
}


# hgtest <- function(p = NULL, sets = NULL, threshold = 0.05) {
#   N <- length(p)
#   Na <- sum(p < threshold)
#   Nna <- N - Na
#   Nf <- sapply(sets, length)
#   Naf <- sapply(sets, function(x) {
#     sum(p[x] < threshold)
#   })
#   Nnaf <- Nf - Naf
#   Nanf <- Na - Naf
#   Nnanf <- Nna - Nnaf
#   phyperg <- 1 - phyper(Naf - 1, Nf, N - Nf, Na)
#   phyperg
# }

#' @export
magma <- function(stat = NULL, sets = NULL, 
                  method="magma", type = "joint", test = "one-sided",
                  pi=0.001, nit=5000, nburn=1000) {
  
  # Check if stat and sets are provided
  if (is.null(stat) || is.null(sets)) {
    stop("Both 'stat' and 'sets' must be provided.")
  }
  
  if(is.vector(stat)) stat <- as.matrix(stat)
  if(is.null(rownames(stat))) stop("Please provide names or rownames to stat object")
  y <- scale(stat, center=TRUE, scale=TRUE)

  sets <- mapSets(sets=sets,rsids=rownames(stat), index=FALSE)
  
    
  # Compute X for feature sets (sparse format)
  X <- designMatrix(sets = sets, rowids = rownames(y))

  m <- sapply(sets, length)
  
  if(method=="magma") {
    
    # Compute summary stat for feature sets
    stat_summary <- computeStat(X = X, y = y[rownames(X), ], scale = TRUE)
    
    bMAGMA <- solve(stat_summary$XX + diag(0.001, nrow(stat_summary$XX))) %*% stat_summary$Xy
    bMARG <- (1 / diag(stat_summary$XX)) * stat_summary$Xy
    sebMAGMA <- sqrt(diag(solve(stat_summary$XX + diag(0.001, nrow(stat_summary$XX)))))
    sebMARG <- sqrt((1 / diag(stat_summary$XX)))
    zMAGMA <- bMAGMA / sebMAGMA
    zMARG <- bMARG / sebMARG
    
    # Two-sided
    if(test=="two-sided") {
      pMAGMA <- pnorm(abs(zMAGMA), mean = 0, sd = 1, lower.tail = FALSE)
      pMARG <- pnorm(abs(zMARG), mean = 0, sd = 1, lower.tail = FALSE)
    }
    
    # One-sided
    if(test=="one-sided") {
      pMAGMA <- pnorm(zMAGMA, mean = 0, sd = 1, lower.tail = FALSE)
      pMARG <- pnorm(zMARG, mean = 0, sd = 1, lower.tail = FALSE)
    }
    if (type == "marginal") df <- data.frame(ID=names(sets), m = m, 
                                             b = bMARG, seb = sebMARG, z = zMARG, p = pMARG)
    if (type == "joint") df <- data.frame(ID=names(sets), m = m, 
                                          b = bMAGMA, seb = sebMAGMA, z = zMAGMA, p = pMAGMA)
    o <- order(df$p, decreasing=FALSE)
    df[,3:5] <- round(df[,3:5],4)
    rownames(df) <- NULL
    return(df[o,])
  }
  if(method%in%c("blr","bayesC","bayesR")) {
    
    if(method=="blr") method <- "bayesC" 
    
    # Compute summary stat for feature sets
    stat <- computeStat(X = X, y = y[rownames(X), ], scale = TRUE)

    if(ncol(y)==1) {
      # Fit BLR model
      fit <- qgg:::blr(yy=stat$yy, XX=stat$XX, Xy=stat$Xy, n=stat$n,
                       method=method, pi=pi,
                       nit=nit, nburn=nburn)
      
      
      # Dataframe with BLR results
      df <- data.frame(ID=names(sets), 
                       m = m, 
                       b = fit$bm, 
                       PIP=fit$dm)
      o <- order(df$PIP, decreasing=TRUE)
      df[,3:4] <- round(df[,3:4],4)
      rownames(df) <- NULL
      return(df[o,])
    }  
    if(ncol(y)>1) {
      stat$XX <- rep(list(stat$XX),ncol(y))
      # Fit BLR model
      fit <- qgg:::mtblr(yy=stat$yy, XX=stat$XX, Xy=stat$Xy, n=stat$n,
                       method=method, pi=pi,
                       nit=nit, nburn=nburn)
      
      o <- order(rowSums(fit$dm), decreasing=TRUE)
      
      # Dataframe with BLR results
      df <- list( feature=data.frame(ID=names(sets),
                        m = m)[o,],
                        b=round(fit$bm[o,], 4),
                        PIP=round(fit$dm[o,], 4))
      return(df)
    }  
    
  }
}

#' @export
vegas <- function(Glist=NULL, sets=NULL, stat=NULL, p=NULL, threshold=1e-10, tol=1e-7, minsize=2, verbose=FALSE) {
  
  if(is.null(Glist)) stop("Please provide Glist object")
  #if(is.null(stat)) stop("Please provide stat object")
  if(is.null(sets)) stop("Please provide sets object")
  
  if(!is.null(stat)) {
    if(verbose) message("Map stat to markers in Glist")
    sets <- mapSets(sets=sets, rsids=stat$rsids, index=FALSE)
    sets <- mapSets(sets=sets, rsids=unlist(Glist$rsids), index=FALSE)
    
    isets <- mapSets(sets=sets, rsids=stat$rsids, index=TRUE)
    
    # Compute marker statistics
    p <- as.numeric(stat$p)
    p[p<threshold] <- threshold
    chisq <- qchisq(p, df = 1, lower.tail = FALSE)
    
    chistat <- sapply(isets,function(x){sum(chisq[x])})
    chr <- sapply(isets,function(x){stat$chr[x][1]})
    # This is just a preliminary fix
    if(length(Glist$bedfiles)==1) chr <- rep(1,length(chr))
    m <- sapply(sets,function(x){length(x)})
    
    pg <- rep(1,length(sets))
    names(pg) <- names(chr) <- names(m) <- names(sets)
    for(i in 1:length(sets)) {
      if(length(sets[[i]])>1) {
        B <- getG(Glist=Glist, chr[i], rsids=sets[[i]], scale=TRUE)
        ev <- eigen(cor(B))$values
        ev[ev < tol] <- tol
        try(pg[i] <- pchisqsum(chistat[i], df = rep(1, length(ev)), a = ev, lower.tail = FALSE))
        if(verbose) message(paste("Finished processing gene" ,i))
      } 
    }
    zstat <- -qnorm(pg/2,TRUE)
    df <- data.frame(Gene=names(pg),Chr=chr,m=m,x=chistat,z=zstat,p=pg)
    colnames(df) <- c("EnsemblID", "Chr", "m", "X2","z","p")
    return(df)
  }
  if(!is.null(p)) {
    if(is.vector(p)) p <- as.matrix(p)
    rsids <- rownames(p)
    nstudy <- ncol(p)
    if(is.null(rsids)) stop("Please provide names/rownames in for your p object")
    #p <- apply(p,2,as.numeric)
    p[p<threshold] <- threshold
    #rownames(p) <- rsids
    
    sets <- mapSets(sets=sets, rsids=rsids, index=FALSE)
    if(!is.null(Glist$rsidsLD)) sets <- mapSets(sets=sets, rsids=unlist(Glist$rsidsLD), index=FALSE)
    if(is.null(Glist$rsidsLD)) sets <- mapSets(sets=sets, rsids=unlist(Glist$rsids), index=FALSE)
    
    # Compute some relevant statistics
    msets <- sapply(sets,function(x){length(x)})
    sets <- sets[msets>minsize]
    msets <- sapply(sets,function(x){length(x)})
    chrSets <- mapSets(sets=sets, Glist=Glist, index=TRUE)
    chr <- unlist(Glist$chr)
    chr <- sapply(chrSets,function(x){as.numeric(unique(chr[x]))[1]})
    # This is just a preliminary fix
    if(length(Glist$bedfiles)==1) chr <- rep(1,length(chr))
      
    # set indices
    isets <- mapSets(sets=sets, rsids=rsids, index=TRUE)
    
    # Compute marker statistics
    chisq <- qchisq(p, df = 1, lower.tail = FALSE)
    
    # Compute gene statistics
    chistat <- sapply(isets,function(x){colSums(chisq[x,])})
    chistat <- t(chistat)
    
    pg <- matrix(1,ncol=nstudy, nrow=length(sets))
    rownames(pg) <- names(sets)
    colnames(pg) <- colnames(p)
    for(i in 1:length(sets)) {
      B <- getG(Glist=Glist, chr[i], rsids=sets[[i]], scale=TRUE)
      ev <- eigen(cor(B))$values
      ev[ev < tol] <- tol
      for (j in 1:nstudy) {
        try(pg[i,j] <- qgg:::pchisqsum(chistat[i,j], df = rep(1, length(ev)), a = ev, lower.tail = FALSE))
      }
      if(verbose) message(paste("Finished processing gene" ,i))
    }
    zstat <- -qnorm(pg/2,TRUE)
    
    if(ncol(pg)==1) {
      res <- data.frame("EnsemblID"=rownames(pg),chr=chr,m=msets,X=chistat,Z=zstat,p=pg)
      return(res)
    }
    if(ncol(pg)>1) {
      res <- list( genes=data.frame(EnsemblID=rownames(pg),Chr=chr,m=msets),
                   X2=chistat,z=zstat,p=pg)
      return(res)
    }
  }
}

pchisqsum <- function (x, df, a, lower.tail = TRUE) {
  sat <- satterthwaite(a, df)
  guess <- pchisq(x / sat$scale, sat$df, lower.tail = lower.tail)
  for (i in seq(length = length(x))) {
    lambda <- rep(a, df)
    sad <- sapply(x, saddle, lambda = lambda)
    if (lower.tail) sad <- 1 - sad
    guess <- ifelse(is.na(sad), guess, sad)
  }
  return(guess)
}

satterthwaite <- function(a, df) {
  if (any(df > 1)) {
    a <- rep(a, df)
  }
  tr <- mean(a)
  tr2 <- mean(a^2) / (tr^2)
  list(scale = tr * tr2, df = length(a) / tr2)
}

saddle <- function(x, lambda) {
  d <- max(lambda)
  lambda <- lambda / d
  x <- x / d
  k0 <- function(zeta) {
    -sum(log(1 - 2 * zeta * lambda)) / 2
  }
  kprime0 <- function(zeta) {
    sapply(zeta, function(zz) sum(lambda / (1 - 2 * zz * lambda)))
  }
  kpprime0 <- function(zeta) {
    2 * sum(lambda^2 / (1 - 2 * zeta * lambda)^2)
  }
  if (any(lambda < 0)) {
    lmin <- max(1 / (2 * lambda[lambda < 0])) * 0.99999
  } else if (x > sum(lambda)) {
    lmin <- -0.01
  } else {
    lmin <- -length(lambda) / (2 * x)
  }
  lmax <- min(1 / (2 * lambda[lambda > 0])) * 0.99999
  hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin,
                     upper = lmax, tol = 1e-08)$root
  w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
  v <- hatzeta * sqrt(kpprime0(hatzeta))
  if (abs(hatzeta) < 1e-04) {
    NA
  }  else {
    pnorm(w + log(v / w) / w, lower.tail = FALSE)
  }
}



#' @export
pops <- function(stat = NULL, sets = NULL, validate=NULL, threshold=NULL,
                 method="bayesC", pi=0.001, nit=5000, nburn=1000,
                 updateB=TRUE, updateE=TRUE, updatePi=TRUE, updateG=TRUE) {
  
  if(!is.null(validate)) {
    fit <- cvpops
    return(fit)
  }
  # Check if stat and sets are provided
  if (is.null(stat) || is.null(sets)) {
    stop("Both 'stat' and 'sets' must be provided.")
  }
  
  if(is.vector(stat)) stat <- as.matrix(stat)
  if(is.null(rownames(stat))) stop("Please provide names or rownames to stat object")
  
  # Map sets to rownames in stat
  sets <- mapSets(sets=sets,rsids=rownames(stat), index=FALSE)
  
  # Center and scale y
  y <- scale(stat, center=TRUE, scale=TRUE)
  orig_stat <- stat
  
  # Compute X for feature sets (sparse format)
  X <- designMatrix(sets = sets, rowids = rownames(y))
  X <- X[rownames(X)%in%rownames(y),]
  y <- as.matrix(y[rownames(y)%in%rownames(X),])
  X <- X[rownames(y),]
  
  if(!is.null(threshold)) {
    fit <- qgg:::magma(stat=orig_stat, sets=sets, type = "marginal", test = "one-sided",
                       method="magma")
    selected <- fit$p<threshold
    if(sum(selected)<2) stop("Number of selected features less than 2")
    cls <- fit$ID[selected]
    X <- X[,cls]
  }
  
  # Compute summary stat for feature sets
  stat <- computeStat(X = X, y = y, scale = TRUE)
  isNA <- is.na(stat$Xy)
  stat$Xy <- stat$Xy[!isNA]
  stat$XX <- stat$XX[!isNA,!isNA]
  
  if (method=="rr") {
    b <- solve(stat$XX + diag(0.001, nrow(stat$XX))) %*% stat$Xy
  }
  
  if (method%in%c("bayesC","bayesR")) {
    # Fit BLR model
    fit <- qgg:::blr(yy=stat$yy, XX=stat$XX, Xy=stat$Xy, n=stat$n,
                     method=method, pi=pi,
                     nit=nit, nburn=nburn,
                     updateB=updateB, updateE=updateE, updatePi=updatePi)
    b <- fit$bm
  }
  
  ypred <- as.matrix(X%*%b)
  colnames(ypred) <- colnames(orig_stat)
  return(ypred)
}


cvpops <- function(stat = NULL, sets = NULL, validate=NULL, threshold=NULL,
                 method="bayesC", pi=0.001, nit=5000, nburn=1000) {
  
  # Check if stat and sets are provided
  if (is.null(stat) || is.null(sets)) {
    stop("Both 'stat' and 'sets' must be provided.")
  }
  
  if(is.vector(stat)) stat <- as.matrix(stat)
  if(is.null(rownames(stat))) stop("Please provide names or rownames to stat object")
  
  # Map sets to rownames in stat
  sets <- mapSets(sets=sets,rsids=rownames(stat), index=FALSE)
  
  # Center and scale y
  y <- scale(stat, center=TRUE, scale=TRUE)
  orig_stat <- stat
  
  # Compute X for feature sets (sparse format)
  X <- designMatrix(sets = sets, rowids = rownames(y))
  X <- X[rownames(X)%in%rownames(y),]
  y <- as.matrix(y[rownames(y)%in%rownames(X),])
  X <- X[rownames(y),]
  
  
  if (is.matrix(validate)) {
    cvnames <- colnames(validate)
    validate <- as.data.frame(validate, stringsAsFactors=FALSE)
    names(validate) <- cvnames
  }
  nv <- length(validate)
  cvnames <- names(validate)
  if(is.null(cvnames)) {
    cvnames <- paste0("CV",1:nv)
    names(validate) <- cvnames
  }
  
  validate <- lapply(validate, function(x){x[x%in%rownames(y)]})
  
  nv <- length(validate)
  pa <- NULL
  
  for(v in 1:nv) {
    
    ensg <- validate[[v]]
    train <- !rownames(y)%in%ensg
    
    Xt <- X
    
    
    if(!is.null(threshold)) {
      fit <- qgg:::magma(stat=orig_stat[train,], sets=sets, type = "marginal", test = "one-sided",
                         method="magma")
      selected <- fit$p<threshold
      if(sum(selected)<2) stop("Number of selected fetaures less than 2")
      cls <- fit$ID[selected]
      Xt <- X[,cls]
    }
    
    # Compute summary stat for feature sets
    #stat <- computeStat(X = X[train,], y = y[train, ], scale = TRUE)
    stat <- computeStat(X = Xt[train,], y = y[train, ], scale = TRUE)
    isNA <- is.na(stat$Xy)
    stat$Xy <- stat$Xy[!isNA]
    stat$XX <- stat$XX[!isNA,!isNA]
    
    if (method=="rr") {
      b <- solve(stat$XX + diag(0.001, nrow(stat$XX))) %*% stat$Xy
    }
    
    if (method%in%c("bayesC","bayesR")) {
      # Fit BLR model
      fit <- qgg:::blr(yy=stat$yy, XX=stat$XX, Xy=stat$Xy, n=stat$n,
                       method=method, pi=pi,
                       nit=nit, nburn=nburn)
      b <- fit$bm
    }
    
    ypred <- Xt[,!isNA]%*%b
    yobs <- y[rownames(y)%in%ensg,]
    pa <- rbind(pa,acc(yobs=yobs, ypred=ypred[names(yobs),]))
  }
  return(pa)
}
