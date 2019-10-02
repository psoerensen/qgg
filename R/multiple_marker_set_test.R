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
#' summary statistics can be obtained from linear model analyses (from PLINK or using the qgg lma approximation),
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
#'  ma <- lma(y=y,X=X,W=W)
#'
#'  # Create marker sets
#'  f <- factor(rep(1:100,each=10), levels=1:100)
#'  sets <- split(as.character(1:1000),f=f)
#'
#'  # Set test based on sums
#'  mma <- gsea(stat = ma[,"stat"]**2, sets = sets, method = "sum", nperm = 10000)
#'  head(mma)
#'
#'  # Set test based on hyperG
#'  mma <- gsea(stat = ma[,"p"], sets = sets, method = "hyperg", threshold = 0.05)
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
  if (method == "sum") {
    m <- length(stat)
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
    return(res)
  }
}



gsets <- function(stat = NULL, sets = NULL, ncores = 1, np = 1000, method = "sum") {
  m <- length(stat)
  nsets <- length(sets)
  msets <- sapply(sets, length)
  setstat <- sapply(sets, function(x) {
    sum(stat[x])
  })

  res <- .Fortran("psets",
    m = as.integer(m),
    stat = as.double(stat),
    nsets = as.integer(nsets),
    setstat = as.double(setstat),
    msets = as.integer(msets),
    p = as.integer(rep(0, nsets)),
    np = as.integer(np),
    ncores = as.integer(ncores),
    PACKAGE = "qgg"
  )

  p <- res$p / np
}


mapSets <- function(sets = NULL, rsids = NULL, Glist = NULL, index = TRUE) {
  if (!is.null(Glist)) rsids <- unlist(Glist$rsids)
  nsets <- sapply(sets, length)
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
    W <- .Fortran("readbed",
      n = as.integer(n),
      nr = as.integer(nr),
      rws = as.integer(rws),
      nc = as.integer(nc),
      cls = as.integer(cls[[j]]),
      impute = as.integer(impute),
      scale = as.integer(scale),
      direction = as.integer(direction),
      W = matrix(as.double(0), nrow = nr, ncol = nc),
      nbytes = as.integer(nbytes),
      fnRAWCHAR = as.integer(unlist(sapply(as.character(fnRAW),charToRaw),use.names=FALSE)),
      nchars = nchar(as.character(fnRAW)),
      PACKAGE = "qgg"
    )$W
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


#' LD pruning of summary statistics
#'
#' @description
#' Perform LD pruning of summary statistics before they are used in gene set enrichment analyses.
#' @param stat vector or matrix of single marker statistics (e.g. coefficients, t-statistics, p-values)
#' @param statistics is the type of statistics used in stat (e.g. statistics="p-value")
#' @param ldSets list of marker sets - names corresponds to row names in stat
#' @param r2 threshold for r2 used in LD pruning
#' @param threshold p-value threshold used in LD pruning
#' @param Glist list providing information about genotypes stored on disk
#' @param method used including method="pruning" which is default or "clumping"


#' @export

adjLD <- function(stat = NULL, statistics = "p-value", Glist = NULL, r2 = 0.9, ldSets = NULL, threshold = 1,
                  method = "pruning") {
  rsidsStat <- rownames(stat)
  if (statistics == "p-value") pstat <- stat[, "p"]
  pstat <- as.matrix(pstat)
  rownames(pstat) <- rsidsStat
  colnames(pstat) <- "p"
  if (method %in% c("pruning", "clumping")) {
    if (!is.null(ldSets)) nchr <- length(ldSets)
    if (!is.null(Glist)) nchr <- Glist$nchr
    for (i in 1:ncol(pstat)) {
      m <- length(rsidsStat)
      indx1 <- rep(T, m)
      indx2 <- rep(F, m)
      for (chr in 1:nchr) {
        if (!is.null(Glist)) {
          setsChr <- getLDsets(Glist = Glist, r2 = r2, chr = chr)
          setsChr <- mapSets(sets = setsChr, rsids = rsidsStat)
        }
        if (!is.null(ldSets)) setsChr <- ldSets[[chr]]
        rsidsChr <- names(setsChr)
        rwsChr <- match(rsidsChr, rsidsStat)
        p <- pstat[rwsChr, i]
        o <- order(p, decreasing = FALSE)
        for (j in o) {
          if (p[j] <= threshold) {
            if (indx1[rwsChr[j]]) {
              rws <- setsChr[[j]]
              indx1[rws] <- F
              indx2[rwsChr[j]] <- T
            }
          }
        }
        message(paste("Finished pruning chromosome:", chr, "for stat column:", colnames(pstat)[i]))
      }
      if (method == "clumping") {
        pstat[indx1, i] <- 0
        p <- pstat[, i]
        pstat[p > threshold, i] <- 0
      }
      if (method == "pruning") pstat[!indx2, i] <- 0
    }
  }
  stat <- stat[!rowSums(pstat == 0) == ncol(pstat), ]
  return(stat)
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
  N <- length(p)
  Na <- sum(p < threshold)
  Nna <- N - Na
  Nf <- sapply(sets, length)
  Naf <- sapply(sets, function(x) {
    sum(p[x] < threshold)
  })
  Nnaf <- Nf - Naf
  Nanf <- Na - Naf
  Nnanf <- Nna - Nnaf
  phyperg <- 1 - phyper(Naf - 1, Nf, N - Nf, Na)
  phyperg
}
