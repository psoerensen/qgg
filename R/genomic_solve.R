####################################################################################################################
#    Module 5: GSOLVE
####################################################################################################################
#'
#' Genomic prediction based on a linear mixed model
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
#' @param Glist list of information about centered and scaled genotype matrix
#' @param rsids vector of marker rsids used in the analysis
#' @param ids vector of individuals used in the analysis
#' @param lambda vector of single marker weights used in BLUP
#' @param maxit maximum number of iterations used in the Gauss-Seidel procedure
#' @param tol tolerance, i.e. the maximum allowed difference between two consecutive iterations of reml to declare convergence
#' @param sets	list containing marker sets rsids
#' @param validate	dataframe or list of individuals used in cross-validation (one column for each set)
#' @param scale logical if TRUE the genotypes in Glist has been scaled to mean zero and variance one
#' @param ncores number of cores used in the analysis


#' @author Peter Soerensen

#' @examples
#'
#' # Simulate data
#' W <- matrix(rnorm(20000000), ncol = 10000)
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


gsolve <- function(y = NULL, X = NULL, Glist = NULL, W = NULL, ids = NULL, rsids = NULL,
                   sets = NULL, validate = NULL, scale = TRUE, lambda = NULL, weights = FALSE,
                   maxit = 500, tol = 0.00001, method = "gsru", ncores = 1) {
  if (!is.null(W)) {
    if (method == "gsru") {
      fit <- qgg:::gsru(
        y = y, W = W, X = X, sets = sets, lambda = lambda,
        weights = weights, maxit = maxit, tol = tol
      )
    }
    if (method == "gsqr") {
      fit <- qgg:::gsqr(
        y = y, W = W, X = X, sets = sets, msets = msets,
        lambda = lambda, weights = weights, maxit = maxit, tol = tol
      )
    }
    return(fit)
  }
  if (!is.null(Glist)) {
    ids <- names(y)
    fnRAW <- Glist$fnRAW
    n <- Glist$n
    nbytes <- ceiling(n / 4)
    cls <- match(rsids, Glist$rsids)
    nc <- length(cls)
    rws <- match(ids, Glist$ids)
    nr <- length(rws)

    maf <- unlist(Glist$maf)[cls]
    meanw <- 2 * maf
    sdw <- sqrt(2 * maf * (1 - maf))

    b <- bold <- bnew <- NULL
    # if (!is.null(X)) {
    #     b <- (solve(t(Xt)%*%Xt)%*%t(Xt))%*%yt     # initialize b
    #     bold <- rep(0,ncol(Xt))              # initialize b
    # }

    # if (!is.null(Xt)) yt <- yt-Xt%*%b         # initialize e

    if (length(lambda) == 1) {
      lambda <- rep(lambda, nc)
    }
    s <- rep(0, nc) # initialize diagonal elements of the W'W matrix

    yn <- gn <- en <- rep(0, n)
    yn[rws] <- as.vector(y)

    fit <- qgg:::fsolve(
      n = n, nr = nr, rws = rws, nc = nc, cls = cls, scale = scale,
      nbytes = nbytes, fnRAW = fnRAW, ncores = ncores, nit = maxit,
      lambda = lambda, tol = tol, y = yn, g = gn, e = en, s = s, meanw = meanw, sdw = sdw
    )
    # fit <- fsolve(n=n,nc=nc,cls=cls,nr=nr,rws=rws,fnRAW=fnRAW,nit=maxit,lambda=lambda,tol=tol,y=y,g=g,e=e,s=s,meanw=meanw,sdw=sdw)

    # if (!is.null(X)) yhat <- ghat[names(y)] + X[names(y),]%*%b
    yhat <- NULL
    names(fit$g) <- Glist$ids
    # e <- y - yhat
    delta <- NULL
    return(list(s = fit$s, b = b, nit = fit$nit, delta = delta, e = fit$e, yhat = yhat, g = fit$g))
  }
}



gsru <- function(y = NULL, X = NULL, W = NULL, sets = NULL, lambda = NULL, weights = FALSE, maxit = 500, tol = 0.0000001) {
  n <- length(y) # number of observations
  m <- ncol(W) # number of markers
  dww <- rep(0, m) # initialize diagonal elements of the W'W matrix
  for (i in 1:m) {
    dww[i] <- sum(W[, i]**2)
  }
  b <- bold <- bnew <- NULL
  if (!is.null(X)) {
    b <- (solve(t(X) %*% X) %*% t(X)) %*% y # initialize b
    bold <- rep(0, ncol(X)) # initialize b
  }
  if (length(lambda) == 1) {
    lambda <- rep(lambda, m)
  }
  e <- y
  if (!is.null(X)) e <- y - X %*% b # initialize e
  s <- (crossprod(W, e) / dww) / m # initialize s
  sold <- rep(0, m) # initialize s
  if (weights) {
    p <- apply(W, 2, function(x) {
      cor.test(x, e)$p.value
    })
    logP <- -log10(p)
    lambda <- lambda / ((logP / sum(logP)) * ncol(W)) # 1/logP
  }
  if (is.null(sets)) {
    sets <- as.list(1:m)
  }
  nsets <- length(sets)
  nit <- 0
  delta <- 1
  while (delta > tol) {
    nit <- nit + 1
    for (i in 1:nsets) {
      rws <- sets[[i]]
      lhs <- dww[rws] + lambda[rws] # form lhs
      rhs <- crossprod(W[, rws], e) + dww[rws] * s[rws] # form rhs with y corrected by other effects
      snew <- rhs / lhs
      e <- e - tcrossprod(W[, rws], matrix((snew - s[rws]), nrow = 1)) # update e with current estimate of b
      s[rws] <- snew # update estimates
    }
    gc()
    # if (!is.null(X)) {
    #  bnew <- solve(t(X)%*%X)%*%t(X)%*%e
    #  e  <- e - X%*%(bnew-bold)
    # }
    delta <- sum((s - sold)**2)
    # if (!is.null(X)) delta <- sum((s-sold)**2) + sum((b-bold)**2)
    delta <- delta / sqrt(m)
    sold <- s
    bold <- bnew
    if (nit == maxit) break
    print(paste("Iteration", nit, "delta", delta))
  }
  ghat <- W %*% s
  if (!is.null(X)) b <- (solve(t(X) %*% X) %*% t(X)) %*% (y - ghat) # initialize b
  if (!is.null(X)) yhat <- ghat + X %*% b
  e <- y - yhat
  return(list(s = s, b = b, nit = nit, delta = delta, e = e, yhat = yhat, g = ghat))
}


rsolve <- function(y = NULL, X = NULL, Glist = NULL, ids = NULL, rsids = NULL, sets = NULL, lambda = NULL, weights = FALSE, maxit = 500, tol = 0.0000001) {
  cls <- match(rsids, Glist$rsids)

  maf <- unlist(Glist$maf)[cls]
  meanW <- 2 * maf
  sdW <- sqrt(2 * maf * (1 - maf))
  gc()


  if (!is.null(ids)) yt <- y[ids]
  if (!is.null(ids)) Xt <- as.matrix(X[ids, ])

  n <- Glist$n # number of observations
  m <- Glist$m # number of markers
  if (!is.null(rsids)) m <- length(rsids)
  rwsW <- 1:Glist$n
  if (!is.null(ids)) rwsW <- match(ids, Glist$ids)
  b <- bold <- bnew <- NULL
  if (!is.null(X)) {
    b <- (solve(t(Xt) %*% Xt) %*% t(Xt)) %*% yt # initialize b
    bold <- rep(0, ncol(Xt)) # initialize b
  }
  e <- yt
  if (!is.null(Xt)) e <- yt - Xt %*% b # initialize e

  if (length(lambda) == 1) {
    lambda <- rep(lambda, m)
  }
  dww <- rep(0, m) # initialize diagonal elements of the W'W matrix
  s <- rep(0, m) # initialize diagonal elements of the W'W matrix

  current <- 0
  bfW <- file(Glist$fnRAW, "rb")
  for (i in 1:m) {
    where <- (cls[i] - current - 1) * Glist$n
    current <- cls[i]
    seek(bfW, where = where, origin = "current", rw = "read")
    w <- as.double(readBin(bfW, "raw", n = Glist$n, size = 1, endian = "little"))
    # w[w>0] <- as.vector(scale(w[w>0]))
    w[w > 0] <- (w[w > 0] - 1 - meanW[i]) / sdW[i]
    dww[i] <- sum(w[rwsW]**2)
    if (!dww[i] == 0) s[i] <- (sum(w[rwsW] * e) / dww[i]) / m # initialize s
  }
  close(bfW)
  gc()
  s[dww == 0] <- 0
  sold <- rep(0, m) # initialize s
  nit <- 0
  delta <- 1
  while (delta > tol) {
    nit <- nit + 1
    bfW <- file(Glist$fnRAW, "rb")
    current <- 0
    for (i in 1:m) {
      where <- (cls[i] - current - 1) * Glist$n
      current <- cls[i]
      seek(bfW, where = where, origin = "current", rw = "read")
      w <- as.double(readBin(bfW, "raw", n = Glist$n, size = 1, endian = "little"))
      # w[w>0] <- as.vector(scale(w[w>0]))
      w[w > 0] <- (w[w > 0] - 1 - meanW[i]) / sdW[i]
      lhs <- dww[i] + lambda[i] # form lhs
      rhs <- crossprod(w[rwsW], e) + dww[i] * s[i] # form rhs with y corrected by other effects
      snew <- rhs / lhs
      if (dww[i] == 0) snew <- 0
      e <- e - w[rwsW] * (snew - s[i]) # update e with current estimate of b
      s[i] <- snew # update estimates
    }
    close(bfW)
    gc()
    delta <- sum((s - sold)**2) / sqrt(m)
    sold <- s
    bold <- bnew
    if (nit == maxit) break
    print(paste("Iteration", nit, "delta", delta))
  }
  names(s) <- rsids
  ghat <- rep(0, Glist$n)
  bfW <- file(Glist$fnRAW, "rb")
  current <- 0
  for (i in 1:m) {
    where <- (cls[i] - current - 1) * Glist$n
    current <- cls[i]
    seek(bfW, where = where, origin = "current", rw = "read")
    w <- as.double(readBin(bfW, "raw", n = Glist$n, size = 1, endian = "little"))
    # w[w>0] <- as.vector(scale(w[w>0]))
    w[w > 0] <- (w[w > 0] - meanW[i]) / sdW[i]
    ghat <- ghat + w * s[i]
  }
  close(bfW)
  names(ghat) <- Glist$ids
  if (!is.null(X)) yhat <- ghat[names(y)] + X[names(y), ] %*% b
  e <- y - yhat
  return(list(s = s, b = b, nit = nit, delta = delta, e = e, yhat = yhat, g = ghat[names(y)]))
}

qrSets <- function(W = NULL, sets = NULL, msets = 100, return.level = "Q") {
  m <- ncol(W)
  if (is.null(sets)) sets <- split(1:m, ceiling(seq_along(1:m) / msets))
  qrR <- list()
  for (i in 1:length(sets)) {
    qrW <- qr(W[, sets[[i]]])
    W[, sets[[i]]] <- qr.Q(qrW)
    qrR[[i]] <- qr.R(qrW)
    gc()
  }
  QRlist <- W
  if (return.level == "QR") QRlist <- list(Q = W, R = qrR, sets = sets)
  return(QRlist)
}

plotgs <- function(fit = NULL, s = NULL, sets = NULL) {
  if (is.null(s)) s <- fit$s
  m <- length(s)
  plot(
    y = s, x = 1:m, ylab = "Coefficients", xlab = "Position", col = 1,
    pch = ".", frame.plot = FALSE
  )
  points(y = s[sets], x = (1:m)[sets], col = 2)
}

gsqr <- function(y = NULL, X = NULL, W = NULL, sets = NULL, msets = 100,
                 lambda = NULL, weights = FALSE, maxit = 500, tol = 0.0000001) {
  QRlist <- qgg:::qrSets(W = W, msets = msets, return.level = "QR")
  # lambdaR <- sapply(QRlist$R,function(x){(1/diag(x))**2})
  # lambda <- lambda*lambdaR
  fit <- qgg:::gsru(y = y, X = X, W = QRlist$Q, sets = QRlist$sets, lambda = lambda, weights = weights)
  nsets <- length(QRlist$sets)
  for (i in 1:nsets) {
    rws <- QRlist$sets[[i]]
    # fit$s[rws] <- solve(QRlist$R[[i]])%*%fit$s[rws,1]
    fit$s[rws] <- backsolve(QRlist$R[[i]], fit$s[rws, 1])
  }
  return(fit)
}

fsolve <- function(n = NULL, nr = NULL, rws = NULL, nc = NULL, cls = NULL,
                   scale = TRUE, nbytes = NULL, fnRAW = NULL, ncores = 1, nit = NULL,
                   lambda = NULL, tol = NULL, y = NULL, g = NULL, e = NULL, s = NULL,
                   meanw = NULL, sdw = NULL) {
  OS <- .Platform$OS.type
  if (OS == "windows") fnRAW <- tolower(gsub("/", "\\", fnRAW, fixed = T))
  fit <- .Fortran("solvebed",
    n = as.integer(n),
    nr = as.integer(nr),
    rws = as.integer(rws),
    nc = as.integer(nc),
    cls = as.integer(cls),
    scale = as.integer(scale),
    nbytes = as.integer(nbytes),
    fnRAW = as.character(fnRAW),
    nchars = nchar(as.character(fnRAW)),
    ncores = as.integer(ncores),
    nit = as.integer(nit),
    lambda = as.double(lambda),
    tol = as.double(tol),
    y = as.double(y),
    g = as.double(g),
    e = as.double(e),
    s = as.double(s),
    mean = as.double(meanw),
    sd = as.double(sdw),
    PACKAGE = "qgg"
  )
  return(fit)
}


#'
#' Genomic prediction based on single marker summary statistics
#'
#'
#' @description
#' The gscore function is used for genomic predictions based on single marker summary statistics
#' (coefficients, log-odds ratios, z-scores) and observed genotypes.
#'

#' @param stat matrix of single marker effects
#' @param Glist list of information about genotype matrix
#' @param ids vector of individuals used in the analysis
#' @param rws rows in genotype matrix used in the analysis
#' @param scale logical if TRUE the genotype markers have been scale to mean zero and variance one
#' @param msize number of genotype markers used for batch processing
#' @param ncores number of cores used in the analysis



#' @export
#'


gscore <- function(Glist = NULL, stat = NULL, ids = NULL, impute = TRUE, msize = 100, ncores = 1) {

  # Prepase summary stat
  if (!sum(colnames(stat)[1:3] == c("rsids", "alleles", "af")) == 3) {
    stop("First three columns in data frame stat should be: rsids, alleles, af ")
  }
  rsidsOK <- stat$rsids %in% Glist$rsids
  if (any(!rsidsOK)) {
    warning("Some variants not found in genotype files")
    print(paste("Number of variants missing;", sum(!rsidsOK)))
    print(paste("Number of variants used;", sum(!rsidsOK)))
    stat <- stat[rsidsOK, ]
    stat$rsids <- as.character(stat$rsids)
    stat$alleles <- as.character(stat$alleles)
  }
  S <- stat[, -c(1:3)]
  if (is.vector(S)) S <- as.matrix(S)
  S <- apply(S, 2, as.numeric)
  colnames(S) <- colnames(stat)[-c(1:3)]
  rsids <- as.character(stat$rsids)
  af <- stat$af

  # Prepare input data for mpgrs
  n <- Glist$n
  nbytes <- ceiling(n / 4)
  rws <- 1:n
  if (!is.null(ids)) rws <- match(ids, Glist$ids)
  nr <- length(rws)

  cls <- match(rsids, Glist$rsids)
  nc <- length(cls)

  direction <- as.integer(stat$alleles == Glist$a2[cls])
  # S[direction==0,] <- -S[direction==0,]

  fnRAW <- as.character(Glist$fnRAW)

  nprs <- ncol(S)
  prs <- matrix(0, nrow = nr, ncol = nprs)
  rownames(prs) <- Glist$ids[rws]
  colnames(prs) <- colnames(S)

  m <- nrow(S)
  #sets <- split(1:m, ceiling(seq_along(1:m) / msize))
  #nsets <- length(sets)

  #for (i in 1:nsets) {
  #  nc <- length(sets[[i]])
  #  prsSet <- .Fortran("mpgrs",
  #    n = as.integer(n),
  #    nr = as.integer(nr),
  #    rws = as.integer(rws),
  #    nc = as.integer(nc),
  #    cls = as.integer(cls[ sets[[i]] ]),
  #    nbytes = as.integer(nbytes),
  #    fnRAW = as.character(fnRAW),
  #    nchars = nchar(as.character(fnRAW)),
  #    nprs = as.integer(nprs),
  #    s = matrix(as.double(S[ sets[[i]], ]), nrow = nc, ncol = nprs),
  #    prs = matrix(as.double(0), nrow = nr, ncol = nprs),
  #    af = as.double(af[ sets[[i]] ]),
  #    impute = as.integer(impute),
  #    direction = as.integer(direction[ sets[[i]] ]),
  #    ncores = as.integer(ncores),
  #    PACKAGE = "qgg"
  #  )$prs
  #
  #  prs <- prs + prsSet
  #}
  
   prsSet <- .Fortran("mpgrs",
                          n = as.integer(n),
                          nr = as.integer(nr),
                          rws = as.integer(rws),
                          nc = as.integer(nc),
                          cls = as.integer(cls),
                          nbytes = as.integer(nbytes),
                          fnRAW = as.character(fnRAW),
                          nchars = nchar(as.character(fnRAW)),
                          nprs = as.integer(nprs),
                          s = matrix(as.double(S), nrow = nc, ncol = nprs),
                          prs = matrix(as.double(0), nrow = nr, ncol = nprs),
                          af = as.double(af),
                          impute = as.integer(impute),
                          direction = as.integer(direction),
                          ncores = as.integer(ncores),
                          PACKAGE = "qgg"
       )$prs

  return(prs)
}
