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
#' @param Glist list of information about genotype matrix stored on disk
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



#' @export


gsolve <- function(y = NULL, X = NULL, Glist = NULL, W = NULL, ids = NULL, rsids = NULL,
                   sets = NULL, validate = NULL, scale = TRUE, lambda = NULL, weights = FALSE,
                   maxit = 500, tol = 0.00001, method = "gsru", ncores = 1) {
  if (!is.null(W)) {
    if (method == "gsru") {
      fit <- gsru(
        y = y, W = W, X = X, sets = sets, lambda = lambda,
        weights = weights, maxit = maxit, tol = tol
      )
    }
    if (method == "gsqr") {
      fit <- gsqr(
        y = y, W = W, X = X, sets = sets, msets = 100,
        lambda = lambda, weights = weights, maxit = maxit, tol = tol
      )
    }
    return(fit)
  }
  if (!is.null(Glist)) {

    nt <- 1
    if(!is.list(y)){
      ids <- names(y)
      if(!is.null(ids)) rws <- match(ids, Glist$ids)
    }
    if(is.list(y)){
      nt <- length(y)
      ids <- names(y[[1]])
      if(!is.null(ids)) rws <- match(ids, Glist$ids)
    }
    
    for ( chr in 1:Glist$nchr) {
      cls <- 1:Glist$mchr[[chr]]
      cls <- cls[ Glist$rsids[[chr]]%in%rsids ]
      af <- Glist$af[[chr]][cls]
      if(nt==1) {
        b <- rep(0.0, length(cls))
        lambda <- rep(lambda, length(cls))
        fit <- .Call("_qgg_solvebed", Glist$bedfiles[chr], Glist$n, cls, maxit, af, b, lambda, y)
      }
      if(nt>1) {
        b <- rep(0.0, length(cls))
        b <- rep(list(b),nt)
        if(length(lambda)!=nt) stop("lambda not of length nt= number of traits")
        lambda <- lapply(lambda,function(x){rep(x,length(cls))})
        fit <- .Call("_qgg_mtsolvebed", Glist$bedfiles[chr], Glist$n, cls, maxit, af, b, lambda, y)
      }
    }
    return(fit)
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
    delta <- sum((s - sold)**2)
    delta <- delta / sqrt(m)
    sold <- s
    bold <- bnew
    if (nit == maxit) break
    message(paste("Iteration", nit, "delta", delta))
  }
  ghat <- W %*% s
  if (!is.null(X)) b <- (solve(t(X) %*% X) %*% t(X)) %*% (y - ghat) # initialize b
  if (!is.null(X)) yhat <- ghat + X %*% b
  e <- y - yhat
  return(list(s = s, b = b, nit = nit, delta = delta, e = e, yhat = yhat, g = ghat))
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


gsqr <- function(y = NULL, X = NULL, W = NULL, sets = NULL, msets = 100,
                 lambda = NULL, weights = FALSE, maxit = 500, tol = 0.0000001) {
  QRlist <- qrSets(W = W, msets = msets, return.level = "QR")
  fit <- gsru(y = y, X = X, W = QRlist$Q, sets = QRlist$sets, lambda = lambda, weights = weights)
  nsets <- length(QRlist$sets)
  for (i in 1:nsets) {
    rws <- QRlist$sets[[i]]
    fit$s[rws] <- backsolve(QRlist$R[[i]], fit$s[rws, 1])
  }
  return(fit)
}

