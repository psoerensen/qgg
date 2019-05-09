####################################################################################################################
#    Module 3: LMM, Bayesian
####################################################################################################################
#'
#' Genomic Feature Model analyses implemented using Bayesian Methods (small data)
#'
#' @description
#' Bayesian Genomic Feature models implemented using Bayesian Methods.
#' Multiple features and multiple traits models can be fitted.
#' Different genetic models (e.g. additive, dominance, gene by gene and gene by environment interactions) can be specified.
#'
#' @details
#' The models are implemented using empirical Bayesian methods. The hyperparameters of the dispersion parameters of the Bayesian model can
#' be obtained from prior information or estimated by maximum likelihood, and conditional on these, the model is fitted using
#' Markov chain Monte Carlo. Furthermore, a spectral decomposition of genomic feature relationship matrices plays an
#' important computational role in the Markov chain Monte Carlo strategy implemented.
#'
#' @param y is a matrix of phenotypes
#' @param g is a list of G matrices (g$G), prior (co)variances (g$sigma), prior degrees of freedom (g$sigma),  prior distribution (g$ds)
#' @param nsamp is the number of samples after burnin
#' @param nburn is the number of burnin samples
#' @param nsave is the number of samples to save
#' @param tol is the tolerance
#' @return Returns a list including
#' \item{sigmas}{posterior samples variance components}
#' \item{alphas}{posterior samples alphas}
#' \item{mus}{posterior samples mu}
#' \item{logCPO}{log CPOs}
#' \item{g}{g is a list of G matrices (g$G), prior (co)variances (g$sigma), prior degrees of freedom (g$sigma),  prior distribution (g$ds)}
#' @author Peter Sørensen
#' @references Sørensen, P., de los Campos, G., Morgante, F., Mackay, T. F., & Sorensen, D. (2015). Genetic control of environmental variation of two quantitative traits of Drosophila melanogaster revealed by whole-genome sequencing. Genetics, 201(2), 487-497.
#'
#' @export
#'

gbayes <- function(y = NULL, g = NULL, nsamp = 50, nburn = 10, nsave = 10000, tol = 0.001) {

  # nsamp is the number of samples
  y <- as.matrix(y)
  n <- nrow(y) # number of observation
  nt <- ncol(y) # number of traits
  mu <- colMeans(y)
  samples <- as.integer(seq(nburn + 1, nsamp, length.out = min(nsave, nsamp)))

  psum <- matrix(0, nrow = n, ncol = nt)

  nset <- length(g$G) # number of sets
  for (i in 1:nt) { # genetic values (n*nset)
    g$values[[i]] <- matrix(0, nrow = n, ncol = nset)
  }
  if (nt > 1) sigmae <- g$sigma[[nset + 1]][1, 1]
  if (nt == 1) sigmae <- g$sigma[[nset + 1]]
  sigma <- g$sigma
  mus <- matrix(NA, nrow = nsamp, ncol = nt)

  alphas <- vector(length = nt, mode = "list")

  U <- D <- sigmas <- vector(length = nset, mode = "list")
  for (i in 1:nset) {
    e <- eigen(g$G[[i]]) # eigen value decomposition of the matrix G
    ev <- e$values
    U[[i]] <- e$vectors[, ev > tol] # keep eigen vector if ev > tol
    D[[i]] <- e$values[ev > tol] # keep eigen value if ev > tol
    sigmas[[i]] <- matrix(NA, nrow = nsamp, ncol = (nt * (nt + 1)) / 2)
    for (j in 1:nt) {
      as <- matrix(NA, nrow = sum(ev > tol), ncol = length(samples))
      colnames(as) <- as.character(samples)
      alphas[[j]][[i]] <- as
    }
  }

  for (i in 1:nsamp) { # loop over the number of samples (nsamp)

    for (j in 1:nset) { # loop over the number of sets
      rhs <- NULL

      for (t in 1:nt) {
        yadj <- y[, t] - mu[t] - rowSums(as.matrix(g$values[[t]][, -j]))
        rhst <- crossprod(U[[j]], yadj)
        rhs <- cbind(rhs, rhst)
      }

      Vi <- solve(sigma[[j]])
      a <- matrix(0, nrow = nrow(rhs), ncol = nt)

      for (k in 1:nrow(rhs)) {
        iC <- solve(diag(1, nt) + (sigmae / D[[j]][k]) * Vi)
        ahat <- iC %*% rhs[k, ]
        a[k, ] <- mvrnorm(n = 1, mu = ahat, Sigma = iC * sigmae)
      }

      for (t in 1:nt) {
        g$values[[t]][, j] <- U[[j]] %*% a[, t]
      }
      df <- nrow(a) + g$df[j]

      if (any(i %in% samples)) {
        for (t in 1:nt) {
          alphas[[t]][[j]][, as.character(i)] <- a[, t]
        }
      }

      # inverse chisquare
      if (nt == 1) {
        # scg <- sum((1 / D[[j]]) * a**2) + g$sigma[[j]] * g$df[j]
        scg <- sum((1 / D[[j]]) * a**2) + (g$sigma[[j]] * (g$df[j] + 2)) / g$df[j] # => S = (mode * (df + 2)) / df
        sigma[j] <- scg / rchisq(n = 1, df = df, ncp = 0)
        sigmas[[j]][i, ] <- sigma[j]
      }

      # inverse wishart
      if (nt > 1) {
        if (g$ds[j] == "invWish") {
          S <- t(a * (1 / D[[j]])) %*% a + g$sigma[[j]] * (g$df[j] + nt + 1) # => S = mode * (df + nt + 1)
          S <- riwish(df, S)
          sigma[[j]] <- S
          sigmas[[j]][i, ] <- S[as.vector(lower.tri(S, diag = TRUE))]
          # print(c(j, g$ds[j]))
        }

        if (g$ds[j] == "invChisq") {
          scg <- t(a * (1 / D[[j]])) %*% a + (g$sigma[[j]] * (g$df[j] + 2)) / g$df[j] # => S = (mode * (df + 2)) / df

          for (t in 1:nt) {
            sigma[[j]][t, t] <- scg[t, t] / rchisq(n = 1, df = df, ncp = 0)
          }

          sigmas[[j]][i, ] <- sigma[[j]][as.vector(lower.tri(sigma[[j]], diag = TRUE))]
        }
      }
    }

    for (t in 1:nt) {
      yadj <- y[, t] - rowSums(as.matrix(g$values[[t]]))
      rhs <- sum(yadj)
      lhs <- (n + sigmae / 100000)
      mu[t] <- rnorm(1, mean = rhs / lhs, sd = 1 / sqrt(lhs))
    }
    mus[i, ] <- mu
    print(i)
    for (t in 1:nt) {
      e <- y[, t] - mu[t] - rowSums(as.matrix(g$values[[t]]))
      p <- (1 / sqrt(2 * pi)) * exp(-0.5 * (e**2))
      if (i > nburn) psum[, t] <- psum[, t] + 1 / p
    }
  }

  # summary of model fit
  logCPO <- NULL
  for (t in 1:nt) {
    logCPO[t] <- sum(log((nsamp - nburn) * (1 / psum[, t])))
  }

  return(list(sigmas = sigmas, mus = mus, logCPO = logCPO, g = g, alphas = alphas, nsamp = nsamp, nburn = nburn)) # return posterior samples of sigma
}



################################################################################
# Bayesian lasso function
################################################################################

blasso <- function(y = NULL, X = NULL, lambda0 = NULL, sigma0 = NULL, nsamp = 100) {
  n <- length(y) # number of observations
  p <- ncol(X) # number of regression variables
  dxx <- colSums(X**2) # diagonal elements of the X'X matrix
  b <- as.numeric(rep(0, p)) # initialize b
  mu <- 0 # initilaize mu
  invtau2 <- rep(0, p) # initialize tau2
  lambda2 <- p # initialize lambda2
  lrate0 <- 0.1 # lambda2 hyperparameter rate
  lshape0 <- 1 # lambda2 hyperparameter shape

  Xs <- as.vector(X)

  e <- as.vector(y - mu - X %*% b) # initialize residuals
  sigma2 <- var(e) / 2 # initialize sigma2

  for (i in 1:nsamp) {

    # sample mu and update residuals
    muC <- mu
    mu <- rnorm(1, mean = mean(e + muC), sd = sigma2 / n)
    e <- as.vector(e - mu + muC)

    # sample b and update residuals
    bC <- b
    xxs <- as.numeric(dxx + invtau2) # Park & Casella
    for (j in 1:p) {
      rhs <- as.numeric(sum(X[, j] * e) + dxx[j] * bc[j])
      b[j] <- rnorm(1, mean = rhs / xxs[j], sd = sigma2 / xxs[j])
      e <- as.numeric(e - X[, j] * (b[j] - bC[j]))
    }
    e <- as.numeric(y - mu - X %*% b) # update residuals

    # sample invtau2
    mutau <- sqrt(lambda2 * sigma2) / abs(b) # Park & Casella
    lambdatau <- rep(lambda2, p) # Park & Casella
    invtau2 <- rinvgauss(p, mu = mutau, lambda = lambdatau)

    # update lambda
    shl2 <- p + lshape0 # Park & Casella
    ratel2 <- sum(1 / (invtau2)) / 2 + lrate0 # Park & Casella
    lambda2 <- rgamma(n = 1, shape = shl2, rate = ratel2) # Park & Casella
    print(lambda2)
    if (!is.null(lambda0)) lambda2 <- lambda0

    # sample sigma
    she <- (n + p) / 2 # Park & Casella ((n-1+p)/2)
    sce <- sum(e**2) / 2 + sum((b * invtau2 * b)) / 2 # Park & Casella
    sigma2 <- 1 / rgamma(1, shape = she, rate = sce) # Park & Casella
    she <- n + p # Gustavo
    sce <- sum(e**2) + sum((b * invtau2 * b)) # Gustavo
    sigma2 <- sce / rchisq(n = 1, df = she) # Gustavo
    if (!is.null(sigma0)) sigma2 <- sigma0
    print(c(sum(e**2), sum((b * (1 / invtau2) * b)), sum((1 / invtau2))))
    plot(b)
    print(c(i, mu, lambda2, sigma2, cor(y, X %*% b)))
  }
}


#' # Simulate data and test function
#'
#' X <- matrix(rnorm(10000000),nrow=1000)
#' set <- sample(1:ncol(X),10)
#' y <- rowSums(X[,set]) + rnorm(nrow(X),mean=0,sd=1)
#' n <- length(y)
#' h2 <- var(rowSums(X[,set]))/var(y)
#'
#' lambda0 <- (sum(colSums(X**2))/n)*(1-h2)/h2
#'
#' blasso(y=y, X=X, lambda0=lambda0)

################################################################################
# Multi Class Bayesian Regression function
################################################################################



mcbr <- function(y = NULL, X = NULL, nc = NULL, l1 = NULL, l2 = NULL, phi = NULL, nsamp = 100) {
  n <- length(y) # number of observations
  p <- ncol(X) # number of regression variables

  dxx <- colSums(X**2) # diagonal elements of the X'X matrix

  b <- as.numeric(rep(0, p)) # initialize b
  mu <- 0 # initilaize mu

  if (is.null(l1)) l1 <- 1 / (10**((1:nc) - 4)) # prior shape parameter lambda
  if (is.null(l2)) l2 <- rep(1 / 10**2, nc) # prior rate parameter lambda
  if (is.null(phi)) phi0 <- c(0.5, 0.25, 0.15, 0.05, 0.025, 0.01, 0.001) # prior phi
  phi <- phi0
  logphi <- log(phi) # prior phi on the log scale

  lambdak <- rep(1000, nc) # initialize lambdak one for each class
  lambda <- rep(1000, p) # initialize lambda one for each regressor

  g <- rep(1, p) # initialize class indicator variable
  cors <- abs(cor(y, X))
  quants <- quantile(cors, probs = 1 - phi)
  for (i in 1:nc) {
    #  g[cors>quants[i]] <- i
  }

  e <- as.vector(y - mu - X %*% b) # initialize residuals

  sigma2 <- var(e) / 2 # initialize sigma2


  for (i in 1:nsamp) {

    # sample mu and update residuals
    muC <- mu
    mu <- rnorm(1, mean = mean(e + muC), sd = sigma2 / n)
    e <- as.vector(e - mu + muC)

    # sample b and update residuals
    bC <- b
    xxs <- as.numeric(dxx + lambda * sigma2)
    for (j in 1:p) {
      rhs <- as.numeric(sum(X[, j] * e) + dxx[j] * bC[j])
      b[j] <- rnorm(1, mean = rhs / xxs[j], sd = sigma2 / (xxs[j]))
      e <- as.numeric(e - X[, j] * (b[j] - bC[j]))
    }
    e <- as.numeric(y - mu - X %*% b) # update residuals

    # sample lambda
    b2 <- b**2
    l1k <- l1
    l2k <- l2
    for (k in 1:nc) {
      if (sum(g == k) > 0) l1k[k] <- l1k[k] + sum(g == k) / 2
      if (sum(g == k) > 0) l2k[k] <- l2k[k] + sum(b2[g == k]) / 2
      lambdak[k] <- rgamma(n = 1, shape = l1k[k], rate = l2k[k])
      if (sum(g == k) > 0) lambda[g == k] <- lambdak[k]
    }
    loglambdak <- log(lambdak)

    # sample gamma
    for (j in 1:p) {
      probs <- -0.5 * b2[j] * lambdak + logphi + 0.5 * loglambdak
      g[j] <- (1:nc)[rmultinom(1, 1, prob = exp(probs)) == 1]
    }

    # sample phi
    # phi <- rep(0.0001,nc)
    rws <- as.numeric(names(table(g)))
    # phi[rws] <- table(g)/p
    print(c(i, phi))

    # sample sigma
    she <- n
    sce <- sum(e**2)
    sigma2 <- sce / rchisq(n = 1, df = she)

    plot(b)
    print(table(g))
    print(c(mu, sigma2))
  }
}


#' # Simulate data and test function
#'
#' X <- matrix(rnorm(10000000),nrow=1000)
#' set <- sample(1:ncol(X),10)
#' y <- rowSums(X[,set]) + rnorm(nrow(X),mean=0,sd=1)
#' h2 <- var(rowSums(X[,set]))/var(y)
#'
#' mcbr(y=y, X=X, nc=7)


################################################################################
# SSVS function
################################################################################



hgprior <- list(sce0 = 0.001, scg0 = 0.001, dfe0 = 4, dfg0 = 4)
p1 <- 0.001
g0 <- NULL

hssvs <- function(y = NULL, X = NULL, set = NULL, p1 = NULL, g0 = NULL, nsamp = 100, hgprior = list(sce0 = 0.001, scg0 = 0.001, dfe0 = 4, dfg0 = 4)) {
  n <- length(y) # number of observations
  p <- ncol(X) # number of regression variables
  dxx <- colSums(X**2) # diagonal elements of the X'X matrix
  p0 <- 1 - p1 # prior proportion of elements in gamma to be set to one

  b <- rep(0, p) # initialize b
  g <- rep(0, p) # intialize g

  nset <- length(set) # number of sets
  pset <- sapply(set, length)
  for (i in 1:nset) {
    n1 <- as.integer(pset[i] * p1)
    cset <- abs(cor(y, X[, set[[i]]]))
    o <- order(cset, decreasing = TRUE)[1:n1]
    o <- set[[i]][o]
    g[o] <- 1
  }

  bset <- NULL # beta's within sets
  gset <- NULL # g's within sets
  for (i in 1:nset) {
    bset[[i]] <- b[set[[i]]]
    gset[[i]] <- g[set[[i]]]
  }

  mu <- 0 # initilaize mu
  e <- as.vector(y) # initialize residuals
  sigma2 <- var(e) / 2 # initialize sigma2
  g1 <- 1 # initialize slab variance
  if (is.null(g0)) g0 <- 1e-09 # initialize spike variance

  g0 <- rep(g0, nset) # spike variance for each set
  g1 <- rep(g1, nset) # slab variance for each set

  Xs <- matrix(0, nrow = n, ncol = nset) # matrix (n*nset) of genomic values one for each set


  for (i in 1:nsamp) {
    for (j in 1:nset) {
      cls <- set[[j]]
      samp <- ssvs2(e = e, X = X[, cls], b = bset[[j]], dxx = dxx[cls], mu = mu, g = gset[[j]], sigma2 = sigma2, p0 = p0, p1 = p1, g0 = g0[j], g1 = g1[j], hgprior = hgprior)
      e <- samp$e
      bset[[j]] <- samp$b
      mu <- samp$mu
      sigma2 <- samp$sigma2
      gset[[j]] <- samp$g
      g0[j] <- samp$g0
      g1[j] <- samp$g1
      Xs[, j] <- X[, cls] %*% bset[[j]]
      print(c("round", i, j, "g", set[[j]][gset[[j]] == 1]))
    }
    h2 <- apply(Xs, 2, var)
    print(c("h2", h2))
    h2 <- h2 / sum(h2)
    print(c("h2", h2))

    barplot(h2)
    print("")
    print(c("round", i, "gset", sapply(gset, sum)))
    print(c("round", i, "mu, sigma2", c(mu, sigma2)))
    print(c("round", i, "g0", g0))
    print(c("round", i, "g1", g1))
    print(c("round", i, "cor", c(cor(y, Xs), cor(y, rowSums(Xs)), cor(y, rowSums(Xs[, 1:2])))))
  }
}

ssvs2 <- function(e = e, X = X, b = b, dxx = dxx, mu = mu, g = g, sigma2 = sigma2, p0 = p0, p1 = p1, g0 = g0, g1 = g1, hgprior = hgprior) {
  n <- length(e) # number of observations
  p <- ncol(X) # number of regression variables

  sce0 <- hgprior$sce0 # prior residual sums of squares
  scg0 <- hgprior$scg0 # prior slab residual sums of squares
  dfe0 <- hgprior$dfe0 # prior residual degrees of freedom
  dfg0 <- hgprior$dfg0 # prior slab residual degrees of freedom


  # sample b and update residuals
  bC <- b
  shrinkage <- rep(sigma2 / g1, times = p)
  shrinkage[g == 0] <- sigma2 / g0
  xxs <- dxx + shrinkage
  for (j in 1:p) {
    rhs <- sum(X[, j] * e) + dxx[j] * bC[j]
    b[j] <- rnorm(1, mean = rhs / xxs[j], sd = sigma2 / xxs[j])
    e <- e - X[, j] * (b[j] - bC[j])
  }

  # sample g0
  scg00 <- sum((b**2)[g == 0]) + scg0
  dfg00 <- sum(1 - g) + dfg0 - 2
  g0 <- scg00 / rchisq(n = 1, df = dfg00, ncp = 0)

  # sample g1
  scg1 <- sum((b**2)[g == 1]) + scg0
  dfg1 <- sum(g) + dfg0 - 2
  for (k in 1:1000) { # scaling needs to be revised
    g1 <- scg1 / rchisq(n = 1, df = dfg1, ncp = 0)
    if (g1 > g0) break # scaling needs to be revised
    if (k == 1000) g1 <- g0 # scaling needs to be revised
  } # scaling needs to be revised


  # sample g
  ratio1 <- p1 * (exp(-(b**2) / g1) / sqrt(g1))
  ratio0 <- p0 * (exp(-(b**2) / g0) / sqrt(g0))
  ratio <- ratio1 / ratio0
  u <- 1 / (1 - runif(p, 0, 1))
  g[1:p] <- 0
  g[ratio > u] <- 1


  # sample mu and update residuals
  muC <- mu
  mu <- rnorm(1, mean = mean(e + muC), sd = sigma2 / n)
  e <- e - mu + muC

  # sample sigma
  sce <- sum(e**2) + sce0
  dfe <- n + dfe0 - 2
  sigma2 <- sce / rchisq(n = 1, df = dfe, ncp = 0)


  return(list(e = e, mu = mu, b = b, sigma2 = sigma2, g = g, g0 = g0, g1 = g1))
}

#' # Simulate data and test function
#'
#' X <- matrix(rnorm(10000000),nrow=1000)
#' set1 <- sample(1:ncol(X),5)
#' set2 <- sample(1:ncol(X),5)
#' set1a <- unique(c(set1,sample(1:ncol(X),100)))
#' set2a <- unique(c(set2,sample(1:ncol(X),100)))
#' set3 <- (1:ncol(X))[-unique(c(set1a,set2a))]
#'
#' y <- rowSums(X[,set1]) + rowSums(X[,set2]) + rnorm(nrow(X),mean=0,sd=1)
#'
#' varE <- var( y - rowSums(X[,set1]) - rowSums(X[,set2]) )
#' varE
#'
#' h2 <- ( var(rowSums(X[,set1])) + var(rowSums(X[,set2])) )/var(y)
#' h2
#'
#' set <- list(set1=set1a,set2=set2a,set3=set3)
#'
#' res <- hssvs(y=y,X=X,set=set, p1=0.01,g0=0.0000001, nsamp=100)



ssvs <- function(y = NULL, X = NULL, p1 = NULL, g0 = NULL, nsamp = 100, hgprior = list(sce0 = 0.001, scg0 = 0.001, dfe0 = 4, dfg0 = 4)) {
  n <- length(y) # number of observations
  p <- ncol(X) # number of regression variables
  dxx <- colSums(X**2) # diagonal elements of the X'X matrix
  p0 <- 1 - p1 # prior proportion of elements in gamma to be set to one

  b <- rep(0, p) # initialize b
  g <- rep(0, p) # intialize g
  Xy <- colSums(X * y)
  n1 <- as.integer(p * p1)
  o <- order(Xy, decreasing = TRUE)[1:n1]
  g[o] <- 1
  mu <- 0 # initilaize mu
  sigma2 <- var(y) / 2 # initialize sigma2
  g1 <- 1 # initialize slab variance
  if (is.null(g0)) g0 <- 1e-09 # initialize spike variance

  sce0 <- hgprior$sce0 # prior residual sums of squares
  scg0 <- hgprior$scg0 # prior slab residual sums of squares
  dfe0 <- hgprior$dfe0 # prior residual degrees of freedom
  dfg0 <- hgprior$dfg0 # prior slab residual degrees of freedom

  e <- as.vector(y - mu - X %*% b) # initialize residuals

  Xs <- as.vector(X)




  for (i in 1:nsamp) {

    # sample mu and update residuals
    muC <- mu
    mu <- rnorm(1, mean = mean(e + muC), sd = sigma2 / n)
    e <- e - mu + muC

    # sample b and update residuals
    bC <- b
    shrinkage <- rep(sigma2 / g1, times = p)
    shrinkage[g == 0] <- sigma2 / g0
    xxs <- dxx + shrinkage
    for (j in 1:p) {
      rhs <- sum(X[, j] * e) + dxx[j] * bC[j]
      b[j] <- rnorm(1, mean = rhs / xxs[j], sd = sigma2 / xxs[j])
      e <- e - X[, j] * (b[j] - bC[j])
    }


    # sample sigma
    sce <- sum(e**2) + sce0
    dfe <- n + dfe0 - 2
    sigma2 <- sce / rchisq(n = 1, df = dfe, ncp = 0)
    # sigma2 <- sum(e**2)/rchisq(n=1, df=n-2, ncp=0)


    # sample g1
    scg1 <- sum((b**2)[g == 1]) + scg0
    dfg1 <- sum(g) + dfg0 - 2
    g1 <- scg1 / rchisq(n = 1, df = dfg1, ncp = 0)
    # dfg <- max(1,sum(g)-2)
    # sc1 <- max(sum((b**2)[g==1]),g0)
    # g1 <- sc1/rchisq(n=1, df=dfg, ncp=0)

    # sample g0
    scg00 <- sum((b**2)[g == 0]) + scg0
    dfg00 <- sum(1 - g) + dfg0 - 2
    g00 <- scg00 / rchisq(n = 1, df = dfg00, ncp = 0)
    print(c(g1, g00))
    g0 <- g00

    # sample g
    ratio1 <- p1 * (exp(-(b**2) / g1) / sqrt(g1))
    ratio0 <- p0 * (exp(-(b**2) / g0) / sqrt(g0))
    ratio <- ratio1 / ratio0
    u <- 1 / (1 - runif(p, 0, 1))
    g[1:p] <- 0
    g[ratio > u] <- 1

    print("")
    print(c("round", i, (1:p)[g == 1]))
    print(c("round", i, c(mu, sigma2, g1)))
    print(c("round", i, b[g == 1]))
  }
}


#' # Simulate data and test function
#'
#' X <- matrix(rnorm(10000000),nrow=1000)
#' set <- sample(1:ncol(X),10)
#' y <- rowSums(2*X[,set]) + rnorm(nrow(X),mean=0,sd=1)
#'
#' res <- ssvs(y=y,X=X,p1=0.001,g0=0.0000001, nsamp=100)
#' sort(set)

################################################################################
# Bayesian group lasso
################################################################################

bglasso <- function(y = NULL, X = NULL, g = NULL, nsamp = 100, hgprior = list(sce0 = 0.001, scg0 = 0.001, dfe0 = 4, dfg0 = 4)) {
  n <- length(y) # number of observations
  nb <- ncol(X) # number of regression variables
  ng <- length(unique(g)) # number of groups
  dxx <- colSums(X**2) # diagonal elements of the X'X matrix

  b <- rep(0, nb) # initialize b
  mu <- 0 # initilaize mu
  sigma <- 1 # initialize sigma
  g1 <- rep(1 / nb, ng) # initialize slab variance
  g1 <- rep(0, ng) # initialize slab variance
  for (j in 1:ng) {
    g1[j] <- 1 / sum(g == j)
  }
  print(g1)
  # g1 <- c(1,0.000001)                     # initialize slab variance

  sce0 <- hgprior$sce0 # prior residual sums of squares
  scg0 <- hgprior$scg0 # prior slab residual sums of squares
  dfe0 <- hgprior$dfe0 # prior residual degrees of freedom
  dfg0 <- hgprior$dfg0 # prior slab residual degrees of freedom

  e <- as.vector(y - mu - X %*% b) # initialize residuals

  for (i in 1:nsamp) {

    # sample mu and update residuals
    mu0 <- mu
    mu <- rnorm(1, mean = mean(e + mu0), sd = sigma / n)
    e <- e - mu + mu0

    # sample b and update residuals
    b0 <- b
    shrinkage <- rep(0, times = nb)
    for (j in 1:ng) {
      shrinkage[g == j] <- sigma / g1[j]
    }
    xxs <- dxx + shrinkage
    for (j in 1:nb) {
      rhs <- sum(X[, j] * e) + dxx[j] * b0[j]
      b[j] <- rnorm(1, mean = rhs / xxs[j], sd = sigma / xxs[j])
      e <- e - X[, j] * (b[j] - b0[j])
    }

    # sample sigma
    sce <- sum(e**2) + sce0
    dfe <- n + dfe0 - 2
    sigma <- sce / rchisq(n = 1, df = dfe, ncp = 0)
    # sigma <- sum(e**2)/rchisq(n=1, df=n-2, ncp=0)


    # sample g1
    for (j in 1:ng) {
      scg1 <- sum((b**2)[g == j]) + scg0
      dfg1 <- sum(g == j) + dfg0 - 2
      g1[j] <- scg1 / rchisq(n = 1, df = dfg1, ncp = 0)
    }

    # sample g

    print("")
    print(c("round", i, (1:nb)[g == 1]))
    print(c("round", i, c(mu, sigma, g1)))
    print(c("round", i, b[g == 1]))
  }
}

#' # Simulate data and test function
#'
#' X <- matrix(rnorm(10000000),nrow=1000)
#' set1 <- sample(1:ncol(X),50)
#' set2 <- sample(1:ncol(X),50)
#' set3 <- sample(1:ncol(X),50)
#' set4 <- sample(1:ncol(X),50)
#' set5 <- sample(1:ncol(X),50)
#' g <- rep(6,times=ncol(X))
#' g[set1] <- 1
#' g[set2] <- 2
#' g[set3] <- 3
#' g[set4] <- 4
#' g[set5] <- 5
#'
#' set <- unique(c(set1,set2,set3,set4,set5))
#'
#' y <- rowSums(X[,set]) + rnorm(nrow(X),mean=0,sd=20)
#'
#' res <- bglasso(y=y,X=X,g=g, nsamp=10)
