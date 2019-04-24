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
# ’ @references Sørensen, P., de los Campos, G., Morgante, F., Mackay, T. F., & Sorensen, D. (2015). Genetic control of environmental variation of two quantitative traits of Drosophila melanogaster revealed by whole-genome sequencing. Genetics, 201(2), 487-497.
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
