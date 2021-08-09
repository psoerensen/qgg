auc <- function(yobs = NULL, ypred = NULL) {
  n0 <- length(yobs[yobs == 0])
  n1 <- length(yobs[yobs == 1])
  y <- cbind(yobs, ypred)
  y <- y[order(y[, 2], decreasing = TRUE), ]
  y <- cbind(y, seq(from = nrow(y), to = 1))
  rd <- mean(y[y[, 1] == 1, ][, 3])
  auc <- (1 / n0) * (rd - (n1 / 2) - (1 / 2))
  auc
}

rnag <- function(yobs = NULL, ypred = NULL) {
  fit0 <- glm(yobs ~ 1, family = binomial(link = "logit"))
  fit1 <- glm(yobs ~ 1 + ypred, family = binomial(link = "logit"))
  n <- length(yobs)
  LR <- anova(fit1)$Deviance[2]
  L0 <- as.numeric(logLik(fit0))
  r2nag <- (1 - exp(-LR / n)) / (1 - exp(-(-2 * L0) / n))
  return(r2nag)
}

acc <- function(yobs = NULL, ypred = NULL, typeoftrait = "quantitative") {
  fit <- lm(yobs ~ ypred)
  r2 <- summary(fit)$r.squared
  pa <- cor(ypred, yobs)
  mspe <- sum((ypred - yobs)^2) / length(yobs)
  intercept <- fit$coef[1]
  slope <- fit$coef[2]
  aurc <- r2nag <- NA
  if (typeoftrait == "binary") aurc <- auc(yobs = yobs, ypred = ypred)
  if (typeoftrait == "binary") r2nag <- rnag(yobs = yobs, ypred = ypred)
  res <- round(c(pa, r2, r2nag, aurc, intercept, slope, mspe), 3)
  names(res) <- c("Corr", "R2", "Nagel R2", "AUC", "intercept", "slope", "MSPE")
  if(typeoftrait == "quantitative") res <- res[-c(3:4)]
  return(res)
}

fastlm <- function(y = NULL, X = NULL, sets = NULL) {
  XX <- crossprod(X)
  XXi <- chol2inv(chol(XX))
  Xy <- crossprod(X, y)
  coef <- crossprod(XXi, Xy)
  rownames(coef) <- colnames(X)
  yhat <- crossprod(t(X), coef)

  sse <- sum((y - yhat)**2)
  dfe <- length(y) - ncol(X)

  se <- sqrt(sse / dfe) * sqrt(diag(XXi))
  stat <- coef / se
  p <- 2 * pt(-abs(stat), df = dfe)
  names(se) <- colnames(X)

  sigma_e <- sse / dfe
  ftest <- NULL
  if (!is.null(sets)) {
    nsets <- length(sets)
    for (i in 1:nsets) {
      rws <- sets[[i]]
      dfq <- length(rws)
      q <- crossprod(coef[rws, ], crossprod(solve(XXi[rws, rws] * sigma_e), coef[rws, ]))
      pq <- pchisq(q, df = dfq, lower.tail = FALSE)
      pfstat <- pf(q / dfq, dfq, dfe, lower.tail = FALSE)
      ftest <- rbind(ftest, c(q / dfq, dfq, dfe, pfstat))
    }
    colnames(ftest) <- c("F-stat", "dfq", "dfe", "p")
    rownames(ftest) <- names(sets)
  }

  fit <- list(coef = coef, se = se, stat = stat, p = p, ftest = ftest, yhat = yhat)

  return(fit)
}

hwe <- function(Glist=NULL){
  HWE <- vector(length(Glist$n0),mode="list")
  for(CHR in 1:length(Glist$n0)){
    obs <- cbind(Glist$n0[[CHR]],Glist$n1[[CHR]],Glist$n2[[CHR]])
    freq <- (2*obs[,1] + obs[,2]) / (2*rowSums(obs))
    exp <- cbind(rowSums(obs)*freq**2,2*rowSums(obs)*freq*(1-freq) ,rowSums(obs)*(1-freq)**2)
    chi <- rowSums((obs - exp)^2/exp)
    p <- pchisq(chi,1,lower.tail=F)
    names(p) <- Glist$rsids[[CHR]]
    HWE[[CHR]] <- p
   }
  return(HWE)
}
