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
  if (any(is.na(yobs))) stop(paste("NAs found in yobs"))
  if (any(is.na(ypred))) stop(paste("NAs found in ypred"))
  if(is.vector((ypred))) ypred <- as.matrix(ypred) 
  if (is.null(rownames(ypred))) stop(paste("rownames/names missing for ypred"))
  if (is.null(names(yobs))) stop(paste("rownames/names missing for yobs"))
  if (sum(names(yobs)%in%rownames(ypred))==0) stop(paste("No id overlap between yobs and ypred"))
  yobs <- yobs[names(yobs)%in%rownames(ypred)]
  ypred <- ypred[names(yobs),,drop=FALSE]
  result <- NULL
  for(i in 1:ncol(ypred)) {
    fit <- lm(yobs ~ ypred[,i])
    r2 <- summary(fit)$r.squared
    pa <- cor(ypred[,i], yobs)
    mspe <- sum((ypred[,i] - yobs)^2) / length(yobs)
    intercept <- fit$coef[1]
    slope <- fit$coef[2]
    aurc <- r2nag <- NA
    if (typeoftrait == "binary") aurc <- auc(yobs = yobs, ypred = ypred[,i])
    if (typeoftrait == "binary") r2nag <- rnag(yobs = yobs, ypred = ypred[,i])
    res <- round(c(pa, r2, r2nag, aurc, intercept, slope, mspe), 3)
    names(res) <- c("Corr", "R2", "Nagel R2", "AUC", "intercept", "slope", "MSPE")
    if(typeoftrait == "quantitative") res <- res[-c(3:4)]
    result <- rbind(result,res)
  }
  if(ncol(ypred)>1) rownames(result) <- colnames(ypred)
  return(result)
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


###########################################################################################
# Define global parameters
###########################################################################################

# h2x is heritability of target trait
# h2y is heritability of correlated trait
# Nx number of samples for target trait
# Ny number of samples for correlated trait
# rg genetic correlation between target and correlated trait
# Kx is prevalence of target trait
# Ky is prevalence of correlated trait
# Px is case-control proportion of target trait
# Py is case-control proportion of correlated trait


# Single trait prediction of R2 for continous trait:
predict_r2_st = function(h2x, Nx, Me, M){
  b = M / (M + Me)
  r_uu = sqrt(h2x)*sqrt((b * h2x ) / (b * h2x  +  Me / (Nx)))
  return(r_uu**2)
}

# Multi-trait (2-traits) prediction of R2 for continous trait
predict_r2_mt = function(h2x, Nx, h2y, Ny, rg, Me, M){
  b = M / (M + Me)
  r2_uu1 = h2x * (b * h2x ) / (b * h2x  +  Me / ( Nx))
  r2_uu2 = h2y * (b * h2y ) / (b * h2y  +  Me / ( Ny))
  g = matrix(c(r2_uu1, rg^2 * r2_uu2), nrow=2)
  P = matrix(c(r2_uu1, r2_uu1 * r2_uu2 * rg * rg,
               r2_uu2 * r2_uu1 * rg * rg, rg^2 * r2_uu2), nrow=2, byrow=T)
  r_uu = sqrt(t(g) %*% solve(P) %*% g)
  return(r_uu[1,1]**2)
}


# Multi-trait (2-traits) prediction of AUC for binary trait using information on correlated continous trait:
predict_auc_mt_continous = function(h2x, Nx, Kx, Px, h2y, Ny, rg, Me, M){
  b = M / (M + Me)
  zx = dnorm(-qnorm(Kx))
  r2_uu1 = (b * h2x * zx^2) / (b * h2x * zx^2 + (Kx * (1-Kx))^2 * Me / (Px * (1-Px) * Nx))
  r2_uu2 = (b * h2y) / (b * h2y +  Me / Ny)
  g = matrix(c(r2_uu1, rg^2 * r2_uu2), nrow=2)
  P = matrix(c(r2_uu1, r2_uu1 * r2_uu2 * rg * rg,
               r2_uu2 * r2_uu1 * rg * rg, rg^2 * r2_uu2), nrow=2, byrow=T)
  r_uu = sqrt(t(g) %*% solve(P) %*% g)
  i_case  = zx / Kx
  i_con   = -i_case*Kx / (1-Kx)
  auc  = pnorm(((i_case - i_con) * b*h2x * r_uu^2) / sqrt(b*h2x * r_uu^2 * ((1- b *  h2x*r_uu^2*i_case*(i_case-(-qnorm(Kx)))) + (1 - b*h2x*r_uu^2*i_con*(i_con - -qnorm(Kx))))))
  return(auc[1,1])
}

# Single trait prediction of AUC for binary trait:
predict_auc_st = function(h2x, Nx, Kx, Px, Me, M){
  b = M / (M + Me)
  zx = dnorm(-qnorm(Kx))
  r_uu = sqrt((b * h2x * zx^2) / (b * h2x * zx^2 + (Kx * (1-Kx))^2 * Me / (Px * (1-Px) * Nx)))
  i_case  = zx / Kx
  i_con   = -i_case*Kx / (1-Kx)
  auc  = pnorm(((i_case - i_con) * b*h2x * r_uu^2) / sqrt(b*h2x * r_uu^2 * ((1- b*h2x*r_uu^2*i_case*(i_case-(-qnorm(Kx)))) + (1 - b*h2x*r_uu^2*i_con*(i_con - -qnorm(Kx))))))
  return(auc)
}

# Multi-trait (2-traits) prediction of AUC for binary trait
predict_auc_mt_cc = function(h2x, Nx, Kx, Px, h2y, Ny, Ky, Py, rg, Me, M){
  b = M / (M + Me)
  zx = dnorm(-qnorm(Kx))
  zy = dnorm(-qnorm(Ky))
  r2_uu1 = (b * h2x * zx^2) / (b * h2x * zx^2 + (Kx * (1-Kx))^2 * Me / (Px * (1-Px) * Nx))
  r2_uu2 = (b * h2y * zy^2) / (b * h2y * zy^2 + (Ky * (1-Ky))^2 * Me / (Py * (1-Py) * Ny))
  g = matrix(c(r2_uu1, rg^2 * r2_uu2), nrow=2)
  P = matrix(c(r2_uu1, r2_uu1 * r2_uu2 * rg * rg,
               r2_uu2 * r2_uu1 * rg * rg, rg^2 * r2_uu2), nrow=2, byrow=T)
  r_uu = sqrt(t(g) %*% solve(P) %*% g)
  i_case  = zx / Kx
  i_con   = -i_case*Kx / (1-Kx)
  auc  = pnorm(((i_case - i_con) * b*h2x * r_uu^2) / sqrt(b*h2x * r_uu^2 * ((1- b *  h2x*r_uu^2*i_case*(i_case-(-qnorm(Kx)))) + (1 - b*h2x*r_uu^2*i_con*(i_con - -qnorm(Kx))))))
  return(auc)
}


