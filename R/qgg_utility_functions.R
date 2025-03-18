#' Compute AUC
#' @description
#' Compute Area Under the Curve (AUC)
#' @keywords internal
#' @param yobs is a vector of observed phenotypes
#' @param ypred is a vector of predicted phenotypes
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

#' Compute Nagelkerke R2
#' @description
#' Compute Nagelkerke R2 ££ perhaps add: r2nag?
#' @keywords internal
#' @param yobs is a vector of observed phenotypes
#' @param ypred is a vector of predicted phenotypes
rnag <- function(yobs = NULL, ypred = NULL) {
  fit0 <- glm(yobs ~ 1, family = binomial(link = "logit"))
  fit1 <- glm(yobs ~ 1 + ypred, family = binomial(link = "logit"))
  n <- length(yobs)
  LR <- anova(fit1)$Deviance[2]
  L0 <- as.numeric(logLik(fit0))
  r2nag <- (1 - exp(-LR / n)) / (1 - exp(-(-2 * L0) / n))
  return(r2nag)
}

#' Compute prediction accuracy for a quantitative or binary trait
#' @description
#' Compute prediction accuracy for a quantitative or binary trait
#' @param yobs is a vector of observed phenotypes
#' @param ypred is a vector of predicted phenotypes
#' @param typeoftrait is a character with possible values "binary" or "quantitative" (default)
#' @export
acc <- function(yobs = NULL, ypred = NULL, typeoftrait = "quantitative") {
  if (any(is.na(yobs))) stop(paste("NAs found in yobs"))
  if (any(is.na(ypred))) stop(paste("NAs found in ypred"))
  if(is.vector((ypred))) ypred <- as.matrix(ypred) 
  if (is.null(rownames(ypred))) stop(paste("rownames/names missing for ypred"))
  if (is.null(names(yobs))) stop(paste("rownames/names missing for yobs"))
  if (sum(names(yobs)%in%rownames(ypred))==0) stop(paste("No id overlap between yobs and ypred"))
  if (any(duplicated(names(yobs)))) stop(paste("Duplicated names in yobs"))
  if (any(duplicated(rownames(ypred)))) stop(paste("Duplicated rownames in ypred"))
  
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

#' Perform Hardy Weinberg Equilibrium Test
#' @description
#' Perform Hardy Weinberg Equilibrium Test ££(HWE)
#' @keywords internal
#' @param Glist is a list structure with information about genotypes stored on disk
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

#' Expected R2 for single trait prediction of a continuous trait
#'
#' Computes the expected R2 value for the single trait prediction of a continuous trait.
#'
#' @param h2x Heritability of the target trait.
#' @param Nx Number of samples for the target trait.
#' @param M Number of markers.
#' @param Me Number of independent chromosome segments.
#'
#' @return A numeric value representing the expected R2 for the single trait prediction.
#'
#' @keywords internal
#' @export
#' 
predict_r2_st = function(h2x, Nx, Me, M){
  b = M / (M + Me)
  r_uu = sqrt(h2x)*sqrt((b * h2x ) / (b * h2x  +  Me / (Nx)))
  return(r_uu**2)
}

#' Expected R2 for multiple trait prediction of continuous traits
#'
#' Computes the expected R2 value for the multiple trait prediction of continuous traits.
#'
#' @param h2x Heritability of the target trait.
#' @param h2y Heritability of the correlated trait.
#' @param Nx Number of samples for the target trait.
#' @param Ny Number of samples for the correlated trait.
#' @param rg Genetic correlation between the target and correlated trait.
#' @param M Number of markers.
#' @param Me Number of independent chromosome segments.
#'
#' @return A numeric value representing the expected R2 for the multiple trait prediction.
#' 
#' @keywords internal
#' @export
#' 
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


#' Expected AUC for prediction of a binary trait using information on a correlated continuous trait
#'
#' Computes the expected Area Under the Curve (AUC) for predicting a binary trait using 
#' information from a correlated continuous trait.
#'
#' @param h2x Heritability of the target trait.
#' @param h2y Heritability of the correlated trait.
#' @param Nx Number of samples for the target trait.
#' @param Ny Number of samples for the correlated trait.
#' @param rg Genetic correlation between the target and correlated trait.
#' @param Kx Prevalence of the target trait.
#' @param Px Case-control proportion of the target trait.
#' @param M Number of markers.
#' @param Me Number of independent chromosome segments.
#'
#' @return A numeric value representing the expected AUC.
#' 
#' @keywords internal
#' @export
#' 
predict_auc_mt_continuous = function(h2x, Nx, Kx, Px, h2y, Ny, rg, Me, M){
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


#' Expected AUC for prediction of a binary trait
#'
#' Computes the expected Area Under the Curve (AUC) for predicting a binary trait.
#'
#' @param h2x Heritability of the target trait.
#' @param Nx Number of samples for the target trait.
#' @param Kx Prevalence of the target trait.
#' @param Px Case-control proportion of the target trait.
#' @param M Number of markers.
#' @param Me Number of independent chromosome segments.
#'
#' @return A numeric value representing the expected AUC.
#' 
#' @keywords internal
#' @export
#' 
predict_auc_st = function(h2x, Nx, Kx, Px, Me, M){
  b = M / (M + Me)
  zx = dnorm(-qnorm(Kx))
  r_uu = sqrt((b * h2x * zx^2) / (b * h2x * zx^2 + (Kx * (1-Kx))^2 * Me / (Px * (1-Px) * Nx)))
  i_case  = zx / Kx
  i_con   = -i_case*Kx / (1-Kx)
  auc  = pnorm(((i_case - i_con) * b*h2x * r_uu^2) / sqrt(b*h2x * r_uu^2 * ((1- b*h2x*r_uu^2*i_case*(i_case-(-qnorm(Kx)))) + (1 - b*h2x*r_uu^2*i_con*(i_con - -qnorm(Kx))))))
  return(auc)
}

#' Expected AUC for prediction of a binary trait using information on correlated binary trait
#'
#' Computes the expected Area Under the Curve (AUC) for predicting a binary trait using 
#' information on a correlated binary trait.
#'
#' @param h2x Heritability of the target trait.
#' @param h2y Heritability of the correlated trait.
#' @param Nx Number of samples for the target trait.
#' @param Ny Number of samples for the correlated trait.
#' @param rg Genetic correlation between the target and the correlated trait.
#' @param Kx Prevalence of the target trait.
#' @param Ky Prevalence of the correlated trait.
#' @param Px Case-control proportion of the target trait.
#' @param Py Case-control proportion of the correlated trait.
#' @param M Number of markers.
#' @param Me Number of independent chromosome segments.
#'
#' @return A numeric value representing the expected AUC.
#' 
#' @keywords internal
#' @export
#' 
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


#' Forest plot
#'
#' This function generates a forest plot, which is commonly used to visualize
#' effect sizes and their confidence intervals.
#'
#' @param x A vector of point estimates or effect sizes.
#' @param sd A vector of standard deviations corresponding to the values in `x`.
#' @param cex A numerical value indicating the amount by which plotting text and symbols should be scaled. Default is 1.
#' @param mar A numerical vector of the form `c(bottom, left, top, right)` which gives the number of lines of margin to be specified on the four sides of the plot. Default is `c(5,12,3,1)`.
#' @param mai A numerical vector indicating the margins in inches.
#' @param xlim The x limits (x1, x2) of the plot.
#' @param pos Position of y-axis labels.
#' @param reorder A logical value. If `TRUE`, data points are reordered based on the values in `x`. Default is `TRUE`.
#' @param xaxis A logical value. If `TRUE`, x-axis is drawn. Default is `TRUE`.
#' @param main An overall title for the plot.
#' @param xlab A label for the x-axis. Default is "x".
#'
#' @keywords internal
#' @export
#' 
plotForest <- function(x=NULL,sd=NULL,cex=1, mar=NULL, mai=NULL, xlim=NULL, pos=NULL, reorder=TRUE, xaxis=TRUE, main=NULL, xlab="x") {
  if(is.null(mar)) par(mar=c(5,12,3,1))
  # mai <- c(0.5,5.2,0.3,0.1)
  if(!is.null(mai)) par(mai=mai)
  o <- 1:length(x)
  if(reorder) o <- order(x,decreasing=TRUE)
  x <- x[o]
  if(is.null(sd)) sd <- rep(0,length(x)) 
  sd <- sd[o]
  #if(is.null(xlim)) xlim <- c(min(x)*0.9,min(max(x)*1.1,1)) 
  lower <- x - sd*1.96
  upper <- x + sd*1.96
  labels <- names(x)
  if(is.null(pos)) {
    if(min(lower)<0) pos <- min(lower)*1.01
    if(min(lower)>0) pos <- min(lower)*0.99
  }
  if(min(lower)<0) xlim_lower <- min(lower)*1.01
  if(min(lower)>0) xlim_lower <- min(lower)*0.99
  if(max(upper)<0) xlim_upper <- max(upper)*0.99
  if(max(upper)>0) xlim_upper <- max(upper)*1.01
  if(is.null(xlim)) xlim <- c(xlim_lower,xlim_upper)
  plot(x=x, y=1:length(x), xlim=xlim, pch = 20, xlab=xlab, bty='n', ylab='', yaxt='n', xaxt='n', cex.lab = cex, main=main)
  for(i in 1:length(x)){
    lines(x=c(lower[i], upper[i]), y = rep(i, each=2), lwd=2)
  }
  axis(1, cex.axis=cex)
  if(xaxis) axis(2, at=1:length(x), labels=labels, las=1, lwd=0, pos=pos, outer=TRUE, cex.axis=cex) # pos is x-axis location of labels
  abline(v=pos, lty=2)
}
