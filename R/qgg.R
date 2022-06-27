
.onLoad <- function(libname = find.package("qgg"), pkgname = "qgg") {
  OS <- .Platform$OS.type
  if (OS == "windows") dll <- paste(find.package("qgg"), "/libs/",.Platform$r_arch,"/qgg.dll", sep = "")
  if (!OS == "windows") dll <- paste(find.package("qgg"), "/libs/qgg.so", sep = "")
  if (file.exists(dll)) dyn.load(dll)
}

.onUnload <- function(libname = find.package("qgg"), pkgname = "qgg") {
  OS <- .Platform$OS.type
  if (OS == "windows") dll <- paste(find.package("qgg"), "/libs/",.Platform$r_arch,"/qgg.dll", sep = "")
  if (!OS == "windows") dll <- paste(find.package("qgg"), "/libs/qgg.so", sep = "")
  if (file.exists(dll)) dyn.unload(dll)
}

#' @import Rcpp
#' @importFrom grDevices dev.off gray tiff
#' @importFrom graphics abline barplot boxplot layout par plot points text
#' @importFrom stats anova as.dist binomial cor cor.test glm hclust lm logLik pchisq pf phyper
#' @importFrom stats pt quantile rchisq residuals rgamma rmultinom rnorm runif sd var
#' @importFrom utils write.table askYesNo
#' @importFrom statmod rinvgauss
#' @importFrom data.table fread
#' @importFrom parallel mclapply
#' @importFrom MASS mvrnorm
#' @importFrom MCMCpack riwish
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis hist image legend lines
#' @importFrom stats cov cov2cor dnorm median model.matrix na.omit pnorm qnorm
NULL

#' @useDynLib qgg, .registration = TRUE
NULL
