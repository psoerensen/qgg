####################################################################################################################
#
#    Module 1: LMM, REML, BLUP
#    Module 2: Set Test
#    Module 3: LMM, Bayesian
#    Module 4: GREML
#    Module 5: LMMA Test
#
####################################################################################################################
#
# Modules needs to be separated into individual files
# change u -> g 
#

####################################################################################################################
#    Module 1: LMM, REML, BLUP
####################################################################################################################
#' 
#' Genomic Feature Model analyses implemented using Likelihood Methods
#'
#' @description
#' Genomic Feature Best Linear Unbiased Prediction models implemented using REML. 
#' Multiple features and multiple traits models can be fitted.
#' Different genetic models (e.g. additive, dominance, gene by gene and gene by environment interactions) can be specified.
#'
#' @details 
#' The models are implemented using restricted maximum likelihood methods. 
#' Variance components estimated using REML and predictions are based on BLUP. 
#' First and second derivatives of log-likehood can be obtained in addition to asymtotic standard deviation of parameter estimates.
#' Predicted random effects and single marker effects and statistics can be obtained.
#' Covariance decomposition procedures.
#' Cross validation procedures for assessment of prediction accuracy and model selection. 
#' Currently depend on R library regress. 
#' Plug-in avaliable for DMU. Future release may use customized REML procedures.    
#' 
#' @param fm formula with model statement for the fixed factors in the linear mixed model
#' @param weights vector of weights for the residual variance
#' @param W matrix centered and scaled genotypes or other types of molecular data
#' @param sets list of marker sets corresponding to column names in W 
#' @param G list of relationship / correlation matrices
#' @param data data frame containing the phenotypic observations and fixed factors specified in the model statements
#' @param validate matrix or a list with the ids of validation sets corresponding to the rows in data
#' @param mkplots logical indicating whether or not to make plots
#' @return Returns results in a list structure including 
#' \item{f}{matrix of predicted random effects (training set)} 
#' \item{fv}{matrix of predicted random effects validation set} 
#' \item{Vf}{covariance of predicted random effects} 
#' \item{s}{single marker effects} 
#' \item{vs}{variance of single marker effects} 
#' \item{sigmas}{estimated variance components} 
#' \item{pa}{predictive ability based on validation sets} 
#' \item{mspe}{mean squared prediction error based on validation sets} 
#' \item{ypred}{predicted phenotypes validation set} 
#' \item{yobs}{observed phenotypes validation set} 
#' \item{fit}{fit object from regress} 
#' \item{validate}{matrix of validation sets as provided in arguments} 
#' @author Peter Sørensen
#’ @references Edwards, S. M., Sørensen, I. F., Sarup, P., Mackay, T. F., & Sørensen, P. (2016). Genomic prediction for quantitative traits is improved by mapping variants to gene ontology categories in Drosophila melanogaster. Genetics, 203(4), 1871-1883.
#’ @references Rohde, P. D., Demontis, D., Cuyabano, B. C. D., Børglum, A. D., & Sørensen, P. (2016). Covariance Association Test (CVAT) Identifies Genetic Markers Associated with Schizophrenia in Functionally Associated Biological Processes. Genetics, 203(4), 1901-1913.
#’ @references Edwards, S. M., Thomsen, B., Madsen, P., & Sørensen, P. (2015). Partitioning of genomic variance reveals biological pathways associated with udder health and milk production traits in dairy cattle. Genetics Selection Evolution, 47(1), 60.
#’ @references Sarup, P., Jensen, J., Ostersen, T., Henryon, M., & Sørensen, P. (2016). Increased prediction accuracy using a genomic feature model including prior information on quantitative trait locus regions in purebred Danish Duroc pigs. BMC genetics, 17(1), 11.
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(2000000), ncol = 10000)
#'   colnames(W) <- as.character(1:ncol(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#' 
#'
#' # Ex. model 1: simple 
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ mu
#' 
#' # Create framework for lists
#' setsGB <- list(A = colnames(W)) # gblup model
#' setsGF <- list(C1 = colnames(W)[1:1000], C2 = colnames(W)[1001:2000], C3 = colnames(W)[2000:10000])   # gfblup model
#' setsGT <- list(C1 = colnames(W)[1:10], C2 = colnames(W)[1001:1010], C3 = colnames(W)[1:10000]) # true model
#' 
#' # REML analyses and cross validation
#' n <- length(y)
#' fold <- 10
#' nsets <- 50
#' 
#' # Validation sets may be inputted as matrices or lists
#' validate <- replicate(nsets, sample(1:n, as.integer(n / fold))) # matrix input
#'   fitGB <- gfm(fm = fm, W = W, sets = setsGB, data = data, validate = validate)
#'   fitGF <- gfm(fm = fm, W = W, sets = setsGF, data = data, validate = validate)
#'   fitGT <- gfm(fm = fm, W = W, sets = setsGT, data = data, validate = validate)
#' 
#' validClust <- list(C1 = 1:40, C2 = 41:77) # list input
#'   fitGB <- gfm(fm = fm, W = W, sets = setsGB, data = data, validate = validClust)
#'   fitGF <- gfm(fm = fm, W = W, sets = setsGF, data = data, validate = validClust)
#'   fitGT <- gfm(fm = fm, W = W, sets = setsGT, data = data, validate = validClust)
#' 
#' 
#' # Ex. model 2: more complex
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ mu
#' fit <- gfm(fm = fm, W = W, sets = list(colnames(W)), data = data)
#' 
#' # REML analyses and cross validation
#' pc1 <- rnorm(200)
#' pc2 <- rnorm(200)
#' f1 <- factor(rep(1:2, 100))
#' 
#' fmf <- y ~ pc1 + pc2 + f1
#' 
#' fitGB <- gfm(fm = fmf, W = W, data = data.frame(y, pc1, pc2, f1), validate = validate)
#' fitGF <- gfm(fm = fmf, W = W, sets = setsGF, data = data.frame(y, pc1, pc2, f1), validate = validate)
#' fitGT <- gfm(fm = fmf, W = W, sets = setsGT, data = data.frame(y, pc1, pc2, f1), validate = validate)
#' 
#' @export
#' 

gfm <- function(fm = NULL, weights = NULL, W = NULL, sets = NULL, G = NULL, data = NULL, validate = NULL, mkplots = TRUE) {
     
     mf <- model.frame(fm, data = data, na.action = na.pass)
     mf <- eval(mf, parent.frame())
     y <- model.response(mf)
     X <- model.matrix(mf, data = data)
     if (ncol(X) == 2) {if (sum(colnames(X) == c("(Intercept)", "mu")) == 2) X <- matrix(X[, 1], ncol = 1)}
     n <- length(y)
     
     # Validation sets
     if (is.null(validate)) {validate <- matrix(-c(1:n), nrow = n, ncol = 1)}
     if (is.matrix(validate)) nvsets <- ncol(validate)
     if (is.list(validate)) nvsets <- length(validate)
     
     # Check input data
     if (!is.null(sets)) {if (any(!sapply(sets, is.character)))
          stop("Sets not character variables")}
     if (!is.null(W)) {if (any(!is.character(colnames(W))))
          stop("Column names of W is not character variables")}
     if (!is.null(sets)) {if (any(!unlist(sapply(sets, function(x) {x %in% colnames(W)})))) 
          stop("Sets does not match column names of W")}
     if (!is.null(sets)) {if (any(!sapply(sets, is.character))) 
          stop("Sets not character variables")}
     #if (any(!apply(validate, 2, is.integer))) stop("Validation sets not integer variables")
     if (is.matrix(validate)) {if (max(validate) > n)
          stop("Validation sets contains integer values larger than n = number of observations")}
     if (is.list(validate)) {if (max(unlist(validate)) > n)
          stop("Validation sets contains integer values larger than n = number of observations")}
     if (!is.null(weights)) {if (!length(weights) == nrow(W))
          stop("nrow(W) not equal to length(weights)")}
     
     # Compute G matrices
     if (is.null(sets) & is.null(G)) {G[[1]] <- (W %*% t(W)) / ncol(W)}
     nk <- length(G)
     nsets <- NULL
     if (!is.null(sets)) {
          nsets <- length(sets)
          m <- sapply(sets, length)
          for (i in 1:nsets) {G[[i + nk]] <- (W[, sets[[i]]] %*% t(W[, sets[[i]]])) / m[i]}
          names(G)[(nk + 1):(nk + nsets)] <- names(sets)
     }
     if (is.null(nsets)) {nsets <- nk}
     if (!is.null(sets)) {nsets <- nsets + nk}
     if (is.null(names(G))) names(G) <- paste("C", 1:nsets, sep = "")
     G[[nsets + 1]] <- diag(1, n)
     identity <- TRUE
     if (!is.null(weights)) {
          G[[nsets + 1]] <- diag(1 / weights)
          identity <- FALSE
     }
     names(G)[nsets + 1] <- "e"
     
     # Model statements random effects
     vfm <- paste("~", paste(paste("G[[", 1:nsets, "]][t, t]", sep = ""), collapse = "+"))
     if (!is.null(weights)) {vfm <- paste("~", paste(paste("G[[", 1:(nsets + 1), "]][t, t]", sep = ""), collapse = "+"))}
     vfm <- as.formula(vfm)
     
     # Fit mixed model and predict random effect (with or without observation)
     pa <- mspe <- sigmas <- ypred <- yobs <- NULL
     s <- vs <- Vf <- NULL
     for (i in 1:nvsets) {
          if (is.matrix(validate)) v <- validate[, i]
          if (is.list(validate)) v <- validate[[i]]
          t <- (1:n)[-v]
          fit <- regress(fm, vfm, identity = identity, verbose = 1, pos = rep(TRUE, nsets + 1),
                         data = droplevels(data[t, ]))
          V <- matrix(0, ncol = n, nrow = n)
          for (j in 1:(nsets + 1)) {V <- V + G[[j]] * fit$sigma[j]}
          f <- fv <- NULL
          for (j in 1:nsets) {f <- cbind(f, (G[[j]][t, t] * fit$sigma[j]) %*% fit$W %*% (y[t] - fit$fitted))}
          for (j in 1:nsets) {fv <- cbind(fv, (G[[j]][v, t] * fit$sigma[j]) %*% fit$W %*% (y[t] - fit$fitted))}
          if (nsets == 1 & nk == 0) {
               Vf <- NULL
               P <- fit$W %*% fit$Q
               for (j in 1:nsets) { 
                    Vf[[j]] <- (fit$sigma[j] * G[[j]][t, t]) %*% P %*% V[t, t] %*% P %*% (G[[j]][t, t] * fit$sigma[j])
                    bvb <- f2b(W = W[t, sets[[j]]], f = f[, j], Vf = Vf[[j]])
                    if (nsets == 1) {
                         s <- bvb$b
                         vs <- bvb$vb
                    }
                    if (nsets > 1) {
                         s[[j]] <- bvb$b
                         vs[[j]] <- bvb$vb
                    }
               }
          }
          if (length(t) < n) {
               if (any(!fit$X == X[t, ])) stop("!fit$X[t, ] == X[t, ] problem with design matrices for fixed effects")
               #yhat <- fit$fitted[1] + V[v, t] %*% fit$W %*% (y[t] - fit$fitted)
               yhat <- X[v, ] %*% fit$beta + V[v, t] %*% fit$W %*% (y[t] - fit$fitted)
               pa <- c(pa, cor(yhat, y[v])) 
               mspe <- c(mspe, sum((yhat - y[v])**2) / length(v))
               sigmas <- rbind(sigmas, fit$sigma)
               ypred <- c(ypred, yhat)
               yobs <- c(yobs, y[v])
          }
     }
     
     # Make plots
     if (nvsets > 1 & mkplots) {
          colnames(sigmas) <- names(G)
          layout(matrix(1:4, ncol = 2))
          boxplot(pa, main = "Predictive Ability", ylab = "Correlation")
          boxplot(mspe, main = "Prediction Error", ylab = "MSPE")
          boxplot(sigmas, main = "Estimates", ylab = "Variance")
          plot(yobs, ypred, ylab = "Predicted", xlab = "Observed") 
          coef  <- lm(ypred ~ yobs)$coef
          abline(a = coef[1], b = coef[2], lwd = 2, col = 2, lty = 2)
     }
     
     return(list(f = f, fv = fv, Vf = Vf, s = s, vs = vs, sigmas = sigmas, pa = pa, mspe = mspe, ypred = ypred, 
                 yobs = yobs, fit = fit, validate = validate))

}



####################################################################################################################

     f2b <- function(W = NULL, f = NULL, Vf = NULL) {
          
          if (!nrow(W) == length(f)) {stop("nrow(W) == length(f) is false")}
          WW <- W %*% t(W)
          WWi <- (MASS:::ginv)(WW)
          b <- t(W) %*% WWi %*% f
          
          if (!is.null(Vf)) {
               WWWi <- crossprod(W, WWi)
               remove(WWi)
               remove(W)
               WWWiVf <- tcrossprod(WWWi, Vf)
               vb <- rowSums(WWWiVf * WWWi)
               remove(WWWi)
               remove(Vf)
          }
          
          return(list(b = b[, 1], vb = vb))
     
     }

#' @export

     computeG <- function(W = NULL, miss = 0, pdf = TRUE) {
          
          SS <- tcrossprod(W)                              # compute crossproduct, all SNPs
          N <- tcrossprod(!W == miss)                      # compute number of observations, all SNPs
          G <- SS / N
          if (pdf) G <- makepdf(G)
          #return(list(G = G, SS = SS, N = N))
          return(G)
          
     }  
     
     
     makepdf <- function(G = NULL, tol = 0.0001) {
          
          rn <- rownames(G)
          e <- eigen(G)                                    # eigen value decomposition, matrix G
          U <- e$vectors                                   # eigen vectors
          e <- e$values                                    # eigen values
          e[e < tol] <- tol
          D <- diag(e)                                   
          G <- U %*% D %*% t(U)                            # compute pdf
          colnames(G) <- rownames(G) <- rn
          
          return(G)                      
     
     } 
     
     qggginv <- function(G = NULL, tol = NULL) {
          
          rn <- rownames(G)
          e <- eigen(G)                                    # eigen value decomposition, matrix G
          U <- e$vectors                                   # eigen vectors
          e <- e$values                                    # eigen values
          ie <- e
          ie[e > tol] <- 1 / e[e > tol]
          ie[e < tol] <- 0
          D <- diag(ie)                                    # set inverse D to 1 / e
          G <- U %*% D %*% t(U)                            # compute inverse
          ldet <- sum(log(e[e > tol])) 
          colnames(G) <- rownames(G) <- rn
          
          return(list(G = G, ldet = ldet))                 # log determinant 
     
     } 



####################################################################################################################
#    Module 2: SetTests
####################################################################################################################
#' 
#' Genetic marker set tests based on sum statistics
#'
#' @description
#' Set test based on summing the single genetic marker test statistics.
#' The sum test is powerful if the genomic feature harbors many genetic markers having small to moderate effects. 
#'                       
#' @details
#' The singler marker test statistics can be obtained from GBLUP and GFBLUP model fits or from standard GWAS. 
#' The distribution of this test statistic under the null hypothesis (associated markers are picked at random from the total 
#' number of tested genetic markers) is difficult to describe in terms of exact or approximate 
#' distributions, and an empirical distribution is required.
#'                        
#' @param stat vector of single marker statistics (e.g. marker effects, t-stat, p-values)
#' @param sets list of marker sets - names corresponds to rownames in stat
#' @param nperm number of permutations
#' @param W matrix of centered and scaled genotypes (used if method = cvat or score)
#' @param method including sum, cvat, hyperG, score
#' @param threshold used if method = hyperG
#' @return Returns a dataframe including 
#' \item{setT}{marker set test statistics} 
#' \item{nset}{number of markers in the set}
#' \item{p}{p-value for marker set}
#' @author Peter Sørensen
#’ @references Rohde, P. D., Demontis, D., Cuyabano, B. C. D., Børglum, A. D., & Sørensen, P. (2016). Covariance Association Test (CVAT) Identifies Genetic Markers Associated with Schizophrenia in Functionally Associated Biological Processes. Genetics, 203(4), 1901-1913.
#’ @references Rohde, P. D., Edwards S. M., Sarup P., Sørensen, P. (August, 2014). Gene-based Association Approach Identify Genes Across Stress Traits in Fruit Flies. Poster presented at the 10th World Congress of Genetics Applied to Livestock Production (WCGALP), Vancouver, Canada.
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(2000000), ncol = 10000)
#'   colnames(W) <- as.character(1:ncol(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#' 
#' # REML analyses 
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ mu
#' fit <- gfm(fm = fm, W = W, sets = list(colnames(W)), data = data)
#' fit$df <- 10
#' fit$p <- pt(fit$s / sqrt(fit$vs), df = fit$df, lower.tail = FALSE) 
#' 
#' sets <- list(A = as.character(1:100), B = as.character(101:1000), C = as.character(1001:5000), D = as.character(5001:10000))
#'
#' # Set test based on sums 
#' res <- setTest(stat = fit$s**2, sets = sets, method = "sum", nperm = 100)
#' 
#' # Set test based on cvat 
#' res <- setTest(stat = fit$s, W = W, sets = sets, method = "cvat", nperm = 100)
#' 
#' # Set test based on hyperG 
#' res <- setTest(stat = fit$p, sets = sets, method = "hyperG", threshold = 0.05)
#' 
#' @export
#'

setTest <- function(stat = NULL, W = NULL, sets = NULL, nperm = NULL, method = "sum", threshold = 0.05) {
     
  if (method == "sum") setT <- sumTest(stat = stat, sets = sets, nperm = nperm) 
  if (method == "cvat") setT <- cvat(s = stat, W = W, sets = sets, nperm = nperm) 
  if (method == "hyperG") setT <- hgTest(p = stat, sets = sets, threshold = threshold) 
  if (method == "score") setT <- scoreTest(e = e, W = W, sets = sets, nperm = nperm)
     
  return(setT)

}

sumTest <- function(stat = NULL, sets = NULL, nperm = NULL, method = "sum") {
     
  if (method == "mean") setT <- sapply(sets, function(x) {mean(stat[x])})
  if (method == "sum") setT <- sapply(sets, function(x) {sum(stat[x])})
  if (method == "max") setT <- sapply(sets, function(x) {max(stat[x])})
  if (!is.null(nperm)) {
    p <- rep(0, length(sets)) 
    n <- length(stat)
    nset <- sapply(sets, length)
    rws <- 1:n
    names(rws) <- names(stat)
    sets <- lapply(sets, function(x) {rws[x]}) 
    for (i in 1:nperm) {
      rws <- sample(1:n, 1)
      o <- c(rws:n, 1:(rws - 1))
      pstat <- stat[o]
      if (method == "mean") setTP <- sapply(sets, function(x) {mean(pstat[x])})
      if (method == "sum") setTP <- sapply(sets, function(x) {sum(pstat[x])})
      if (method == "max") setTP <- sapply(sets, function(x) {max(pstat[x])})
      p <- p + as.numeric(setT > setTP) 
    }  
    p <- 1 - p / nperm
    setT <- data.frame(setT, nset, p)
  }
     
  return(setT)

}

msetTest <- function(stat = NULL, sets = NULL, nperm = NULL, method = "sum") {
     
     setT <- apply(stat, 2, function(x) {setTest(stat = x, sets = sets, nperm = nperm, method = method)})
     names(setT) <- colnames(stat)
     setT  

} 

####################################################################################################################
#' 
#' Genetic marker set tests based on the covariance test
#'
#' @description
#' Genetic marker set tests based on the covariance statistics for a set of genetic markers.
#' 
#' @details
#' The covariance test statistic is derived from a GBLUP (or GFBLUP) model fit. It is a measure of covariance between the total genomic effect for all markers 
#' and the genomic effect for the genetic markers in the genomic feature. It also relates to the explained sums of
#' squares for the genetic markers. 
#' The distribution of this test statistic under the null hypothesis is difficult to describe in terms of exact or approximate 
#' distributions, and an empirical distribution is required.
#' 
#' @param fit is the fit object obtained from a linear mixed model fit using the greml function
#' @param g vector (or list) of genetic effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param s vector (or list) of single marker effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param W matrix of centered and scaled genotypes (n x m)
#' @param sets list of marker sets corresponding to columns in W
#' @param nperm number of permutations
#' @return Returns a dataframe including 
#' \item{setT}{covariance test statistics} 
#' \item{nset}{number of markers in the set}
#' \item{p}{p-value}
#' @author Peter Sørensen
#’ @references Rohde, P. D., Demontis, D., Cuyabano, B. C. D., Børglum, A. D., & Sørensen, P. (2016). Covariance Association Test (CVAT) Identifies Genetic Markers Associated with Schizophrenia in Functionally Associated Biological Processes. Genetics, 203(4), 1901-1913.
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(20000000), ncol = 10000)
#' 	colnames(W) <- as.character(1:ncol(W))
#' 	rownames(W) <- as.character(1:nrow(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#'
#' # Create model
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ 0 + mu
#' X <- model.matrix(fm, data = data)
#'
#' # Create framework for lists
#' setsGB <- list(A = colnames(W)) # gblup model
#' setsGF <- list(C1 = colnames(W)[1:1000], C2 = colnames(W)[1001:2000], C3 = colnames(W)[2000:10000]) # gfblup model
#' setsGT <- list(C1 = colnames(W)[1:10], C2 = colnames(W)[1001:1010], C3 = colnames(W)[1:10000]) # true model
#'
#' # Compute G
#' G <- computeG(W = W)
#' GB <- lapply(setsGB, function(x) {computeG(W = W[, x])})
#' GF <- lapply(setsGF, function(x) {computeG(W = W[, x])})
#' GT <- lapply(setsGT, function(x) {computeG(W = W[, x])})
#'
#' # REML analyses and multi marker association (set) test
#' fitGB <- greml(y = y, X = X, G = GB, verbose = TRUE)
#'
#' # Use fit object as input
#' cvat(fit = fitGB, W = W, sets = setsGF, nperm = 1000)
#' cvat(fit = fitGB, W = W, sets = setsGT, nperm = 1000)
#'
#' # Use single coefficients as input 
#' s <- crossprod(W / ncol(W), fitGB$Py) * fitGB$theta[1]
#' cvat(s = s, W = W, sets = setsGF, nperm = 1000)
#' cvat(s = s, W = W, sets = setsGT, nperm = 1000)
#' 
#' @export
#'

cvat <- function(fit = NULL, s = NULL, g = NULL, W = NULL, sets = NULL, nperm = 100) {
     if (!is.null(fit)) {s <- crossprod(W/ncol(W),fit$Py)*fit$theta[1]}
     Ws <- t(t(W) * as.vector(s))
     if (is.null(g)) g <- W %*% s   
     cvs <- colSums(as.vector(g) * Ws)
     #setT <- setTest(stat = cvs, sets = sets, nperm = nperm, method = "sum")$p
     #names(setT) <- names(sets)
     setT <- setTest(stat = cvs, sets = sets, nperm = nperm, method = "sum")
     if (!is.null(names(sets))) rownames(setT) <- names(sets)
     return(setT)

}

scoreTest <- function(e = NULL, W = NULL, sets = NULL, nperm = 100) {
     
     we2 <- as.vector((t(W) %*% e)**2)   
     names(we2) <- colnames(W)       
     setT <- setTest(stat = we2, sets = sets, nperm = nperm, method = "sum")$p
     return(setT)

}



####################################################################################################################
#' 
#' Genetic marker set tests based on the hyperG test
#' 
#' @description
#' Genetic marker set tests based on the hyperG test statistics for a set of genetic markers.
#'
#' @details
#' The hyperG marker set test tests a predefined set of markers (i.e. those within a particular genomic feature)
#' for an association with the trait phenotype.
#' Under the null hypothesis (associated markers are picked at random from the total number of tested 
#' genetic markers) it is assumed that the observed count statistic is a realization from a hypergeometric 
#' distribution.
#' This hypothesis can be formulated and tested in a number of ways. Here we consider a test statistic based 
#' on counting the number of genetic markers in the feature that are associated to trait phenotype. 
#' A test based on the count test statistic is likely to have high power to detect association if the genomic 
#' feature harbours genetic markers with large effects. 
#' 
#' @param p vector of single marker p-values (e.g. based on a t-stat)
#' @param sets list of marker sets - names corresponds to rownames in stat
#' @param threshold single marker p-value cut-off
#' @return Returns vector of p values with length equal to the number of sets 
#' @author Peter Sørensen
#’ @references Rohde, P. D., Edwards S. M., Sarup P., Sørensen, P. (August, 2014). Gene-based Association Approach Identify Genes Across Stress Traits in Fruit Flies. Poster presented at the 10th World Congress of Genetics Applied to Livestock Production (WCGALP), Vancouver, Canada.
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(2000000), ncol = 10000)
#'   colnames(W) <- as.character(1:ncol(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#' 
#' # REML analyses 
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ mu
#' fit <- gfm(fm = fm, W = W, sets = list(colnames(W)), data = data)
#' 
#' # hyperG set test 
#' sets <- list(A = as.character(1:100), B = as.character(101:1000), C = as.character(1001:5000), D = as.character(5001:10000))
#' p <- pt(fit$s / fit$vs)
#' res <- hgTest(p = p, sets = sets, threshold = 0.05)
#' 
#' @export
#'

hgTest <- function(p = NULL, sets = NULL, threshold = 0.05) {
     
     N <- length(p)
     Na <- sum(p < threshold)
     Nna <- N - Na
     Nf <- sapply(sets, length)
     Naf <- sapply(sets, function(x) {sum(p[x] < threshold)})
     Nnaf <- Nf - Naf
     Nanf <- Na - Naf
     Nnanf <- Nna - Nnaf
     phyperg <- 1 - phyper(Naf - 1, Nf, N - Nf, Na)
     phyperg

}



####################################################################################################################
#    Module 3: LMM, Bayesian
####################################################################################################################
#' 
#' Genomic Feature Model analyses implemented using Bayesian Methods
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
#’ @references Sørensen, P., de los Campos, G., Morgante, F., Mackay, T. F., & Sorensen, D. (2015). Genetic control of environmental variation of two quantitative traits of Drosophila melanogaster revealed by whole-genome sequencing. Genetics, 201(2), 487-497.
#'
#' @export
#'

bgfm <- function(y = NULL, g = NULL, nsamp = 50, nburn = 10, nsave = 10000, tol = 0.001) {
     
  # nsamp is the number of samples
  y <- as.matrix(y)
  n <- nrow(y)                            # number of observation
  nt <- ncol(y)                           # number of traits
  mu <- colMeans(y)
  samples <- as.integer(seq(nburn + 1, nsamp, length.out = min(nsave, nsamp)))   
  
  psum <- matrix(0, nrow = n, ncol = nt)
  
  nset <- length(g$G)                      # number of sets
  for (i in 1:nt) {                        # genetic values (n*nset)
    g$values[[i]] <- matrix(0, nrow = n, ncol = nset)         
  }
  if (nt > 1) sigmae <- g$sigma[[nset + 1]][1, 1]
  if (nt == 1) sigmae <- g$sigma[[nset + 1]]
  sigma <- g$sigma
  mus <- matrix(NA, nrow = nsamp, ncol = nt)
  
  alphas <- vector(length = nt, mode = "list")
  
  U <- D <- sigmas <- vector(length = nset, mode = "list")
  for (i in 1:nset) {                    
    e <- eigen(g$G[[i]])                   # eigen value decomposition of the matrix G
    ev <- e$values
    U[[i]] <- e$vectors[, ev > tol]        # keep eigen vector if ev > tol
    D[[i]] <- e$values[ev > tol]           # keep eigen value if ev > tol
    sigmas[[i]] <- matrix(NA, nrow = nsamp, ncol = (nt * (nt + 1)) / 2)
    for (j in 1:nt) {
      as <-  matrix(NA, nrow = sum(ev > tol), ncol = length(samples))
      colnames(as) <- as.character(samples)
      alphas[[j]][[i]] <- as
    }    
  }
  
  for (i in 1:nsamp) {                     # loop over the number of samples (nsamp)
       
    for (j in 1:nset) {                    # loop over the number of sets
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
         
      for (t in 1:nt) {g$values[[t]][, j] <- U[[j]] %*% a[, t]}
      df <- nrow(a) + g$df[j]
      
      if (any(i %in% samples)) {
        for (t in 1:nt) {alphas[[t]][[j]][, as.character(i)] <- a[, t]}
      }
      
      # inverse chisquare
      if (nt == 1) {
        #scg <- sum((1 / D[[j]]) * a**2) + g$sigma[[j]] * g$df[j]         
        scg <- sum((1 / D[[j]]) * a**2) + (g$sigma[[j]] * (g$df[j] + 2)) / g$df[j]	# => S = (mode * (df + 2)) / df         
        sigma[j] <- scg / rchisq(n = 1, df = df, ncp = 0)    
        sigmas[[j]][i, ] <- sigma[j]
      }
      
      # inverse wishart
      if (nt > 1) {
           
        if (g$ds[j] == "invWish") {
          S <- t(a * (1 / D[[j]])) %*% a + g$sigma[[j]] * (g$df[j] + nt + 1)		# => S = mode * (df + nt + 1)
          S <- riwish(df, S)
          sigma[[j]] <- S
          sigmas[[j]][i, ] <- S[as.vector(lower.tri(S, diag = TRUE))]
          #print(c(j, g$ds[j]))
        }
           
        if (g$ds[j] == "invChisq") {
          scg <- t(a * (1 / D[[j]])) %*% a + (g$sigma[[j]] * (g$df[j] + 2)) / g$df[j]	# => S = (mode * (df + 2)) / df
             
          for (t in 1:nt) {
            sigma[[j]][t, t] <- scg[t, t] / rchisq(n = 1, df = df, ncp = 0)    
          }
             
          sigmas[[j]][i, ] <- sigma[[j]][as.vector(lower.tri(sigma[[j]], diag = TRUE))]
             
        }
        
      }
         
    }
       
    for (t in 1:nt) {
      yadj <- y[,t] - rowSums(as.matrix(g$values[[t]])) 
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



####################################################################################################################
#    Module 4: GREML analysis
####################################################################################################################
#'
#' Genomic REML analysis
#'
#' @description
#' Genomic restricted maximum likelihood estimation (REML) is an analysis used to estimate genomic and residual variance.
#' Genomic variance is the variance associated with the genomic relationship matrix.
#'
#' @details
#' Linear mixed model (LMM) that models covariance among individuals using realized relationships at genotyped loci.
#' This modeling scheme is achieved through the genomic relationship matrix (G).
#' This matrix can be inputted 'as is' or with a more efficient list structure, Glist, that contains information about G.
#' The model can accomodate fixed effects.
#' Individuals may be subsetted for additional analyses such as cross validation.
#' 
#' @param y vector of phenotypes
#' @param X design matrix of fixed effects
#' @param Glist list of information about G matrix
#' @param G genomic relationship matrix
#' @param ids vector of subsetted individuals to retain for analysis, e.g. cross validation
#' @param theta initial values for reml estimation
#' @param maxit maximum number of iterations of reml analysis
#' @param tol tolerance, i.e. the maximum allowed difference between two consecutive iterations of reml to declare convergence
#' @param bin executable file in fortran
#' @param nthreads number of threads
#' @param wkdir working directory
#' @return Returns a list structure, fit, including
#' \item{llik}{log-likelihood at convergence}
#' \item{theta}{initial values for reml estimation}
#' \item{asd}{asymptotic standard deviation}
#' \item{b}{vector of fixed effect estimates}
#' \item{varb}{vector of variances of fixed effect estimates}
#' \item{u}{vector of random effect estimates}
#' \item{e}{vector of residual effects}
#' \item{Vy}{product of variance-covariance matrix of y at convergence and y}
#' \item{Py}{product of projection matrix of y and y}
#' \item{trPG}{trace of product of projection matrix of y and G}
#' \item{trVG}{trace of product of variance-covariance matrix of y at convergence and G}
#' \item{y}{vector of phenotypes}
#' \item{X}{design matrix of fixed effects}
#' \item{ids}{vector of subsetted individuals retained for analysis}
#' \item{yVy}{product of y, variance-covariance matrix of y at convergence, and y}
#' \item{fnamesG}{filename(s) and locations of of G}
#' \item{wd}{working directory}
#' \item{Glist}{list of information about G matrix}
#' @author Peter Sørensen
#' @references Lee, S. H., & van Der Werf, J. H. (2006). An efficient variance component approach implementing an average information REML suitable for combined LD and linkage mapping with a general complex pedigree. Genetics Selection Evolution, 38(1), 25.
#' @examples
#'
#' # Simulate data
#' W <- matrix(rnorm(20000000), ncol = 10000)
#' 	colnames(W) <- as.character(1:ncol(W))
#' 	rownames(W) <- as.character(1:nrow(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#'
#' # Create model
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ 0 + mu
#' X <- model.matrix(fm, data = data)
#'
#' # Create framework for lists
#' setsGB <- list(A = colnames(W)) # gblup model
#' setsGF <- list(C1 = colnames(W)[1:1000], C2 = colnames(W)[1001:2000], C3 = colnames(W)[2000:10000]) # gfblup model
#' setsGT <- list(C1 = colnames(W)[1:10], C2 = colnames(W)[1001:1010], C3 = colnames(W)[1:10000]) # true model
#'
#' # Compute G
#' G <- computeG(W = W)
#' GB <- lapply(setsGB, function(x) {computeG(W = W[, x])})
#' GF <- lapply(setsGF, function(x) {computeG(W = W[, x])})
#' GT <- lapply(setsGT, function(x) {computeG(W = W[, x])})
#'
#' # REML analyses
#' fitGB <- greml(y = y, X = X, G = GB, verbose = TRUE)
#' fitGF <- greml(y = y, X = X, G = GF, verbose = TRUE)
#' fitGT <- greml(y = y, X = X, G = GT, verbose = TRUE)
#'
#' # REML analyses and cross validation
#' n <- length(y)
#' fold <- 10
#' nsets <- 5
#' 
#' validate <- replicate(nsets, sample(1:n, as.integer(n / fold)))
#' 
#' cvGB <- greml(y = y, X = X, G = GB, validate = validate)
#' cvGF <- greml(y = y, X = X, G = GF, validate = validate)
#' cvGT <- greml(y = y, X = X, G = GT, validate = validate)
#'
#' cvGB
#' cvGF
#' cvGT
#' 
#' boxplot(cbind(cvGB[,1:4],cvGF[,1:4],cvGT[,1:4]))
#' 
#' @export
#'

greml <- function(y = NULL, X = NULL, Glist=NULL, G=NULL, theta=NULL, ids=NULL, validate=NULL, maxit=100, tol=0.00001,bin=NULL,nthreads=1,wkdir=getwd(), verbose=FALSE)
{
  if(is.null(bin)) { 
    if (is.null(validate)) fit <- remlR(y=y, X=X, Glist=Glist, G=G, theta=theta, ids=ids, maxit=maxit, tol=tol, bin=bin, nthreads=nthreads, verbose=verbose, wkdir=wkdir)
    if (!is.null(validate)) fit <- cvreml(y=y, X=X, Glist=Glist, G=G, theta=theta, ids=ids, validate=validate, maxit=maxit, tol=tol, bin=bin, nthreads=nthreads, verbose=verbose, wkdir=wkdir)
  }
  if(!is.null(bin)) { 
    fit <- remlF(y=y, X=X, Glist=Glist, G=G, ids=ids, theta=theta, maxit=maxit, tol=tol, bin=bin, nthreads=nthreads, verbose=verbose, wkdir=wkdir)
  }
  return(fit)  
}  


####################################################################################################################

# REML interface functions for fortran

remlF <- function(y = NULL, X = NULL, Glist = NULL, G = NULL, theta = NULL, ids = NULL, maxit = 100, tol = 0.00001, bin = NULL, nthreads = 1, wkdir = getwd(), verbose = FALSE ) {
#greml <- function(y = NULL, X = NULL, Glist = NULL, G = NULL, ids = NULL, theta = NULL, maxit = 100, tol = 0.00001, bin = NULL, nthreads = 1, wkdir = getwd()) {
    
	write.reml(y = as.numeric(y), X = X, G = G)
	n <- length(y)
	nf <- ncol(X)
	if (!is.null(G)) fnamesG <- paste("G", 1:length(G), sep = "")
	if (!is.null(Glist$fnG)) fnamesG <- Glist$fnG
	nr <- length(fnamesG) + 1
 	if (is.null(ids)) {indxG <- c(n, 1:n)} 
	if (!is.null(ids)) {indxG <- c(Glist$n, match(ids, Glist$idsG))} 
	write.table(indxG, file = "indxg.txt", quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)

	write.table(paste(n, nf, nr, maxit, nthreads), file = "param.txt", quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
	if (is.null(theta)) theta <- rep(sd(y) / nr, nr)
	#if (is.null(theta)) theta <- rep(var(y) / nr, nr)
	write.table(t(theta), file = "param.txt", quote = FALSE, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE)
	write.table(tol, file = "param.txt", quote = FALSE, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE)
	write.table(fnamesG, file = "param.txt", quote = TRUE, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE)

	execute.reml(bin = bin,  nthreads = nthreads)
	fit <- read.reml(wkdir = wkdir)
	fit$y <- y
	fit$X <- X
	fit$ids <- names(y)
	fit$yVy <- sum(y * fit$Vy)
	fit$fnamesG <- fnamesG
	fit$wd <- getwd()
	fit$Glist <- Glist

	clean.reml(wkdir = wkdir)
      
	return(fit)
      
}


write.reml <- function(y = NULL, X = NULL, G = NULL) {
    
	fileout <- file("y", "wb")
	writeBin(y, fileout)
	close(fileout)
      
	filename <- "X"
	fileout <- file(filename, "wb")
	for (i in 1:nrow(X)) {writeBin(X[i, ], fileout)}
	close(fileout)
  
	if (!is.null(G)) {
		for (i in 1:length(G)) {
			fileout <- file(paste("G", i, sep = ""), "wb")
			#writeBin(G[[i]][upper.tri(G[[i]], diag = TRUE)], fileout)
			nr <- nrow(G[[i]])
			for (j in 1:nr) {
				writeBin(G[[i]][j, j:nr], fileout)
			}
			close(fileout)
		}
	}
          
}

execute.reml  <- function (bin = NULL, nthreads = nthreads) {

	HW <- Sys.info()["machine"]
	OS <- .Platform$OS.type
	if (OS == "windows") {
		"my.system" <- function(cmd) {return(system(paste(Sys.getenv("COMSPEC"), "/c", cmd)))}
        
		#my.system(paste("set MKL_NUM_THREADS = ", nthreads))
		test <- my.system(paste(shQuote(bin), " < param.txt > reml.lst", sep = ""))
	}
	if (!OS == "windows") {
		system(paste("cp ", bin, " reml.exe", sep = ""))
		#system(paste("export MKL_NUM_THREADS=", nthreads))
		system("time ./reml.exe < param.txt > reml.lst")
	}
      
} 
   
read.reml <- function (wkdir = NULL) {
    
	llik <- read.table(file = "llik.qgg", header = FALSE, colClasses = "numeric")
	names(llik) <- "logLikelihood" 
	theta <- read.table(file = "theta.qgg", header = FALSE, colClasses = "numeric")
	colnames(theta) <- "Estimate"
	rownames(theta) <- 1:nrow(theta)
	asd <- read.table(file = "thetaASD.qgg", header = FALSE, colClasses = "numeric") 
	colnames(asd) <- rownames(asd) <- 1:ncol(asd)
	b <- read.table(file = "beta.qgg", header = FALSE, colClasses = "numeric")    
	colnames(b) <- "Estimate"
	rownames(b) <- 1:nrow(b)
	varb <- read.table(file = "betaASD.qgg", header = FALSE, colClasses = "numeric")    
	colnames(varb) <- rownames(varb) <- 1:nrow(b)
	u <- read.table(file = "uhat.qgg", header = FALSE, colClasses = "numeric")
	colnames(u) <- 1:(nrow(theta) - 1)
	e <- read.table(file = "residuals.qgg", header = FALSE, colClasses = "numeric")    
	colnames(e) <- "residuals"
	rownames(e) <- rownames(u) <- 1:nrow(u)
	Vy <- read.table(file = "Vy.qgg", header = FALSE, colClasses = "numeric")    
	rownames(Vy) <- 1:nrow(u)
	Py <- read.table(file = "Py.qgg", header = FALSE, colClasses = "numeric")    
	rownames(Py) <- 1:nrow(u)
	trPG <- as.vector(unlist(read.table(file = "trPG.qgg", header = FALSE, colClasses = "numeric")[, 1]))    
	names(trPG) <- 1:nrow(theta)
	trVG <- as.vector(unlist(read.table(file = "trVG.qgg", header = FALSE, colClasses = "numeric")[, 1]))    
	names(trVG) <- 1:nrow(theta)
	fit <- list(llik = llik, theta = theta, asd = asd, b = b, varb = varb, g = u, e = e, Vy = Vy, Py = Py, trPG = trPG, trVG = trVG)
	fit <- lapply(fit, as.matrix)
      
	return(fit)
      
}

clean.reml <- function(wkdir = NULL) {
    
	fnames <- c("llik.qgg", "theta.qgg", "thetaASD.qgg", "beta.qgg", "betaASD.qgg", 
				"uhat.qgg", "residuals.qgg", "Vy.qgg", "Py.qgg", "trPG.qgg", "trVG.qgg") 
	file.remove(fnames)
      
}

####################################################################################################################

# REML R functions 

remlR <- function(y=NULL, X=NULL, Glist=NULL, G=NULL, theta=NULL, ids=NULL, maxit=100, tol=0.00001, bin=NULL,nthreads=1,wkdir=getwd(), verbose=FALSE )
  
  #reml <- function( y=NULL, X=NULL, Glist=NULL, G=NULL,theta=NULL, ids=NULL, maxit=100, verbose=FALSE)
{
  
  np <- length(G) + 1
  if (is.null(theta)) theta <- rep(sd(y)/np**2,np)
  n <- length(y)
  ai <- matrix(0, ncol=np, nrow=np)
  s <- matrix(0, ncol=1, nrow=np)
  tol <- 0.00001
  delta <- 100
  it <- 0
  G[[np]] <- diag(1,length(y))
  
  while ( max(delta)>tol ) {
    V <- matrix(0,n,n)
    u <- Pu <- matrix(0,nrow=n,ncol=np)
    it <- it + 1
    for ( i in 1:np) { V <- V + G[[i]]*theta[i] }
    Vi <- chol2inv(chol(V))
    remove(V)
    XViXi <- chol2inv(chol(crossprod(X,crossprod(Vi,X) ) ) )
    ViX <- crossprod(Vi,X) 
    ViXXViXi <- tcrossprod(ViX,XViXi)
    remove(XViXi)
    P <- Vi - tcrossprod(ViXXViXi,ViX)
    remove(Vi)
    Py <- crossprod(P,y)
    for ( i in 1:np) {
      u[,i] <- crossprod(G[[i]],Py)
      Pu[,i] <- crossprod(P,u[,i])
    }
    for ( i in 1:np) {
      for ( j in i:np) {
        ai[i,j] <- 0.5*sum(u[,i]*Pu[,j])
        ai[j,i] <- ai[i,j]
      }
      if (i<np) s[i,1] <- -0.5*(sum(G[[i]]*P)-sum(u[,i]*Py))
      if (i==np) s[i,1] <- -0.5*(sum(diag(P))-sum(Py*Py))
    }
    theta.cov <- solve(ai)
    theta0 <- theta + solve(ai)%*%s
    theta0[theta0<0] <- 0.000000001
    delta <- abs(theta - theta0)
    theta <- theta0
    if (verbose) print(paste(c("Iteration:",it,"Theta:",round(theta,5)), sep=""))
    if (it==maxit) break
  }
  V <- matrix(0,n,n)
  for ( i in 1:np) { V <- V + G[[i]]*theta[i] }
  chlV <- chol(V)
  remove(V)
  ldV <- log(sum(diag(chlV)))
  Vi <- chol2inv(chlV)
  remove(chlV)
  chlXVX <- chol(crossprod(X,crossprod(Vi,X) ))
  ldXVX <- log(sum(diag(chlXVX)))
  XViXi <- chol2inv(chlXVX)
  ViX <- crossprod(Vi,X)
  ViXXViXi <- tcrossprod(ViX,XViXi)
  b <- crossprod(ViXXViXi,y)
  vb <- XViXi
  P <- Vi - tcrossprod(ViXXViXi,ViX)
  trPG <- trVG <- rep(0,length(theta))
  for (i in 1:np) {
    trVG[i] <- sum(Vi*G[[i]])
    trPG[i] <- sum(P*G[[i]])
  } 
  Vy <- crossprod(Vi,y)
  remove(Vi)
  Py <- crossprod(P,y)
  yPy <- sum(y*Py)
  yVy <- sum(y*Vy)
  llik <- -0.5*(ldV+ldXVX+yPy)
  
  u <- NULL
  for (i in 1:(length(theta)-1)) {
    u <- cbind(u, crossprod(G[[i]]*theta[i],Py) )       
  }
  fitted <- X%*%b
  predicted <- rowSums(u)+fitted
  
  return(list( y=y, X=X, b=b, vb=vb, g=u, fitted=fitted, predicted=predicted, Py=Py, Vy=Vy, theta=theta, asd=theta.cov, llik=llik, niter=it,trPG=trPG, trVG=trVG,ids=names(y),yVy=yVy   ))
}


cvreml <- function(y=NULL, X=NULL, Glist=NULL, G=NULL, theta=NULL, ids=NULL, validate=NULL, maxit=100, tol=0.00001,bin=NULL,nthreads=1,wkdir=getwd(), verbose=FALSE)
{
  n <- length(y)     
  theta <- pa <- mspe <- yobs <- ypred <- r2 <- llik <- slope <- intercept <- NULL
  for (i in 1:ncol(validate)) {
    v <- validate[,i]
    t <- (1:n)[-v]
    fit <- remlR( y=y[t], X=X[t,], G=lapply(G,function(x){x[t,t]}), verbose=verbose)
    theta <- rbind(theta, as.vector(fit$theta))
    np <- length(fit$theta)
    yhat <- X[v, ] %*% fit$b
    for (j in 1:(np-1)) {
      yhat <- yhat + G[[j]][v,t]%*%fit$Py*fit$theta[j]
    }
    pa <- c(pa, cor(yhat, y[v]))
    mspe <- c(mspe, sum((yhat - y[v])^2)/length(v))
    intercept <- c(intercept,lm( y[v] ~ yhat )$coef[1])
    slope <- c(slope,lm( y[v] ~ yhat )$coef[2])
    r2 <- c(r2,summary(lm( y[v] ~ yhat ))$r.squared)
    llik <- c(llik,fit$llik)
    yobs <- c(yobs, y[v])
    ypred <- c(ypred, yhat)
    if (i > 1) {
      colnames(theta) <- c(names(G),"E")
      layout(matrix(1:4, ncol = 2))
      boxplot(pa, main = "Predictive Ability", ylab = "Correlation")
      boxplot(mspe, main = "Prediction Error", ylab = "MSPE")
      boxplot(theta, main = "Estimates", ylab = "Variance")
      plot(yobs, ypred, ylab = "Predicted", xlab = "Observed")
      coef <- lm(yobs ~ ypred)$coef
      abline(a = coef[1], b = coef[2], lwd = 2, col = 2, lty = 2)
    }
    
  }    
  res <- data.frame(pa,r2,intercept,slope,mspe,llik,theta)
  #return(list(pa=pa,mspe=mspe,theta=theta,ypred=ypred,yobs=yobs))
}





####################################################################################################################
#    Module 5: LMM marker association test
####################################################################################################################
#'
#' Linear mixed model marker association test
#'
#' @description
#' Marker test using linear mixed model (LMM) to test for an association of single markers with a phenotype.
#'
#' @details
#' Linear mixed model single marker association (LMMA) statistics are based on either exact or approximate methods.
#' Exact methods estimate variance components and effects of single markers jointly.
#' Approximate methods estimate single marker effects conditionally.
#'
#' @param fit list of information about linear model fit (output from greml)
#' @param W matrix of centered and scaled genotypes (n x m)
#' @param m is the total number of markers in W 
#' @param statistic is the single marker test statistic used. Default is the "mastor", alternatives include "gblup", "bolt-lmm" and "grammar-gamma"
#' @return Returns a dataframe including 
#' \item{coef}{single marker coefficients} 
#' \item{se}{standard error of coefficients}
#' \item{stat}{single marker test statistic}
#' \item{p}{p-value}
#' @author Peter Sørensen
#' @references Chen, W. M., & Abecasis, G. R. (2007). Family-based association tests for genomewide association scans. The American Journal of Human Genetics, 81(5), 913-926.
#' @references Loh, P. R., Tucker, G., Bulik-Sullivan, B. K., Vilhjalmsson, B. J., Finucane, H. K., Salem, R. M., ... & Patterson, N. (2015). Efficient Bayesian mixed-model analysis increases association power in large cohorts. Nature genetics, 47(3), 284-290.
#' @references Kang, H. M., Sul, J. H., Zaitlen, N. A., Kong, S. Y., Freimer, N. B., Sabatti, C., & Eskin, E. (2010). Variance component model to account for sample structure in genome-wide association studies. Nature genetics, 42(4), 348-354.
#' @references Lippert, C., Listgarten, J., Liu, Y., Kadie, C. M., Davidson, R. I., & Heckerman, D. (2011). FaST linear mixed models for genome-wide association studies. Nature methods, 8(10), 833-835.
#' @references Listgarten, J., Lippert, C., Kadie, C. M., Davidson, R. I., Eskin, E., & Heckerman, D. (2012). Improved linear mixed models for genome-wide association studies. Nature methods, 9(6), 525-526.
#' @references Listgarten, J., Lippert, C., & Heckerman, D. (2013). FaST-LMM-Select for addressing confounding from spatial structure and rare variants. Nature Genetics, 45(5), 470-471.
#' @references Lippert, C., Quon, G., Kang, E. Y., Kadie, C. M., Listgarten, J., & Heckerman, D. (2013). The benefits of selecting phenotype-specific variants for applications of mixed models in genomics. Scientific reports, 3.
#' @references Zhou, X., & Stephens, M. (2012). Genome-wide efficient mixed-model analysis for association studies. Nature genetics, 44(7), 821-824.
#' @references Svishcheva, G. R., Axenovich, T. I., Belonogova, N. M., van Duijn, C. M., & Aulchenko, Y. S. (2012). Rapid variance components-based method for whole-genome association analysis. Nature genetics, 44(10), 1166-1170.
#' @references Yang, J., Zaitlen, N. A., Goddard, M. E., Visscher, P. M., & Price, A. L. (2014). Advantages and pitfalls in the application of mixed-model association methods. Nature genetics, 46(2), 100-106.
#' @references Bulik-Sullivan, B. K., Loh, P. R., Finucane, H. K., Ripke, S., Yang, J., Patterson, N., ... & Schizophrenia Working Group of the Psychiatric Genomics Consortium. (2015). LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nature genetics, 47(3), 291-295.
#' @examples
#'
#' # Simulate data
#' W <- matrix(rnorm(20000000), ncol = 10000)
#' 	colnames(W) <- as.character(1:ncol(W))
#' 	rownames(W) <- as.character(1:nrow(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))
#'
#' # Create model
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ 0 + mu
#' X <- model.matrix(fm, data = data)
#'
#' # Create framework for lists
#' setsGB <- list(A = colnames(W)) # gblup model
#' setsGF <- list(C1 = colnames(W)[1:1000], C2 = colnames(W)[1001:2000], C3 = colnames(W)[2000:10000]) # gfblup model
#' setsGT <- list(C1 = colnames(W)[1:10], C2 = colnames(W)[1001:1010], C3 = colnames(W)[1:10000]) # true model
#'
#' # Compute G
#' G <- computeG(W = W)
#' GB <- lapply(setsGB, function(x) {computeG(W = W[, x])})
#' GF <- lapply(setsGF, function(x) {computeG(W = W[, x])})
#' GT <- lapply(setsGT, function(x) {computeG(W = W[, x])})
#'
#' # REML analyses and single marker association test
#' fitGB <- greml(y = y, X = X, G = GB, verbose = TRUE)
#' maGB <- lmma(fit = fitGB, W = W)
#'
#' @export
#'

lmma <- function( fit=NULL, W=NULL, m=NULL, statistic="mastor") {
  
  if (is.null(m)) m <- ncol(W)
  n <- nrow(W)
  # get linear model fit results
  Sg <- fit$theta[1]
  Py <- fit$Py
  Vy <- fit$Vy
  yVy <- fit$yVy
  trPG <- fit$trPG[1]
  trVG <- fit$trVG[1]
  
  # compute single marker, coefficients, test statistics and p-values
  nW <- colSums(!W==0)
  WPy <- crossprod(W,Py)
  WVy <- crossprod(W,Vy)
  ww <- colSums(W**2)
  
  if (statistic=="mastor") {
    coef <- (WPy/(trPG/m))/m
    se <- (1/(trPG/m))/m
    stat <- WPy**2/sum(Py**2)
    p <- pchisq(stat,df=1,ncp=0, lower.tail=FALSE)
  }
  
  mma <- data.frame(coef=coef,se=se,stat=stat,p=p)
  rownames(mma) <- colnames(W)
  
  return(mma)
  
}

#' @export

plotma <- function(ma=NULL,chr=NULL,rsids=NULL,thresh=5) {
  
  mlogObs <- -log10(ma$p)
  m <- length(mlogObs)
  rwsSNP <- rsids
  if(!is.null(thresh)) { rwsSNP <- mlogObs > thresh }  
  if(is.null(chr)) { chr <- rep(1,m) }  
  
  layout(matrix(1:2,ncol=1))
  o <- order(mlogObs,decreasing=TRUE)
  mlogExp <- -log10( (1:m)/m ) 
  plot(y=mlogObs[o],x=mlogExp, col=2, pch="+",
       frame.plot=FALSE, main="", xlab="Expected -log10(p)", ylab="Observed -log10(p)")
  abline(a=0,b=1)
  
  colSNP <- "red"
  
  is.even <- function(x) x %% 2 == 0
  colCHR <- rep(gray(0.3),m)
  colCHR[is.even(chr)] <- gray(0.9)
  
  plot(y=mlogObs,x=1:m,ylab="Observed -log10(p)",xlab="Position",col=colCHR,   
       pch=".",frame.plot=FALSE)
  points(y=mlogObs[rwsSNP],x=(1:m)[rwsSNP],col=colSNP)  
  abline(h=thresh,col=2,lty=2)
}   


####################################################################################################################
#    Module 5: GSOLVE 
####################################################################################################################

#' @export

gsolve <- function( y=NULL, X=NULL, W=NULL, sets=NULL, msets=100, lambda=NULL, validate=NULL, weights=FALSE, maxit=500, tol=0.0000001) { 
     if(is.null(validate)) fit <- gsqr(y=y, W=W, X=X, sets=sets,msets=msets,lambda=lambda,weights=weights, maxit=maxit, tol=tol)
     if(!is.null(validate)) { 
          n <- length(y)     
          pa <- mspe <- intercept <- slope <- r2 <- NULL
          for ( k in 1:ncol(validate)) {
               v <- validate[, k]
               t <- (1:n)[-v]
               fit <- gsqr(y=y[t], X=as.matrix(X[t,]), W=W[t,], sets=sets,msets=msets,lambda=lambda,weights=weights, maxit=maxit, tol=tol)
               yv <- y[v]
               yvhat <- W[v,]%*%fit$s
               if(!is.null(X)) yvhat <- yvhat + X[v,]%*%fit$b
               r2 <- c(r2, summary(lm(yv ~ yvhat))$r.squared)
               pa <- c(pa, cor(yvhat, yv))
               mspe <- c(mspe, sum((yvhat - yv)^2)/length(yv))
               intercept <- c(intercept, lm(yv ~ yvhat )$coef[1])
               slope <- c(slope, lm(yv ~ yvhat)$coef[2])
          }
          res <- data.frame(pa, r2, intercept, slope, mspe)
          fit <- res
     }   
     return(fit)         
}

gsru <- function( y=NULL, X=NULL, W=NULL, sets=NULL, lambda=NULL, weights=FALSE, maxit=500, tol=0.0000001) { 
     n <- length(y)                        # number of observations
     m <- ncol(W)                          # number of markers
     dww <- rep(0,m)                       # initialize diagonal elements of the W'W matrix
     for( i in 1:m) { dww[i] <- sum(W[,i]**2) }                 
     b <- bold <- bnew <- NULL
     if (!is.null(X)) {
          b <- (solve(t(X)%*%X)%*%t(X))%*%y     # initialize b
          bold <- rep(0,ncol(X))              # initialize b
     }
     if(length(lambda)==1) { lambda <- rep(lambda,m)}
     e <- y
     if (!is.null(X)) e <- y-X%*%b                                # initialize e
     s <- (crossprod(W,e)/dww)/m      # initialize s
     sold <- rep(0,m)                      # initialize s
     if (weights) {
     p<- apply(W,2,function(x){cor.test(x,e)$p.value})
     logP <- -log10(p)
     lambda <- 1/logP
     }
     if(is.null(sets)) { sets <- as.list(1:m)} 
     nsets <- length(sets)
     nit <- 0
     delta <- 1
     while ( delta>tol ) {
          nit <- nit + 1
          for( i in 1:nsets) {
               rws <- sets[[i]] 
               lhs <- dww[rws] + lambda[rws]          # form lhs
               rhs <- crossprod(W[,rws],e) + dww[rws]*s[rws]  # form rhs with y corrected by other effects
               snew <- rhs/lhs
               e  <- e - tcrossprod(W[,rws],matrix((snew-s[rws]),nrow=1))          # update e with current estimate of b
               s[rws] <- snew                         # update estimates
          }
          gc()
          #if (!is.null(X)) {
          #  bnew <- solve(t(X)%*%X)%*%t(X)%*%e 
          #  e  <- e - X%*%(bnew-bold)            
          #}
          delta <- sum((s-sold)**2)
          #if (!is.null(X)) delta <- sum((s-sold)**2) + sum((b-bold)**2)
          delta <- delta/sqrt(m)
          sold <- s
          bold <- bnew 
          if (nit==maxit) break
          print(paste("Iteration",nit,"delta",delta))
     }
     ghat <- W%*%s
     if (!is.null(X)) yhat <- ghat + X%*%b
     e <- y - yhat
     return(list(s=s,b=b,nit=nit,delta=delta, e=e, yhat=yhat, g=ghat))
}

qrSets <- function( W=NULL, sets=NULL, msets=100, return.level="Q") {
     m <- ncol(W)
     if(is.null(sets)) sets <- split(1:m, ceiling(seq_along(1:m)/msets))
     qrR <- list() 
     for ( i in 1:length(sets) ) {
          qrW <- qr(W[,sets[[i]]])
          W[,sets[[i]]] <- qr.Q(qrW)
          qrR[[i]] <- qr.R(qrW) 
          gc()
     }
     QRlist <- W
     if(return.level=="QR") QRlist <- list(Q=W,R=qrR,sets=sets)
     return(QRlist)
}  

#' @export

plotGS <- function( fit=NULL, s=NULL, sets=NULL ) {
     if(is.null(s)) s <- fit$s
     m <- length(s) 
     plot(y=s,x=1:m,ylab="Coefficients",xlab="Position",col=1,   
          pch=".",frame.plot=FALSE)
     points(y=s[sets],x=(1:m)[sets],col=2)  
}  

gsqr <- function( y=NULL, X=NULL, W=NULL, sets=NULL, msets=100, lambda=NULL, weights=FALSE, maxit=500, tol=0.0000001) { 
     QRlist <- qrSets(W=W,msets=msets,return.level="QR")
     fit <- gsru(y=y, X=X, W=QRlist$Q, sets=QRlist$sets, lambda=lambda, weights=weights) 
     nsets <- length(QRlist$sets)
     for ( i in 1:nsets) {
          rws <- QRlist$sets[[i]]
          fit$s[rws] <- solve(QRlist$R[[i]])%*%fit$s[rws,1]
     }
     return(fit)
}  

