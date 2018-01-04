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

     computeGRM <- function(W = NULL, miss = 0, pdf = TRUE) {
          
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


