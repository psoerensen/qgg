####################################################################################################################
#    Module 1: LMM, REML, BLUP
####################################################################################################################


#' @export
#' 

gfm <- function(fm = NULL, weights = NULL, W = NULL, sets = NULL, G = NULL, data = NULL, validate = NULL, mkplots = FALSE) {
     
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


