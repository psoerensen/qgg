####################################################################################################################
#    Module 4: GREML analysis
####################################################################################################################
#'
#' Genomic REML analysis
#'
#' @description
#' The greml function is used for estimation of genomic parameters (co-variance, heritability and correlation)
#' for linear mixed models using restricted maximum likelihood estimation (REML) and genomic prediction using
#' best linear unbiased prediction (BLUP).
#'
#' The linear mixed model can account for multiple genetic factors (fixed and random genetic marker effects),
#' adjust for complex family relationships or population stratification, and adjust for other non-genetic factors
#' including lifestyle characteristics. Different genetic architectures (infinitesimal, few large and many
#' small effects) is accounted for by modeling genetic markers in different sets as fixed or random effects
#' and by specifying individual genetic marker weights. Different genetic models (e.g. additive and non-additive)
#' can be specified by providing additive and non-additive genomic relationship matrices (GRMs) (constructed using grm).
#' The GRMs can be accessed from the R environment or from binary files stored on disk facilitating analyses of
#' large-scale genetic data.
#'
#' The output contains estimates of variance components, fixed and random effects, first and second derivatives of
#' log-likelihood, and the asymptotic standard deviation of parameter estimates.
#'
#' Assessment of predictive accuracy (including correlation and R2, and AUC for binary phenotypes) can be obtained
#' by providing greml with a dataframe or list containing sample IDs used in the validation, see examples for details.
#'
#' Genomic parameters can also be estimated with DMU (http://www.dmu.agrsci.dk/DMU/) if interface =”DMU”.
#' This option requires DMU to be installed locally, and the path to the DMU binary files has to be specified
#' (see examples below for details).

#' @param y vector or matrix of phenotypes
#' @param X design matrix for factors modeled as fixed effects
#' @param GRM list of one or more genomic relationship matrices
#' @param GRMlist list providing information about GRM matrix stored in binary files on disk
#' @param theta vector of initial values of co-variance for REML estimation
#' @param ids vector of individuals used in the analysis
#' @param validate dataframe or list of individuals used in cross-validation (one column/row for each validation set)
#' @param maxit maximum number of iterations used in REML analysis
#' @param tol tolerance, i.e. convergence criteria used in REML
#' @param ncores number of cores used for the analysis
#' @param fm formula with model statement for the linear mixed model
#' @param data data frame containing the phenotypic observations and fixed factors specified in the model statements
#' @param interface used for specifying whether to use R or Fortran implementations of REML
#' @param wkdir is the working directory used for REML
#' @param makeplots logical if TRUE makes some plots or parameter estimates and prediction accuracy during cross validation
#' @param verbose logical if TRUE print more details during optimization
#' @param bin directory for fortran binaries (e.g. DMU binaries dmu1 and dmuai)

#'
#' @return Returns a list structure including
#' \item{llik}{log-likelihood at convergence}
#' \item{theta}{covariance estimates from REML}
#' \item{asd}{asymptotic standard deviation}
#' \item{b}{vector of fixed effect estimates}
#' \item{varb}{vector of variances of fixed effect estimates}
#' \item{g}{vector or matrix of random effect estimates}
#' \item{e}{vector or matrix of residual effects}
#' \item{accuracy}{matrix of prediction accuracies (only returned if validate is provided)}


#' @author Peter Soerensen


#' @references Lee, S. H., & van Der Werf, J. H. (2006). An efficient variance component approach implementing an average information REML suitable for combined LD and linkage mapping with a general complex pedigree. Genetics Selection Evolution, 38(1), 25.

#' @examples
#'
#' \donttest{
#'
#' # Simulate data
#' W <- matrix(rnorm(1000000), ncol = 1000)
#' 	colnames(W) <- as.character(1:ncol(W))
#' 	rownames(W) <- as.character(1:nrow(W))
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 501:510]) + rnorm(nrow(W))
#'
#' # Create model
#' data <- data.frame(y = y, mu = 1)
#' fm <- y ~ 0 + mu
#' X <- model.matrix(fm, data = data)
#'
#' # Compute GRM
#' GRM <- grm(W = W)
#'
#' # REML analyses
#' fitG <- greml(y = y, X = X, GRM = list(GRM))
#'
#'
#' # REML analyses and cross validation
#'
#' # Create marker sets
#' setsGB <- list(A = colnames(W)) # gblup model
#' setsGF <- list(C1 = colnames(W)[1:500], C2 = colnames(W)[501:1000]) # gfblup model
#' setsGT <- list(C1 = colnames(W)[1:10], C2 = colnames(W)[501:510]) # true model
#'
#' GB <- lapply(setsGB, function(x) {grm(W = W[, x])})
#' GF <- lapply(setsGF, function(x) {grm(W = W[, x])})
#' GT <- lapply(setsGT, function(x) {grm(W = W[, x])})
#'
#' n <- length(y)
#' fold <- 10
#' nvalid <- 5
#'
#' validate <- replicate(nvalid, sample(1:n, as.integer(n / fold)))
#' cvGB <- greml(y = y, X = X, GRM = GB, validate = validate)
#' cvGF <- greml(y = y, X = X, GRM = GF, validate = validate)
#' cvGT <- greml(y = y, X = X, GRM = GT, validate = validate)
#'
#' cvGB$accuracy
#' cvGF$accuracy
#' cvGT$accuracy
#'
#' }

#'
#' @export
#'

greml <- function(y = NULL, X = NULL, GRMlist = NULL, GRM = NULL, theta = NULL,
                  ids = NULL, validate = NULL, maxit = 100, tol = 0.00001, bin = NULL,
                  ncores = 1, wkdir = getwd(), verbose = FALSE, makeplots = FALSE,
                  interface = "R", fm = NULL, data = NULL) {
  if (!is.null(GRM)) {
    if (is.null(validate)) {
      fit <- remlr(
        y = y, X = X, GRMlist = GRMlist, G = GRM, theta = theta, ids = ids,
        maxit = maxit, tol = tol, bin = bin, ncores = ncores, verbose = verbose,
        wkdir = wkdir
      )
    }
    if (!is.null(validate)) {
      fit <- cvreml(
        y = y, X = X, GRMlist = GRMlist, G = GRM, theta = theta, ids = ids,
        validate = validate, maxit = maxit, tol = tol, bin = bin,
        ncores = ncores, verbose = verbose, wkdir = wkdir,
        makeplots = makeplots
      )
    }
  }
  if (!is.null(GRMlist)) {
    fit <- remlf(
      y = y, X = X, GRMlist = GRMlist, G = GRM, theta = theta, ids = ids, maxit = maxit,
      tol = tol, ncores = ncores, verbose = verbose
    )
  }
  return(fit)
}





####################################################################################################################
# REML R functions



remlr <- function(y = NULL, X = NULL, GRMlist = NULL, G = NULL, theta = NULL, ids = NULL, maxit = 100,
                  tol = 0.00001, bin = NULL, ncores = 1, wkdir = getwd(), verbose = FALSE) {
  np <- length(G) + 1
  if (is.null(theta)) theta <- rep(sd(y) / np**2, np)
  n <- length(y)
  if(is.null(X)) X <- model.matrix(y~1)
  ai <- matrix(0, ncol = np, nrow = np)
  s <- matrix(0, ncol = 1, nrow = np)
  tol <- 0.00001
  delta <- 100
  it <- 0
  G[[np]] <- diag(1, length(y))

  while (max(delta) > tol) {
    V <- matrix(0, n, n)
    u <- Pu <- matrix(0, nrow = n, ncol = np)
    it <- it + 1
    for (i in 1:np) {
      V <- V + G[[i]] * theta[i]
    }
    Vi <- chol2inv(chol(V))
    remove(V)
    XViXi <- chol2inv(chol(crossprod(X, crossprod(Vi, X))))
    ViX <- crossprod(Vi, X)
    ViXXViXi <- tcrossprod(ViX, XViXi)
    remove(XViXi)
    P <- Vi - tcrossprod(ViXXViXi, ViX)
    remove(Vi)
    Py <- crossprod(P, y)
    for (i in 1:np) {
      u[, i] <- crossprod(G[[i]], Py)
      Pu[, i] <- crossprod(P, u[, i])
    }
    for (i in 1:np) {
      for (j in i:np) {
        ai[i, j] <- 0.5 * sum(u[, i] * Pu[, j])
        ai[j, i] <- ai[i, j]
      }
      if (i < np) s[i, 1] <- -0.5 * (sum(G[[i]] * P) - sum(u[, i] * Py))
      if (i == np) s[i, 1] <- -0.5 * (sum(diag(P)) - sum(Py * Py))
    }
    theta.cov <- solve(ai)
    theta0 <- theta + solve(ai) %*% s
    theta0[theta0 < 0] <- 0.000000001
    delta <- abs(theta - theta0)
    theta <- theta0
    output <- c(1:10, seq(11, maxit, 5))
    #if (verbose & it < 10 | it %in% output) print(paste(c("Iteration:", it, "Theta:", round(theta, 2)), sep = ""))
    if (verbose) print(paste(c("Iteration:", it, "Theta:", round(theta, 2)), sep = ""))
    if (it == maxit) {
      warning("Maximum number of iterations reached. Try increasing maxit.")
      break
    }
  }
  if (it < maxit & verbose) print(paste(c("Converged at Iteration:", it, "Theta:", round(theta, 2)), sep = ""))
  V <- matrix(0, n, n)
  for (i in 1:np) {
    V <- V + G[[i]] * theta[i]
  }
  chlV <- chol(V)
  remove(V)
  ldV <- log(sum(diag(chlV)))
  Vi <- chol2inv(chlV)
  remove(chlV)
  chlXVX <- chol(crossprod(X, crossprod(Vi, X)))
  ldXVX <- log(sum(diag(chlXVX)))
  XViXi <- chol2inv(chlXVX)
  ViX <- crossprod(Vi, X)
  ViXXViXi <- tcrossprod(ViX, XViXi)
  b <- crossprod(ViXXViXi, y)
  vb <- XViXi
  P <- Vi - tcrossprod(ViXXViXi, ViX)
  trPG <- trVG <- rep(0, length(theta))
  for (i in 1:np) {
    trVG[i] <- sum(Vi * G[[i]])
    trPG[i] <- sum(P * G[[i]])
  }
  Vy <- crossprod(Vi, y)
  remove(Vi)
  Py <- crossprod(P, y)
  yPy <- sum(y * Py)
  yVy <- sum(y * Vy)
  llik <- -0.5 * (ldV + ldXVX + yPy)

  u <- NULL
  for (i in 1:(length(theta) - 1)) {
    u <- cbind(u, crossprod(G[[i]] * theta[i], Py))
  }
  fitted <- X %*% b
  predicted <- rowSums(u) + fitted
  e <- y - predicted
  theta <- as.vector(theta)
  if (is.null(names(G))) names(theta) <- c(paste("G", 1:(np - 1), sep = ""), "E")
  if (!is.null(names(G))) names(theta) <- c(names(G)[-np], "E")
  if (is.null(names(G))) colnames(u) <- c(paste("G", 1:(np - 1), sep = ""))
  if (!is.null(names(G))) colnames(u) <- names(G)[-np]

  return(list(
    y = y, X = X, b = b, vb = vb, g = u, e = e, fitted = fitted, predicted = predicted,
    Py = Py, Vy = Vy, theta = theta, asd = theta.cov, llik = llik, niter = it, trPG = trPG,
    trVG = trVG, ids = names(y), yVy = yVy
  ))
}

cvreml <- function(y = NULL, X = NULL, GRMlist = NULL, G = NULL, theta = NULL, ids = NULL, validate = NULL,
                   maxit = 100, tol = 0.00001, bin = NULL, ncores = 1, wkdir = getwd(), verbose = FALSE,
                   makeplots = FALSE) {
  
  n <- length(y)
  if(is.null(X)) X <- model.matrix(y~1)
  theta <- NULL
  #theta <- yobst <- yobsv <- ypredt <- ypredv <- NULL
  #yot <- ypt <- yft <- yat <- yrt <- yov <- ypv <- yfv <- yav <- yrv <- NULL
  res <- NULL
  
  if (is.matrix(validate)) {
    cvnames <- colnames(validate)
    validate <- as.data.frame(validate, stringsAsFactors=FALSE)
    names(validate) <- cvnames
  }
  nv <- length(validate)
  cvnames <- names(validate)
  if(is.null(cvnames)) {
    cvnames <- paste0("CV",1:nv)
    names(validate) <- cvnames
  }
  
  typeoftrait <- "quantitative"
  
  if (nlevels(factor(y)) == 2) typeoftrait <- "binary"

  #training <- validation <- NULL
  training <- validation <- vector(mode="list",length=nv)
  #ghatt <- ghatv <- vector(mode="list",length=nv)
  
  for (i in 1:nv) {
    v <- validate[[i]] # index for validation
    t <- (1:n)[-v] # index for training
    try( fit <- remlr(y = y[t], X = X[t, ], G = lapply(G, function(x) {
      x[t, t]
    }), verbose = verbose))

    
    if(!class(fit)=="try-error") {
      theta <- rbind(theta, as.vector(fit$theta))
      np <- length(fit$theta)
      ypredt <- X[t, ] %*% fit$b # fixed for training
      ypredv <- X[v, ] %*% fit$b # fixed for validation
      ghattrain <- ghatval <- NULL # random for training and validation
      for (j in 1:(np - 1)) {
        ypredt <- ypredt + G[[j]][t, t] %*% fit$Py * fit$theta[j]  # fixed + random for training
        ghattrain <- cbind(ghattrain, G[[j]][t, t] %*% fit$Py * fit$theta[j]) # random for validation
        ypredv <- ypredv + G[[j]][v, t] %*% fit$Py * fit$theta[j] # fixed + random for validation
        ghatval <- cbind(ghatval, G[[j]][v, t] %*% fit$Py * fit$theta[j]) # random for validation
      }
      yobst <- y[t] # observation for training
      yobsv <- y[v] # observation for validation
      
      if (!is.atomic(validate)) res <- rbind(res, acc(yobs = yobsv, ypred = ypredv, typeoftrait = typeoftrait))
      #yot <- c(yot, yobst)
      #ypt <- c(ypt, ypredt)
      #yov <- c(yov, yobsv)
      #ypv <- c(ypv, ypredv)

      yft <- X[t, ] %*% fit$b # compute fixed effects for training
      yfv <- X[v, ] %*% fit$b # compute fixed effects for validation
      yat <- yobst - yft # compute phenotype adjusted for fixed effects for training
      yav <- yobsv - yfv # compute phenotype adjusted for fixed effects for validation
      yrt <- yobst - ypredt # compute residuals for training
      yrv <- yobsv - ypredv # compute residuals for validation
      
      training[[i]]  <- cbind(yobst, ypredt, yft, yat, yrt, ghattrain)
      validation[[i]]  <- cbind(yobsv, ypredv, yfv, yav, yrv, ghatval)

      tnames <- as.character(1:length(y))      
      if(!is.null(names(y))) tnames <- names(y) 
      rownames(training[[i]]) <- tnames[t]
      rownames(validation[[i]]) <- tnames[v]

      colnames(training[[i]]) <- colnames(validation[[i]]) <- c("yobs", "ypred", "yfix", "yadj", "yres", names(G))
    }
    
    #colnames(ghattrain) <- colnames(ghatval) <- names(G)
    #colnames(ghatval) <- names(G)
    #ghatt[[i]] <- ghattrain
    #ghatv[[i]] <- ghatval
    
  }
  
  
  # PSO this should be left out  
  #  if (is.atomic(validate)) res <- matrix(acc(yobs = yov, ypred = ypv, typeoftrait = typeoftrait), nrow = 1)
  # END PSO this should be left out  
  
  # if(is.atomic(validate)) res <- matrix(qcpred(yobs=yo,ypred=yp,typeoftrait=typeoftrait),nrow=1)
  
  res <- as.data.frame(res)
  rownames(res) <- cvnames
  names(training) <- names(validation) <- cvnames
  
  #names(res) <- c("Corr", "R2", "Nagel R2", "AUC", "intercept", "slope", "MSPE")
  if (is.null(names(G))) names(G) <- paste("G", 1:(np - 1), sep = "")
  colnames(theta) <- c(names(G), "E")
  theta <- as.data.frame(round(theta, 3))
  rownames(theta) <- cvnames
  if (makeplots) {
    layout(matrix(1:4, ncol = 2))
    boxplot(res$Corr, main = "Predictive Ability", ylab = "Correlation")
    boxplot(res$MSPE, main = "Prediction Error", ylab = "MSPE")
    boxplot(theta, main = "Estimates", ylab = "Variance")
    plot(y = yo, x = yp, xlab = "Predicted", ylab = "Observed")
    coef <- lm(yo ~ yp)$coef
    abline(a = coef[1], b = coef[2], lwd = 2, col = 2, lty = 2)
  }
  
  # return(list(accuracy = res, theta = theta, yobst = yot, ypredt = ypt, yrest = yrest, ghatt=ghatt, yobsv = yov, ypredv = ypv, yresv = yresv, ghatv=ghatv))
  #return(list(accuracy = res, theta = theta, training = training, validation = validation, ghatt=ghatt, ghatv=ghatv))
  return(list(accuracy = res, theta = theta, training = training, validation = validation))
  
}  

# cvreml <- function(y = NULL, X = NULL, GRMlist = NULL, G = NULL, theta = NULL, ids = NULL, validate = NULL,
#                    maxit = 100, tol = 0.00001, bin = NULL, ncores = 1, wkdir = getwd(), verbose = FALSE,
#                    makeplots = FALSE) {
#   n <- length(y)
#   theta <- yobs <- ypred <- yo <- yp <- NULL
#   res <- NULL
#   if (is.matrix(validate)) validate <- as.data.frame(validate)
#   nv <- length(validate)
#   typeoftrait <- "quantitative"
#   if (nlevels(factor(y)) == 2) typeoftrait <- "binary"
#   
#   ghat <- vector(mode="list",length=nv)
#   
#   for (i in 1:nv) {
#     v <- validate[[i]]
#     t <- (1:n)[-v]
#     try( fit <- remlr(y = y[t], X = X[t, ], G = lapply(G, function(x) {
#       x[t, t]
#     }), verbose = verbose))
#     if(!class(fit)=="try-error") {
#     theta <- rbind(theta, as.vector(fit$theta))
#     np <- length(fit$theta)
#     ypred <- X[v, ] %*% fit$b
#     ghatv <- NULL
#     for (j in 1:(np - 1)) {
#       ypred <- ypred + G[[j]][v, t] %*% fit$Py * fit$theta[j]
#       ghatv <- cbind(ghatv, G[[j]][v, t] %*% fit$Py * fit$theta[j])
#     }
#     yobs <- y[v]
#     if (!is.atomic(validate)) res <- rbind(res, acc(yobs = yobs, ypred = ypred, typeoftrait = typeoftrait))
#     yo <- c(yo, yobs)
#     yp <- c(yp, ypred)
#     }
#     colnames(ghatv) <- names(G)
#     ghat[[i]] <- ghatv
#   }
# 
#   if (is.atomic(validate)) res <- matrix(acc(yobs = yo, ypred = yp, typeoftrait = typeoftrait), nrow = 1)
#   # if(is.atomic(validate)) res <- matrix(qcpred(yobs=yo,ypred=yp,typeoftrait=typeoftrait),nrow=1)
#   res <- as.data.frame(res)
#   names(res) <- c("Corr", "R2", "Nagel R2", "AUC", "intercept", "slope", "MSPE")
#   if (is.null(names(G))) names(G) <- paste("G", 1:(np - 1), sep = "")
#   colnames(theta) <- c(names(G), "E")
#   theta <- as.data.frame(round(theta, 3))
#   if (makeplots) {
#     layout(matrix(1:4, ncol = 2))
#     boxplot(res$Corr, main = "Predictive Ability", ylab = "Correlation")
#     boxplot(res$MSPE, main = "Prediction Error", ylab = "MSPE")
#     boxplot(theta, main = "Estimates", ylab = "Variance")
#     plot(y = yo, x = yp, xlab = "Predicted", ylab = "Observed")
#     coef <- lm(yo ~ yp)$coef
#     abline(a = coef[1], b = coef[2], lwd = 2, col = 2, lty = 2)
#   }
#   return(list(accuracy = res, theta = theta, yobs = yo, ypred = yp, ghat=ghat))
# }


####################################################################################################################
# REML interface functions for fortran linked library




remlf <- function(y = NULL, X = NULL, GRMlist = NULL, G = NULL, theta = NULL, ids = NULL, maxit = 100,
                  tol = 0.00001, ncores = 1, verbose = FALSE) {
  if (!is.null(G)) writeGRM(GRM = G)

  ids <- names(y)

  n <- length(y)
  nf <- ncol(X)
  #if (!is.null(G)) rfnames <- paste("G", 1:length(G), sep = "")
  #if (!is.null(G)) rfnames <- paste(getwd(), rfnames, sep = "/")
  if (!is.null(GRMlist)) rfnames <- GRMlist$fnG
  nr <- length(rfnames) + 1
  if (!is.null(G)) ngr <- nrow(G[[1]])
  if (!is.null(G)) indx <- match(ids, rownames(G[[1]]))
  if (!is.null(GRMlist)) ngr <- GRMlist$n
  if (!is.null(GRMlist)) indx <- match(ids, GRMlist$idsG)

  fnr <- paste(paste(sample(letters, 10, replace = TRUE), collapse = ""), ".qgg", sep = "")
  #write.table(as.character(rfnames), file = "param.qgg", quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)

  fnGCHAR = matrix(as.integer(0), nrow = nr-1, ncol = 1000)
  ncharsg = rep(as.integer(0), length = nr-1)
  for (i in 1:(nr-1)) {
    nchars <- as.integer(nchar(as.character(rfnames[i])))       
    fnGCHAR[i,1:nchars] <- as.integer(unlist(sapply(as.character(rfnames[i]),charToRaw),use.names=FALSE))
    ncharsg[i] = nchars
  }

  # reml(n,nf,nr,tol,maxit,ncores,ngr,indx,y,X,theta,ai,b,varb,u,Vy,Py,llik,trPG,trVG,ncharsg,fnGCHAR)

  fit <- .Fortran("reml",
    n = as.integer(n),
    nf = as.integer(nf),
    nr = as.integer(nr),
    tol = as.double(tol),
    maxit = as.integer(maxit),
    ncores = as.integer(ncores),
    #fnr = as.character(fnr),
    # rfnames = as.character(rfnames),
    ngr = as.integer(ngr),
    indx = as.integer(indx),
    y = as.double(y),
    X = matrix(as.double(X), nrow = nrow(X)),
    theta = as.double(theta),
    ai = matrix(as.double(0), nrow = nr, ncol = nr),
    b = as.double(rep(0, nf)),
    varb = matrix(as.double(0), nrow = nf, ncol = nf),
    u = matrix(as.double(0), nrow = n, ncol = nr),
    Vy = as.double(rep(0, n)),
    Py = as.double(rep(0, n)),
    llik = as.double(0),
    trPG = as.double(rep(0, nr)),
    trVG = as.double(rep(0, nr)),
    ncharsg = as.integer(ncharsg),
    fnGCHAR = matrix(as.integer(fnGCHAR), nrow = nr-1, ncol = 1000),
    PACKAGE = "qgg"
  )
  #file.remove("param.qgg")
  
  fit$ids <- names(y)
  fit$yVy <- sum(y * fit$Vy)
  fit$wd <- getwd()
  fit$GRMlist <- GRMlist
  rownames(fit$u) <- names(y)
  colnames(fit$u) <- c(paste("G", 1:(ncol(fit$u) - 1), sep = ""), "E1")
  fit$g <- fit$u[, 1:(nr - 1)]
  fit$e <- fit$u[, nr]
  fit$u <- NULL
  np <- length(fit$theta)
  names(fit$theta) <- c(paste("G", 1:(np - 1), sep = ""), "E")

  return(fit)
}

#' Compute Genomic BLUP values

#' @description
#' Compute Genomic BLUP values based on linear mixed model fit output from greml

#' @param GRM list of one or more genomic relationship matrices
#' @param GRMlist list providing information about GRM matrix stored in binary files on disk
#' @param fit list object output from greml function
#' @param ids vector of ids for which BLUP values is computed
#' @param idsRWS vector of row ids in GRM for which BLUP values is computed
#' @param idsCLS vector of column ids in GRM for which BLUP values is computed

#' @export
#'

gblup <- function(GRMlist = NULL, GRM = NULL, fit = NULL, ids = NULL, idsCLS = NULL, idsRWS = NULL) {
  GRMlist <- fit$GRMlist
  fnG <- GRMlist$fnG
  Py <- fit$Py
  names(Py) <- fit$ids
  g <- NULL
  nr <- length(fnG)
  if (is.null(idsCLS)) idsCLS <- fit$ids
  if (is.null(idsRWS)) idsRWS <- fit$ids
  if (!is.null(ids)) idsRWS <- ids
  if (sum(!idsRWS %in% GRMlist$idsG) > 0) stop("Error some ids not found in idsG")
  for (i in 1:nr) {
    GRMlist$fnG <- fnG[i]
    G <- getGRM(GRMlist = GRMlist, idsCLS = idsCLS, idsRWS = idsRWS)
    g <- cbind(g, G %*% Py[idsCLS] * fit$theta[i])
  }
  colnames(g) <- paste("G", 1:nr, sep = "")
  return(g)
}
