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
#' @param ncores number of cores
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

greml <- function(y = NULL, X = NULL, Glist=NULL, G=NULL, theta=NULL, ids=NULL, validate=NULL, maxit=100, tol=0.00001,bin=NULL,ncores=1,wkdir=getwd(), verbose=FALSE, makeplots=FALSE,interface="R")
{
  if(interface=="R") { 
    if (is.null(validate)) fit <- remlR(y=y, X=X, Glist=Glist, G=G, theta=theta, ids=ids, maxit=maxit, tol=tol, bin=bin, ncores=ncores, verbose=verbose, wkdir=wkdir)
    if (!is.null(validate)) fit <- cvreml(y=y, X=X, Glist=Glist, G=G, theta=theta, ids=ids, validate=validate, maxit=maxit, tol=tol, bin=bin, ncores=ncores, verbose=verbose, wkdir=wkdir, makeplots=makeplots)
  }
  if(!is.null(bin)) { 
    fit <- remlF(y=y, X=X, Glist=Glist, G=G, ids=ids, theta=theta, maxit=maxit, tol=tol, bin=bin, ncores=ncores, verbose=verbose, wkdir=wkdir)
  }
  if(interface=="fortran") { 
    fit <- freml(y=y, X=X, Glist=Glist, G=G, theta=theta, ids=ids, maxit=maxit, tol=tol, ncores=ncores, verbose=verbose) 
  }        
  return(fit)  
}  


####################################################################################################################
# REML interface functions for fortran executable

remlF <- function(y = NULL, X = NULL, Glist = NULL, G = NULL, theta = NULL, ids = NULL, maxit = 100, tol = 0.00001, bin = NULL, ncores = 1, wkdir = getwd(), verbose = FALSE ) {
#greml <- function(y = NULL, X = NULL, Glist = NULL, G = NULL, ids = NULL, theta = NULL, maxit = 100, tol = 0.00001, bin = NULL, ncores = 1, wkdir = getwd()) {
    
	write.reml(y = as.numeric(y), X = X, G = G)
	n <- length(y)
	nf <- ncol(X)
	if (!is.null(G)) fnamesG <- paste("G", 1:length(G), sep = "")
	if (!is.null(Glist$fnG)) fnamesG <- Glist$fnG
	nr <- length(fnamesG) + 1
 	if (is.null(ids)) {indxG <- c(n, 1:n)} 
	if (!is.null(ids)) {indxG <- c(Glist$n, match(ids, Glist$idsG))} 
	write.table(indxG, file = "indxg.txt", quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)

	write.table(paste(n, nf, nr, maxit, ncores), file = "param.txt", quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
	if (is.null(theta)) theta <- rep(sd(y) / nr, nr)
	#if (is.null(theta)) theta <- rep(var(y) / nr, nr)
	write.table(t(theta), file = "param.txt", quote = FALSE, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE)
	write.table(tol, file = "param.txt", quote = FALSE, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE)
	write.table(fnamesG, file = "param.txt", quote = TRUE, sep = " ", append = TRUE, col.names = FALSE, row.names = FALSE)

	execute.reml(bin = bin,  ncores = ncores)
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
	   nr <- nrow(G[[i]])
	   for (j in 1:nr) {
		writeBin(as.double(G[[i]][1:nr, j]), fileout,size=8,endian = "little")
	   }
	   close(fileout)
	  }
	}
          
}

execute.reml  <- function (bin = NULL, ncores = ncores) {

	HW <- Sys.info()["machine"]
	OS <- .Platform$OS.type
	if (OS == "windows") {
		"my.system" <- function(cmd) {return(system(paste(Sys.getenv("COMSPEC"), "/c", cmd)))}
        
		#my.system(paste("set MKL_NUM_THREADS = ", ncores))
		test <- my.system(paste(shQuote(bin), " < param.txt > reml.lst", sep = ""))
	}
	if (!OS == "windows") {
		system(paste("cp ", bin, " reml.exe", sep = ""))
		#system(paste("export MKL_NUM_THREADS=", ncores))
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

remlR <- function(y=NULL, X=NULL, Glist=NULL, G=NULL, theta=NULL, ids=NULL, maxit=100, tol=0.00001, bin=NULL,ncores=1,wkdir=getwd(), verbose=FALSE )
  
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
    output <- c(1:10,seq(11,maxit,5))     
    if (verbose & it%in%output) print(paste(c("Iteration:",it,"Theta:",round(theta,2)), sep=""))
    if (it==maxit) break
  }
  if (verbose) print(paste(c("Converged at Iteration:",it,"Theta:",round(theta,2)), sep=""))
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
  e <- y-predicted
  theta <- as.vector(theta)
  if(is.null(names(G))) names(theta) <- c(paste("G",1:(np-1),sep=""),"E")
  if(!is.null(names(G))) names(theta) <- c(names(G)[-np],"E")
  
  return(list( y=y, X=X, b=b, vb=vb, g=u, e=e, fitted=fitted, predicted=predicted, Py=Py, Vy=Vy, theta=theta, asd=theta.cov, llik=llik, niter=it,trPG=trPG, trVG=trVG,ids=names(y),yVy=yVy   ))
}


cvreml <- function(y=NULL, X=NULL, Glist=NULL, G=NULL, theta=NULL, ids=NULL, validate=NULL, maxit=100, tol=0.00001,bin=NULL,ncores=1,wkdir=getwd(), verbose=FALSE, makeplots=FALSE)
{
  n <- length(y)     
  theta <- yobs <- ypred <- yo <- yp <- NULL
  res <- NULL
  if(is.matrix(validate)) validate <- as.data.frame(validate)
  nv <- length(validate)
  for (i in 1:nv) {
    v <- validate[[i]]  
    t <- (1:n)[-v]
    fit <- remlR( y=y[t], X=X[t,], G=lapply(G,function(x){x[t,t]}), verbose=verbose)
    theta <- rbind(theta, as.vector(fit$theta))
    np <- length(fit$theta)
    ypred <- X[v, ] %*% fit$b
    for (j in 1:(np-1)) {
      ypred <- ypred + G[[j]][v,t]%*%fit$Py*fit$theta[j]
    }
    yobs <- y[v]
    if(!is.atomic(validate)) res <- rbind(res,qcpred(yobs=yobs,ypred=ypred,typeoftrait=typeoftrait))
    yo <- c(yo, yobs)
    yp <- c(yp, ypred)
  }
  typeoftrait <- "quantitative"
  if(nlevels(factor(y))==2) typeoftrait <- "binary" 
 
  if(is.atomic(validate)) res <- matrix(qcpred(yobs=yo,ypred=yp,typeoftrait=typeoftrait),nrow=1)
  res <- as.data.frame(res)
  names(res) <- c("Corr","R2","Nagel R2", "AUC", "intercept", "slope", "MSPE")
  if(is.null(names(G))) names(G) <- paste("G",1:(np-1),sep="")
  colnames(theta) <- c(names(G),"E")
  theta <- as.data.frame(round(theta,3))
  if (makeplots) {
   layout(matrix(1:4, ncol = 2))
   boxplot(res$Corr, main = "Predictive Ability", ylab = "Correlation")
   boxplot(res$MSPE, main = "Prediction Error", ylab = "MSPE")
   boxplot(theta, main = "Estimates", ylab = "Variance")
   plot(y=yo, x=yp, xlab = "Predicted", ylab = "Observed")
   coef <- lm(yo ~ yp)$coef
   abline(a = coef[1], b = coef[2], lwd = 2, col = 2, lty = 2)
  }
  return(list(pred=res,theta=theta,yobs=yo,ypred=yp))
}


####################################################################################################################
# REML interface functions for fortran linked library

#' @export
#'

freml <- function(y = NULL, X = NULL, Glist = NULL, G = NULL, theta = NULL, ids = NULL, maxit = 100, tol = 0.00001, ncores = 1, verbose = FALSE ) 

{

   if(!is.null(G)) writeG(G = G)

   ids <- names(y)

   n <- length(y)
   nf <- ncol(X)
   if (!is.null(G)) rfnames <- paste("G", 1:length(G), sep = "")
   if (!is.null(G)) rfnames <- paste(getwd(),rfnames,sep="/")	 
   if (!is.null(Glist)) rfnames <- Glist$fnG
   nr <- length(rfnames) + 1
   if (!is.null(G)) ngr <- nrow(G[[1]])
   if (!is.null(G)) indx <- match(ids, rownames(G[[1]]))
   if (!is.null(Glist)) ngr <- Glist$n
   if (!is.null(Glist)) indx <- match(ids, Glist$idsG)
   
   fit <- .Fortran("reml", 
          n = as.integer(n),
          nf = as.integer(nf),
          nr = as.integer(nr),
          tol = as.double(tol),
          maxit = as.integer(maxit),
          ncores = as.integer(ncores),
          rfnames = as.character(rfnames),
          ngr = as.integer(ngr),
          indx = as.integer(indx),
          y = as.double(y),
          X = matrix(as.double(X),nrow=nrow(X)),
          theta = as.double(theta),
          ai = matrix(as.double(0),nrow=nr,ncol=nr),
          b = as.double(rep(0,nf)),
          varb = matrix(as.double(0),nrow=nf,ncol=nf),
          u = matrix(as.double(0),nrow=n,ncol=nr),
          Vy = as.double(rep(0,n)),
          Py = as.double(rep(0,n)),
          llik = as.double(0),
          trPG = as.double(rep(0,nr)),
          trVG = as.double(rep(0,nr)),
          PACKAGE = 'qgg'
   )
   
   fit$ids <- names(y)
   fit$yVy <- sum(y * fit$Vy)
   fit$wd <- getwd()
   fit$Glist <- Glist
   fit$g <- fit$u[,1:(nr-1)]
   fit$e <- fit$u[,nr]
   fit$u <- NULL 

   return(fit)
}

#' @export
#'


writeG <- function(G = NULL) {
   if (!is.null(G)) {
     for (i in 1:length(G)) {
      fileout <- file(paste("G", i, sep = ""), "wb")
      nr <- nrow(G[[i]])
      for (j in 1:nr) {
        writeBin(as.double(G[[i]][1:nr, j]), fileout,size=8,endian = "little")
      }
      close(fileout)
     }
   }
}
