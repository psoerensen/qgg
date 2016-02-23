####################################################################################################################
#
#    Module 1: LMM, REML, BLUP
#    Module 2: Set Test
#    Module 3: LMM, Bayesian
#
####################################################################################################################



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
#' @param fm a formula with model statement for the fixed factors in the linear mixed model
#' @param weights a vector of weights for the residual variance
#' @param W a matrix centered and scaled genotypes or other types of molecular data
#' @param sets a list of marker sets corresponding to column names in W 
#' @param K a list of correlation matrices
#' @param data a data frame containing the phenotypic observations and fixed factors specified in the model statements
#' @param validate a matrix or a list with the ids of validation sets corresponding to the rows in data
#' @param mkplots a logical indicating whether or not to make plots
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
#' @references Mapping Variants to Gene Ontology Categories Improves Genomic Prediction for Quantitative Traits in Drosophila melanogaster. Under review Genetics (2016). Edwards SM, Sørensen IF, Sarup P, Mackay TF, Sørensen P. 
#' @references Genomic BLUP Derived Set Tests Identify Genetic Variants Associated with Schizophrenia in Functionally Associated Biological Processes. Under review, Genetics (2015). Rohde PD, Demontis D, The GEMS Group, Børglum AD, Sørensen P.
#' @references Partitioning of genomic variance reveals biologic pathways associated with udder health and milk production traits in dairy cattle. GSE (2015) 47:60. Edwards SM, Thomsen B, Madsen P, Sørensen P.
#' @references Increased prediction accuracy using a genomic feature model including prior information on quantitative trait locus regions in purebred Danish Duroc pigs. BMC Genetics (2016) 17:11. Sarup P, Jensen J, OstersenT, Henryon M, Sørensen P.
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(2000000),ncol=10000)
#'   colnames(W) <- as.character(1:ncol(W))
#' y <- rowSums(W[,1:10]) + rowSums(W[,1001:1010]) + rnorm(nrow(W))
#' 
#' # REML analyses and cross validation 
#' data <- data.frame(y=y,mu=1)
#' fm <- y ~ mu
#' 
#' setsGB <- list(A=colnames(W))
#' setsGF <- list(C1=colnames(W)[1:1000],C2=colnames(W)[1001:2000],C2=colnames(W)[2000:10000])
#' setsGT <- list(C1=colnames(W)[1:10],C2=colnames(W)[1001:1010],C2=colnames(W)[1:10000])
#' 
#' n <- length(y)
#' fold <- 10
#' nsets <- 50
#' validate <- replicate(nsets,sample(1:n,as.integer(n/fold)))
#' 
#' fitGB <- gfm(fm=fm,W=W,sets=setsGB,data=data,validate=validate)
#' fitGF <- gfm(fm=fm,W=W,sets=setsGF,data=data,validate=validate)
#' fitGT <- gfm(fm=fm,W=W,sets=setsGT,data=data,validate=validate)
#' 
#' 
#' validClust <- list(C1=1:40,C2=41:77)
#' 
#' fitGB <- gfm(fm=fm,W=W,sets=setsGB,data=data,validate=validClust)
#' fitGF <- gfm(fm=fm,W=W,sets=setsGF,data=data,validate=validClust)
#' fitGT <- gfm(fm=fm,W=W,sets=setsGT,data=data,validate=validClust)
#' 
#' 
#' # REML analyses and variance decomposition
#' 
#' fitGB <- gfm(fm=fm,W=W,data=data)
#' covSets(W=W,g=fitGB$f,sets=setsGF)$cvf
#' covSets(W=W,g=fitGB$f,sets=setsGT)$cvf
#' 
#' 
#' fitGT <- gfm(fm=fm,W=W,sets=setsGT,data=data)
#' covSets(W=W,g=fitGT$f[,1],sets=setsGT)$cvf
#' covSets(W=W,g=fitGT$f[,2],sets=setsGT)$cvf
#' 
#' cvgf1 <- covSets(W=W,g=fitGT$f[,1],sets=setsGT, level2=TRUE)
#' cvgf2 <- covSets(W=W,g=fitGT$f[,2],sets=setsGT, level2=TRUE)
#' 
#' bars1 <- cvgf1$cvf[,2]
#' bars11 <- cvgf1$cvfSet[[1]][,5]
#' bars12 <- cvgf1$cvfSet[[2]][,5]
#' 
#' bars2 <- cvgf2$cvf[,2]
#' bars21 <- cvgf2$cvfSet[[1]][,5]
#' bars22 <- cvgf2$cvfSet[[2]][,5]
#' 
#' layout(matrix(1:6,ncol=2,byrow=TRUE))
#' barplot(bars1)
#' barplot(bars2)
#' barplot(bars11)
#' barplot(bars12)
#' barplot(bars21)
#' barplot(bars22)
#' 
#' 
#' 
#' # REML analyses 
#' data <- data.frame(y=y,mu=1)
#' fm <- y ~ mu
#' fit <- gfm(fm=fm,W=W,sets=list(colnames(W)),data=data)
#' 
#' # REML analyses 
#' pc1 <- rnorm(200)
#' pc2 <- rnorm(200)
#' f1 <- factor(rep(1:2,100))
#' fmf <- y~pc1+pc1+f1
#' fitGB <- gfm(fm=fmf,W=W,data=data.frame(y,pc1,pc1,f1),validate=validate)
#' fitGF <- gfm(fm=fmf,W=W,sets=setsGF,data=data.frame(y,pc1,pc1,f1),validate=validate)
#' fitGT <- gfm(fm=fmf,W=W,sets=setsGT,data=data.frame(y,pc1,pc1,f1),validate=validate)
#' 
#' @export
#' 

gfm <- function(fm=NULL, weights=NULL, W=NULL,sets=NULL,K=NULL, data=NULL,validate=NULL, mkplots=TRUE) {
     
     mf <- model.frame(fm, data = data, na.action = na.pass)
     mf <- eval(mf, parent.frame())
     y <- model.response(mf)
     X <- model.matrix(mf, data=data)
     if ( ncol(X)==2) {
          if ( sum(colnames(X)==c("(Intercept)", "mu"))==2) X <- matrix(X[,1],ncol=1)
     }
     n <- length(y)
     
     # Validation sets
     if(is.null(validate)) { validate <- matrix(-c(1:n),nrow=n,ncol=1)}
     if(is.matrix(validate)) nvsets <- ncol(validate)
     if(is.list(validate)) nvsets <- length(validate)
     
     # Check input data
     if (!is.null(sets)) { if (any(!sapply(sets,is.character))) stop("Sets not character variables") }
     if (!is.null(W)) { if (any(!is.character(colnames(W)))) stop("Column names of W is not character variables") }
     if (!is.null(sets)) { if (any(!unlist(sapply(sets,function(x){x%in%colnames(W)})) )) 
          stop("Sets does not match column names of W") }
     if (!is.null(sets)) { if (any(!sapply(sets,is.character))) 
          stop("Sets not character variables") }
     #if (any(!apply(validate,2,is.integer))) stop("Validation sets not integer variables")
     if(is.matrix(validate)) {if(max(validate)>n) stop("Validation sets contains integer values larger than n=number of observations") }
     if(is.list(validate)) {if(max(unlist(validate))>n) stop("Validation sets contains integer values larger than n=number of observations") }
     if (!is.null(weights)) { if (!length(weights)==nrow(W)) stop("nrow(W) not equal to length(weights)") }
     
     # Compute K matrices
     if(is.null(sets)&is.null(K)) { K[[1]] <- (W%*%t(W))/ncol(W) }
     nk <- length(K)
     nsets <- NULL
     if(!is.null(sets)) {
          nsets <- length(sets)
          m <- sapply(sets,length)
          for (i in 1:nsets) { K[[i+nk]] <- (W[,sets[[i]]]%*%t(W[,sets[[i]]]))/m[i] }
          names(K)[(nk+1):(nk+nsets)] <- names(sets)
     }
     if(is.null(nsets)) {nsets <- nk}
     if(!is.null(sets)) {nsets <- nsets + nk}
     if (is.null(names(K))) names(K) <- paste("C",1:nsets,sep="")
     K[[nsets+1]] <- diag(1,n)
     identity <- TRUE
     if (!is.null(weights)) {
          K[[nsets+1]] <- diag(1/weights)
          identity <- FALSE
     }
     names(K)[nsets+1] <- "e"
     
     # Model statements random effects
     vfm <- paste("~",paste(paste("K[[",1:nsets,"]][t,t]",sep=""),collapse="+"))
     if (!is.null(weights)) { vfm <- paste("~",paste(paste("K[[",1:(nsets+1),"]][t,t]",sep=""),collapse="+"))}
     vfm <- as.formula(vfm)
     
     # Fit mixed model and predict random effect (with or without observation)
     pa <- mspe <- sigmas <- ypred <- yobs <- NULL
     s <- vs <- Vf <- NULL
     for (i in 1:nvsets) {
          if(is.matrix(validate)) v <- validate[,i]
          if(is.list(validate)) v <- validate[[i]]
          t <- (1:n)[-v]
          fit <- regress(fm, vfm ,identity=identity,verbose=1, pos=rep(TRUE,nsets+1),
                         data=droplevels(data[t,]))
          V <- matrix(0,ncol=n,nrow=n)
          for (j in 1:(nsets+1)) { V <- V + K[[j]]*fit$sigma[j] }
          f <- fv <- NULL
          for (j in 1:nsets) { f <- cbind(f,(K[[j]][t,t]*fit$sigma[j])%*%fit$W%*%(y[t]-fit$fitted)) }
          for (j in 1:nsets) { fv <- cbind(fv,(K[[j]][v,t]*fit$sigma[j])%*%fit$W%*%(y[t]-fit$fitted)) }
          if(nsets==1&nk==0) {
               Vf <- NULL
               P <- fit$W%*%fit$Q
               for (j in 1:nsets) { 
                    Vf[[j]] <- (fit$sigma[j]*K[[j]][t,t])%*%P%*%V[t,t]%*%P%*%(K[[j]][t,t]*fit$sigma[j])
                    bvb <- f2b(W=W[t,sets[[j]] ],f=f[,j], Vf=Vf[[j]])
                    if (nsets==1) {
                         s <- bvb$b
                         vs <- bvb$vb
                    }
                    if (nsets>1) {
                         s[[j]] <- bvb$b
                         vs[[j]] <- bvb$vb
                    }
               }
          }
          
          if (length(t)<n) {
               if(any(!fit$X==X[t,])) stop("!fit$X[t,]==X[t,] problem with design matrices for fixed effects")
               #yhat <- fit$fitted[1] + V[v,t]%*%fit$W%*%(y[t]-fit$fitted)
               yhat <- X[v,]%*%fit$beta + V[v,t]%*%fit$W%*%(y[t]-fit$fitted)
               pa <- c(pa, cor(yhat,y[v])) 
               mspe <- c(mspe,sum((yhat-y[v])**2)/length(v))
               sigmas <- rbind(sigmas,fit$sigma)
               ypred <- c(ypred,yhat)
               yobs <- c(yobs,y[v])
          }
     }
     # Make plots
     if (nvsets>1&mkplots) {
          colnames(sigmas) <- names(K)
          layout(matrix(1:4,ncol=2))
          boxplot(pa,main="Predictive Ability",ylab="Correlation")
          boxplot(mspe,main="Prediction Error",ylab="MSPE")
          boxplot(sigmas,main="Estimates",ylab="Variance")
          plot(yobs,ypred,ylab="Predicted",xlab="Observed") 
          coef  <- lm(ypred~yobs)$coef
          abline(a=coef[1],b=coef[2],lwd=2,col=2,lty=2)
     }
     return(list(f=f,fv=fv,Vf=Vf,s=s,vs=vs,sigmas=sigmas,pa=pa,mspe=mspe,ypred=ypred,yobs=yobs,fit=fit,validate=validate))
}



####################################################################################################################
     f2b <- function(W=NULL,f=NULL, Vf=NULL) {
          if (!nrow(W)==length(f)) {  stop("nrow(W)==length(f) is false") }
          WW <- W%*%t(W)
          WWi <- (MASS:::ginv)(WW)
          b <- t(W)%*%WWi%*%f
          if (!is.null(Vf)) {
               WWWi <- crossprod(W,WWi)
               remove(WWi)
               remove(W)
               WWWiVf <- tcrossprod(WWWi,Vf)
               vb <- rowSums(WWWiVf*WWWi)
               remove(WWWi)
               remove(Vf)
          }
          return(list(b=b[,1],vb=vb))
     }

     computeG <- function(W=NULL,pdf=TRUE) {
          SS <- tcrossprod(W)        # compute crossproduct all SNPs
          N <- tcrossprod(!W==0)     # compute number of observation all SNPs
          G <- SS/N
          if (pdf) G <- makepdf(G)
          return(list(G=G,SS=SS,N=N))
     }  
     
     makepdf <- function(G=NULL, tol=0.0001) {
          rn <- rownames(G)
          e <- eigen(G)                                    # eigen value decomposition of the matrix G
          U <- e$vectors                                   # eigen vectors
          e <- e$values                                    # eigen values
          e[e<tol] <- tol
          D <- diag(e)                                   
          G <- U%*%D%*%t(U)          	                        # compute pdf
          colnames(G) <- rownames(G) <- rn
          return(G)                      
     } 
     
     ginv <- function(G=NULL, tol=NULL) {
          rn <- rownames(G)
          e <- eigen(G)                                    # eigen value decomposition of the matrix G
          U <- e$vectors                                   # eigen vectors
          e <- e$values                                    # eigen values
          #e[e<tol] <- tol
          #D <- diag(1/e)                                   # set inverse D to 1/e
          ie <- e
          ie[e>tol] <- 1/e[e>tol]
          #ie[e<tol] <- tol
          ie[e<tol] <- 0
          D <- diag(ie)                                   # set inverse D to 1/e
          G <- U%*%D%*%t(U)           		                        # compute inverse
          ldet <- sum(log(e[e>tol])) 
          #ldet <- sum(log(e)) 
          colnames(G) <- rownames(G) <- rn
          return(list(G=G,ldet=ldet))                      # log determinant 
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
#'                        
#' @details
#' The singler marker test statistics can be obtained from GBLUP and GFBLUP model fits or from standard GWAs. 
#' The distribution of this test statistic under the null hypothesis (associated markers are picked at random from the total 
#' number of tested genetic markers) is difficult to describe in terms of exact or approximate 
#' distributions, and an empirical distribution is required.
#'                        
#' @param stat is a vector of single marker statistics (e.g. marker effects, t-stat, p-values)
#' @param sets is a list of marker sets  - names corresponds to rownames in stat
#' @param nperm is the number of permutations
#' @param W matrix of centered and scaled genotypes (used if method=cvs or score)
#' @param method including sum, cvs, hyperG, score
#' @param threshold used if method=hyperG
#' @return Returns a dataframe including 
#' \item{setT}{marker set test statistics} 
#' \item{nset}{number of markers in the set }
#' \item{p}{p-value for marker set}
#' @author Peter Sørensen
#' @references Genomic BLUP Derived Set Tests Identify Genetic Variants Associated with Schizophrenia in Functionally Associated Biological Processes. Under review, Genetics (2015). Rohde PD, Demontis D, The GEMS Group, Børglum AD, Sørensen P.
#' @references Gene-based Association Approach Identify Genes Across Stress Traits in Fruit Flies. In: Proceedings, 10th World Congress of Genetics Applied to Livestock Production (WCGALP), Vancouver, Canada, 2014. Dunn Rohde P, Edwards SM, Sarup PM, Sørensen P. 
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(2000000),ncol=10000)
#'   colnames(W) <- as.character(1:ncol(W))
#' y <- rowSums(W[,1:10]) + rowSums(W[,1001:1010]) + rnorm(nrow(W))
#' 
#' # REML analyses 
#' data <- data.frame(y=y,mu=1)
#' fm <- y ~ mu
#' fit <- gfm(fm=fm,W=W,sets=list(colnames(W)),data=data)
#' fit$df <- 10
#' fit$p <- pt(fit$s/sqrt(fit$vs),df=fit$df, lower.tail=FALSE) 
#' 
#' sets <- list(A=as.character(1:100),B=as.character(101:1000),C=as.character(1001:5000),D=as.character(5001:10000))

#' # Set test based on sums 
#' res <- setTest(stat=fit$s**2,sets=sets, method="sum", nperm=100)
#' 
#' # Set test based on cvs 
#' res <- setTest(stat=fit$s,W=W, sets=sets,method="cvs", nperm=100)
#' 
#' # Set test based on hyperG 
#' res <- setTest(stat=fit$p,sets=sets, method="hyperG", threshold=0.05)
#' 
#' @export
setTest <- function(stat=NULL, W=NULL, sets=NULL, nperm=NULL, method="sum", threshold=0.05) {
  if(method=="sum") setT <- sumTest(stat=stat, sets=sets, nperm=nperm) 
  if(method=="cvs") setT <- cvsTest(s=stat, W=W, sets=sets, nperm=nperm) 
  if(method=="hyperG") setT <- hgTest(p=stat, sets=sets, threshold=threshold) 
  if(method=="score") setT <- scoreTest(e=e, W=W, sets=sets, nperm=nperm) 
  return(setT)
}


sumTest <- function(stat=NULL, sets=NULL, nperm=NULL, method="sum") {
  if(method=="mean") setT <- sapply( sets,function(x){ mean(stat[x]) })
  if(method=="sum") setT <- sapply( sets,function(x){ sum(stat[x]) })
  if(method=="max") setT <- sapply( sets,function(x){ max(stat[x]) })
  if (!is.null(nperm)) {
    p <- rep(0,length(sets)) 
    n <- length(stat)
    nset <- sapply(sets,length)
    rws <- 1:n
    names(rws) <- names(stat)
    sets <- lapply(sets, function(x) { rws[x] })
    for ( i in 1:nperm ) {
      rws <- sample(1:n,1)
      o <- c(rws:n,1:(rws-1))
      pstat <- stat[o]
      if(method=="mean") setTP <- sapply(sets,function(x){ mean(pstat[x]) })
      if(method=="sum") setTP <- sapply(sets,function(x){ sum(pstat[x]) })
      if(method=="max") setTP <- sapply(sets,function(x){ max(pstat[x]) })
      p <- p + as.numeric(setT>setTP) 
    }
    p <- 1-p/nperm
    setT <- data.frame(setT,nset,p)
  }
  return(setT)
}
msetTest <- function(stat=NULL, sets=NULL, nperm=NULL, method="sum") {
     setT <- apply(stat,2,function(x) { setTest(stat=x,sets=sets,nperm=nperm,method=method) })
     names(setT) <- colnames(stat)
     setT  
} 

####################################################################################################################
#' 
#' Genetic marker set tests based on the covariance test
#'
#' @description
#' Genetic marker set tests based on the covariance statistics for a set of genetic markers
#' 
#' @details
#' The covariance test statistic is derived from a GBLUP( or GFBLUP) model fit. It is a measure of covariance between the total genomic effect for all markers 
#' and the genomic effect for the genetic markers in the genomic feature. It also relates to the explained sums of
#' squares for the genetic markers. 
#' The distribution of this test statistic under the null hypothesis is difficult to describe in terms of exact or approximate 
#' distributions, and an empirical distribution is required.
#' 
#' @param g vector (or list) of genetic effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param s vector (or list) of single marker effects obtained from a linear mixed model fit (GBLUP of GFBLUP)
#' @param W matrix of centered and scaled genotypes (nxm)
#' @param sets list of marker sets corresponding to columns in W
#' @param nperm number of permutations
#' @return Returns a dataframe including 
#' \item{setT}{covariance test statistics} 
#' \item{nset}{number of markers in the set }
#' \item{p}{p-value }
#' @author Peter Sørensen
#' @references Genomic BLUP Derived Set Tests Identify Genetic Variants Associated with Schizophrenia in Functionally Associated Biological Processes. Under review, Genetics (2015). Rohde PD, Demontis D, The GEMS Group, Børglum AD, Sørensen P.
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(2000000),ncol=10000)
#'   colnames(W) <- as.character(1:ncol(W))
#' y <- rowSums(W[,1:10]) + rowSums(W[,1001:1010]) + rnorm(nrow(W))
#' 
#' # REML analyses 
#' data <- data.frame(y=y,mu=1)
#' fm <- y ~ mu
#' fit <- gfm(fm=fm,W=W,sets=list(colnames(W)),data=data)
#' 
#' # CVS set test 
#' sets <- list(A=as.character(1:100),B=as.character(101:1000),C=as.character(1001:5000),D=as.character(5001:10000))
#' res <- cvsTest(s=fit$s,W=W,sets=sets,nperm=100)
#' res
#' 
#' 
#' @export
cvsTest <- function(s=NULL,g=NULL, W=NULL, sets=NULL, nperm=100){
     Ws <-   t(t(W)*as.vector(s))
     if (is.null(g)) g <- W%*%s   
     cvs <- colSums(as.vector(g)*Ws)
     #setT <- setTest( stat=cvs, sets=sets, nperm=nperm, method="sum")$p
     #names(setT) <- names(sets)
     setT <- setTest( stat=cvs, sets=sets, nperm=nperm, method="sum")
     if(!is.null(names(sets)))  rownames(setT) <- names(sets)
     return(setT)
}

scoreTest <- function(e=NULL,W=NULL, sets=NULL, nperm=100){
     we2 <- as.vector((t(W)%*%e)**2)   
     names(we2) <- colnames(W)       
     setT <- setTest( stat=we2, sets=sets, nperm=nperm, method="sum")$p
     return(setT)
}



####################################################################################################################
#' 
#' Genetic marker set tests based on the hyperG test
#' 
#' 
#' @description
#' Genetic marker set tests based on the hyperG test statistics for a set of genetic markers.
#'
#'
#' @details
#' The hyperG marker set test is based on the idea to identify the association between two types of classification 
#' of the markers: 1) classified as being in a predefined set of markers (i.e. genomic feature), 
#' and 2) classified as being associated to the trait phenotype. 
#' This can be formulated and tested in a number of ways. Here we consider a test statistic based 
#' on counting the number of genetic markers in the feature that are associated to trait phenotype. 
#' Test based on the count test statistic is likely to have high power to detect association if the genomic 
#' feature harbours genetic markers with large effects. 
#' Under the null hypothesis (associated markers are picked at random from the total number of tested 
#' genetic markers) it is assumed that the observed count statistics is a realization from a hypergeometric 
#' distribution
#' Test based on the count test statistic is likely to have high power to detect association if the genomic 
#' feature harbours genetic markers with large effects.

#' 
#' @param p is a vector of single marker p-values (e.g. based on a t-stat)
#' @param sets is a list of marker sets  - names corresponds to rownames in stat
#' @param threshold is the single marker p-value cut-off
#' @return Returns a vector of p values with length equal to the number of sets 
#' @author Peter Sørensen
#' @references Gene-based Association Approach Identify Genes Across Stress Traits in Fruit Flies. In: Proceedings, 10th World Congress of Genetics Applied to Livestock Production (WCGALP), Vancouver, Canada, 2014. Dunn Rohde P, Edwards SM, Sarup PM, Sørensen P. 
#' @examples
#' 
#' # Simulate data
#' W <- matrix(rnorm(2000000),ncol=10000)
#'   colnames(W) <- as.character(1:ncol(W))
#' y <- rowSums(W[,1:10]) + rowSums(W[,1001:1010]) + rnorm(nrow(W))
#' 
#' # REML analyses 
#' data <- data.frame(y=y,mu=1)
#' fm <- y ~ mu
#' fit <- gfm(fm=fm,W=W,sets=list(colnames(W)),data=data)
#' 
#' # hyperG set test 
#' sets <- list(A=as.character(1:100),B=as.character(101:1000),C=as.character(1001:5000),D=as.character(5001:10000))
#' p <- pt(fit$s/fit$vs)
#' res <- hgTest(p=p,sets=sets,threshold=0.05)
#' 
#' @export
hgTest <- function(p=NULL,sets=NULL,threshold=0.05) {
     N <- length(p)
     Na <- sum(p<threshold)
     Nna <- N-Na
     Nf <- sapply(sets,length)
     Naf <- sapply(sets, function(x) { sum(p[x]<threshold) })
     Nnaf <- Nf-Naf
     Nanf <- Na-Naf
     Nnanf <- Nna-Nnaf
     phyperg <- 1-phyper(Naf-1, Nf, N-Nf, Na)
     phyperg
}

####################################################################################################################
#' 
#' covSets 
#' 
#'
#' @description
#' Partitioning of covariance using genetic marker sets.
#' 
#' @export
covSets <- function(W=NULL,g=NULL,sets=NULL, level2=FALSE) {
     
     W <- W[,unlist(sets)]
     nsets <- sapply(sets,length)
     sets <- sets[nsets>0]
     nsets <- nsets[nsets>0]
     
     WWi <- ginv(W%*%t(W),tol=0.0001)$G
     b <- as.vector(crossprod(W,WWi%*%g))
     names(b) <- colnames(W)
     Wb <-   t(t(W)*as.vector(b))
     gWb <- colSums(as.vector(g)*Wb)
     
     cv <- cvSet <- NULL
     
     ##########################################################
     # level 1 
     ##########################################################
     
     f <- sapply(sets,function(x){ rowSums(as.matrix(Wb[,x]))})
     cvs <- apply(f,2,function(x){x%*%g})
     Vf <- var(f)
     vf <- diag(Vf)
     covf <- rowSums(Vf) - vf
     hvf <- vf/var(g)
     hcovf <- (vf+covf)/var(g)
     cv <- cbind( nsets, vf,hvf,covf,hcovf, cvs,
                  vf/nsets,hvf/nsets,covf/nsets,hcovf/nsets,cvs/nsets)
     colnames(cv) <- c("nset","varf","hf","covf","hcovf","cvs","varf pr snp","hf pr snp","covf pr snp","hcovf pr snp","cvs pr snp")
     
     ##########################################################
     # level 2 
     ##########################################################
     if (level2) {
          cvSet <- lapply(sets,function(x){
               Wf <- as.matrix(Wb[,x])
               gf <- rowSums(Wf)
               Vf <- var(Wf)
               vf <- diag(Vf)
               covf <- rowSums(Vf) - vf
               hvf <- vf/var(gf)
               hcovf <- (vf+covf)/var(gf)
               cvs <- apply(Wf,2,function(x){x%*%gf})
               cf <- cbind(vf,hvf,covf,hcovf,cvs)
               colnames(cf) <- c("varf","hf","covf","hcovf","cvs")
               cf
          })
     }
     list(cvf=cv,cvfSet=cvSet, f=f)
}



####################################################################################################################
#    Module 3: LMM, Bayesian
####################################################################################################################

#' 
#' Genomic Feature Model analyses implemented using Bayesian Methods.
#'
#' @description
#' Bayesian Genomic Feature models implemented using Bayesian Methods. 
#' Multiple features and multiple traits models can be fitted.
#' Different genetic models (e.g. additive, dominance, gene by gene and gene by environment interactions) can be specified.
#' 
#' @details
#' The models are implemented using empirical Bayesian methods.The hyperparameters of the dispersion parameters of the Bayesian model can
#' be obtained from prior information or estimated by maximum likelihood, and conditional on these, the model is fitted using
#' Markov chain Monte Carlo. Furthermore a spectral decomposition of genomic feature relationship matrices plays an 
#' important computational role in the Markov chain Monte Carlo strategy implemented,

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
#' @references Genetic control of environmental variation of two quantitative traits of Drosophila melanogaster revealed by whole-genome sequencing. Genetics (2015) 201(2):487-97. Sørensen P, de los Campos G, Morgante F, Mackay TFC and Sorensen D.
#' @export
bgfm <- function(y=NULL, g=NULL, nsamp=50, nburn=10, nsave=10000, tol=0.001) {
  # nsamp is the number of samples
  y <- as.matrix(y)
  n <- nrow(y)                            # number of observation
  nt <- ncol(y)                           # number of traits
  mu <- colMeans(y)
  samples <- as.integer(seq(nburn+1,nsamp,length.out=min(nsave,nsamp)))   
  
  psum <- matrix(0,nrow=n,ncol=nt)
  
  nset <- length(g$G)                       # number of sets
  for ( i in 1:nt) {                        # genetic values (n*nset)
    g$values[[i]] <- matrix(0,nrow=n,ncol=nset)         
  }
  if(nt>1) sigmae <- g$sigma[[nset+1]][1,1]
  if(nt==1) sigmae <- g$sigma[[nset+1]]
  sigma <- g$sigma
  mus <- matrix(NA,nrow=nsamp,ncol=nt)
  
  alphas <- vector(length=nt,mode="list")
  
  U <- D <- sigmas <- vector(length=nset,mode="list")
  for ( i in 1:nset ){                    
    e <- eigen(g$G[[i]])                    # eigen value decomposition of the matrix G
    ev <- e$values
    U[[i]] <- e$vectors[,ev>tol]            # keep eigen vector if ev>tol
    D[[i]] <- e$values[ev>tol]              # keep eigen value if ev>tol
    sigmas[[i]] <- matrix(NA,nrow=nsamp,ncol=(nt*(nt+1))/2)
    for (j in 1:nt) {
      as <-  matrix(NA,nrow=sum(ev>tol),ncol=length(samples))
      colnames(as) <- as.character(samples)
      alphas[[j]][[i]] <- as
    }
  }
  
  for ( i in 1:nsamp) {                     # loop over the number of samples (nsamp)
    for ( j in 1:nset ) {                   # loop over the number of sets
      rhs <- NULL
      for (t in 1:nt) {
        yadj <- y[,t]-mu[t]-rowSums(as.matrix(g$values[[t]][,-j]))
        rhst <- crossprod( U[[j]],yadj ) 
        rhs <- cbind(rhs,rhst)
      }
      Vi <- solve(sigma[[j]])
      a <- matrix(0,nrow=nrow(rhs),ncol=nt)
      for (k in 1:nrow(rhs)) {
        iC <- solve( diag(1,nt) + (sigmae/D[[j]][k])*Vi )      
        ahat <- iC%*%rhs[k,]
        a[k,] <- mvrnorm(n=1,mu=ahat,Sigma=iC*sigmae) 
      }
      for (t in 1:nt) { g$values[[t]][,j] <- U[[j]]%*%a[,t] }
      df <- nrow(a) + g$df[j]
      
      if (any(i%in%samples)) {
        for (t in 1:nt) { alphas[[t]][[j]][,as.character(i)] <- a[,t] }
      }
      
      # inverse chisquare
      if (nt==1) {
        #scg <- sum((1/D[[j]])*a**2) + g$sigma[[j]]*g$df[j]         
        scg <- sum((1/D[[j]])*a**2) + (g$sigma[[j]]*(g$df[j]+2))/g$df[j]	# => S = (mode*(df+2))/df         
        sigma[j] <- scg/rchisq(n=1, df=df, ncp=0)    
        sigmas[[j]][i,] <- sigma[j]
      }
      
      # inverse wishart
      if (nt>1) {
        if (g$ds[j]=="invWish") {
          S <- t(a*(1/D[[j]]))%*%a + g$sigma[[j]]*(g$df[j]+nt+1)		# => S = mode*(df+nt+1)
          S <- riwish(df,S)
          sigma[[j]] <- S
          sigmas[[j]][i,] <- S[as.vector(lower.tri(S, diag=TRUE))]
          #print(c(j,g$ds[j]))
        }
        if (g$ds[j]=="invChisq") {
          scg <- t(a*(1/D[[j]]))%*%a + (g$sigma[[j]]*(g$df[j]+2))/g$df[j]	# => S = (mode*(df+2))/df
          for (t in 1:nt) {
            sigma[[j]][t,t] <- scg[t,t]/rchisq(n=1, df=df, ncp=0)    
          }
          sigmas[[j]][i,] <- sigma[[j]][as.vector(lower.tri(sigma[[j]], diag=TRUE))]
        }
        
      }
    }
    for (t in 1:nt) {
      yadj <- y[,t]-rowSums(as.matrix(g$values[[t]])) 
      rhs <- sum(yadj)
      lhs <- (n+sigmae/100000)
      mu[t] <- rnorm(1,mean=rhs/lhs,sd=1/sqrt(lhs))
    }  
    mus[i,] <- mu
    print(i)
    for (t in 1:nt) {
      e <- y[,t]-mu[t]-rowSums(as.matrix(g$values[[t]]))
      p <-  (1/sqrt(2*pi))*exp(-0.5*(e**2)) 
      if(i>nburn) psum[,t] <- psum[,t] + 1/p
    }
    
  }
  # summary of model fit
  logCPO <- NULL
  for (t in 1:nt) {
    logCPO[t] <- sum(log((nsamp-nburn)*(1/psum[,t])))
  }
  return(list(sigmas=sigmas,mus=mus, logCPO=logCPO, g=g, alphas=alphas, nsamp=nsamp, nburn=nburn) ) # return posterior samples of sigma 
}


# #' 
# #' Genomic feature model analyses using Bayesian Mixture Models
# #'
# #' @description
# #' Genomic Feature Linear Mixed Models for predicting quantitative trait phenotypes from high 
# #' resolution genomic polymorphism data. Genomic features are regions on the genome that are 
# #' hypothesized to be enriched for causal variants affecting the trait. Several genomic feature 
# #' classes can be formed based on previous studies and different sources of information including 
# #' genes, chromosomes, biological pathways, gene ontologies, sequence annotation, prior QTL regions, 
# #' or other types of external evidence. This is important because prediction is difficult for 
# #' populations of unrelated individuals when the number of causal variants is low relative to the 
# #' total number of polymorphisms, and causal variants individually have small effects on the traits. 
# #' The models can implemented using Bayesian Methods. Single trait and multiple random effects. 
# #' 
# #' hssvs
# #'
# #' @param y is a vector of phenotypes 
# #' @param X is a matrix of centered and scaled genotypes
# #' @param set is a list of marker sets  - names corresponds to rownames in stat
# #' @param hgprior is a list of prior scale parameters (sce0 and scg0), prior degrees of freedom (dfe0 and dgf0),  prior distribution (not yet implemented)
# #' @param nsamp is the number of samples after burnin
# #' @param p1 prior proportion of elements in gamma to be set to one
# #' @param g0 spike variance
# #' @return Returns a list including 
# #' \item{sigmas}{posterior samples variance components} 
# #' \item{mus}{posterior samples mu}
# #' \item{logCPO}{log CPOs}
# #' \item{g}{posterior samples g}
# #' \item{alphas}{posterior samples alphas}
# #' \item{mus}{number of markers in the set }
# #' @author Peter Sørensen
# #' @references Genetic control of environmental variation of two quantitative traits of Drosophila melanogaster revealed by whole-genome sequencing. Genetics (2015) 201(2):487-97. Sørensen P, de los Campos G, Morgante F, Mackay TFC and Sorensen D.
# #' @export
# hssvs <- function( y=NULL, X=NULL, set=NULL, p1=NULL, g0=NULL, nsamp=100, hgprior=list(sce0=0.001, scg0=0.001, dfe0=4, dfg0=4)) {
#      
#      n <- length(y)                           # number of observations
#      p <- ncol(X)                             # number of regression variables
#      dxx <- colSums(X**2)                     # diagonal elements of the X'X matrix
#      p0 <- 1-p1                               # prior proportion of elements in gamma to be set to one
#      
#      b <- rep(0,p)                            # initialize b
#      g <- rep(0,p)                            # intialize g
#      
#      nset <- length(set)                      # number of sets
#      pset <- sapply(set,length)
#      for ( i in 1:nset) {
#           n1 <- as.integer(pset[i]*p1) 
#           cset <- abs(cor(y, X[,set[[i]]]))
#           o <- order(cset,decreasing=TRUE)[1:n1]
#           o <- set[[i]][o]
#           g[o] <- 1
#      }
#      
#      bset <- NULL                             # beta's within sets
#      gset <- NULL                             # g's within sets
#      for ( i in 1:nset ) {
#           bset[[i]] <- b[set[[i]]]
#           gset[[i]] <- g[set[[i]]]
#      }
#      
#      mu <- 0                                  # initilaize mu
#      sigma2 <- var(y)/2                       # initialize sigma2
#      g1 <- 1                                  # initialize slab variance
#      if (is.null(g0)) g0 <- 1e-09             # initialize spike variance
#      
#      g0 <- rep(g0,nset)                       # spike variance for each set
#      g1 <- rep(g1,nset)                       # slab variance for each set
#      
#      #Xs <- matrix(0,nrow=n,ncol=nset)         # matrix (n*nset) of genomic values one for each set
#      
#      e <- as.vector(y)                        # initialize residuals
#      
#      for ( i in 1:nsamp ) {
#           
#           for ( j in 1:nset) {
#                cls <- set[[j]]
#                samp <- ssvs(e=e, X=X[,cls], b=bset[[j]], dxx=dxx[cls], mu=mu, g=gset[[j]], sigma2=sigma2, p0=p0, p1=p1, g0=g0[j], g1=g1[j], hgprior=hgprior)
#                e <- samp$e
#                bset[[j]] <- samp$b
#                mu <- samp$mu
#                sigma2 <- samp$sigma2
#                gset[[j]] <- samp$g
#                g0[j] <- samp$g0
#                g1[j] <- samp$g1
#                # Xs[,j] <- X[,cls]%*%bset[[j]]
#           }
#           #  h2 <- apply(Xs,2,var)
#           #  h2 <- h2/var(y)
#           #  barplot(h2)
#           print("")
#           print(c("round",i, "gset", sapply(gset,sum)))
#           #  print(c("round",i, "mu, sigma2, cor(y,Xs)", c(mu,sigma2,cor(y,rowSums(Xs)))))
#           #  print(c("round",i, "g0",g0 ))
#           #  print(c("round",i, "g1",g1 ))
#      }
#      
#      
# }
# 
# 
# ssvs <- function(e=e, X=X, b=b, dxx=dxx, mu=mu, g=g, sigma2=sigma2, p0=p0, p1=p1, g0=g0, g1=g1, hgprior=hgprior) {
#      
#      n <- length(e)                           # number of observations
#      p <- ncol(X)                             # number of regression variables
#      
#      sce0 <- hgprior$sce0                     # prior residual sums of squares
#      scg0 <- hgprior$scg0                     # prior slab residual sums of squares
#      dfe0 <- hgprior$dfe0                     # prior residual degrees of freedom
#      dfg0 <- hgprior$dfg0                     # prior slab residual degrees of freedom
#      
#      
#      # sample b and update residuals 
#      bC <- b
#      shrinkage <- rep(sigma2/g1,times=p)
#      shrinkage[g==0] <- sigma2/g0 
#      xxs <- dxx + shrinkage
#      for ( j in 1:p){
#           rhs <- sum(X[,j]*e) + dxx[j]*bC[j]
#           b[j] <- rnorm(1,mean=rhs/xxs[j],sd=sigma2/xxs[j])
#           e  <- e - X[,j]*(b[j] - bC[j]) 
#      }
#      
#      # sample g0 
#      scg00 <- sum((b**2)[g==0]) + scg0
#      dfg00 <- sum(1-g) + dfg0 - 2 
#      g0 <- scg00/rchisq(n=1, df=dfg00, ncp=0) 
#      
#      # sample g1 
#      scg1 <- sum((b**2)[g==1]) + scg0
#      dfg1 <- sum(g) + dfg0 - 2 
#      for (k in 1:1000) {                      # scaling needs to be revised
#           g1 <- scg1/rchisq(n=1, df=dfg1, ncp=0) 
#           if(g1>g0) break                         # scaling needs to be revised
#           if(k==1000) g1 <- g0                    # scaling needs to be revised
#      }                                        # scaling needs to be revised
#      
#      
#      # sample g 
#      ratio1 <- p1*(exp(-(b**2)/g1)/sqrt(g1))
#      ratio0 <- p0*(exp(-(b**2)/g0)/sqrt(g0))
#      ratio <- ratio1/ratio0
#      u <- 1/(1-runif(p,0,1))
#      g[1:p] <- 0
#      g[ratio>u] <- 1
#      
#      
#      # sample mu and update residuals 
#      muC <- mu
#      mu <- rnorm(1,mean=mean(e+muC),sd=sigma2/n)
#      e <- e - mu + muC 
#      
#      # sample sigma 
#      sce <- sum(e**2) + sce0
#      dfe <- n + dfe0 - 2 
#      sigma2 <- sce/rchisq(n=1, df=dfe, ncp=0) 
#      
#      
#      return(list(e=e, mu=mu,b=b,sigma2=sigma2,g=g, g0=g0, g1=g1))
#      
# }
# 
# ################################################################################
# # Simulate data and test SSVS function
# ################################################################################
# 
# 
# # # simulated data set
# # X <- matrix(rnorm(10000000),nrow=1000)
# # set1 <- sample(1:ncol(X),5)
# # set2 <- sample(1:ncol(X),5)
# # set1a <- unique(c(set1,sample(1:ncol(X),100)))
# # set2a <- unique(c(set2,sample(1:ncol(X),100)))
# # set3 <- (1:ncol(X))[-unique(c(set1a,set2a))]
# # 
# # y <- rowSums(X[,set1]) + rowSums(X[,set2]) + rnorm(nrow(X),mean=0,sd=1)       
# # 
# # varE <- var( y - rowSums(X[,set1]) - rowSums(X[,set2]) ) 
# # varE
# # 
# # h2 <- ( var(rowSums(X[,set1])) + var(rowSums(X[,set2])) )/var(y)
# # h2
# # 
# # set <- list(set1=set1a,set2=set2a,set3=set3)
# # 
# # res <- hssvs(y=y,X=X,set=set, p1=0.01,g0=0.0000001, nsamp=100)
# 
# 
# 
# 
# ################################################################################
# # Multi Class Bayesian Regression function 
# ################################################################################
# 
# 
# mcbr <- function(y=NULL, X=NULL, nc=NULL, l1=NULL, l2=NULL, phi=NULL, nsamp=100 ) {
#      
#      n <- length(y)             # number of observations
#      p <- ncol(X)               # number of regression variables
#      
#      dxx <- colSums(X**2)       # diagonal elements of the X'X matrix
#      
#      b <- as.numeric(rep(0,p))  # initialize b
#      mu <- 0                    # initilaize mu
#      
#      if (is.null(l1))  l1 <- 1/(10**((1:nc)-4))                       # prior shape parameter lambda
#      if (is.null(l2))  l2 <- rep(1/10**2,nc)                          # prior rate parameter lambda
#      if (is.null(phi)) phi0 <- c(0.5,0.25,0.15,0.05,0.025,0.01,0.001)  # prior phi
#      phi <- phi0
#      logphi <- log(phi)                                               # prior phi on the log scale
#      
#      lambdak <- rep(1000,nc)                             # initialize lambdak one for each class
#      lambda <- rep(1000,p)                               # initialize lambda one for each regressor
#      
#      g <- rep(1,p)                                       # initialize class indicator variable
#      cors <- abs(cor(y, X))
#      quants <- quantile(cors, probs=1-phi)
#      for ( i in 1:nc) {
#           #  g[cors>quants[i]] <- i 
#      }
#      
#      e <- as.vector(y - mu - X%*%b)     # initialize residuals
#      
#      sigma2 <- var(e)/2                 # initialize sigma2
#      
#      
#      for ( i in 1:nsamp ) {
#           
#           # sample mu and update residuals 
#           muC <- mu
#           mu <- rnorm(1,mean=mean(e+muC),sd=sigma2/n)
#           e <- as.vector(e - mu + muC)
#           
#           # sample b and update residuals                      
#           bC <- b
#           xxs <- as.numeric(dxx + lambda*sigma2)                    
#           for ( j in 1:p){
#                rhs <- as.numeric(sum(X[,j]*e) + dxx[j]*bC[j])
#                b[j] <- rnorm(1,mean=rhs/xxs[j],sd=sigma2/(xxs[j]))
#                e  <- as.numeric(e - X[,j]*(b[j] - bC[j])) 
#           }
#           e <- as.numeric(y - mu - X%*%b)                     # update residuals
#           
#           # sample lambda
#           b2 <- b**2
#           l1k <- l1
#           l2k <- l2
#           for ( k in 1:nc ) {
#                if(sum(g==k)>0) l1k[k] <- l1k[k] + sum(g==k)/2
#                if(sum(g==k)>0) l2k[k] <- l2k[k] + sum(b2[g==k])/2
#                lambdak[k] <- rgamma( n=1, shape=l1k[k], rate=l2k[k])      
#                if(sum(g==k)>0) lambda[g==k] <- lambdak[k]
#           }
#           loglambdak <- log(lambdak)
#           
#           # sample gamma
#           for ( j in 1:p ) {
#                probs <- -0.5*b2[j]*lambdak + logphi + 0.5*loglambdak
#                g[j] <- (1:nc)[rmultinom(1,1,prob=exp(probs))==1]
#           }
#           
#           # sample phi
#           #phi <- rep(0.0001,nc)
#           rws <- as.numeric(names(table(g)))
#           #phi[rws] <- table(g)/p
#           print(c(i,phi))
#           
#           # sample sigma                                       
#           she <- n                                       
#           sce <- sum(e**2)
#           sigma2 <- sce/rchisq(n=1,df=she)                   
#           
#           plot(b)
#           print(table(g))
#           print(c(mu,sigma2))
#           
#      }
#      
# }
# 
# 
# # Simulate data and test function
# #X <- matrix(rnorm(10000000),nrow=1000)
# #set <- sample(1:ncol(X),10)
# #y <- rowSums(X[,set]) + rnorm(nrow(X),mean=0,sd=1)       
# #h2 <- var(rowSums(X[,set]))/var(y)
# #
# #mcbr(y=y, X=X, nc=7)
# 
# 
# 
# 
# bgL <- function(y=NULL, X=NULL, g=NULL, nsamp=100, hgprior=list(sce0=0.001, scg0=0.001, dfe0=4, dfg0=4)) {
#      n <- length(y)           # number of observations
#      nb <- ncol(X)            # number of regression variables
#      ng <- length(unique(g))  # number of groups
#      dxx <- colSums(X**2)     # diagonal elements of the X'X matrix
#      
#      b <- rep(0,nb)                           # initialize b
#      mu <- 0                                  # initilaize mu
#      sigma <- 1                               # initialize sigma
#      g1 <- rep(1/nb,ng)                     # initialize slab variance
#      g1 <- rep(0,ng)                     # initialize slab variance
#      for ( j in 1:ng ){
#           g1[j] <- 1/sum(g==j)
#      }
#      print(g1)
#      #g1 <- c(1,0.000001)                     # initialize slab variance
#      
#      sce0 <- hgprior$sce0                     # prior residual sums of squares
#      scg0 <- hgprior$scg0                     # prior slab residual sums of squares
#      dfe0 <- hgprior$dfe0                     # prior residual degrees of freedom
#      dfg0 <- hgprior$dfg0                     # prior slab residual degrees of freedom
#      
#      e <- as.vector(y - mu - X%*%b)           # initialize residuals
#      
#      for ( i in 1:nsamp ) {
#           
#           # sample mu and update residuals 
#           mu0 <- mu
#           mu <- rnorm(1,mean=mean(e+mu0),sd=sigma/n)
#           e <- e - mu + mu0 
#           
#           # sample b and update residuals 
#           b0 <- b
#           shrinkage <- rep(0,times=nb)
#           for ( j in 1:ng ){
#                shrinkage[g==j] <- sigma/g1[j]
#           }
#           xxs <- dxx + shrinkage
#           for ( j in 1:nb){
#                rhs <- sum(X[,j]*e) + dxx[j]*b0[j]
#                b[j] <- rnorm(1,mean=rhs/xxs[j],sd=sigma/xxs[j])
#                e  <- e - X[,j]*(b[j] - b0[j]) 
#           }
#           
#           # sample sigma 
#           sce <- sum(e**2) + sce0
#           dfe <- n + dfe0 - 2 
#           sigma <- sce/rchisq(n=1, df=dfe, ncp=0) 
#           #sigma <- sum(e**2)/rchisq(n=1, df=n-2, ncp=0) 
#           
#           
#           # sample g1 
#           for ( j in 1:ng){
#                scg1 <- sum((b**2)[g==j]) + scg0
#                dfg1 <- sum(g==j) + dfg0 - 2 
#                g1[j] <- scg1/rchisq(n=1, df=dfg1, ncp=0) 
#           }
#           
#           # sample g 
#           
#           print("")
#           print(c("round",i, (1:nb)[g==1]))
#           print(c("round",i, c(mu,sigma,g1)))
#           print(c("round",i, b[g==1]))
#      }
# }
# 
# 
# 
# 
# ################################################################################
# # Simulate data and test SSVS functions
# ################################################################################
# 
# #library(BLR)
# #data(wheat)
# #
# #res1 <- ssvs(y=Y[,1],X=X,p1=0.01,g0=0.0000001, nsamp=100)     
# #res2 <- ssvs(y=Y[,2],X=X,p1=0.01,g0=0.0000001, nsamp=100)     
# #res4 <- ssvs(y=Y[,3],X=X,p1=0.01,g0=0.0000001, nsamp=100)     
# #res5 <- ssvs(y=Y[,4],X=X,p1=0.01,g0=0.0000001, nsamp=100)     
# 
# 
# # # simulated data set
# # X <- matrix(rnorm(10000000),nrow=1000)
# # set1 <- sample(1:ncol(X),50)
# # set2 <- sample(1:ncol(X),50)
# # set3 <- sample(1:ncol(X),50)
# # set4 <- sample(1:ncol(X),50)
# # set5 <- sample(1:ncol(X),50)
# # g <- rep(6,times=ncol(X))
# # g[set1] <- 1
# # g[set2] <- 2
# # g[set3] <- 3
# # g[set4] <- 4
# # g[set5] <- 5
# # 
# # set <- unique(c(set1,set2,set3,set4,set5))
# # 
# # y <- rowSums(X[,set]) + rnorm(nrow(X),mean=0,sd=20)       
# # 
# # system.time( res <- bgL(y=y,X=X,g=g, nsamp=10) ) 
# # 
# # 
