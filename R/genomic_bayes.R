####################################################################################################################
#    Module 6: Bayesian models
####################################################################################################################
#'
#' Genomic prediction models implemented using Bayesian Methods (small data)
#'
#' @description
#' Genomic prediction models implemented using Bayesian Methods (small data).
#' The models are implemented using empirical Bayesian methods. The hyperparameters of the dispersion parameters of the Bayesian model can
#' be obtained from prior information or estimated by maximum likelihood, and conditional on these, the model is fitted using
#' Markov chain Monte Carlo. These functions are currently under development and future release will be able to handle large data sets.
#'
#'
#' @param y is a matrix of phenotypes
#' @param W is a matrix of centered and scaled genotypes
#' @param nsamp is the number of samples after burnin
#' @param sets is a list of markers defining a group
#' @param nsets is a list of number of marker groups
#' @param phi is the proportion of markers in each marker variance class (phi=c(0.999,0.001),used if method="ssvs")
#' @param h2 is the trait heritability
#' @param method specifies the methods used (method="mixed","bayesC",blasso","blr")
#' @param nburn is the number of burnin samples
#' @param nsave is the number of samples to save
#' @param tol is the tolerance
#'

#' @author Peter Sørensen


#' @examples
#'
#'
#' # Simulate data and test functions
#'
#' W <- matrix(rnorm(100000),nrow=1000)
#' set1 <- sample(1:ncol(W),5)
#' set2 <- sample(1:ncol(W),5)
#' sets <- list(set1,set2)
#' g <- rowSums(W[,c(set1,set2)])
#' e <- rnorm(nrow(W),mean=0,sd=1)
#' y <- g + e
#'
#'
#' fitM <- gbayes(y=y, W=W, method="mixed")
#' fitA <- gbayes(y=y, W=W, method="bayesA")
#' fitL <- gbayes(y=y, W=W, method="blasso")
#' fitC <- gbayes(y=y, W=W, method="bayesC")
#'
#' fitM <- gbayes(y=y, W=W, method="mixed", algorithm="fastbayes")
#' fitA <- gbayes(y=y, W=W, method="bayesA", algorithm="fastbayes")
#' fitL <- gbayes(y=y, W=W, method="blasso", algorithm="fastbayes")
#' fitC <- gbayes(y=y, W=W, method="bayesC", algorithm="fastbayes")


#'
#' @export
#'

gbayes <- function(y=NULL, X=NULL, W=NULL, b=NULL, badj=NULL, seb=NULL, LD=NULL, n=NULL,
                   vara=NULL, varb=NULL, vare=NULL, lambda=NULL, scaleY=TRUE,
                   h2=NULL, pi=NULL, updateB=FALSE, updateE=FALSE, updatePi=FALSE, models=NULL,
                   nub=4, nue=4, nit=100, method="mixed", algorithm="default") {
     
     method <- match(method, c("blup","mixed","bayesA","blasso","bayesC")) - 1
     if( !sum(method%in%c(0:4))== 1 ) stop("Method specified not valid") 
     
     nt <- 1
     if(is.list(y)) nt <- length(y)
     
     if(nt==1 && !algorithm=="sbayes") {
          
          if(scaleY) y <- as.vector(scale(y)) 
          
          if(is.null(pi)) pi <- 0.01
          
          if(is.null(h2)) h2 <- 0.5
          
          n <- nrow(W)
          m <- ncol(W)
          
          if(is.null(b)) b <- rep(0,m)
          e=y-mean(y)
          if(is.null(vare)) vare <- var(e)
          if(is.null(varb)) varb <- (vare/m)*h2
          if(is.null(lambda)) lambda <- rep(vare/varb,m)
          if(is.null(vara)) vara <- vare*h2
          
          if(algorithm=="default") {
               fit <- .Call("_qgg_bayes",
                            y=y, 
                            W=split(W, rep(1:ncol(W), each = nrow(W))), 
                            b=b,
                            lambda = lambda,
                            pi = pi,
                            vara = vara,
                            varb = varb,
                            vare = vare,
                            nub=nub,
                            nue=nue,
                            updateB = updateB,
                            updateE = updateE,
                            updatePi = updatePi,
                            nit=nit,
                            method=as.integer(method)) 
               names(fit[[1]]) <- colnames(W)
               names(fit) <- c("b","p","mu","B","E","Pi","g","e")
          } 
          

          if(algorithm=="fastbayes") {
               LD <- crossprod(W)
               
               fit <- .Call("_qgg_fbayes",
                              y=y, 
                             W=split(W, rep(1:ncol(W), each = nrow(W))), 
                             LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
                             b=b,
                             lambda = lambda,
                             pi = pi,
                             vara = vara,
                             varb = varb,
                             vare = vare,
                             nub=nub,
                             nue=nue,
                             updateB = updateB,
                             updateE = updateE,
                             updatePi = updatePi,
                             nit=nit,
                             method=as.integer(method)) 
               names(fit[[1]]) <- colnames(W)
               names(fit) <- c("b","p","mu","B","E","Pi")
          }
          
          
          
          
     } 
     
     if(algorithm=="sbayes") {
          
          if(!is.null(W) && is.null(LD)) {
               n <- nrow(W)
               LD <- crossprod(W)/(n-1)
          }    
          
          m <- ncol(LD)
          
          if(is.null(pi)) pi <- 0.01
          
          if(is.null(h2)) h2 <- 0.5
          
          if(is.null(vare)) vare <- 1
          if(is.null(varb)) varb <- (vare/m)*h2
          if(is.null(lambda)) lambda <- rep(vare/varb,m)
          if(is.null(vara)) vara <- vare*h2
          
          wy <- b*n
          b <- rep(0,m)
          if(!is.null(badj)) b <- badj 
          
          fit <- sbayes(wy=wy, 
                        LD=split(LD, rep(1:ncol(LD), each = nrow(LD))), 
                        b = b,
                        lambda = lambda,
                        pi = pi,
                        vara = vara,
                        varb = varb,
                        vare = vare,
                        nub=nub,
                        nue=nue,
                        updateB = updateB,
                        updateE = updateE,
                        updatePi = updatePi,
                        n=n,
                        nit=nit,
                        method=as.integer(method))
          names(fit[[1]]) <- rownames(LD)
          names(fit) <- c("b","p","mu","B","E","Pi")
     }
     
     
     if(nt>1) {
          
          if(method==2) stop("Multiple trait not yet implemented for Bayes A") 
          if(method==3) stop("Multiple trait not yet implemented for Bayesian Lasso")
          
          n <- nrow(W)
          m <- ncol(W)
          
          if(scaleY) y <- lapply(y,function(x){as.vector(scale(x))})
          
          
          if(is.null(b)) b <- lapply(1:nt,function(x){rep(0,m)}) 
          
          if(is.null(models)) {
               models <- rep(list(0:1), nt)
               models <- t(do.call(expand.grid, models))
               models <- split(models, rep(1:ncol(models), each = nrow(models)))
          } 
          if(is.character(models)) {
               if(models=="restrictive") {
                    models <- list(rep(0,nt),rep(1,nt))
                    pi <- c(0.999,0.001)
               }
          }
          
          if(is.null(pi)) {
               pi <- c(0.99,rep(0.01,length(models)-1)) 
          }
          
          if(is.null(h2)) h2 <- 0.005
          if(is.null(vare)) {
               vare <- diag(sapply(y,var))
          }
          if(is.null(varb)) varb <- diag(sapply(y,var)/(m*pi[length(models)]))*h2
          
          
          fit <- mtbayes(y=y, 
                         W=split(W, rep(1:ncol(W), each = nrow(W))), 
                         b=b,
                         B = varb,
                         E = vare,
                         models=models,
                         pi=pi,
                         nub=nub,
                         nue=nue,
                         updateB=updateB,
                         updateE=updateE,
                         nit=nit,
                         method=as.integer(method))
          
     }
     
     return(fit)
     
}

#' Plot fit from gbayes
#'
#' @description
#' Summary plots from gbayes fit

#' @param fit object from gbayes

#'
#' @export
#'

plotBayes <- function(fit=NULL, causal=NULL) {
     if(!is.list(fit[[1]])) {
          layout(matrix(1:6,nrow=3,ncol=2))
          plot(fit[[1]],ylab="Posterior", xlab="Marker", main="Marker effect", frame.plot=FALSE)  
          if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch=4, cex=2, lwd=3 )
          plot(fit[[2]],ylab="Posterior mean", xlab="Marker", main="Marker indicator", frame.plot=FALSE)  
          if(!is.null(causal)) points(x=causal,y=rep(0,length(causal)),col="red", pch=4, cex=2, lwd=3 )
          hist(fit[[6]],xlab="Posterior", main="Pi")  
          hist(fit[[4]],xlab="Posterior", main="Marker variance")  
          hist(fit[[5]],xlab="Posterior", main="Residual variance")  
          plot(fit[[6]],xlab="Sample", ylab="Pi", frame.plot=FALSE)  
     } 
     if(is.list(fit[[1]])) {
          layout(matrix(1:4,2,2))
          matplot(as.data.frame(fit[[1]]),ylab="Marker effect", frame.plot=FALSE)  
          matplot(as.data.frame(fit[[2]]),ylab="Marker indicator", frame.plot=FALSE)  
          matplot(as.data.frame(fit[[4]]),ylab="Marker variance", frame.plot=FALSE)  
          matplot(as.data.frame(fit[[5]]),ylab="Residual variance", frame.plot=FALSE)  
     } 
     
}

splitWithOverlap <- function(vec, seg.length, overlap) {
     starts = seq(1, length(vec), by=seg.length-overlap)
     ends   = starts + seg.length - 1
     ends[ends > length(vec)] = length(vec)
     lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

neff <- function(seb=NULL,af=NULL,Vy=1) {
     seb2 <- seb**2
     vaf <- 2*af*(1-af)
     neff <- round(median(Vy/(vaf*seb2)))
     return(neff)
}

