####################################################################################################################
#    Module 5: GSOLVE 
####################################################################################################################

#' @export

gsolve <- function( y=NULL, X=NULL, W=NULL, sets=NULL, msets=100, lambda=NULL, validate=NULL, weights=FALSE, method="gsru", maxit=500, tol=0.0000001) { 
     if(is.null(validate)) {
        if(method=="gsru")  fit <- gsru(y=y, W=W, X=X, sets=sets, lambda=lambda, weights=weights, maxit=maxit, tol=tol)
        if(method=="gsqr")  fit <- gsqr(y=y, W=W, X=X, sets=sets,msets=msets,lambda=lambda,weights=weights, maxit=maxit, tol=tol)
     }
     if(!is.null(validate)) { 
          n <- length(y)     
          pa <- mspe <- intercept <- slope <- r2 <- NULL
          for ( k in 1:ncol(validate)) {
               v <- validate[, k]
               t <- (1:n)[-v]
               if(method=="gsru")  fit <- gsru(y=y[t], X=as.matrix(X[t,]), W=W[t,], sets=sets, lambda=lambda, weights=weights, maxit=maxit, tol=tol)
               if(method=="gsqr")  fit <- gsqr(y=y[t], X=as.matrix(X[t,]), W=W[t,], sets=sets, msets=msets, lambda=lambda, weights=weights, maxit=maxit, tol=tol)
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


#' @export

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

#' @export

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
          #fit$s[rws] <- solve(QRlist$R[[i]])%*%fit$s[rws,1]
          fit$s[rws] <- backsolve(QRlist$R[[i]],fit$s[rws,1])
     }
     return(fit)
}  

