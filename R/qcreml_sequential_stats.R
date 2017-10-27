
#######################################################################
# seq.sum(x)
#######################################################################
# obtains a sequential sum for a sequence of values
#######################################################################
#---------------------------------------------------------------------#
# ARGUMENTS
#---------------------------------------------------------------------#
# x            vector containing a sequence: numeric
#---------------------------------------------------------------------#
# VALUE
#---------------------------------------------------------------------#
# sequential sum
#######################################################################

seq.sum <- function(x)
{
   return(rev(cumsum(rev(x))))
}

#######################################################################
# seq.mean(x)
#######################################################################
# obtains a sequential mean for a sequence of values
#######################################################################
#---------------------------------------------------------------------#
# ARGUMENTS
#---------------------------------------------------------------------#
# x            vector containing a sequence: numeric
#---------------------------------------------------------------------#
# VALUE
#---------------------------------------------------------------------#
# sequential mean
#######################################################################

seq.mean <- function(x)
{
   return(rev(cumsum(rev(x)))/seq(from=length(x),to=1))
}

#######################################################################
# seq.var(x)
#######################################################################
# obtains a sequential variance for a sequence of values
#######################################################################
#---------------------------------------------------------------------#
# ARGUMENTS
#---------------------------------------------------------------------#
# x            vector containing a sequence: numeric
#---------------------------------------------------------------------#
# VALUE
#---------------------------------------------------------------------#
# sequential variance
#######################################################################

seq.var <- function(x)
{
   mx <- seq.mean(x)
   tmp <- numeric(0)
   for(i in 1:length(x))
   {
      tmp <- c(tmp,(1/(length(x)-i+1))*sum((x[i:length(x)]-mx[i])^2))
   }
   rm(i)
   return(tmp)
}

#######################################################################
# seq.cov(x,y)
#######################################################################
# obtains a sequential covariance between two sequences of values
#######################################################################
#---------------------------------------------------------------------#
# ARGUMENTS
#---------------------------------------------------------------------#
# x,y          vectors containing the sequences: numeric
# var.return   return sequential variances of the sequences: logical
#---------------------------------------------------------------------#
# VALUE
#---------------------------------------------------------------------#
# sequential covariance
#######################################################################

seq.cov <- function(x,y,var.return=FALSE)
{
   if(length(x) != length(y))
   {
      stop("ERROR: incompatible dimensions")
   }
   mx <- seq.mean(x)
   my <- seq.mean(y)
   s2x <- seq.var(x)
   s2y <- seq.var(y)
   tmp <- numeric(0)
   for(i in 1:length(x))
   {
      tmp <- c(tmp,(1/(length(x)-i+1))*sum((x[i:length(x)]-mx[i])*(y[i:length(y)]-my[i])))
   }
   rm(i)
   return(tmp)
}

#######################################################################
# seq.cov.weighted(x,y,weight=2)
#######################################################################
# obtains a sequential covariance between two sequences of values
#######################################################################
#---------------------------------------------------------------------#
# ARGUMENTS
#---------------------------------------------------------------------#
# x,y          vectors containing the sequences: numeric
# weight       indicates which variable is the weighting variable
#              weight=1, x is the weighting variable
#              weight=2, y is the weighting variable (default)
#---------------------------------------------------------------------#
# VALUE
#---------------------------------------------------------------------#
# weighted sequential covariances
#######################################################################

seq.cov.weighted <- function(x,y,weight=2)
{
   if(!(weight == 1 | weight == 2))
   {
      stop("ERROR: 'weight' must be 1 or 2")
   }
   tmp <- seq.cov(x,y)
   if(weight == 1)
   {
      w <- x
   } else w <- y
   tmp <- tmp*seq.sum(w)/sum(w)
}

#######################################################################
# seq.cor(x,y)
#######################################################################
# obtains a sequential correlation between two sequences of values
#######################################################################
#---------------------------------------------------------------------#
# ARGUMENTS
#---------------------------------------------------------------------#
# x,y          vectors containing the sequences: numeric
#---------------------------------------------------------------------#
# VALUE
#---------------------------------------------------------------------#
# sequential correlation
#######################################################################

seq.cor <- function(x,y)
{
   if(length(x) != length(y))
   {
      stop("ERROR: incompatible dimensions")
   }
   s2x <- seq.var(x)
   s2y <- seq.var(y)
   tmp <- seq.cov(x,y)
   tmp[-length(tmp)] <- tmp[-length(tmp)]/sqrt(s2x[-length(s2x)]*s2y[-length(s2y)])
   return(tmp)
}

#######################################################################
# seq.cor.weighted(x,y,weight=2)
#######################################################################
# obtains a sequential correlation between two sequences of values
#######################################################################
#---------------------------------------------------------------------#
# ARGUMENTS
#---------------------------------------------------------------------#
# x,y          vectors containing the sequences: numeric
# weight       indicates which variable is the weighting variable
#              weight=1, x is the weighting variable
#              weight=2, y is the weighting variable (default)
#---------------------------------------------------------------------#
# VALUE
#---------------------------------------------------------------------#
# weighted sequential correlation
#######################################################################

seq.cor.weighted <- function(x,y,weight=2)
{
   if(!(weight == 1 | weight == 2))
   {
      stop("ERROR: 'weight' must be 1 or 2")
   }
   tmp <- seq.cor(x,y)
   if(weight == 1)
   {
      w <- x
   } else w <- y
   tmp <- tmp*seq.sum(w)/sum(w)
}

#####################################################################################################################
# REML.seq.eval(y,g,W=NULL,d=NULL,G=NULL,U=NULL,L=NULL,draw.plot=TRUE,save.plot=FALSE,path.plot=getwd(),pch=1,col=1)
#####################################################################################################################
# obtains the sequential statistics and plot for evaluation of h2.REML
#####################################################################################################################
#-------------------------------------------------------------------------------------------------------------------#
# ARGUMENTS
#-------------------------------------------------------------------------------------------------------------------#
# y            vector containing the phenotypes
# g            vector containing the fitted genetic values from GBLUP
# W            matrix containing the standardized SNP genotypes
# D            vector of the diagonal elements of matrix of weights for SNPs
# G            G-matrix
# U            eigen-vectors of G-matrix
# L            eigen-values of G-matrix
# draw.plot    if TRUE, a plot of Upsilon1 x zeta.y will be on the output
# save.plot    if TRUE, the prot will be saved in a file
# path.plot    the path that the plot will be saved, if save.plot=TRUE
# pch          pch to values in plot
# col          color to values in plot
#-------------------------------------------------------------------------------------------------------------------#
# VALUE
#-------------------------------------------------------------------------------------------------------------------#
# 
#####################################################################################################################

#' @export
#'

qcreml <- function(y,g,W=NULL,d=NULL,G=NULL,U=NULL,L=NULL,draw.plot=TRUE,save.plot=FALSE,path.plot=getwd(),pch=1,col=1)
{
   if(is.null(U) & is.null(L))
   {
      if(is.null(G))
      {
         if(is.null(d))
         {
            G <- tcrossprod(W)/ncol(W)
         } else G <- tcrossprod(W,tcrossprod(W,diag(d)))/sum(d)
      }
      tmp <- eigen(G)
      U <- tmp$vectors
      L <- tmp$values
      rm(tmp)
   }
   
   uy <- as.numeric(crossprod(y,U))^2
   ug <- as.numeric(crossprod(g,U))^2
   zeta.y <- seq.cor.weighted(L,uy)
   zeta.g <- seq.cor.weighted(L,ug)
   tmp <- eigen(tcrossprod(cbind(zeta.y,zeta.g)))
   Upsilon1 <- tmp$vectors[,1]
   
   if(draw.plot)
   {
      if(save.plot)
      {
         png(path.plot,width=600,height=600)
            par(mar=c(4.5,4.5,1,1))
            plot(Upsilon1,zeta.y,xlab=expression(Upsilon[1]),ylab=expression(u[y]),pch=pch,col=col,cex.lab=1.5,cex.axis=1.2)
            abline(lm(zeta.y ~ Upsilon1)$coef,lty=2)
         dev.off()
      } else
      {
         plot(Upsilon1,zeta.y,xlab=expression(Upsilon[1]),ylab=expression(u[y]),pch=pch,col=col,cex.lab=1.5,cex.axis=1.2)
         abline(lm(zeta.y ~ Upsilon1)$coef,lty=2)
      }
   }
   
   return(list(uy=uy,ug=ug,zeta.y=zeta.y,zeta.g=zeta.g,Upsilon1=Upsilon1,zeta.eigen.values=tmp$values[1:2]))
}

















