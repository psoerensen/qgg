####################################################################################################################
#    Module 6: Genomic Scoring
####################################################################################################################
#'
#' Genomic scoring based on single marker summary statistics
#'
#' @description
#' The gscore function is used for genomic predictions based on single marker summary statistics
#' (coefficients, log-odds ratios, z-scores) and observed genotypes.
#'

#' @param stat matrix of single marker effects
#' @param Glist list of information about genotype matrix
#' @param bedfiles name of the PLINK bed-files
#' @param famfiles name of the PLINK fam-files
#' @param bimfiles name of the PLINK bim-files
#' @param ids vector of individuals used in the analysis
#' @param scale logical if TRUE the genotype markers have been scale to mean zero and variance one
#' @param impute logical if TRUE missing genotypes are set to its expected value (2*af where af is allele frequency)
#' @param msize number of genotype markers used for batch processing
#' @param ncores number of cores used in the analysis
#' @param verbose is a logical; if TRUE it prints more details during optimization
#' @param chr chromosome for which genomic scores is computed
#' @param fit fit object output from gbayes

#' @author Peter Soerensen

#' @examples
#'

#'  ## Plink bed/bim/fam files
#'  bedfiles <- system.file("extdata", paste0("sample_chr",1:2,".bed"), package = "qgg")
#'  bimfiles <- system.file("extdata", paste0("sample_chr",1:2,".bim"), package = "qgg")
#'  famfiles <- system.file("extdata", paste0("sample_chr",1:2,".fam"), package = "qgg")
#'  
#'  # Summarize bed/bim/fam files
#'  Glist <- gprep(study="Example", bedfiles=bedfiles, bimfiles=bimfiles, famfiles=famfiles)
#'  
#'  # Simulate phenotype
#'  sim <- gsim(Glist=Glist, chr=1, nt=1)
#'  
#'  # Compute single marker summary statistics
#'  stat <- glma(y=sim$y, Glist=Glist, scale=FALSE)
#'  
#'  # Compute genomic scores
#'  gsc <- gscore(Glist = Glist, stat = stat)
#'  
#' @export
#'

gscore <- function(Glist = NULL, chr = NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, stat = NULL, fit = NULL, ids = NULL, scale = TRUE, impute = TRUE, msize = 100, ncores = 1, verbose=FALSE) {
     
     if ( !is.null(Glist))  {
          prs <- NULL
          if (!is.null(chr)) chromosomes <- chr
          if (is.null(chr)) chromosomes <- 1:length(Glist$bedfiles)
          if(sum(colnames(stat)%in%c("rsids","chr","pos", "ea","nea", "eaf","bm","dm")) == 8) {
            # Output from gbayes detected
            stat <- stat[,1:7]
          }
          for (chr in chromosomes) {
               if( any(stat$rsids %in% Glist$rsids[[chr]]) ) {
                 prschr <- run_gscore(Glist=Glist, chr=chr, stat = stat, 
                                      ids = ids, scale = scale, ncores = ncores, msize = msize, verbose=verbose)
                 if (chr==chromosomes[1]) prs <- prschr
                 if (!chr==chromosomes[1]) prs <- prs + prschr
                 
               }
          }
          prs <- scale(prs[,1:ncol(prs),drop=FALSE])
     }
     if ( !is.null(bedfiles))  {
          prs <- run_gscore(bedfiles=bedfiles, bimfiles=bimfiles, famfiles=famfiles, stat = stat, 
                            ids = ids, scale = scale, impute = impute, msize = msize, ncores = ncores, verbose=verbose)
     }   
     return(prs)
}


run_gscore <- function(Glist = NULL, chr=NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, stat = NULL, ids = NULL, scale = NULL, impute = TRUE, msize = 100, ncores = 1, verbose=FALSE) {
     
     #if(sum(is.na(stat))>0) stop("stat object contains NAs")
     if(sum(is.na(stat))>0) {
       warning("stat object contains NAs")
       stat <- na.omit(stat)
     }
     #if(is.null(Glist) & is.null(bedfiles)) stop("Please provide Glist or bedfile")
     
     if (!is.null(bedfiles)) {
          if(!file.exists(bedfiles)) stop(paste("bedfiles does not exists:"),bedfiles) 
          Glist <- NULL
          Glist$bedfiles <- bedfiles[1]
          if (!is.null(bimfiles)) Glist$bimfiles <- bimfiles[1]
          if (!is.null(famfiles)) Glist$famfiles <- famfiles[1]
          if (is.null(bimfiles)) Glist$bimfiles <- gsub(".bed", ".bim", bedfiles[1])
          if (is.null(famfiles)) Glist$famfiles <- gsub(".bed", ".fam", bedfiles[1])
          if(!file.exists(Glist$bedfiles)) stop(paste("bedfiles does not exists:"),Glist$famfiles) 
          if(!file.exists(Glist$famfiles)) stop(paste("famfiles does not exists:"),Glist$famfiles) 
          if(!file.exists(Glist$bimfiles)) stop(paste("bimfiles does not exists:"),Glist$bimfiles) 
          
          # Read fam information
          fam <- data.table::fread(input = Glist$famfiles[1], header = FALSE, data.table = FALSE, colClasses = "character")
          Glist$ids <- as.character(fam[, 2])
          if (any(duplicated(Glist$ids))) stop("Duplicated ids found in famfiles")
          message(paste("Finished processing fam file", famfiles[1]))
          
          # Read bim information
          bim <- data.table::fread(input = Glist$bimfiles[1], header = FALSE, data.table = FALSE, colClasses = "character")
          Glist$rsids <- list(as.character(bim[, 2]))
          Glist$a1 <- list(as.character(bim[, 5]))
          Glist$a2 <- list(as.character(bim[, 6]))
          Glist$pos <- list(as.numeric(bim[, 4]))
          Glist$chr <- as.character(bim[, 1])
          message(paste("Finished processing bim file", bimfiles[1]))
          
          Glist$n <- length(Glist$ids)
          Glist$m <- length(Glist$rsids)
          chr <- 1
     }
     
     
     # Prepase summary stat
     if (!sum(colnames(stat)[1:6] == c("rsids","chr","pos", "ea","nea", "eaf")) == 6) {
       stop("First six columns in data frame stat should be: rsids, chr, pos, ea, nea, eaf")
     }
       
     rsidsOK <- stat$rsids %in% Glist$rsids[[chr]]
     
     if (any(!rsidsOK)) {
       message(paste("Number of variants used:", sum(rsidsOK)))
       stat <- stat[rsidsOK, ]
     }
     S <- stat[, -c(1:6),drop=FALSE]
     if (is.vector(S)) S <- as.matrix(S)
     S <- apply(S, 2, as.numeric)
     colnames(S) <- colnames(stat)[-c(1:6)]
     rownames(S) <- rownames(stat)
     rsids <- as.character(stat$rsids)
     af <- stat$eaf
     
     # Prepare input data for mpgrs
     rws <- 1:Glist$n
     if (!is.null(ids)) rws <- match(ids, Glist$ids)
     cls <- match(rsids, Glist$rsids[[chr]])
     if(any( !stat$ea == Glist$a1[[chr]][cls] )) {
       warning("Some variants appear to be flipped => changing sign of variant effect for those variants ")
       flipped <- !stat$ea == Glist$a1[[chr]][cls]
       S[flipped,] <- -S[flipped,]  
     }
     if(!file.exists(Glist$bedfiles[chr])) stop(paste("bed file does not exists:"),Glist$bedfiles[chr]) 
     
     
     # multiple core using openblas
     if(ncores>1) {
          message(paste("Processing bed file", Glist$bedfiles[chr]))
          nt <- ncol(S)
          grs <- matrix(0,nrow=nt,ncol=Glist$n)
          rsids <- splitWithOverlap(rsids,msize,0)
          nsets <- length(rsids)
          print(paste("Processing chromosome", chr))
          for (set in 1:nsets) {
               cls <- match(rsids[[set]],Glist$rsids[[chr]])
               W <- getG(Glist=Glist, cls=cls, chr=chr, scale=scale)
               grs <- grs + tcrossprod(t(S[rsids[[set]],]),W)
               if(verbose) print(paste("Processing segment",set, "of", nsets,"on chromosome", chr))
          }
          grs <- t(grs)
          rownames(grs) <- Glist$ids
          colnames(grs) <- colnames(S)
          gc()
     }
     
     # single core
     if(ncores==1) {
          message(paste("Processing bed file", Glist$bedfiles[chr]))
       
          Slist <- vector(ncol(S),mode="list")
          for (j in 1:ncol(S)) {
               Slist[[j]] <- S[,j]
          }
          cls <- match(stat$rsids, Glist$rsids[[chr]])
          #af <- stat$eaf
          af <- Glist$a1[[chr]][cls]
          if(scale) grs <- .Call("_qgg_mtgrsbed", Glist$bedfiles[chr], Glist$n, cls, af, scale, Slist)
          if(!scale) grs <- .Call("_qgg_mtgrsbed", Glist$bedfiles[chr], n=Glist$n, cls=cls, af=af, scale=FALSE, Slist)
          grs <- do.call(cbind, grs)
          #grs <- as.matrix(as.data.frame(grs))
          rownames(grs) <- Glist$ids
          colnames(grs) <- colnames(S)
     }

          
     # # multicore     
     # if(ncores>1) {
     #      size <- ceiling(length(cls)/ncores)
     #      af <- 1-Glist$af[cls]
     #      af <- split(af, ceiling(seq_along(af) / size))
     #      cls <- split(cls, ceiling(seq_along(cls) / size))
     #      rwsS <- nrow(S)
     #      rwsS <- split(rwsS, ceiling(seq_along(rwsS) / size))
     #      Slist <- vector(length(rwsS),mode="list")
     #      for (j in 1:length(rwsS)) {
     #           Slist[[j]] <- lapply(seq_len(ncol(S)), function(i) S[rwsS[[j]],i])
     #      }
     #      
     #      grslist <- mclapply(1:length(cls), function(set) { .Call("_qgg_mtgrsbed", Glist$bedfiles, Glist$n, cls[[set]], af[[set]], scale, Slist[[set]]  ) }, mc.cores = ncores)
     #      grs <- as.matrix(as.data.frame(grslist[[1]]))
     #      for (j in 2:length(grslist)) {
     #           grs <- grs + as.matrix(as.data.frame(grslist[[j]]))      
     #      }
     #      rownames(grs) <- Glist$ids
     #      colnames(grs) <- colnames(S)
     #      # end multicore
     # }
     
     return(grs)
}


# compute effective number of observations
neff <- function(seb=NULL,af=NULL,Vy=1) {
  seb2 <- seb**2
  vaf <- 2*af*(1-af)
  neff <- round(median(Vy/(vaf*seb2)))
  return(neff)
}


#' Adjust marker effects based on correlated information
#'
#'
#' @description
#' The mtadj function use selection index theory to find the optimal weights across n traits, which is used to adjust marker effects by n correlated traits.
#' (https://www.nature.com/articles/s41467-017-02769-6)
#'    
#' @param h2 vector of heritability estimates
#' @param rg n-by-n matrix of genetic correlations
#' @param n vector of sample size used to estimate marker effects for each trait
#' @param b matrix of marker effects
#' @param z matrix of z-scores
#' @param meff effective number of uncorrelated genomic segments (default=60,000)
#' @param mtotal total number of markers 
#' @param statistics which kind of statistics ("b" or "z") used in the analysis
#' @param method method used to estimate marker effects; OLS: ordinary least square (default), or BLUP: best linear unbiased prediction
#' @param stat dataframe with marker summary statistics
#' 
#' @return Matrix of adjusted marker effects for each trait
#' 
#' @author Palle Duun Rohde and Peter Soerensen
#' 
#' @examples
#' 
#'  #bedfiles <- system.file("extdata", "sample_22.bed", package = "qgg")
#'  #bimfiles <- system.file("extdata", "sample_22.bim", package = "qgg")
#'  #famfiles <- system.file("extdata", "sample_22.fam", package = "qgg")
#'  #Glist <- gprep(study="1000G", bedfiles=bedfiles, bimfiles=bimfiles,famfiles=famfiles)
#'  #Glist <- gprep(Glist, task="sparseld",  msize=200)
#'  #
#'  ##Simulate data
#'  #set.seed(23)
#'  #
#'  #W <- getG(Glist, chr=1, scale=TRUE)
#'  #causal <- sample(1:ncol(W),50)
#'  #set1 <- c(causal, sample(c(1:ncol(W))[-causal],10))
#'  #set2 <- c(causal, sample(c(1:ncol(W))[-set1],10))
#'  #
#'  #b1 <- rnorm(length(set1))
#'  #b2 <- rnorm(length(set2))
#'  #y1 <- W[, set1]%*%b1 + rnorm(nrow(W))
#'  #y2 <- W[, set2]%*%b2 + rnorm(nrow(W))
#'  #
#'  ## Create model
#'  #data1 <- data.frame(y = y1, mu = 1)
#'  #data2 <- data.frame(y = y2, mu = 1)
#'  #X1 <- model.matrix(y ~ 0 + mu, data = data1)
#'  #X2 <- model.matrix(y ~ 0 + mu, data = data2)
#'  #
#'  ## Linear model analyses and single marker association test
#'  #maLM1 <- glma(y=y1, X=X1,W = W)
#'  #maLM2 <- glma(y=y2,X=X2,W = W)
#'  #
#'  ## Compute genetic parameters
#'  #z1 <- maLM1[,"stat"]
#'  #z2 <- maLM2[,"stat"]
#'  #
#'  #z <- cbind(z1=z1,z2=z2)
#'  #
#'  #h2 <- ldsc(Glist, z=z, n=c(500,500), what="h2")
#'  #rg <- ldsc(Glist, z=z, n=c(500,500), what="rg")
#'  #
#'  ## Adjust summary statistics using estimated genetic parameters
#'  #b <- cbind(b1=maLM1[,"b"],b2=maLM2[,"b"])
#'  #bm <- mtadj( h2=h2, rg=rg, b=b, n=c(500,500), method="ols")
#'  
#'  
#' @export
#' 

mtadj <- function(h2=NULL, rg=NULL, stat=NULL, b=NULL, z=NULL, n=NULL, mtotal=NULL, 
                  meff=60000, method="ols", statistics="z") {
  
  if(!is.null(z)) b <- z
  if(!is.null(stat)) {
    b <- stat$b
    if(statistics=="z") b <- stat$b/stat$seb     
    n <- colMeans(stat$n)
  }
  m <- mtotal
  if(is.null(mtotal)) m <- nrow(b)
  
  if(is.null(b)) stop("Marker effect matrix b is missing")
  if(is.null(h2)) stop("Heritability vector h2 is missing")
  if(is.null(rg)) stop("Correlation matrix rg is missing")
  if(is.null(n)) stop("n missing")
  
  cnames <- colnames(b)
  
  # compute r2
  r2 <- rsq(h2=h2,meff=meff,n=n)
  
  # V_sblup
  if (method=="blup") VS <- diag(r2/m)       # eq. 18
  if (method=="ols") VS <- diag(h2/m + 1/n)  # eq. 25 
  nt <- ncol(VS)
  for(i in 1:nt) {
    for(j in i:nt) {
      if(!i==j) {
        # eq. 19 and 26     
        if (method=="blup") VS[i,j] <- (rg[i,j]*r2[i]*r2[j])/(sqrt(h2[i])*sqrt(h2[j])*m)
        if (method=="ols") VS[i,j] <- (rg[i,j]*sqrt(h2[i])*sqrt(h2[j]))/m
        VS[j,i] <- VS[i,j]
      }
    }
  }
  # C_sblup
  CS <- matrix(0,nt,nt)
  for(i in 1:nt) {
    for(j in 1:nt) {
      # eq. 20 and 27      
      if (method=="blup") CS[j,i] <- rg[i,j]*(r2[j]/m)*(sqrt(h2[i])/sqrt(h2[j]))
      if (method=="ols") CS[j,i] <- (rg[i,j]*sqrt(h2[j])*sqrt(h2[i]))/m
    }
  }
  invVS <- solve(VS)
  weights <- invVS%*%CS
  b <- b%*%weights
  colnames(b) <- cnames
  if(!is.null(stat)) {
    if(statistics=="z") stat$z <- b
    if(statistics=="b") stat$b <- b
    stat$weights <- weights
    stat$cvblup <- CS
    stat$vblup <- VS
    return(stat)
  }
  if(is.null(stat)) {
    if(!is.null(z)) return(list(z=b, weights=weights, CS=CS, VS=VS))
    if(is.null(z)) return(list(b=b, weights=weights, CS=CS, VS=VS))
  }
}

# compute expected r-squared
rsq <- function(h2=NULL,meff=NULL,n=NULL) {
  phi <- meff/n
  rsq <- (phi + h2 - sqrt((phi+h2)**2 - 4*phi*h2**2))/(2*phi)
  rsq
}


#' Compute Receiver Operating Curve statistics
#' 
#' @description
#' Compute ROC
#'
#' @author Palle Duun Rohde
#' 
#' @param yobs vector of observed phenotype
#' @param ypred vector of predicted phenotype
#' 
#' @keywords internal
#' 
#' @export
#' 

computeROC <- function(yobs=NULL, ypred=NULL){
  if(!is.matrix(ypred)){
    y <- data.frame(yobs=yobs==1, ypred=ypred)
    y <- y[order(y[, 2], decreasing = TRUE), ]
    roc.df <- cbind(TPR=cumsum(y[,1])/sum(y[,1]), FPR=cumsum(!y[,1])/sum(!y[,1]), yobs=as.numeric(y[,1]), ypred=y[,2])
  }
  
  if(is.matrix(ypred)){
    roc.df <- vector(ncol(ypred),mode="list")
    names(roc.df) <- colnames(ypred)
    for(i in 1:ncol(ypred)){
      y <- data.frame(yobs=yobs==1, ypred=ypred[,i])
      y <- y[order(y[, 2], decreasing = TRUE), ]
      roc.df[[i]] <- cbind(TPR=cumsum(y[,1])/sum(y[,1]), FPR=cumsum(!y[,1])/sum(!y[,1]), yobs=as.numeric(y[,1]), ypred=y[,2])
    }
  }
  return(roc.df)
}

#' Plot Receiver Operating Curves
#'
#' @description
#' Plot ROC
#'
#' @param roc.data data frame with ROC information (from computeROC)
#' @param cols which columns should be used in the ROC plot
#' 
#' @keywords internal

#' @author Palle Duun Rohde
#' 
#' @export
#' 

plotROC <- function(roc.data=NULL, cols = NULL){
    if(is.matrix(roc.data)){
       if(is.null(cols)){cols <- "black"}
       plot(x=roc.data[,2], y=roc.data[,1], bty="n", type="l", las=1, cex.axis=.8,
            xlab="FRP (1-Specificity)", ylab="TPR (Sensitivity)", col=cols)
       abline(a=0, b=1, lty=2)
       text(x=.95, y=.8,label=paste("AUC = ", round(auc(ypred=roc.data[,4], yobs=roc.data[,3]),3),sep=""),cex=.9)
    }
    if(is.list(roc.data)){
        if(is.null(cols)){cols <- 1:length(roc.data)}
        plot(x=0, y=0, bty="n", type="n", las=1, cex.axis=.8, xlab="FRP (1-Specificity)", ylab="TPR (Sensitivity)",
             xlim=c(0,1), ylim=c(0,1))
        abline(a=0, b=1, lty=2)
        print.auc <- NULL
        for(i in 1:length(roc.data)){
            points(roc.data[[i]][,2], y=roc.data[[i]][,1], type="l",col=cols[i])
            print.auc <- c(print.auc, round(auc(ypred=roc.data[[i]][,4], yobs=roc.data[[i]][,3]),3))
        }
        legend("bottomright", legend=paste(names(roc.data), " (AUC = ", print.auc, ")",sep=""), lty=1, col=cols, bty="n", cex=.8)
    }
}
