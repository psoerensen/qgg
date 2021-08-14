####################################################################################################################
#    Module 6: GSCORE
####################################################################################################################
#'
#' Genomic prediction based on single marker summary statistics
#'
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

#' @author Peter Soerensen

#' @examples
#'

#' bedfiles <- system.file("extdata", "sample_22.bed", package = "qgg")
#' bimfiles <- system.file("extdata", "sample_22.bim", package = "qgg")
#' famfiles <- system.file("extdata", "sample_22.fam", package = "qgg")
#' 
#' Glist <- gprep(study="1000G", bedfiles=bedfiles, bimfiles=bimfiles,
#'                famfiles=famfiles, overwrite=TRUE)
#' 
#' rsids <- Glist$rsids
#' stat <- data.frame(rsids=Glist$rsids,alleles=Glist$a2, af=Glist$af, effect=rnorm(Glist$m))
#' 
#' W <- getW(Glist=Glist,rsids=Glist$rsids)
#' pgs1 <- W%*%stat[,4]
#' 
#' pgs2 <- gscore(Glist = Glist, stat = stat) 
#' 
#' pgs3 <- gscore(bedfiles=bedfiles, stat = stat) 
#' 
#' pgs4 <- gscore(bedfiles=bedfiles,bimfiles=bimfiles,famfiles=famfiles, stat = stat) 
#' 
#' 
#' cor(cbind(pgs1,pgs2,pgs3,pgs4))
#'

#' @export
#'

gscore <- function(Glist = NULL, chr = NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, stat = NULL, fit = NULL, ids = NULL, scale = TRUE, impute = TRUE, msize = 100, ncores = 1) {
     
     if ( !is.null(Glist))  {
          prs <- NULL
          # if(is.null(stat) & !is.null(fit)) {
          #   rsids <- unlist(Glist$rsids)
          #   af <- unlist(Glist$af)
          #   alleles <- unlist(Glist$a2)
          #   cls <- match(names(fit$bm),rsids)
          #   if(any(is.na(cls))) stop("Missing rsids")
          #   stat <- data.frame(rsids=names(fit$bm), alleles=alleles[cls], af=af[cls], effect=fit$b)
          #   rownames(stat) <- names(fit$bm)
          # }
          if (!is.null(chr)) chromosomes <- chr
          if (is.null(chr)) chromosomes <- 1:length(Glist$bedfiles)
          for (chr in chromosomes) {
               if( any(stat$rsids %in% Glist$rsids[[chr]]) ) {
                 prschr <- run_gscore(Glist=Glist, chr=chr, stat = stat, 
                                      ids = ids, scale = scale, ncores = ncores)
                 if (is.null(prs)) prs <- prschr
                 if (!is.null(prs)) prs <- prs + prschr
                 
               }
          }
     }
     if ( !is.null(bedfiles))  {
          prs <- run_gscore(bedfiles=bedfiles, bimfiles=bimfiles, famfiles=famfiles, stat = stat, 
                            ids = ids, scale = scale, impute = impute, msize = msize, ncores = ncores)
     }   
     return(prs)
}


run_gscore <- function(Glist = NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, stat = NULL, ids = NULL, scale = scale, impute = TRUE, msize = 100, ncores = 1) {
     
     if(sum(is.na(stat))>0) stop("stat object contains NAs") 
     if(is.null(Glist) & is.null(bedfiles)) stop("Please provide Glist or bedfile")
     
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
          Glist$position <- list(as.numeric(bim[, 4]))
          Glist$chr <- as.character(bim[, 1])
          message(paste("Finished processing bim file", bimfiles[1]))
          
          Glist$n <- length(Glist$ids)
          Glist$m <- length(Glist$rsids)
          chr <- 1
     }
     
     
     # Prepase summary stat
     if (!sum(colnames(stat)[1:4] == c("chr","rsids", "alleles", "af")) == 4) {
          stop("First three columns in data frame stat should be: chr, rsids, alleles, af ")
     }
     rsidsOK <- stat$rsids %in% Glist$rsids[[chr]]
     
     if (any(!rsidsOK)) {
       #warning("Some variants not found in genotype files")
       message(paste("Number of variants used:", sum(rsidsOK)))
       message(paste("Number of variants missing:", sum(!rsidsOK)))
       stat <- stat[rsidsOK, ]
       stat$rsids <- as.character(stat$rsids)
       stat$alleles <- as.character(stat$alleles)
     }
     S <- stat[, -c(1:4)]
     if (is.vector(S)) S <- as.matrix(S)
     S <- apply(S, 2, as.numeric)
     colnames(S) <- colnames(stat)[-c(1:4)]
     rsids <- as.character(stat$rsids)
     af <- stat$af
     
     # Prepare input data for mpgrs
     rws <- 1:Glist$n
     if (!is.null(ids)) rws <- match(ids, Glist$ids)
     cls <- match(rsids, Glist$rsids[[chr]])
     if(any( !stat$alleles == Glist$a2[[chr]][cls] )) {
       warning("Some variants appear to be flipped => changing sign of variant effect for those variants ")
       flipped <- !stat$alleles == Glist$a2[[chr]][cls]
       S[flipped,] <- -S[flipped,]  
     }
     if(!file.exists(Glist$bedfiles)) stop(paste("bed file does not exists:"),Glist$bedfiles) 
     

     # single core
     if(ncores==1) {
          message(paste("Processing bed file", Glist$bedfiles))
       
          Slist <- vector(ncol(S),mode="list")
          for (j in 1:ncol(S)) {
               Slist[[j]] <- S[,j]
          }
          cls <- match(stat$rsids, Glist$rsids[[chr]])
          af <- stat$af
          if(scale) grs <- .Call("_qgg_mtgrsbed", Glist$bedfiles, Glist$n, cls, af, scale, Slist)
          if(!scale) grs <- .Call("_qgg_mtgrsbed", Glist$bedfiles, n=Glist$n, cls=cls, af=af, scale=FALSE, Slist)
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


neff <- function(seb=NULL,af=NULL,Vy=1) {
  seb2 <- seb**2
  vaf <- 2*af*(1-af)
  neff <- round(median(Vy/(vaf*seb2)))
  return(neff)
}

# compute r-squared
rsq <- function(h2=NULL,me=NULL,n=NULL) {
  phi <- me/n
  rsq <- (phi + h2 - sqrt((phi+h2)**2 - 4*phi*h2**2))/(2*phi)
  rsq
}


#' Adjust marker effects based on correlated information
#'
#'
#' @description
#' The adjustB function use selection index theory to find the optimal weights across n traits, which is used to adjust marker effects by n correlated traits.
#'
#' @param h2 vector of heritability estimates
#' @param rg n-by-n matrix of genetic correlations
#' @param n vector of sample size used to estimate marker effects for each trait
#' @param b matrix of marker effects
#' @param me effective number of uncorrelated genomic segments (default=60,000)
#' @param method method used to estimate marker effects; OLS: ordinary least square (default), or BLUP: best linear unbiased prediction
#' 
#' @return Matrix of adjusted marker effects for each trait
#' 
#' @author Peter Soerensen
#' 
#' @examples
#' bedfiles <- system.file("extdata", "sample_22.bed", package = "qgg")
#' bimfiles <- system.file("extdata", "sample_22.bim", package = "qgg")
#' famfiles <- system.file("extdata", "sample_22.fam", package = "qgg")
#' Glist <- gprep(study="1000G", bedfiles=bedfiles, bimfiles=bimfiles,famfiles=famfiles)
#' Glist <- gprep(Glist, task="sparseld",  msize=200)
#' 
#' #Simulate data
#' set.seed(23)
#' 
#' W <- getG(Glist, chr=1, scale=TRUE)
#' causal <- sample(1:ncol(W),50)
#' set1 <- c(causal, sample(c(1:ncol(W))[-causal],10))
#' set2 <- c(causal, sample(c(1:ncol(W))[-set1],10))
#' 
#' b1 <- rnorm(length(set1))
#' b2 <- rnorm(length(set2))
#' y1 <- W[, set1]%*%b1 + rnorm(nrow(W))
#' y2 <- W[, set2]%*%b2 + rnorm(nrow(W))
#' 
#' # Create model
#' data1 <- data.frame(y = y1, mu = 1)
#' data2 <- data.frame(y = y2, mu = 1)
#' X1 <- model.matrix(y ~ 0 + mu, data = data1)
#' X2 <- model.matrix(y ~ 0 + mu, data = data2)
#' 
#' # Linear model analyses and single marker association test
#' maLM1 <- lma(y=y1, X=X1,W = W)
#' maLM2 <- lma(y=y2,X=X2,W = W)
#' 
#' # Compute genetic parameters
#' z1 <- maLM1[,"stat"]
#' z2 <- maLM2[,"stat"]
#' 
#' z <- cbind(z1=z1,z2=z2)
#' 
#' h2 <- ldsc(Glist, z=z, n=c(500,500), what="h2")
#' rg <- ldsc(Glist, z=z, n=c(500,500), what="rg")
#' 
#' # Adjust summary statistics using estiamted genetic parameters
#' b <- cbind(b1=maLM1[,"b"],b2=maLM2[,"b"])
#' badj <- adjustB( h2=h2[,2], rg=rg, b=b, n=c(500,500), method="ols")
#' 
#' @export
#' 

adjustB <- function(h2=NULL, rg=null, b=NULL,  n=NULL, me=60000, method="ols") {
  
  if(is.null(b)) stop("Marker effect matrix b is missing")
  if(is.null(h2)) stop("Heritability vector h2 is missing")
  if(is.null(rg)) stop("Correlation matrix rg is missing")
  cnames <- colnames(b)
  m <- nrow(b)
  
  # compute r2
  r2 <- rsq(h2=h2,me=me,n=n)
  
  # V_sblup
  if (method=="blup") VS <- diag(r2/m)
  if (method=="ols") VS <- diag(h2/m + 1/n)
  nt <- ncol(VS)
  for(i in 1:nt) {
    for(j in i:nt) {
      if(!i==j) {
        if (method=="blup") VS[i,j] <- (rg[i,j]*r2[i]*r2[j])/(sqrt(h2[i]*h2[j])/m)
        if (method=="ols") VS[i,j] <- rg[i,j]*sqrt(h2[i])*sqrt(h2[j])/m
        VS[j,i] <- VS[i,j]
      }
    }
  }
  # C_sblup
  CS <- matrix(0,nt,nt)
  for(i in 1:nt) {
    for(j in 1:nt) {
      if (method=="blup") CS[i,j] <- (rg[i,j]*r2[j]*sqrt(h2[i])/sqrt(h2[j]))/m
      if (method=="ols") CS[i,j] <- (rg[i,j]*sqrt(h2[i])*sqrt(h2[j]))/m
    }
  }
  invVS <- solve(VS)
  weights <- invVS%*%CS
  b <- t(tcrossprod(weights,b))
  colnames(b) <- cnames
  return(b)
}
