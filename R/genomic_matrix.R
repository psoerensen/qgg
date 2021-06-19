###############################################################################################
#  Prepare (processed genotypes) for target population
###############################################################################################

#' Prepare genotype data for all statistical analyses (initial step)
#'
#' @description
#' All functions in qgg relies on a simple data infrastructure that takes five main input sources;
#' phenotype data (y), covariate data (X), genotype data (G or Glist), a genomic relationship
#' matrix (GRM or GRMlist) and genetic marker sets (sets).
#'
#' The genotypes are stored in a matrix (n x m (individuals x markers)) in memory (G) or in a
#' binary file on disk (Glist).h
#'
#' It is only for small data sets that the genotype matrix (G) can stored in memory. For large data
#' sets the genotype matrix has to stored in a binary file on disk (Glist). Glist is as a list
#' structure that contains information about the genotypes in the binary file.
#'
#' The gprep function prepares the Glist, and is required for downstream analyses of large-scale
#' genetic data. Typically, the Glist is prepared once, and saved as an *.Rdata-file.
#'
#' The gprep function reads genotype information from binary PLINK files, and creates the Glist
#' object that contains general information about the genotypes such as reference alleles,
#' allele frequencies and missing genotypes, and construct a binary file on the disk that contains
#' the genotypes as allele counts of the alternative allele (memory usage = (n x m)/4 bytes).
#'
#' The gprep function can also be used to prepare sparse ld matrices.
#' The r2 metric used is the pairwise correlation between markers (allele count alternative allele)
#' in a specified region of the genome. The marker genotype is allele count of the alternative allele
#' which is assumed to be centered and scaled.
#'
#' The Glist structure is used as input parameter for a number of qgg core functions including:
#' 1) construction of genomic relationship matrices (grm), 2) construction of sparse ld matrices,
#' 3) estimating genomic parameters (greml), 4) single marker association analyses (lma or mlma),
#' 5) gene set enrichment analyses (gsea), and 6) genomic prediction from genotypes
#' and phenotypes (gsolve) or genotypes and summary statistics (gscore).
#'

#'
#' @param Glist only provided if task="summary" or task="sparseld"
#' @param task character specifying which task to perform ("prepare" is default, "summary", or "sparseld")
#' @param study name of the study
#' @param fnBED path and filename of the binary file .bed used for storing genotypes on the disk
#' @param bedfiles vector of names for the PLINK bed-files
#' @param famfiles vector of names for the PLINK fam-files
#' @param bimfiles vector of names for the PLINK bim-files
#' @param ids vector of individuals used in the study
#' @param rsids vector of marker rsids used in the study
#' @param fnLD path and filename of the binary files .ld for storing sparse ld matrix on the disk
#' @param msize number of markers used in compuation of sparseld
#' @param overwrite logical if TRUE overwite binary genotype file
#' @param ncores number of cores used to process the genotypes
#'
#' @return Returns a list structure (Glist) with information about genotypes
#'


#' @author Peter Soerensen

#' @examples
#'
#' bedfiles <- system.file("extdata", "sample_22.bed", package = "qgg")
#' bimfiles <- system.file("extdata", "sample_22.bim", package = "qgg")
#' famfiles <- system.file("extdata", "sample_22.fam", package = "qgg")
#' 
#' if(!grepl("^darwin", R.version$os)) {
#'   fnBED <- tempfile(fileext=".raw")
#' 
#'   Glist <- gprep(study="1000G", fnBED=fnBED, bedfiles=bedfiles, bimfiles=bimfiles,
#'                famfiles=famfiles, overwrite=TRUE)
#' 
#'   file.remove(fnBED)
#' }
#' 


#' @export
#'

gprep <- function(Glist = NULL, task = "prepare", study = NULL, fnBED = NULL, fnLD = NULL,
                  bedfiles = NULL, bimfiles = NULL, famfiles = NULL, ids = NULL, rsids = NULL,
                  overwrite = FALSE, msize = 100, ncores = 1) {

  if (task == "prepare") {
    nfiles <- length(bedfiles)
    Glist <- NULL
    Glist$study <- study
    Glist$fnBED <- fnBED
    if (!is.null(fnBED)) {
      if (file.exists(fnBED)) warning(paste("fnBED allready exist"))
    }
    Glist$bedfiles <- bedfiles
    if (!is.null(bimfiles)) Glist$bimfiles <- bimfiles
    if (!is.null(famfiles)) Glist$famfiles <- famfiles
    if (is.null(bimfiles)) Glist$bimfiles <- gsub(".bed", ".bim", bedfiles)
    if (is.null(famfiles)) Glist$famfiles <- gsub(".bed", ".fam", bedfiles)

    # Read fam information
    fam <- fread(input = famfiles[1], header = FALSE, data.table = FALSE, colClasses = "character")
    Glist$ids <- as.character(fam[, 2])
    #Glist$study_ids <- Glist$ids
    Glist$study_ids <- NULL
    Glist$n <- length(Glist$ids)
    if (!is.null(ids)) {
      if (any(!ids %in% as.character(fam[, 2]))) warning(paste("some ids not found in famfiles"))
      Glist$study_ids <- Glist$ids[Glist$ids %in% as.character(ids)]
    }
    if (any(duplicated(Glist$ids))) stop("Duplicated ids found in famfiles")

    Glist$rsids <- vector(mode = "list", length = nfiles)
    Glist$mchr <- vector( length = nfiles)
    Glist$a1 <- vector(mode = "list", length = nfiles)
    Glist$a2 <- vector(mode = "list", length = nfiles)
    Glist$position <- vector(mode = "list", length = nfiles)
    Glist$chr <- vector(mode = "list", length = nfiles)

    Glist$nmiss <- vector(mode = "list", length = nfiles)
    Glist$af <- vector(mode = "list", length = nfiles)
    Glist$maf <- vector(mode = "list", length = nfiles)
    Glist$hom <- vector(mode = "list", length = nfiles)
    Glist$het <- vector(mode = "list", length = nfiles)
    Glist$n0 <- vector(mode = "list", length = nfiles)
    Glist$n1 <- vector(mode = "list", length = nfiles)
    Glist$n2 <- vector(mode = "list", length = nfiles)
#    Glist$cls <- vector(mode = "list", length = nfiles)

    for (chr in 1:length(bedfiles)) {
      bim <- fread(input = bimfiles[chr], header = FALSE, data.table = FALSE, colClasses = "character")
      rsidsBIM <- as.character(bim[, 2])
      if (!is.null(rsids)) bim <- droplevels(bim[rsidsBIM %in% rsids, ])
      fam <- fread(input = famfiles[chr], header = FALSE, data.table = FALSE, colClasses = "character")
      if (any(!Glist$ids %in% as.character(fam[, 2]))) stop(paste("some ids not found in famfiles"))
      Glist$a1[[chr]] <- as.character(bim[, 5])
      Glist$a2[[chr]] <- as.character(bim[, 6])
      Glist$position[[chr]] <- as.numeric(bim[, 4])
      Glist$rsids[[chr]] <- as.character(bim[, 2])
      Glist$mchr[chr] <- length(Glist$rsids[[chr]]) 
      Glist$chr[[chr]] <- as.character(bim[, 1])
      message(paste("Finished processing bim file", bimfiles[chr]))
      if (is.null(Glist$fnBED)) Glist <- summaryBED(Glist=Glist, chr=chr, ids = Glist$ids, ncores = ncores)
      message(paste("Finished processing bed file", bedfiles[chr]))
      
    }

    Glist$nchr <- length(Glist$bedfiles)
    

  }
     # if (task == "combine") {
     #      writeBED(bedfiles = bedfiles, bimfiles = bimfiles, famfiles = famfiles, 
     #               fnBED=fnBED, ids = ids, rsids = rsids, overwrite = FALSE) {
     #           
     # }

  if (task == "sparseld") {
    message("Computing ld")
    Glist$msize <- msize
    Glist$fnLD <- fnLD
    if (is.null(fnLD)) Glist$fnLD <- gsub(".bed", ".ld", Glist$bedfiles)
    Glist$rsidsLD <- vector(mode = "list", length = length(Glist$fnLD))
    Glist$lscore <- vector(mode = "list", length = length(Glist$fnLD))
    if (is.null(ids)) ids <- Glist$ids
    for( chr in 1:length(Glist$fnLD)) {
      Glist <- sparseLD(Glist = Glist, fnLD = Glist$fnLD[chr], msize = msize, chr = chr, rsids = rsids,
        ids = ids, ncores = 1)
      Glist$lscore[[chr]] <- ldscore( Glist=Glist, chr=chr) 
    }
  }

  return(Glist)
}




summaryBED <- function(Glist = NULL, ids = NULL, rsids = NULL, rws = NULL, cls = NULL, chr = NULL, ncores = 1) {

  n <- Glist$n
  m <- Glist$mchr[chr]
  
  rws <- 1:n
  if (!is.null(ids)) rws <- match(ids, Glist$ids)
  if (!is.null(Glist$study_ids)) warning("Study ids not used anymore")
  nr <- length(rws)

  if (is.null(cls)) cls <- 1:m
  if (!is.null(rsids)) cls <- match(rsids, Glist$rsids)
  nc <- length(cls)

  ######################################################################
  #  00 01 10 11         bit level  corresponds to
  #  0  1  2  3          xij level  corresponds to
  #  2  NA  1  0         number of copies of first allele in bim file
  ######################################################################
  
  freq <- .Call("_qgg_freqbed", Glist$bedfiles[chr], Glist$n, cls)
  nmiss <- freq[2,]
  hom <- (freq[1,] + freq[4,]) / (freq[1,] + freq[3,] + freq[4,])
  het <- (freq[3,]) / (freq[1,] + freq[3,] + freq[4,])
  #af <- (2*freq[1,] + freq[3,])/(2*freq[1,] + 2*freq[3,] + 2*freq[4,])
  nalleles <- 2*(n-nmiss)
  af <- 2*freq[1,] + freq[3,]
  af[nalleles>0] <- af[nalleles>0]/nalleles[nalleles>0]
  af[nalleles==0] <- 0.5
  
  maf <- af
  maf[maf > 0.5] <- 1 - maf[maf > 0.5]
  if(!is.null(chr)) {
    Glist$nmiss[[chr]] <- nmiss
    Glist$af[[chr]] <- af
    Glist$maf[[chr]] <- maf
    Glist$hom[[chr]] <- hom
    Glist$het[[chr]] <- het
    Glist$n0[[chr]] <- freq[4,]
    Glist$n1[[chr]] <- freq[3,]
    Glist$n2[[chr]] <- freq[1,]
  }
  return(Glist)
}


# writeBED <- function(bedfiles = NULL, bimfiles = NULL, famfiles = NULL, 
#                      fnBED=NULL, fnBIM=NULL, fnFAM=NULL, ids = NULL, rsids = NULL, overwrite = FALSE) {
#      if (file.exists(fnBED)) {
#           warning(paste("fnBED file allready exist"))
#           if (!overwrite) stop(paste("fnBED file allready exist"))
#      }
#      if (is.null(bimfiles)) bimfiles <- gsub(".bed", ".bim", bedfiles)
#      if (is.null(famfiles)) famfiles <- gsub(".bed", ".fam", bedfiles)
#      if (is.null(fnBIM)) fnBIM <- gsub(".bed", ".bim", fnBIM)
#      if (is.null(fnFAM)) fnFAM <- gsub(".bed", ".fam", fnFAM)
#      if (file.exists(fnBIM)) {
#           stop(paste("fnBIM file allready exist"))
#      }
#      if (file.exists(fnFAM)) {
#           stop(paste("fnFAM file allready exist"))
#      }
#      
#      bim_combined <- NULL
#      for (chr in 1:length(bedfiles)) {
#           message(paste("Processing bedfile:", bedfiles[chr]))
#           bim <- fread(input = bimfiles[chr], header = FALSE, data.table = FALSE, colClasses = "character")
#           fam <- fread(input = famfiles[chr], header = FALSE, data.table = FALSE)
#           n <- nrow(fam)
#           m <- nrow(bim)
#           rsidsBIM <- as.character(bim[, 2])
#           keep <- rep(TRUE, m)
#           if (!is.null(rsids)) keep <- rsidsBIM %in% rsids
#           cls <- (1:m)[keep]
#           
#           fnBEDCHR <- bedfiles[chr]
#           bfBEDCHR <- file(fnBEDCHR, "rb")
#           magic <- readBin(bfBEDCHR, "raw", n = 3)
#           if (!all(magic[1] == "6c", magic[2] == "1b", magic[3] == "01")) {
#                stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
#           }
#           close(bfBEDCHR)
#           append <- 1
#           if (chr == 1) append <- 0
#           result <- .Call("_qgg_bed2bed", fnBED, bedfiles[chr], n, cls)
#           bim_combined <- rbind(bim_combined,bim[keep,])
#           #res <- .Fortran("bed2raw",
#           #                m = as.integer(m),
#           #                cls = as.integer(keep),
#           #                nbytes = as.integer(nbytes),
#           #                append = as.integer(append),
#           #                fnBEDCHAR = as.integer(unlist(sapply(as.character(fnBED),charToRaw),use.names=FALSE)),
#           #                fnBEDCHAR = as.integer(unlist(sapply(as.character(fnBED),charToRaw),use.names=FALSE)),
#           #                ncharbed = nchar(as.character(fnBED)),
#           #                ncharraw = nchar(as.character(fnBED)),
#           #                PACKAGE = "qgg"
#           #)
#           message(paste("Finished processing bedfile:", bedfiles[chr]))
#      }
#      fwrite(bim_combined, file.name=fnBIM)
#      fwrite(fam_combined, file.name=fnFAM)
# }


#' Extract elements from genotype matrix (W) stored on disk
#'
#' @description
#' Extract elements from genotype matrix W (whole or subset) stored on disk.

#' @param Glist list structure with information about genotypes stored on disk
#' @param bedfiles vector of name for the PLINK bed-file
#' @param ids vector of ids in W to be extracted
#' @param rsids vector of rsids in W to be extracted
#' @param rws vector of rows in W to be extracted
#' @param cls vector of columns in W to be extracted
#' @param scale logical if TRUE the genotype markers have been scale to mean zero and variance one
#' @param impute logical if TRUE missing genotypes are set to its expected value (2*af where af is allele frequency)

#' @export
#'

getG <- function(Glist = NULL, chr = NULL, bedfiles = NULL, bimfiles = NULL, famfiles = NULL, ids = NULL, rsids = NULL,
                 rws = NULL, cls = NULL, impute = TRUE, scale = FALSE) {

if(!is.null(chr)) bedfiles <- Glist$bedfiles[chr]
if(is.null(cls)) cls <-  1:Glist$mchr[chr]
if (!is.null(rsids)) cls <- match(rsids, Glist$rsids[[chr]])
if (any(is.na(cls))) {
  warning(paste("some rsids not found in Glist"))
  message(rsids[is.na(cls)])
}
cls <- cls[!is.na(cls)]

af <- Glist$af[[chr]][cls]
if(scale) W <- .Call("_qgg_readW", bedfiles, Glist$n,cls,af)
if(!scale) W <- .Call("_qgg_readG", bedfiles, Glist$n,cls)
colnames(W) <- Glist$rsids[[chr]][cls]
rownames(W) <- Glist$ids
if(!is.null(rws)) W <- W[rws,]
if(is.integer(impute)) W[W==impute] <- impute
return(W)
}

#' @export
#'

cvs <- function(y=NULL, Glist = NULL, chr = NULL, bedfiles = NULL, bimfiles = NULL, famfiles = NULL, ids = NULL, rsids = NULL,
                 rws = NULL, cls = NULL, impute = TRUE, scale = TRUE) {
  
  if(is.null(cls)) cls <- 1:Glist$mchr[chr]
  af <- Glist$af[[chr]][cls]
  ids <- NULL
  if(is.matrix(y)) ids <- rownames(y)
  if(is.vector(y)) ids <- names(y)
  if(is.vector(y)) y <- as.matrix(y)
  for (i in 1:ncol(y)) {
    y[,i] <- y[,i] - mean(y[,i])
  }
  if(is.null(ids)) warning("No names/rownames provided for y")
  if(!is.null(ids)) {
    if(any(is.na(match(ids,Glist$ids))))  stop("Names/rownames for y does match rownames for W")
    if(any(is.na(Glist$ids%in%ids)))  stop("Names/rownames for y does match rownames for W")
  }
  dfe <- nrow(y)-1
  weights <- list(rep(1.0,nrow(y)))
  rws <- match(ids,Glist$ids)
  if(any(is.na(rws))) stop("Some ids in y does not match individuals in Glist$ids")
  if(!is.null(ids) & any(duplicated(ids)) ) {
    ylist <- apply(y,2,function(x){ split(x,f=factor(ids))})
    weights <- lapply(ylist,function(x){sapply(x,length)}) 
    y <- lapply(ylist,function(x){sapply(x,sum)}) 
  }
  if(!is.list(y)) y <- list(y)
  #if(!length(y)==Glist$n) stop("Length of y does not match number of individuals in Glist$n")
  covs <- .Call("_qgg_summarybed", 
                Glist$bedfiles[chr], 
                n=Glist$n,
                cls=cls,
                af=af,
                weights=weights,
                y=y)
  covs <- c(covs,list(NULL),list(NULL),list(NULL),list(NULL))
  for ( i in 1:length(covs[[1]])) {
    names(covs[[1]][[i]]) <- Glist$rsids[[chr]][cls]
    names(covs[[2]][[i]]) <- Glist$rsids[[chr]][cls]
    covs[[3]][[i]] <- (covs[[2]][[i]]/covs[[1]][[i]])
    covs[[4]][[i]] <- 1/sqrt(covs[[1]][[i]])
    covs[[5]][[i]] <- (covs[[2]][[i]]/covs[[1]][[i]])*sqrt(covs[[1]][[i]])
    ptt <- 2 * pt(-abs(covs[[5]][[i]]), df = dfe)
    covs[[6]][[i]] <- ptt
  }
  names(covs) <- c("XX","Xy","b","seb","tstat","p")
  return(covs)
}



#' Extract elements from genotype matrix (W) stored on disk
#'
#' @description
#' Extract elements from genotype matrix W (whole or subset) stored on disk.

#' @param Glist list structure with information about genotypes stored on disk
#' @param bedfiles vector of name for the PLINK bed-file
#' @param ids vector of ids in W to be extracted
#' @param rsids vector of rsids in W to be extracted
#' @param rws vector of rows in W to be extracted
#' @param cls vector of columns in W to be extracted
#' @param scale logical if TRUE the genotype markers have been scale to mean zero and variance one
#' @param impute logical if TRUE missing genotypes are set to its expected value (2*af where af is allele frequency)
#' @param allele vector of alleles to be extracted

#' @export
#'

getW <- function(Glist = NULL, chr = NULL, bedfiles = NULL, bimfiles = NULL, famfiles = NULL, ids = NULL, rsids = NULL,
                     rws = NULL, cls = NULL, impute = TRUE, scale = FALSE,
                     allele = NULL) {
     

  if (!is.null(Glist)) {
    if (is.null(chr)) stop("please provide chr")
    m <- Glist$mchr[chr]
    if (is.null(cls)) cls <- 1:m
    if (!is.null(rsids)) cls <- match(rsids, Glist$rsids[[chr]])
    if (any(is.na(cls))) {
      warning(paste("some rsids not found in Glist"))
      message(rsids[is.na(cls)])
    }
    if (!is.null(allele)) allele <- allele[!is.na(cls)]
    cls <- cls[!is.na(cls)]
    nc <- length(cls)
    if (is.null(allele)) direction <- rep(1, nc)
    if (!is.null(allele)) direction <- as.integer(allele == Glist$a2[[chr]][cls])
    if (is.null(rws)) rws <- 1:Glist$n
    if (!is.null(ids)) rws <- match(ids, Glist$ids)
    nr <- length(rws)
    ids <- Glist$ids[rws]
    rsids <- Glist$rsids[[chr]][cls]
  }

  if (!is.null(bedfiles)) {
    if (is.null(bimfiles)) bimfiles <- gsub(".bed", ".bim", bedfiles)
    if (is.null(famfiles)) famfiles <- gsub(".bed", ".fam", bedfiles)
    bim <- fread(
      input = bimfiles, header = FALSE, data.table = FALSE, showProgress = FALSE,
      colClasses = "character"
    )
    fam <- fread(
      input = famfiles, header = FALSE, data.table = FALSE, showProgress = FALSE,
      colClasses = "character"
    )
    n <- nrow(fam)
    m <- nrow(bim)
    if (is.null(cls)) cls <- 1:m
    if (!is.null(rsids)) cls <- match(rsids, as.character(bim[, 2]))
    if (sum(is.na(cls))==length(cls)) {
         stop(paste("no rsids found in bimfiles"))
    }
    if (any(is.na(cls))) {
      warning(paste("some rsids not found in bimfiles"))
      message(rsids[is.na(cls)])
    }
    if (!is.null(allele)) allele <- allele[!is.na(cls)]
    cls <- cls[!is.na(cls)]
    nc <- length(cls)
    if (is.null(allele)) direction <- rep(1, nc)
    if (!is.null(allele)) direction <- as.integer(allele == as.character(bim[cls, 6]))
    if (length(cls) == 0) stop("No rsids found in bimfiles")
    if (is.null(rws)) rws <- 1:n
    if (!is.null(ids)) rws <- match(ids, as.character(fam[, 2]))
    nr <- length(rws)
    ids <- as.character(fam[rws, 2])
    rsids <- as.character(bim[cls, 2])
  }

  W <- .Call("_qgg_readW", Glist$bedfiles[chr], Glist$n, cls, Glist$af[[chr]][cls])
  W <- W[rws,]     
  rownames(W) <- ids
  colnames(W) <- rsids
  return(W)
}




sparseLD <- function(Glist = NULL, fnLD = NULL, bedfiles = NULL, bimfiles = NULL, famfiles = NULL, msize = 100, chr = NULL, rsids = NULL, allele = NULL,
                     ids = NULL, ncores = 1) {

  if (file.exists(fnLD)) stop("LD file allready exists - please specify other file names")

  if(!is.null(bedfiles)) {
     if(is.null(bimfiles)) bimfiles <- gsub(".bed",".bim",bedfiles)
     if(is.null(famfiles)) famfiles <- gsub(".bed",".fam",bedfiles)

     Glist <- NULL
     Glist$fnBED <- bedfiles

     bim <- fread(input = bimfiles, header = FALSE, data.table = FALSE, colClasses = "character")
     Glist$a1 <- as.character(bim[, 5])
     Glist$a2 <- as.character(bim[, 6])
     Glist$position <- as.numeric(bim[, 4])
     Glist$rsids <- as.character(bim[, 2])
     Glist$chr <- as.character(bim[, 1])

     fam <- fread(input = famfiles, header = FALSE, data.table = FALSE, colClasses = "character")
     Glist$n <- nrow(fam)
     Glist$ids <- as.character(fam[, 2])
     if (any(!ids %in% as.character(fam[, 2]))) stop(paste("some ids not found in famfiles"))
  }

  n <- Glist$n
  rws <- 1:n
  if (!is.null(ids)) rws <- match(ids, Glist$ids)
  message(paste("Compute LD using individuals listed in ids"))
  nr <- length(rws)

  rsidsLD <- Glist$rsids[[chr]]
  if(!is.null(rsids)) rsidsLD <- rsidsLD[rsidsLD%in%rsids]
  Glist$rsidsLD[[chr]] <- rsidsLD
  cls <- match(rsidsLD, Glist$rsids[[chr]])
  nc <- length(cls)
  af <- rep(0.5, nc)
  af <- Glist$af[[chr]][cls]
  
  cls <- split(cls, ceiling(seq_along(cls) / msize))
  af <- split(af, ceiling(seq_along(af) / msize))
  msets <- sapply(cls, length)
  nsets <- length(msets)
  W1 <- matrix(0, nrow = nr, ncol = msize)
  W2 <- matrix(0, nrow = nr, ncol = msize)
  W3 <- matrix(0, nrow = nr, ncol = msize)
  bfLD <- file(fnLD, "wb")
  for (j in 1:nsets) {
    nc <- length(cls[[j]])
    W1 = W2
    W2 = W3
    W3 <- .Call("_qgg_readW", Glist$bedfiles[chr], Glist$n, cls[[j]], af[[j]])
    W3 <- scale(W3)
    if(j == nsets) W3 <- cbind(W3[rws,],matrix(0, nrow = nr, ncol = msize-nc))     
    LD <- t(crossprod(cbind(W1, W2, W3), W2))
    LD <- LD / (nr - 1) 
    LD[is.na(LD)] <- 0
    if (j > 1) {
      for (k in 1:msize) {
        ld <- as.vector(LD[k, k:(k + 2 * msize)])
        writeBin(ld, bfLD, size = 4, endian = "little")
      }
    }
    if (j == nsets) {
      W1 = W2
      W2 = W3
      W3 = matrix(0, nrow = nr, ncol = msize)
      LD <- t(crossprod(cbind(W1, W2, W3), W2))
      LD <- LD / (nr - 1) 
      LD[is.na(LD)] <- 0
      for (k in 1:msets[j]) {
        ld <- as.vector(LD[k, k:(k + 2 * msize)])
        writeBin(ld, bfLD, size = 4, endian = "little")
      }
    }
    message(paste("Finished block", j, "chromosome", chr))
  }
  close(bfLD)
  return(Glist)
}

getLDsets <- function(Glist = NULL, chr = NULL, r2 = 0.5) {
  msize <- Glist$msize
  rsidsChr <- Glist$rsidsLD[[chr]]
  mchr <- length(rsidsChr)
  rsidsLD <- c(rep("start", msize), rsidsChr, rep("end", msize))
  ldSetsChr <- vector(length = mchr, mode = "list")
  names(ldSetsChr) <- rsidsChr
  
  fnLD <- Glist$fnLD[chr]
  bfLD <- file(fnLD, "rb")
  
  nld <- as.integer(msize * 2 + 1)
  for (i in 1:mchr) {
    ld <- readBin(bfLD, "numeric", n = nld, size = 4, endian = "little")
    ld[msize + 1] <- 1
    cls <- which((ld**2) > r2) + i - 1
    ldSetsChr[[i]] <- rsidsLD[cls]
  }
  close(bfLD)
  
  return(ldSetsChr)
}

getLD <- function(Glist = NULL, chr = NULL, rsids=NULL) {
  msize <- Glist$msize
  rsidsChr <- Glist$rsidsLD[[chr]]
  mchr <- length(rsidsChr)
  ld = matrix(0, ncol = mchr, nrow=(msize * 2 + 1))
  colnames(ld) <- rsidsChr
  rownames(ld) <- c(-(msize:1), 0, 1:msize)
  fnLD <- Glist$fnLD[chr]
  bfLD <- file(fnLD, "rb")
  nld <- as.integer(msize * 2 + 1)
  k = 1
  for (i in 1:mchr) {
    ld[,i] = readBin(bfLD, "numeric", n = nld, size = 4, endian = "little")
    ld[msize + 1,i] <- 1
  }
  close(bfLD)
  return(ld)
}


