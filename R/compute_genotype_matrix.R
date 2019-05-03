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
#' binary file on disk (Glist).
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
#' The Glist structure is used as input parameter for a number of qgg core functions including:
#' 1) construction of genomic relationship matrices (grm), 2) estimating genomic parameters (greml),
#' 3) single marker association analyses (lma or mlma), 4) gene set enrichment analyses (gsea),
#' and 5) genomic prediction from genotypes and phenotypes (gsolve) or genotypes and summary statistics (gscore).
#'

#'
#' @param study name of the study
#' @param fnRAW path and .raw filename of the binary file used for storing genotypes on the disk
#' @param bedfiles vector of names for the PLINK bed-files
#' @param famfiles vector of names for the PLINK fam-files
#' @param bimfiles vector of names for the PLINK bim-files
#' @param ids vector of individuals used in the study
#' @param rsids vector of marker rsids used in the study
#' @param overwrite logical if TRUE overwite binary genotype file
#' @param ncores number of cores used to process the genotypes
#'
#' @return Returns a list structure with information about genotypes
#'


#' @author Peter Soerensen

#' @examples
#

#' #Glist <- gprep( bedfiles, bimfiles, study, path, additional arguments...)
#' #W <- getW( Glist, ids, rsids, additional arguments...)


#' @export
#'

gprep <- function(Glist = NULL, task = "prepare", study = NULL, fnRAW = NULL, fnLD = NULL,
                  bedfiles = NULL, bimfiles = NULL, famfiles = NULL, ids = NULL, rsids = NULL,
                  overwrite = FALSE, msize = 100, ncores = 1) {
  if (task == "prepare") {
    nfiles <- length(bedfiles)
    Glist <- NULL
    Glist$study <- study
    Glist$fnRAW <- fnRAW
    if (file.exists(fnRAW)) warning(paste("fnRAW allready exist"))
    Glist$bedfiles <- bedfiles
    if (!is.null(bimfiles)) Glist$bimfiles <- bimfiles
    if (!is.null(famfiles)) Glist$famfiles <- famfiles
    if (is.null(bimfiles)) Glist$bimfiles <- gsub(".bed", ".bim", bedfiles)
    if (is.null(famfiles)) Glist$famfiles <- gsub(".bed", ".fam", bedfiles)

    # Read fam information
    fam <- fread(input = famfiles[1], header = FALSE, data.table = FALSE, colClasses = "character")
    Glist$ids <- as.character(fam[, 2])
    Glist$study_ids <- Glist$ids
    if (!is.null(ids)) {
      if (any(!ids %in% as.character(fam[, 2]))) warning(paste("some ids not found in famfiles"))
      Glist$study_ids <- Glist$ids[Glist$ids %in% as.character(ids)]
    }
    if (any(duplicated(Glist$ids))) stop("Duplicated ids found in famfiles")

    Glist$rsids <- vector(mode = "list", length = nfiles)
    Glist$a1 <- vector(mode = "list", length = nfiles)
    Glist$a2 <- vector(mode = "list", length = nfiles)
    Glist$position <- vector(mode = "list", length = nfiles)
    Glist$chr <- vector(mode = "list", length = nfiles)

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
      Glist$chr[[chr]] <- as.character(bim[, 1])
      print(paste("Finished processing bim file", bimfiles[chr]))
    }
    Glist$study_rsids <- unlist(Glist$rsids)
    if (!is.null(rsids)) Glist$study_rsids <- as.character(rsids)
    if (!is.null(rsids)) if (any(!rsids %in% unlist(Glist$rsids))) warning(paste("some rsids not found in bimfiles"))
    if (!is.null(rsids)) Glist$study_rsids <- unlist(Glist$rsids)[unlist(Glist$rsids) %in% rsids]

    Glist$mchr <- sapply(Glist$rsids, length)
    Glist$rsids <- unlist(Glist$rsids)
    Glist$a1 <- unlist(Glist$a1)
    Glist$a2 <- unlist(Glist$a2)
    Glist$position <- unlist(Glist$position)
    Glist$chr <- unlist(Glist$chr)
    Glist$nchr <- length(unique(Glist$chr))

    Glist$n <- length(Glist$ids)
    Glist$m <- length(Glist$rsids)

    print("Preparing raw file")
    qgg:::writeRAW(Glist = Glist, ids = ids, rsids = Glist$rsids, overwrite = overwrite) # write genotypes to .raw file

    print("Computing allele frequencies, missingness")
    Glist <- qgg:::summaryRAW(Glist = Glist, ids = ids, rsids = Glist$rsids, ncores = ncores)
    Glist$af1 <- 1 - Glist$af
    Glist$af2 <- Glist$af
  }

  if (task == "summary") {
    print("Computing allele frequencies, missingness")
    Glist <- qgg:::summaryRAW(Glist = Glist, ids = ids, rsids = Glist$rsids, ncores = ncores)
    Glist$af1 <- 1 - Glist$af
    Glist$af2 <- Glist$af
  }

  if (task == "sparseld") {
    print("Computing ld")
    Glist$msize <- msize
    Glist$fnLD <- fnLD
    if (is.null(fnLD)) Glist$fnLD <- gsub(".bed", ".ld", bedfiles)
    if (is.null(ids)) ids <- Glist$ids
    mclapply(1:length(Glist$fnLD), function(x) {
      qgg:::sparseLD(
        Glist = Glist, fnLD = Glist$fnLD[x], msize = msize, chr = x, rsids = NULL,
        impute = TRUE, scale = TRUE, ids = ids, ncores = 1
      )
    }
    ,
    mc.cores = ncores
    )
  }

  return(Glist)
}



writeRAW <- function(Glist = NULL, ids = NULL, rsids = NULL, overwrite = FALSE) {
  qgg:::bed2raw(
    fnRAW = Glist$fnRAW, bedfiles = Glist$bedfiles, bimfiles = Glist$bimfiles,
    famfiles = Glist$famfiles, ids = ids, rsids = rsids, overwrite = overwrite
  )
}


bed2raw <- function(fnRAW = NULL, bedfiles = NULL, bimfiles = NULL, famfiles = NULL,
                    ids = NULL, rsids = NULL, overwrite = FALSE) {
  if (file.exists(fnRAW)) {
    warning(paste("fnRAW file allready exist"))
    if (!overwrite) stop(paste("fnRAW file allready exist"))
  }
  for (chr in 1:length(bedfiles)) {
    print(paste("Processing bedfile:", bedfiles[chr]))
    bim <- fread(input = bimfiles[chr], header = FALSE, data.table = FALSE, colClasses = "character")
    fam <- fread(input = famfiles[chr], header = FALSE, data.table = FALSE)
    n <- nrow(fam)
    m <- nrow(bim)
    rsidsBIM <- as.character(bim[, 2])
    keep <- rep(TRUE, m)
    if (!is.null(rsids)) keep <- rsidsBIM %in% rsids
    nbytes <- ceiling(n / 4)
    fnBED <- bedfiles[chr]
    bfBED <- file(fnBED, "rb")
    magic <- readBin(bfBED, "raw", n = 3)
    if (!all(magic[1] == "6c", magic[2] == "1b", magic[3] == "01")) {
      stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
    }
    close(bfBED)
    append <- 1
    if (chr == 1) append <- 0
    res <- .Fortran("bed2raw",
      m = as.integer(m),
      cls = as.integer(keep),
      nbytes = as.integer(nbytes),
      append = as.integer(append),
      fnBED = as.character(fnBED),
      fnRAW = as.character(fnRAW),
      ncharbed = nchar(as.character(fnBED)),
      ncharraw = nchar(as.character(fnRAW)),
      PACKAGE = "qgg"
    )
    print(paste("Finished processing bedfile:", bedfiles[chr]))
  }
}



summaryRAW <- function(Glist = NULL, ids = NULL, rsids = NULL, rws = NULL, cls = NULL, ncores = 1) {
  n <- Glist$n
  m <- Glist$m
  nbytes <- ceiling(n / 4)
  rws <- 1:n
  if (!is.null(ids)) rws <- match(ids, Glist$ids)
  if (!is.null(Glist$study_ids)) rws <- match(Glist$study_ids, Glist$ids)
  nr <- length(rws)
  fnRAW <- Glist$fnRAW
  OS <- .Platform$OS.type
  if (OS == "windows") fnRAW <- tolower(gsub("/", "\\", fnRAW, fixed = T))

  if (is.null(cls)) cls <- 1:m
  if (!is.null(rsids)) cls <- match(rsids, Glist$rsids)

  nc <- length(cls)
  af <- nmiss <- n0 <- n1 <- n2 <- rep(0, nc)
  
  qc <- .Fortran("summarybed",
    n = as.integer(n),
    nr = as.integer(nr),
    rws = as.integer(rws),
    nc = as.integer(nc),
    cls = as.integer(cls),
    af = as.double(af),
    nmiss = as.double(nmiss),
    n0 = as.double(n0),
    n1 = as.double(n1),
    n2 = as.double(n2),
    nbytes = as.integer(nbytes),
    fnRAW = as.character(fnRAW),
    nchars = nchar(as.character(fnRAW)),
    ncores = as.integer(ncores),
    PACKAGE = "qgg"
  )

  qc$hom <- (qc$n0 + qc$n2) / (qc$nr - qc$nmiss)
  qc$het <- qc$n1 / (qc$nr - qc$nmiss)
  qc$maf <- qc$af
  qc$maf[qc$maf > 0.5] <- 1 - qc$maf[qc$maf > 0.5]
  rsids <- Glist$rsids[cls]

  Glist$nmiss <- qc$nmiss
  Glist$af <- qc$af
  Glist$maf <- qc$maf
  Glist$hom <- qc$hom
  Glist$het <- qc$het
  Glist$n0 <- qc$n0
  Glist$n1 <- qc$n1
  Glist$n2 <- qc$n2

  return(Glist)
}



#' @export
#'

readbed <- function(Glist = NULL, bedfiles = NULL, ids = NULL, rsids = NULL,
                    rws = NULL, cls = NULL, impute = TRUE, scale = FALSE,
                    allele = NULL, ncores = 1) {
  if (!is.null(Glist)) {
    n <- Glist$n
    m <- Glist$m
    nbytes <- ceiling(n / 4)
    if (is.null(cls)) cls <- 1:m
    if (!is.null(rsids)) cls <- match(rsids, Glist$rsids)
    if (any(is.na(cls))) {
      warning(paste("some rsids not found in Glist"))
      print(rsids[is.na(cls)])
    }
    if (!is.null(allele)) allele <- allele[!is.na(cls)]
    cls <- cls[!is.na(cls)]
    nc <- length(cls)
    if (is.null(allele)) direction <- rep(1, nc)
    if (!is.null(allele)) direction <- as.integer(allele == Glist$a2[cls])
    if (is.null(rws)) rws <- 1:n
    if (!is.null(ids)) rws <- match(ids, Glist$ids)
    nr <- length(rws)
    fnRAW <- Glist$fnRAW
    OS <- .Platform$OS.type
    if (OS == "windows") fnRAW <- tolower(gsub("/", "\\", fnRAW, fixed = T))
    ids <- Glist$ids[rws]
    rsids <- Glist$rsids[cls]
  }
  if (!is.null(bedfiles)) {
    fnRAW <- bedfiles
    bimfiles <- gsub(".bed", ".bim", bedfiles)
    famfiles <- gsub(".bed", ".fam", bedfiles)
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
    nbytes <- ceiling(n / 4)
    cls <- match(rsids, as.character(bim[, 2]))
    if (any(is.na(cls))) {
      warning(paste("some rsids not found in bimfiles"))
      print(rsids[is.na(cls)])
    }
    if (!is.null(allele)) allele <- allele[!is.na(cls)]
    cls <- cls[!is.na(cls)]
    nc <- length(cls)
    if (is.null(allele)) direction <- rep(1, nc)
    if (!is.null(allele)) direction <- as.integer(allele == Glist$a2[cls])
    if (length(cls) == 0) stop("No rsids found in bimfiles")
    if (is.null(rws)) rws <- 1:n
    if (!is.null(ids)) rws <- match(ids, as.character(fam[, 2]))
    nr <- length(rws)
    OS <- .Platform$OS.type
    if (OS == "windows") fnRAW <- tolower(gsub("/", "\\", fnRAW, fixed = T))
    ids <- as.character(fam[rws, 2])
    rsids <- as.character(bim[cls, 2])
  }
  W <- .Fortran("readbed",
    n = as.integer(n),
    nr = as.integer(nr),
    rws = as.integer(rws),
    nc = as.integer(nc),
    cls = as.integer(cls),
    impute = as.integer(impute),
    scale = as.integer(scale),
    direction = as.integer(direction),
    W = matrix(as.double(0), nrow = nr, ncol = nc),
    nbytes = as.integer(nbytes),
    fnRAW = as.character(fnRAW),
    nchars = nchar(as.character(fnRAW)),
    PACKAGE = "qgg"
  )$W
  rownames(W) <- ids
  colnames(W) <- rsids
  return(W)
}


#' @export
#'
getW <- function(Glist = NULL, ids = NULL, rsids = NULL, rws = NULL, cls = NULL,
                 impute = FALSE, scale = FALSE, allele = NULL) {
  if (is.null(ids)) ids <- Glist$ids
  if (is.null(cls)) cls <- match(rsids, Glist$rsids)
  W <- qgg::readbed(
    Glist = Glist, ids = ids, rsids = rsids, rws = rws, cls = cls,
    impute = impute, scale = scale, allele = allele
  )
  return(W)
}


#' @export
#'

sparseLD <- function(Glist = NULL, fnLD = NULL, msize = 100, chr = NULL, rsids = NULL,
                     impute = TRUE, scale = TRUE, ids = NULL, ncores = 1) {
  if (file.exists(fnLD)) stop("LD file allready exists - please specify other file names")
  n <- Glist$n
  rws <- 1:n
  if (!is.null(ids)) rws <- match(ids, Glist$ids)
  print(paste("Compute LD using individuals listed in ids"))
  nr <- length(rws)
  nbytes <- ceiling(n / 4)
  if (!is.null(chr)) rsids <- Glist$rsids[Glist$chr == chr]
  cls <- match(rsids, Glist$rsids)
  nc <- length(cls)
  cls <- split(cls, ceiling(seq_along(cls) / msize))
  msets <- sapply(cls, length)
  nsets <- length(msets)
  m <- length(rsids)
  W1 <- matrix(logical(0), nrow = nr, ncol = msize)
  W2 <- matrix(logical(0), nrow = nr, ncol = msize)
  W3 <- matrix(logical(0), nrow = nr, ncol = msize)
  fnRAW <- Glist$fnRAW
  bfLD <- file(fnLD, "wb")
  for (j in 1:nsets) {
    nc <- length(cls[[j]])
    direction <- rep(1, nc)
    W1 <- W2
    W2 <- W3
    W3[, 1:nc] <- .Fortran("readbed",
      n = as.integer(n),
      nr = as.integer(nr),
      rws = as.integer(rws),
      nc = as.integer(nc),
      cls = as.integer(cls[[j]]),
      impute = as.integer(impute),
      scale = as.integer(scale),
      direction = as.integer(direction),
      W = matrix(as.double(0), nrow = nr, ncol = nc),
      nbytes = as.integer(nbytes),
      fnRAW = as.character(fnRAW),
      nchars = nchar(as.character(fnRAW)),
      PACKAGE = "qgg"
    )$W
    LD <- t(crossprod(cbind(W1, W2, W3), W2))
    LD <- LD / (nr - 1) # nr-1 accounts for sample mean is estimated
    LD[is.na(LD)] <- 0
    if (j > 1) {
      for (k in 1:msize) {
        ld <- as.vector(LD[k, k:(k + 2 * msize)])
        writeBin(ld, bfLD, size = 8, endian = "little")
      }
    }
    if (j == nsets) {
      W1 <- W2
      W2 <- W3
      W3 <- matrix(0, nrow = nr, ncol = msize)
      LD <- t(crossprod(cbind(W1, W2, W3), W2))
      LD <- LD / (nr - 1) # nr-1 accounts for sample mean is estimated
      LD[is.na(LD)] <- 0
      for (k in 1:msets[j]) {
        ld <- as.vector(LD[k, k:(k + 2 * msize)])
        writeBin(ld, bfLD, size = 8, endian = "little")
      }
    }
    print(paste("Finished block", j, "chromosome", chr))
  }
  close(bfLD)
}


#' @export
#'

getLDCHRsets <- function(Glist = NULL, chr = NULL, r2 = 0.5) {
  msize <- Glist$msize
  ldSets <- NULL
  n <- Glist$n
  rsidsChr <- Glist$rsids[Glist$chr == chr]
  if (!is.null(Glist$study_rsids)) rsidsChr <- rsidsChr[rsidsChr %in% Glist$study_rsids]
  mchr <- length(rsidsChr)
  rsidsLD <- c(rep("start", msize), rsidsChr, rep("end", msize))
  fnLD <- Glist$fnLD[chr]
  bfLD <- file(fnLD, "rb")
  nld <- as.integer(mchr * (msize * 2 + 1))
  ld <- readBin(bfLD, "double", n = nld, size = 8, endian = "little")
  ld <- matrix(ld, nrow = mchr, byrow = TRUE)
  close(bfLD)
  ld[, msize + 1] <- 1
  ldSetsChr <- vector(length = mchr, mode = "list")
  names(ldSetsChr) <- rsidsChr
  for (i in 1:mchr) {
    cls <- which((ld[i, ]**2) > r2) + i - 1
    ldSetsChr[[i]] <- rsidsLD[cls]
  }
  print(paste("Finished chromosome", chr))
  return(ldSetsChr)
}


#' @export
#'

mapLDsets <- function(ldSets = NULL, rsids = NULL, Glist = NULL, index = TRUE) {
  mpsets <- NULL
  if (!is.null(Glist)) rsids <- unlist(Glist$rsids)
  for (chr in 1:length(ldSets)) {
    rsSets <- ldSets[[chr]]
    rsidsChr <- Glist$rsids[Glist$chr == chr]
    rsidsChr <- rsidsChr[rsidsChr %in% names(rsSets)]
    nsets <- sapply(rsSets, length)
    rsChr <- rep(names(rsSets), times = nsets)
    rsSets <- unlist(rsSets, use.names = FALSE)
    rsSets <- match(rsSets, rsids)
    inW <- !is.na(rsSets)
    rsSets <- rsSets[inW]
    if (!index) rsSets <- rsids[rsSets]
    rsChr <- rsChr[inW]
    rsChr <- factor(rsChr, levels = unique(rsChr))
    rsSets <- split(rsSets, f = rsChr)
    mpsets[[chr]] <- rsSets[rsidsChr]
  }
  return(mpsets)
}



#' @export
#'

getLD <- function(Glist = NULL, chr = NULL) {
  msize <- Glist$msize
  rsidsChr <- Glist$rsids[Glist$chr == chr]
  if (!is.null(Glist$study_rsids)) rsidsChr <- rsidsChr[rsidsChr %in% Glist$study_rsids]
  mchr <- length(rsidsChr)
  fnLD <- Glist$fnLD[chr]
  bfLD <- file(fnLD, "rb")
  nld <- as.integer(mchr * (msize * 2 + 1))
  ld <- readBin(bfLD, "double", n = nld, size = 8, endian = "little")
  ld <- matrix(ld, nrow = mchr, byrow = TRUE)
  close(bfLD)
  ld[, msize + 1] <- 1
  rownames(ld) <- rsidsChr
  colnames(ld) <- c(-(msize:1), 0, 1:msize)
  return(ld)
}
