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
#' @param Glist list of information about genotype matrix stored on disk - only provided if task="summary" or task="sparseld"
#' @param task character specifying which task to perform ("prepare" is default, "summary", or "sparseld")
#' @param study name of the study
#' @param fnBED path and filename of the binary file .bed used for storing genotypes on the disk
#' @param bedfiles vector of names for the PLINK bed-files
#' @param famfiles vector of names for the PLINK fam-files
#' @param bimfiles vector of names for the PLINK bim-files
#' @param mapfiles vector of names for the mapfiles
#' @param ids vector of individuals used in the study
#' @param rsids vector of marker rsids used in the study
#' @param ldfiles path and filename of the binary files .ld for storing sparse ld matrix on the disk
#' @param msize number of markers used in compuation of sparseld
#' @param overwrite logical if TRUE overwite binary genotype/ld file
#' @param ncores number of cores used to process the genotypes
#' @param assembly character name of assembly
#' @param r2 threshold
#' @param kb size of genomic region in kb
#' @param cm size of genomic region in cm
#'
#' @return Returns a list structure (Glist) with information about genotypes
#'


#' @author Peter Soerensen

#' @examples
#'
#' bedfiles <- system.file("extdata", "sample_chr1.bed", package = "qgg")
#' bimfiles <- system.file("extdata", "sample_chr1.bim", package = "qgg")
#' famfiles <- system.file("extdata", "sample_chr1.fam", package = "qgg")
#' 
#' Glist <- gprep(study="Example", bedfiles=bedfiles, bimfiles=bimfiles,
#'              famfiles=famfiles)
#' 


#' @export
#'

gprep <- function(Glist = NULL, task = "prepare", study = NULL, fnBED = NULL, ldfiles = NULL,
                  bedfiles = NULL, bimfiles = NULL, famfiles = NULL, mapfiles=NULL, 
                  ids = NULL, rsids = NULL, assembly=NULL,
                  overwrite = FALSE, msize = 100, r2=NULL, kb=NULL, cm=NULL, ncores = 1) {

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
    fam <- data.table::fread(input = famfiles[1], header = FALSE, data.table = FALSE, colClasses = "character")
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
    Glist$cpra <- vector(mode = "list", length = nfiles)
    Glist$map <- vector(mode = "list", length = nfiles)
    
    Glist$nmiss <- vector(mode = "list", length = nfiles)
    Glist$af <- vector(mode = "list", length = nfiles)
    Glist$af1 <- vector(mode = "list", length = nfiles)
    Glist$af2 <- vector(mode = "list", length = nfiles)
    Glist$maf <- vector(mode = "list", length = nfiles)
    Glist$hom <- vector(mode = "list", length = nfiles)
    Glist$het <- vector(mode = "list", length = nfiles)
    Glist$n0 <- vector(mode = "list", length = nfiles)
    Glist$n1 <- vector(mode = "list", length = nfiles)
    Glist$n2 <- vector(mode = "list", length = nfiles)

    for (chr in 1:length(bedfiles)) {
      bim <- data.table::fread(input = bimfiles[chr], header = FALSE, data.table = FALSE, colClasses = "character")
      rsidsBIM <- as.character(bim[, 2])
      if (!is.null(rsids)) bim <- droplevels(bim[rsidsBIM %in% rsids, ])
      fam <- data.table::fread(input = famfiles[chr], header = FALSE, data.table = FALSE, colClasses = "character")
      if (any(!Glist$ids %in% as.character(fam[, 2]))) stop(paste("some ids not found in famfiles"))
      message(paste("Finished processing fam file", famfiles[chr]))
      Glist$a1[[chr]] <- as.character(bim[, 5])
      Glist$a2[[chr]] <- as.character(bim[, 6])
      Glist$position[[chr]] <- as.numeric(bim[, 4])
      Glist$map[[chr]] <- as.numeric(bim[, 3])
      Glist$rsids[[chr]] <- as.character(bim[, 2])
      Glist$mchr[chr] <- length(Glist$rsids[[chr]]) 
      Glist$chr[[chr]] <- as.character(bim[, 1])
      Glist$cpra[[chr]] <- paste(Glist$chr[[chr]],Glist$pos[[chr]],Glist$a1[[chr]],Glist$a2[[chr]],sep="_")
      message(paste("Finished processing bim file", bimfiles[chr]))
      #if (is.null(Glist$fnBED)) Glist <- summaryBED(Glist=Glist, chr=chr, ids = Glist$ids, ncores = ncores)
      if (is.null(Glist$fnBED)) Glist <- summaryBED(Glist=Glist, chr=chr, ids = Glist$study_ids, ncores = ncores)
      message(paste("Finished processing bed file", bedfiles[chr]))
      names(Glist$nmiss[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$af[[chr]]) <- Glist$rsids[[chr]]
      Glist$af1[[chr]] <- Glist$af[[chr]]
      Glist$af2[[chr]] <- 1-Glist$af[[chr]]
      names(Glist$maf[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$a1[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$a2[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$position[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$map[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$het[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$hom[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$n0[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$n1[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$n2[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$cpra[[chr]]) <- Glist$rsids[[chr]]
      names(Glist$rsids[[chr]]) <- Glist$cpra[[chr]]
      names(Glist$chr[[chr]]) <- Glist$rsids[[chr]]
    }

    Glist$nchr <- length(Glist$bedfiles)
    Glist$assembly <- assembly

  }
     # if (task == "combine") {
     #      writeBED(bedfiles = bedfiles, bimfiles = bimfiles, famfiles = famfiles, 
     #               fnBED=fnBED, ids = ids, rsids = rsids, overwrite = FALSE) {
     #           
     # }

  if (task == "sparseld") {
    message("Computing ld")
    Glist$msize <- msize
    Glist$ldfiles <- ldfiles
    if (is.null(ldfiles)) Glist$ldfiles <- gsub(".bed", ".ld", Glist$bedfiles)
    Glist$rsidsLD <- vector(mode = "list", length = length(Glist$ldfiles))
    Glist$ldscores <- vector(mode = "list", length = length(Glist$ldfiles))
    if (is.null(ids)) ids <- Glist$ids
    Glist$idsLD <- ids
    for( chr in 1:length(Glist$ldfiles)) {
      message(paste("Compute sparse LD matrix for chromosome:",chr))
      Glist <- sparseLD(Glist = Glist, fnLD = Glist$ldfiles[chr], msize = msize, chr = chr, rsids = rsids,
        ids = ids, ncores = 1, overwrite=overwrite)
      Glist$ldscores[[chr]] <- ldscore( Glist=Glist, chr=chr) 
    }
  }
  if (task == "ldsets") {
    if(is.null(r2)) stop("Please specify r2 threshold - can be a vector of values (e.g. r2=c(0.7,0.8,0.9) )")
    Glist$ldSets <- vector(length=length(r2), mode="list")
    names(Glist$ldSets) <- r2
    for (i in 1:length(r2) ) {
      ldSets <- vector(length=Glist$nchr, mode="list")
      for (chr in 1:Glist$nchr) {
        message(paste("Extract LD information for chromosome:", chr))
        ldSets[[chr]] <- getLDsets(Glist = Glist, r2 = 0.5, chr = chr)
      }
      Glist$ldSets[[i]] <- ldSets
    }
  }  
  
  if (task == "geneticmap") {
    message("Add Genetic map to Glist")
    if(is.null(Glist)) stop("Please provide Glist")
    if(is.null(mapfiles)) stop("Please provide mapfiles")
    Glist$map <- NULL
    for (chr in 1:Glist$nchr) {
      Glist$map[[chr]] <- rep(NA,Glist$mchr[chr])
      names(Glist$map[[chr]]) <- Glist$rsids[[chr]]
      map <- data.table::fread(mapfiles[chr],data.table=F)
      dups <- map[duplicated(map[,2]),2]
      map <- map[!map[,2]%in%dups,]
      map <- map[map[,2]%in%Glist$rsids[[chr]],]
      rownames(map) <- map[,2]
      Glist$map[[chr]][rownames(map)] <- map[,3]
      message(paste("Finished processing map file", mapfiles[chr]))
    }
    return(Glist)
  }
  if (task == "ldscores") {
    message("Computing ldscores")
    if(is.null(Glist)) stop("Please provide Glist")
    if(!is.null(cm)) message(paste("Computing ldscores using cm:",cm))
    if(!is.null(kb)) message(paste("Computing ldscores using kb:",kb))
    ldscores <- NULL
    for (chr in 1:22) { 
      ldscores[[chr]] <- ldscore(Glist=Glist, chr=chr, cm=cm) 
    }
    return(ldscores)
  }
  
  return(Glist)
}




summaryBED <- function(Glist = NULL, ids = NULL, rsids = NULL, rws = NULL, cls = NULL, chr = NULL, ncores = 1) {

  n <- Glist$n
  m <- Glist$mchr[chr]
  
  rws <- 1:n
  if (!is.null(ids)) rws <- match(ids, Glist$ids)
  if (!is.null(Glist$study_ids)) {
    warning("Study ids used in calculating genotype frequencies etc.")
  }
  nr <- length(rws)

  if (is.null(cls)) cls <- 1:m
  if (!is.null(rsids)) cls <- match(rsids, Glist$rsids)
  nc <- length(cls)

  ######################################################################
  #  00 01 10 11         bit level  corresponds to
  #  0  1  2  3          xij level  corresponds to
  #  2  NA  1  0         number of copies of first allele in bim file
  ######################################################################
  mask <- rep(0,n)
  mask[rws] <- 1
  freq <- .Call("_qgg_freqbed", Glist$bedfiles[chr], Glist$n, mask, cls)
  nmiss <- freq[2,]
  hom <- (freq[1,] + freq[4,]) / (freq[1,] + freq[3,] + freq[4,])
  het <- (freq[3,]) / (freq[1,] + freq[3,] + freq[4,])
  #af <- (2*freq[1,] + freq[3,])/(2*freq[1,] + 2*freq[3,] + 2*freq[4,])
  #nalleles <- 2*(n-nmiss)
  nalleles <- 2*(nr-nmiss)
  af <- 2*freq[1,] + freq[3,]
  af[nalleles>0] <- af[nalleles>0]/nalleles[nalleles>0]
  af[nalleles==0] <- 0.5
  tol_upper <- 0.99999
  tol_lower <- 0.00001
  af[af>tol_upper] <- tol_upper
  af[af<tol_lower] <- tol_lower

  
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

#'
#' Quality control of marker summary statistics
#'
#' @description
#' Quality control is a critical step for working with summary statistics (in particular
#'                                                                         for external). 
#' Processing and quality control of GWAS summary statistics includes:                                                                      
#'
#' - map marker ids (rsids/cpra (chr, pos, ref, alt)) to LD reference panel data 
#' - check effect allele (flip EA, EAF, Effect)
#' - check effect allele frequency
#' - thresholds for MAF and HWE
#' - exclude INDELS, CG/AT and MHC region
#' - remove duplicated marker ids
#' - check which build version
#' - check for concordance between marker effect and LD data
#'
#' External summary statistics format:
#'  marker, chr, pos, effect_allele, non_effect_allele, effect_allele_freq, effect, effect_se, stat, p, n    
#' 
#' Internal summary statistics format:
#'  rsids, chr, pos, a1, a2, af, b, seb, stat, p, n
#' 
#'
#' @param Glist list of information about genotype matrix stored on disk
#' @param excludeMAF exclude marker if minor allele frequency (MAF) is below threshold (0.01 is default)
#' @param excludeINFO exclude marker if info score (INFO) is below threshold (0.8 is default)
#' @param excludeMISS exclude marker if missingness (MISS) is above threshold (0.05 is default)
#' @param excludeHWE exclude marker if p-value for Hardy Weinberg Equilibrium test is below threshold (0.01 is default)
#' @param excludeCGAT exclude marker if alleles are ambigous (CG or AT)
#' @param excludeMHC exclude marker if located in MHC region 
#' @param excludeINDEL exclude marker if it an insertion/deletion  
#' @param excludeDUPS exclude marker id if duplicated
#' @param assembly character name of assembly


#' @author Peter Soerensen


#'
#' @export
#'



gfilter <- function(Glist = NULL, excludeMAF=0.01, excludeMISS=0.05, excludeINFO=NULL,excludeCGAT=TRUE,
                    excludeINDEL=TRUE, excludeDUPS=TRUE, excludeHWE=1e-12, excludeMHC=FALSE, assembly="GRCh37") {
# excludeINFO is a numeric value of the info score used for filtering
  rsids <- unlist(Glist$rsids)
  if(is.null(Glist$study_ids)) Glist$study_ids <- Glist$ids
  if(!is.null(excludeMAF)) isMAF <- unlist(lapply(Glist$maf, function(x){x<=excludeMAF}))
  if(!is.null(excludeMISS)) isMISS <- unlist(lapply(Glist$nmiss,function(x) {x/length(Glist$study_ids)>excludeMISS}))
  if(!is.null(excludeHWE)) isHWE <- unlist(hwe(Glist)) < excludeHWE 
  
  
  isHWE[is.na(isHWE)] <- TRUE
  if(excludeMHC) {
    if(assembly=="GRCh37"){
        isMHC <-  Glist$position[[6]] > 28477797 & Glist$position[[6]] < 33448354
        rsidsMHC <- names(isMHC)[isMHC]
    }
    if(assembly=="GRCh38"){
        isMHC <-  Glist$position[[6]] > 28510120 & Glist$position[[6]] < 33480577
        rsidsMHC <- names(isMHC)[isMHC]
    }
  }
  a1 <- unlist(Glist$a1)
  a2 <- unlist(Glist$a2)
  isAT <- a1=="A" & a2=="T"
  isTA <- a1=="T" & a2=="A"
  isCG <- a1=="C" & a2=="G"
  isGC <- a1=="G" & a2=="C"
  isCGAT <- isAT | isTA | isCG | isGC
  CGTA <- c("C","G","T","A")
  isINDEL <- !((a1%in%CGTA) & (a2%in%CGTA))

  message(paste("Number of markers excluded by low MAF:", sum(isMAF)))
  message(paste("Number of markers excluded by deviation from HWE:", sum(isHWE)))
  message(paste("Number of markers excluded by missingnes:", sum(isMISS)))

  rsidsQC <- isMAF | isMISS | isHWE

  if(excludeCGAT) {
    rsidsQC <- rsidsQC | isCGAT
    message(paste("Number of markers excluded by ambiguity (CG or AT):", sum(isCGAT)))
  }
  if(excludeINDEL) {
    rsidsQC <- rsidsQC | isINDEL
    message(paste("Number of markers excluded by being INDEL:", sum(isINDEL)))
  }
  
  if(excludeDUPS) {
    rsidsDUPS <- rsids[duplicated(rsids)]
    isDUPS <- rsids%in%rsidsDUPS
    rsidsQC <- rsidsQC | isDUPS
    message(paste("Number of markers excluded by duplicated rsids", sum(isDUPS)))
  }
  
  if(!is.null(excludeINFO)){
    rsidsINFO <- unlist(Glist$info)
    rsidsINFO <- rsidsINFO[rsidsINFO>excludeINFO]
    exINFO <- length(unlist(Glist$info)) - length(rsidsINFO)
    rsidsQC <- rsidsQC[names(rsidsQC)%in%names(rsidsINFO)]
    message(paste("Number of markers excluded by info score", exINFO))
  }

  #if(excludeCG_AT) rsidsQC <- isMAF | isMISS | isHWE | isCGAT
  rsidsQC <- names(rsidsQC)[!rsidsQC]
  if(excludeMHC) {
    rsidsQC <- rsidsQC[!rsidsQC%in%rsidsMHC]
    message(paste("Number of markers excluded in MHC region:", length(rsidsMHC)))
  }
  rsidsQC <- unlist(lapply(Glist$rsids,function(x){x[x%in%rsidsQC]}))
  message(paste("Number of markers excluded:", length(rsids)-length(rsidsQC)))
  message(paste("Number of markers retained:", length(rsidsQC)))
  return(unlist(rsidsQC))  
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
#           bim <- data.table::fread(input = bimfiles[chr], header = FALSE, data.table = FALSE, colClasses = "character")
#           fam <- data.table::fread(input = famfiles[chr], header = FALSE, data.table = FALSE)
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


#' Extract elements from genotype matrix stored on disk
#'
#' @description
#' Extract elements from genotype matrix (rows/columns, ids/rsids) stored on disk.

#' @param Glist list structure with information about genotypes stored on disk
#' @param bedfiles vector of name for the PLINK bed-file
#' @param bimfiles vector of name for the PLINK bim-file
#' @param famfiles vector of name for the PLINK fam-file
#' @param ids vector of ids in W to be extracted
#' @param chr chromosome for which W is to be extracted
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
if(!file.exists(bedfiles)) stop("Glist$bedfiles[chr] does not exist")
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
if (!is.null(ids)) {
  rws <- match(ids,Glist$ids)
  if(any(is.na(rws))) stop("Some ids not found in Glist")
}
if(!is.null(rws)) W <- W[rws,, drop = FALSE]
if(is.integer(impute)) W[W==impute] <- impute
return(W)
}



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
    if (!is.null(allele)) direction <- as.integer(allele == Glist$a1[[chr]][cls])
    if (is.null(rws)) rws <- 1:Glist$n
    if (!is.null(ids)) rws <- match(ids, Glist$ids)
    nr <- length(rws)
    ids <- Glist$ids[rws]
    rsids <- Glist$rsids[[chr]][cls]
  }

  if (!is.null(bedfiles)) {
    if (is.null(bimfiles)) bimfiles <- gsub(".bed", ".bim", bedfiles)
    if (is.null(famfiles)) famfiles <- gsub(".bed", ".fam", bedfiles)
    bim <- data.table::fread(
      input = bimfiles, header = FALSE, data.table = FALSE, showProgress = FALSE,
      colClasses = "character"
    )
    fam <- data.table::fread(
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
                     ids = NULL, ncores = 1, overwrite=FALSE) {

  if (!overwrite && file.exists(fnLD) && !interactive()) stop("LD file allready exists - please specify other file names")
  if (overwrite && file.exists(fnLD)) warning("LD file allready exists - replacing existing file")
  if (!overwrite && file.exists(fnLD) && interactive()) askYesNo(paste("LD file",fnLD,"allready exists - do you want to replace existing file?"))
  if(!is.null(bedfiles)) {
     if(is.null(bimfiles)) bimfiles <- gsub(".bed",".bim",bedfiles)
     if(is.null(famfiles)) famfiles <- gsub(".bed",".fam",bedfiles)

     Glist <- NULL
     Glist$fnBED <- bedfiles

     bim <- data.table::fread(input = bimfiles, header = FALSE, data.table = FALSE, colClasses = "character")
     Glist$a1 <- as.character(bim[, 5])
     Glist$a2 <- as.character(bim[, 6])
     Glist$position <- as.numeric(bim[, 4])
     Glist$rsids <- as.character(bim[, 2])
     Glist$chr <- as.character(bim[, 1])

     fam <- data.table::fread(input = famfiles, header = FALSE, data.table = FALSE, colClasses = "character")
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
  if(length(rsidsLD)<msize) stop("LD marker window size is to big - use smaller number for msize")
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
    W3 <- scale(W3[rws,])
    #if(j == nsets) W3 <- cbind(W3[rws,],matrix(0, nrow = nr, ncol = msize-nc))     
    if(j == nsets) W3 <- cbind(W3,matrix(0, nrow = nr, ncol = msize-nc))     
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

#' Get marker LD sets
#' @description
#' Get marker LD sets based on sparse LD matrix in Glist
#' @keywords internal
#' @param Glist list structure with information about genotypes stored on disk
#' @param chr chromosome for which LD sets are extracted
#' @param r2 threshold used for extracting LD sets
#' @export

getLDsets <- function(Glist = NULL, chr = NULL, r2 = 0.5) {
  if(!is.null(chr)) {
    msize <- Glist$msize
    rsidsChr <- Glist$rsidsLD[[chr]]
    mchr <- length(rsidsChr)
    rsidsLD <- c(rep("start", msize), rsidsChr, rep("end", msize))
    ldSetsChr <- vector(length = mchr, mode = "list")
    names(ldSetsChr) <- rsidsChr
    
    fnLD <- Glist$ldfiles[chr]
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
  
}

#' Get sparse LD matrix
#' @description
#' Get sparse LD matrix in Glist
#' @keywords internal
#' @param Glist list structure with information about genotypes stored on disk
#' @param chr chromosome from which LD sets are extracted
#' @param rsids included in the sparse LD matrix extracted
#' @export

getLD <- function(Glist = NULL, chr = NULL, rsids=NULL) {
  msize <- Glist$msize
  rsidsChr <- Glist$rsidsLD[[chr]]
  mchr <- length(rsidsChr)
  ld = matrix(0, ncol = mchr, nrow=(msize * 2 + 1))
  colnames(ld) <- rsidsChr
  rownames(ld) <- c(-(msize:1), 0, 1:msize)
  fnLD <- Glist$ldfiles[chr]
  bfLD <- file(fnLD, "rb")
  nld <- as.integer(msize * 2 + 1)
  k = 1
  for (i in 1:mchr) {
    ld[,i] = readBin(bfLD, "numeric", n = nld, size = 4, endian = "little")
    ld[msize + 1,i] <- 1
  }
  close(bfLD)
  if(!is.null(rsids)) ld <- ld[,colnames(ld)%in%rsids]
  return(ld)
}

getSparseLD <- function(Glist = NULL, chr = NULL, r2 = 0, onebased=TRUE, rsids=NULL, format="sparse") {
  msize <- Glist$msize
  rsidsChr <- Glist$rsidsLD[[chr]]
  mchr <- length(rsidsChr)
  rsidsLD <- c(rep("start", msize), rsidsChr, rep("end", msize))
  if(onebased) rsids_indices <- c(rep(0, msize), 1:mchr, rep(0, msize))
  if(!onebased) rsids_indices <- c(rep(0, msize), 0:(mchr-1), rep(0, msize))
  ld_indices <- vector(length = mchr, mode = "list")
  ld_values <- vector(length = mchr, mode = "list")
  names(ld_indices) <- names(ld_values) <- rsidsChr
  
  fnLD <- Glist$ldfiles[chr]
  bfLD <- file(fnLD, "rb")
  
  nld <- as.integer(msize * 2 + 1)
  for (i in 1:mchr) {
    ld <- readBin(bfLD, "numeric", n = nld, size = 4, endian = "little")
    ld[msize + 1] <- 1
    cls <- which((ld**2) > r2) + i - 1
    ld_indices[[i]] <- rsids_indices[cls]
    ld_values[[i]] <- ld[(ld**2) > r2]
  }
  close(bfLD)
  if(format=="sparse") return(list(indices=ld_indices,values=ld_values))
  if(format=="dense") {
    if(is.null(rsids)) {
      LD <- matrix(0,nrow=length(ld_values),ncol=length(ld_values), 
                   dimnames=list(names(ld_values),names(ld_values) ) )  
      for(i in 1:length(ld_values)) {
        LD[i,ld_indices[[i]]] <- ld_values[[i]]
      }
    }
    if(!is.null(rsids)) {
      LD <- matrix(0,nrow=length(rsids),ncol=length(rsids), 
                   dimnames=list(rsids,rsids) )
      for(i in 1:length(rsids)) {
        ldrsids <-rsidsChr[ld_indices[[rsids[i]]]]
        ldvals <- ld_values[[rsids[i]]]
        inlist <- ldrsids%in%rsids 
        LD[rsids[i],ldrsids[inlist]] <- ldvals[inlist]
      }
    }
    return(LD)
  }
}


plotLD <- function(LD=NULL, cols=NULL) {
  if(is.null(cols)) cols <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)
  LD <- LD[,ncol(LD):1]
  image(LD, col=cols, axes = FALSE)
  axis(1, at = seq(0, 1, length = nrow(LD)), labels = rownames(LD), las=2, cex.axis=0.5)
  axis(2, at = seq(0, 1, length = ncol(LD)), labels = colnames(LD), las=2, cex.axis=0.5)
}



#' Get marker rsids in a genome region
#' @description
#' Get marker rsids in a specific genome region based on marker in Glist
#' @keywords internal
#' @param Glist list structure with information about genotypes stored on disk
#' @param chr chromosome from which markers are extracted
#' @param region genome region (in base pairs) from which markers are extracted
#' @export

getMarkers <- function(Glist=NULL, chr=NULL, region=NULL) {
  minpos <- min(region)
  maxpos <- max(region)
  select <-  Glist$position[[chr]] > minpos & Glist$position[[chr]] < maxpos
  Glist$rsids[[chr]][select]
}
