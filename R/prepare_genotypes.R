###############################################################################################
#  Compute W (processed genotypes) for target population
###############################################################################################
#
# Description:
#
# The processed genotypes for a particular target population are stored in a matrix W. This 
# matrix is used for computing genomic relationship matrices used for estimation of variance 
# components using REML methods (GREML), linear (mixed) model marker association analyses (LMMA)
# and prediction of polygenic risk scores (PRS).
#   
# By default genotypes are allele counts of alternative allele which are centered and scaled. 
# Below are functions to compute and access the processed genotypes in W as implemented 
# in the qgg package
#
# The genotype matrix is computed using the computeW function based on the raw genotypes. 
# Currently the raw genotype format supported is plink bed/bim/fam files. Other formats such as
# bgen will be supported later. To fascilitate a simple and common interface for using raw genotypes 
# we prepare a Wlist structure which contains information about the raw and processed genotypes. 
# Wlist is obtained using the prepW function based on plink bed/bim/fam files. 
#
# The genotype matrix, W, is too big to fit into memory and is therefore stored in a binary file.
# By default the W matrix is stored in a binary file which is compressed by a ratio
# of 25-30 such that the resulting file size is approximately n*m*8/25 bytes. 
#
# Information about the processed genotype matrix W is provided in a Wlist structure. 
# The W matrix can be accessed using the getW function. 
#
# Examples of simple usage:
#
# 	Wlist <- prepW( bedfiles, bimfiles, study, path, additional arguments...)
# 	Wlist <- computeW( Wlist, chr, additional arguments...)
# 	Wlist <- computeW( bedfiles, bimfiles, study, path,  additional arguments...)
# 	W <- getW( Wlist, ids, rsids, additional arguments...)
#
#
# Required input:
#
# study: the study name
# path: the name of the path where the binary genotype files are stored
# bedfiles: a vector of BED file names
# bimfiles: a vector of BIM file names
#
# Optional input:
#
# ids: a vector of subject IDs used in W (row IDS) 
# rsids: a vector of genetic variant IDs used in W (column IDs) 
# msize: the number of genetic markers used in each chunk
#
# Output:
#
# The output of the prepW and computeW functions is a Wlist structure containing information 
# about the genotype matrix W used in downstream analyses. This should be saved in Rdata file and used 
# for downstream analyses. Most users do not use this directly. Furthermore this structure 
# also allow us to change the underlying data structures without apparent changes in 
# the Wlist interface.  
# 
# The Wlist structure contains the following slots (which are likely to change in future releases):
#
#   Wlist$fnW vector of file names 
#   Wlist$rsids vector of rsids
#   Wlist$alleles vector of allele used as alternative  
#   Wlist$chr vector chromosome
#   Wlist$ids vector of subject ids used in W
#   Wlist$nb number of chunks
#   Wlist$m total number of markers in W
#   Wlist$n total number if subjects
#   Wlist$msize number of markers in each chunk
#   Wlist$study study name
#   Wlist$bedfiles vector bedfile names
#   Wlist$bimfiles vector bimfile names
#
###############################################################################################


#######################################################################################
# compute W functions
#######################################################################################

#' @export
#'

computeW <- function(Wlist=NULL,chr=NULL, scaleW=TRUE) {
     ids <- Wlist$ids
     bedfile <- Wlist$bedfiles[chr]
     fnW <- Wlist$fnW[Wlist$chr==chr]
     rsids <- Wlist$rsids[Wlist$chr==chr]
     alleles <- Wlist$alleles[Wlist$chr==chr]
     bedfile <- Wlist$bedfiles[chr]
     
     for ( i in 1:length(fnW) ) {
          g <- readBed( bedfile=bedfile, rsids=rsids[[i]], alleles=alleles[[i]], ids=ids) 
          if(scaleW) g <- scale(g) 
          g[is.na(g)] <- 0
          bfW <- gzfile(fnW[i],"wb")
          for ( j in 1:ncol(g) ) {
               writeBin( g[,j], bfW, size = 8, endian = "little")
          }
          close(bfW)
          print(paste("Finished block",i,"proportion",i/length(fnW)))
     }
}

#' @export
#'

prepW <- function( study=NULL, path=NULL, bedfiles=NULL, bimfiles=NULL, ids=NULL, rsids=NULL, msize=NULL ){
     
     Wlist <- NULL
     Wlist$fnW <- NULL
     Wlist$rsids <- NULL
     Wlist$alleles <- NULL
     Wlist$chr <- NULL
     
     for ( chr in 1:22 ) {
          snps <- read.table(file=bimfiles[chr], header=FALSE)
          alleles <- as.character(snps[,5])   
          rsids <- as.character(snps[,2])   
          m <- nrow(snps)
          sets <- split(1:m, ceiling(seq_along(1:m)/msize))
          nsets <- length(sets)
          fnWChr <- paste0(path,"/W",chr,"_",1:nsets,"_",study,".gz",sep="")
          if (chr==1) Wlist$fnW <- fnWChr
          if (chr>1) Wlist$fnW <- c(Wlist$fnW,fnWChr)
          if (chr==1) Wlist$rsids <- split(rsids,ceiling(seq_along(1:m)/msize))
          if (chr>1) Wlist$rsids <- c(Wlist$rsids,split(rsids,ceiling(seq_along(1:m)/msize)))
          if (chr==1) Wlist$alleles <- split(alleles,ceiling(seq_along(1:m)/msize))
          if (chr>1) Wlist$alleles <- c(Wlist$alleles,split(alleles,ceiling(seq_along(1:m)/msize)))
          if (chr==1) Wlist$chr <- rep(chr,nsets)
          if (chr>1) Wlist$chr <- c(Wlist$chr,rep(chr,nsets))
     }   
     
     rsids <- unlist(Wlist$rsids)
     Wlist$indx <- lapply( Wlist$rsids, function(x) {match(x,rsids)} )
     
     Wlist$ids <- as.character(ids)
     
     Wlist$nb <- length(Wlist$fnW)
     Wlist$m <- sum(sapply(Wlist$rsids,length))
     Wlist$n <- length(Wlist$ids)
     Wlist$msize <- msize
     Wlist$study <- study
     Wlist$bedfiles <- bedfiles
     Wlist$bimfiles <- bimfiles
     
     return(Wlist)
}


#' @export
#'

readBed <- function( bedfile=NULL, rsids=NULL, alleles=NULL, ids=NULL){
     
     # make temp dir and write snplist to file
     
     temp_dir <- tempdir()
     
     system(paste("mkdir -p", temp_dir) )
     
     fnSNPs <- tempfile(tmpdir = temp_dir)
     fnIDs <- tempfile(tmpdir = temp_dir)
     fnAlleles <- tempfile(tmpdir = temp_dir)
     fnGenotypes <- tempfile(tmpdir = temp_dir)
     
     if (!is.null(alleles)) write.table( cbind(rsids,alleles), row.names = F, col.names = F, quote = F, sep = "\t", file = fnAlleles  )  
     if (!is.null(ids)) write.table( cbind(ids,ids), row.names = F, col.names = F, quote = F, sep = "\t", file = fnIDs  )  
     write.table( rsids, row.names = F, col.names = F, quote = F, sep = "\t", file = fnSNPs )  
     
     # call plink to write out additively coded SNPs
     if (is.null(ids)) {
          system( paste("plink2 --bfile", bedfile,
                        "--extract", fnSNPs,
                        "--keep-allele-order --recode A --recode-allele ",fnAlleles,
                        "--silent",
                        "--out", fnGenotypes ) ) }
     
     # call plink to write out additively coded SNPs
     if (!is.null(ids)) {
          system( paste("plink2 --bfile", bedfile,
                        "--extract", fnSNPs,
                        "--keep", fnIDs,
                        "--keep-allele-order --recode A --recode-allele ",fnAlleles,
                        "--silent",
                        "--out", fnGenotypes ) )}
     
     genotypes <- read.table( paste0( fnGenotypes, ".raw"), head = T, check.names=FALSE )
     ids <- genotypes$IID
     genotypes <- genotypes[ , c(7:ncol(genotypes))]
     rsids_alleles <- colnames(genotypes)
     rsidsG <- t(sapply(rsids_alleles, function(x){ strsplit(x,split="_", fixed=TRUE)[[1]][1]} ))[1,]
     allelG <- t(sapply(rsids_alleles, function(x){ strsplit(x,split="_", fixed=TRUE)[[1]][2]} ))[1,]
     rownames(genotypes) <- ids
     colnames(genotypes) <- rsidsG
     system( paste("rm -r", temp_dir ) )	# clean up
     
     
     if (any(!rsidsG==rsids)) stop("Some variants not found in bedfile")
     
     
     return(genotypes) }

#######################################################################################

