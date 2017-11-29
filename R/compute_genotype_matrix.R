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

computeW <- function(Wlist=NULL,chr=NULL, scaled=TRUE, compressed=TRUE) {
     ids <- Wlist$ids
     bedfile <- Wlist$bedfiles[chr]
     fnW <- Wlist$fnW[Wlist$chr==chr]
     rsids <- Wlist$rsids[Wlist$chr==chr]
     alleles <- Wlist$alleles[Wlist$chr==chr]
     bedfile <- Wlist$bedfiles[chr]
     
     for ( i in 1:length(fnW) ) {
          g <- readBed( bedfile=bedfile, rsids=rsids[[i]], alleles=alleles[[i]], ids=ids) 
          if(scaled) g <- scale(g) 
          g[is.na(g)] <- 0
          if(compressed) bfW <- gzfile(fnW[i],"wb")
          if(!compressed) bfW <- file(fnW[i],"wb")
          for ( j in 1:ncol(g) ) {
               writeBin( g[,j], bfW, size = 8, endian = "little")
          }
          close(bfW)
          print(paste("Finished block",i,"proportion",i/length(fnW)))
     }
}

#' @export
#'

prepW <- function( study=NULL, path=NULL, fnRAW=NULL, bedfiles=NULL, bimfiles=NULL, ids=NULL, rsids=NULL, compressed=FALSE ){

     nfiles <- length(bedfiles)
     
     Wlist <- NULL
     Wlist$fnRAW <- vector(length=nfiles)
     if(!is.null(fnRAW)) Wlist$fnRAW <- fnRAW
     Wlist$rsids <- vector(mode="list",length=nfiles)
     Wlist$alleles <- vector(mode="list",length=nfiles)
     Wlist$chr <- vector(length=nfiles)

     for ( chr in 1:length(bedfiles) ) {
          snps <- read.table(file=bimfiles[chr], header=FALSE)
          Wlist$alleles[[chr]] <- as.character(snps[,5])   
          Wlist$rsids[[chr]] <- as.character(snps[,2])   
          if(is.null(fnRAW)) Wlist$fnRAW[[chr]] <- paste0(path,"/RAW",chr,"_",study,sep="")
          print(paste("Finished processing bim file",bimfiles[chr]))
     }   
     
     rsids <- unlist(Wlist$rsids)
     Wlist$indx <- lapply( Wlist$rsids, function(x) {match(x,rsids)} )
     
     Wlist$ids <- as.character(ids)
     
     Wlist$nb <- length(Wlist$fnRAW)
     Wlist$m <- sum(sapply(Wlist$rsids,length))
     Wlist$n <- length(Wlist$ids)
     Wlist$study <- study
     Wlist$bedfiles <- bedfiles
     Wlist$bimfiles <- bimfiles
     Wlist$compressed <- compressed
     
     return(Wlist)
}

#' @export
#'

writeBED2RAW <- function(Wlist=NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, fnRAW=NULL, ids=NULL, rsids=NULL, alleles=NULL, chr=NULL, msize=100) {
     n <- length(ids)
     m <- length(rsids)
     sets <- split(1:m, ceiling(seq_along(1:m)/msize))
     nsets <- length(sets)
     rsids <- split(rsids,ceiling(seq_along(1:m)/msize))
     alleles <- split(alleles,ceiling(seq_along(1:m)/msize))
     bfW <- file(fnRAW,"wb")
     for ( i in 1:nsets ) {
          g <- readBed( bedfile=bedfiles, rsids=rsids[[i]], alleles=alleles[[i]], ids=ids)
          g <- g+1 
          g[is.na(g)] <- 0
          for ( j in 1:ncol(g) ) {
               writeBin( as.raw(g[,j]), bfW, size = 1, endian = "little")
          }
          gc()
          print(paste("Finished block",i,"proportion",i/nsets))
     }
     close(bfW)
}

#' @export
#'

compressRAW <- function(Wlist=NULL, chr=NULL, type="gzip") {
     rsids <- Wlist$rsids[[chr]]
     n <- Wlist$n
     m <- length(rsids)
     fnRAW <- Wlist$fnRAW[chr]
     wgz <- vector(mode = "list", length = m)
     bfW <- file(fnRAW,"rb")
     if(type=="gzip") {
          for (i in 1:m ) {
               w <- readBin( bfW, "raw", n=n, size = 1, endian = "little")
               wgz[[i]] <- memCompress( w, type=type)
               print(paste("Compressed marker",i,"chromosome",chr))
          }
     }     
     if(type=="lz4") {
          library(lz4)
          for (i in 1:m ) {
               w <- readBin( bfW, "raw", n=n, size = 1, endian = "little")
               wgz[[i]] <- lzCompress( w, level=8)
               print(paste("Compressed marker",i,"chromosome",chr))
          }
     }     
     close(bfW)
     names(wgz) <- rsids
     return(wgz)
}

#' @export
#'

getW <- function(Wlist=NULL, ids=NULL, rsids=NULL, rws=NULL,cls=NULL, scaled=FALSE) {
     if(is.null(ids)) ids <- Wlist$ids
     if(is.null(cls)) cls <- match(rsids,unlist(Wlist$rsids))
     if(is.null(rws)) rws <- match(ids,Wlist$ids)
     maf <- unlist(Wlist$maf)[cls]
     meanW <- 2*maf
     sdW <- sqrt(2*maf*(1-maf))
     n <- Wlist$n
     #W <- matrix(logical(0),nrow=length(rws),ncol=length(cls))
     W <- rep(0,length(rws)*length(cls))
     dim(W) <- c(length(rws),length(cls))
     bfW <- file(Wlist$fnRAW,"rb")
     current <- 0
     for (i in 1:length(cls) ) {
          where <- (cls[i]-current-1)*n
          current <- cls[i]
          seek(bfW, where=where, origin="current", rw="read")
          w <- as.double(readBin( bfW, "raw", n=n, size = 1, endian = "little"))
          if(scaled) w[w>0] <- (w[w>0]-1-meanW[i])/sdW[i]
          W[,i] <- w[rws]
     }
     close(bfW)
     return(W)
}


#' @export
#'

qcW <- function(Wlist=NULL) {
     qcG <- NULL
     for ( chr in 1:Wlist$nb ) {
          qcG[[chr]] <- qcGenotypes(Wlist=Wlist,chr=chr)  
          print(chr) 
     }
     Wlist$nmiss <- sapply(qcG,function(x){ x$nmiss})
     Wlist$af <- sapply(qcG,function(x){ x$af})
     Wlist$maf <- sapply(qcG,function(x){ x$maf})
     return(Wlist)
}

#' @export
#'

qcGenotypes <- function(Wlist=NULL, chr=NULL) {
     rsids <- Wlist$rsids[[chr]]
     n <- Wlist$n
     m <- length(rsids)       
     nmiss <- af <- rep(0,m)      
     bfW <- file(Wlist$fnRAWCHR[chr],"rb")
     for ( i in 1:m ) { 
          g <- as.integer(readBin( bfW, "raw", n=n, size = 1, endian = "little"))
          miss <- g==0
          nmiss[i] <- sum(miss)
          af[i] <- sum(g[!miss]-1)
     }
     close(bfW)
     
     maf <- af <- af/((n-nmiss)*2)    
     maf[maf>0.5] <- 1-maf[maf>0.5]
     names(af) <- names(maf) <- names(nmiss) <- rsids
     return(list(maf=maf,af=af,nmiss=nmiss)) 
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
     #rsidsG <- t(sapply(rsids_alleles, function(x){ strsplit(x,split="_", fixed=TRUE)[[1]][1]} ))[1,]
     #allelG <- t(sapply(rsids_alleles, function(x){ strsplit(x,split="_", fixed=TRUE)[[1]][2]} ))[1,]
     rsidsG <- substr(rsids_alleles,1,nchar(rsids))
     allelG <- substr(rsids_alleles,nchar(rsids)+2,nchar(rsids_alleles))
     rownames(genotypes) <- ids
     colnames(genotypes) <- rsidsG
     system( paste("rm -r", temp_dir ) )	# clean up
     
     
     if (any(!rsidsG==rsids)) stop("Some variants not found in bedfile")
     
     
     return(genotypes) }

#' @export
#'

readBED <- function( bedfiles=NULL, bimfiles=NULL, famfiles=NULL, chr=NULL, rsids=NULL, alleles=NULL, ids=NULL){
     #adapted from https://github.com/andrewparkermorgan/argyle/blob/master/R/plink.R
     bim <- read.table(file=bimfiles[chr], header=FALSE)
     fam <- read.table(file=famfiles[chr], header=FALSE)
     
     n <- nrow(fam)
     m <- nrow(bim)
     bsize <- ceiling(n/4)
     indx <- seq(1,n*2,2)
     rws <- 1:n
     if(!is.null(ids)) rws <- match(ids,fam[,2])
     cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,bim[,2])
     cls <- cls[order(cls)]
     
     W <- matrix(0,nrow=length(rws),ncol=length(cls))
     
     bfBED <- file(bedfiles[chr],"rb")
     magic <- readBin(bfBED, "raw", n=3)
     if (!all(magic[1] == "6c", magic[2] == "1b", magic[3] == "01"))
          stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
     
     current <- 0
     for ( i in 1:length(cls)) {
          where <- (cls[i]-current-1)*bsize
          current <- cls[i]
          if(!current==(i-1)) seek(bfBED, where=where, origin="current", rw="read")
          raw <- as.logical(rawToBits(readBin(bfBED, "raw", bsize)))
          g <- raw[indx] + raw[indx+1]
          g[raw[indx]==1 & raw[indx+1]==0 ] <- NA # 1/0 is missing
          W[,i] <- g[rws]
     }
     close(bfBED)
     rownames(W) <- fam[rws,2]
     colnames(W) <- bim[cls,2]
     return(W)
     }


#' @export
#'

readraw <- function(Wlist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL,scaled=TRUE) { 
     #subroutine readraw(n,m,cls,nc,scaled,W,fnRAW)	
     dll <- paste(find.package("qgg"),"/libs/qgg.so",sep="")    
     #dll <- "/data/home/peters/prs/src/readraw.so"    
     dyn.load(dll)
     n <- Wlist$n
     if(!is.null(rsids)) cls <- match(rsids,unlist(Wlist$rsids))
     #o <- order(cls)
     #cls <- cls[o]
     nc <- length(cls)
     fnRAW <- Wlist$fnRAW
     fit <- .Fortran("readraw", 
                     n = as.integer(n),
                     cls = as.integer(cls),
                     nc = as.integer(nc),
                     scaled = as.integer(scaled),
                     W = matrix(as.double(0),nrow=n,ncol=nc),
                     fnRAW = as.character(fnRAW)
                     
     )
     #dimnames(fit$W) <- list(Wlist$ids,Wlist$rsids[cls])
     #rownames(fit$W) <- Wlist$ids
     #colnames(fit$W) <- Wlist$rsids[cls]
     return(fit$W)
}




#######################################################################################

