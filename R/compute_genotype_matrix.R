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


computeW <- function(Wlist=NULL, ids=NULL, rsids=NULL, overwrite=FALSE) {


     writeBED2RAW( fnRAW=Wlist$fnRAW, bedfiles=Wlist$bedfiles, bimfiles=Wlist$bimfiles, famfiles=Wlist$famfiles, ids=ids, rsids=rsids, overwrite=overwrite)

}


#' @export
#'

prepW <- function( study=NULL, fnRAW=NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, ids=NULL, rsids=NULL){
     
     nfiles <- length(bedfiles)
     
     if(nfiles>1) {
          nchr <- nfiles
          if(file.exists(fnRAW)) warning(paste("fnRAW allready exist"))
          
          Wlist <- NULL
          if(!is.null(fnRAW)) Wlist$fnRAW <- fnRAW
          Wlist$rsids <- vector(mode="list",length=nchr)
          Wlist$alleles <- vector(mode="list",length=nchr)
          Wlist$af <- vector(mode="list",length=nchr)
          Wlist$maf <- vector(mode="list",length=nchr)
          Wlist$nmiss <- vector(mode="list",length=nchr)
          Wlist$het <- vector(mode="list",length=nchr)
          Wlist$n0 <- vector(mode="list",length=nchr)
          Wlist$n1 <- vector(mode="list",length=nchr)
          Wlist$n2 <- vector(mode="list",length=nchr)
          Wlist$chr <- 1:nchr
          Wlist$nchr <- nchr
          
          Wlist$study_ids <- NULL
          fam <- read.table(file=famfiles[1], header=FALSE)
          Wlist$ids <- as.character(fam[,2])
          if(!is.null(ids)) Wlist$study_ids <- as.character(ids) 
          if(!is.null(ids)) if(any(!ids%in%as.character(fam[,2]))) stop(paste("some ids not found in famfiles"))
          
          for ( chr in 1:length(bedfiles) ) {
               bim <- read.table(file=bimfiles[chr], header=FALSE)
               fam <- read.table(file=famfiles[chr], header=FALSE)
               if(any(!Wlist$ids%in%as.character(fam[,2]))) stop(paste("some ids not found in famfiles"))
               Wlist$alleles[[chr]] <- as.character(bim[,6])   
               Wlist$rsids[[chr]] <- as.character(bim[,2])   
               print(paste("Finished processing bim file",bimfiles[chr]))
          }   
          Wlist$study_rsids <- NULL
          if(!is.null(rsids)) Wlist$study_rsids <- as.character(rsids) 
          if(!is.null(rsids)) if( any(!rsids%in%unlist(Wlist$rsids)) ) stop(paste("some rsids not found in bimfiles"))
          
          Wlist$mchr <- sapply(Wlist$rsids,length)
          Wlist$m <- sum(Wlist$mchr)
          Wlist$n <- length(Wlist$ids)
          Wlist$study <- study
          Wlist$bedfiles <- bedfiles
          Wlist$bimfiles <- bimfiles
          Wlist$famfiles <- famfiles

     }

     if(nfiles==1) {

          bim <- read.table(file=bimfiles[1], header=FALSE)
          fam <- read.table(file=famfiles[1], header=FALSE)

          Wlist <- NULL
          
          Wlist$chr <- bim[!duplicated(bim[,1]),1]
          nchr <- length(Wlist$chr)
          Wlist$nchr <- length(Wlist$chr)
          Wlist$rsids <- split(as.character(bim[,2]),f=as.factor(bim[,1]))
          Wlist$alleles <- split(as.character(bim[,6]),f=as.factor(bim[,1]))
          Wlist$position <- split(bim[,4],f=as.factor(bim[,1]))
          
          if(file.exists(fnRAW)) warning(paste("fnRAW allready exist"))
          Wlist$fnRAW <- fnRAW
          
          Wlist$af <- vector(mode="list",length=nchr)
          Wlist$maf <- vector(mode="list",length=nchr)
          Wlist$nmiss <- vector(mode="list",length=nchr)
          Wlist$nhet <- vector(mode="list",length=nchr)
          Wlist$n0 <- vector(mode="list",length=nchr)
          Wlist$n1 <- vector(mode="list",length=nchr)
          Wlist$n2 <- vector(mode="list",length=nchr)
          Wlist$chr <- 1:nchr
          
          Wlist$study_ids <- NULL
          fam <- read.table(file=famfiles[1], header=FALSE)
          Wlist$ids <- as.character(fam[,2])
          if(!is.null(ids)) Wlist$study_ids <- as.character(ids) 
          if(!is.null(ids)) if(any(!ids%in%as.character(fam[,2]))) stop(paste("some ids not found in famfiles"))
          
          if(any(duplicated(Wlist$ids))) stop("Duplicated ids found in famfiles")

          Wlist$study_rsids <- NULL
          if(!is.null(rsids)) Wlist$study_rsids <- as.character(rsids) 
          if(!is.null(rsids)) if( any(!rsids%in%unlist(Wlist$rsids)) ) stop(paste("some rsids not found in bimfiles"))
          
          Wlist$mchr <- sapply(Wlist$rsids,length)
          Wlist$m <- sum(Wlist$mchr)
          Wlist$n <- length(Wlist$ids)
          Wlist$study <- study
          Wlist$bedfiles <- bedfiles
          Wlist$bimfiles <- bimfiles
          Wlist$famfiles <- famfiles

     }
          return(Wlist)
}


#' @export
#'



writeBED2RAW <- function(fnRAW=NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, ids=NULL, rsids=NULL,overwrite=FALSE) {

     if(file.exists(fnRAW)) {
          warning(paste("fnRAW file allready exist"))
          if(!overwrite) stop(paste("fnRAW file allready exist"))
     }
     bfRAW <- file(fnRAW,"wb")
     for ( chr in 1:length(bedfiles)) {
          print(paste("Processing bedfile:",bedfiles[chr]))
          bim <- read.table(file=bimfiles[chr], header=FALSE)
          fam <- read.table(file=famfiles[chr], header=FALSE)
          n <- nrow(fam)
          m <- nrow(bim)
          #indx <- seq(1,n*2,2)
          nbytes <- ceiling(n/4)
          printmarker <- rep(F,m)
          printmarker[seq(1,m,10000)] <- T
          fnBED <- bedfiles[chr]  
          bfBED <- file(fnBED,"rb")
          magic <- readBin(bfBED, "raw", n=3)
          if (!all(magic[1] == "6c", magic[2] == "1b", magic[3] == "01"))
               stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
          for ( j in 1:m) {
               raw <- readBin(bfBED, "raw", n=nbytes)
               writeBin( raw, bfRAW, size = 1, endian = "little")
               #raw <- as.logical(rawToBits(readBin(bfBED, "raw", bsize)))
               #raw1 <- raw[indx]
               #raw2 <- raw[indx+1]
               #isNA <- raw1==1 & raw2==0
               #g <- raw1 + raw2 + 1
               #g[isNA] <- 0
               #writeBin( as.raw(g[rws]), bfRAW, size = 1, endian = "little")
               if(printmarker[j]) print(paste("Finished marker",j))
          }
          close(bfBED)
          print(paste("Finished processing bedfile:",bedfiles[chr]))
     }
     close(bfRAW)
     
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

readbed.plink <- function( bedfile=NULL, rsids=NULL, alleles=NULL, ids=NULL,plink="plink2"){
     
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
          system( paste(plink," --bfile", bedfile,
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
     if(!is.null(ids)) rws <- match(ids,as.character(fam[,2]))
     cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,as.character(bim[,2]))
     cls <- cls[order(cls)]
     
     #W <- matrix(0,nrow=length(rws),ncol=length(cls))
     W <- rep(0,length(rws)*length(cls))
     dim(W) <-c(length(rws),length(cls))
     
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
          raw1 <- raw[indx]
          raw2 <- raw[indx+1]
          isNA <- raw1==1 & raw2==0
          g <- raw1 + raw2
          #g <- raw[indx] + raw[indx+1]
          #g[raw[indx]==1 & raw[indx+1]==0 ] <- NA # 1/0 is missing
          g[isNA] <- NA # 1/0 is missing
          W[,i] <- g[rws]
     }
     close(bfBED)
     #rownames(W) <- as.character(fam[rws,2])
     #colnames(W) <- as.character(bim[cls,2])
     return(W)
     }


#' @export
#'

readbed <- function(Wlist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL,scaled=TRUE) { 
     #readbed(n,nr,rws,nc,cls,scaled,W,nbytes,fnRAW)	
     #dll <- paste(find.package("qgg"),"/libs/qgg.so",sep="")    
     #dyn.load(dll)
     #is.loaded("readbed")
     n <- Wlist$n
     m <- Wlist$m
     nbytes <- ceiling(n/4)
     if(is.null(cls)) cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,unlist(Wlist$rsids))
     nc <- length(cls)
     if (is.null(rws)) {
       rws <- 1:n
       if(!is.null(ids)) rws <- match(ids,Wlist$ids)
       if(!is.null(Wlist$study_ids)) rws <- match(Wlist$study_ids,Wlist$ids)
     }
     nr <- length(rws)
     fnRAW <- Wlist$fnRAW
     res <- .Fortran("readbed", 
                     n = as.integer(n),
                     nr = as.integer(nr),
                     rws = as.integer(rws),
                     nc = as.integer(nc),
                     cls = as.integer(cls),
                     scaled = as.integer(scaled),
                     W = matrix(as.double(0),nrow=nr,ncol=nc),
                     nbytes = as.integer(nbytes),
                     fnRAW = as.character(fnRAW),
                     PACKAGE = 'qgg'
                     
     )
     #dyn.unload(dll)
     #rownames(res$W) <- Wlist$ids[rws]
     #colnames(res$W) <- unlist(Wlist$ids)[cls]
     return(res$W)
}



#' @export
#'

qcraw <- function(Wlist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL) { 
     #subroutine craw(n,nr,rws,nc,cls,af,nmiss,n0,n1,n2,fnRAW)	
     n <- Wlist$n
     m <- Wlist$m
     if(is.null(cls)) cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,unlist(Wlist$rsids))
     nc <- length(cls)
     rws <- 1:n
     if(!is.null(ids)) rws <- match(ids,Wlist$ids)
     nr <- length(rws)
     af <- nmiss <- n0 <- n1 <- n2 <- rep(0,nc)
     fnRAW <- Wlist$fnRAW
     qc <- .Fortran("qcraw", 
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
                     fnRAW = as.character(fnRAW),
                     PACKAGE = 'qgg'
                     
     )
     qc$hom <- (qc$n0+qc$n2)/(qc$n-qc$nmiss)
     qc$het <- qc$n1/(qc$n-qc$nmiss)
     qc$maf <- qc$af
     qc$maf[qc$maf>0.5] <- 1-qc$maf[qc$maf>0.5]
     #names(fit$af) <- unlist(Wlist$rsids)[cls]
     #names(fit$nmiss) <- unlist(Wlist$rsids)[cls]
     #return(list(af=fit$af,nmiss=fit$nmiss))
     return(qc)
}

#' @export
#'

summaryW <- function(Wlist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL) {
     
     #qc <- qcraw(Wlist=Wlist,ids=ids, rsids=rsids, rws=rws, cls=cls)    
     qc <- qcbed(Wlist=Wlist,ids=ids, rsids=rsids, rws=rws, cls=cls)    
     Wlist$nmiss <- qc$nmiss
     Wlist$af <- qc$af
     Wlist$maf <- qc$maf
     Wlist$hom <- qc$hom
     Wlist$het <- qc$het
     Wlist$n0 <- qc$n0
     Wlist$n1 <- qc$n1
     Wlist$n2 <- qc$n2
     
     return(Wlist)
}


#' @export
#'

qcbed <- function(Wlist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL) { 
     #dll <- paste(find.package("qgg"),"/libs/qgg.so",sep="")    
     #dyn.load(dll)
     #is.loaded("qcbed")
     n <- Wlist$n
     m <- Wlist$m
     nbytes <- ceiling(n/4)
     if(is.null(cls)) cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,unlist(Wlist$rsids))
     nc <- length(cls)
     rws <- 1:n
     if(!is.null(ids)) rws <- match(ids,Wlist$ids)
     if(!is.null(Wlist$study_ids)) rws <- match(Wlist$study_ids,Wlist$ids)
     nr <- length(rws)
     af <- nmiss <- n0 <- n1 <- n2 <- rep(0,nc)
     fnRAW <- Wlist$fnRAW
     qc <- .Fortran("qcbed", 
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
                    PACKAGE = 'qgg'
                    
     )
     #dyn.unload(dll)
     qc$hom <- (qc$n0+qc$n2)/(qc$n-qc$nmiss)
     qc$het <- qc$n1/(qc$n-qc$nmiss)
     qc$maf <- qc$af
     qc$maf[qc$maf>0.5] <- 1-qc$maf[qc$maf>0.5]
     rsids <- unlist(Wlist$rsids)[cls]
     names(qc$af) <- rsids
     names(qc$maf) <- rsids
     names(qc$hom) <- rsids
     names(qc$het) <- rsids
     names(qc$n0) <- rsids
     names(qc$n1) <- rsids
     names(qc$n2) <- rsids
     names(qc$nmiss) <- rsids
     return(qc)
}

#' @export
#'

mafbed <- function(Wlist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL) { 
     #dll <- paste(find.package("qgg"),"/libs/qgg.so",sep="")    
     #dyn.load(dll)
     #is.loaded("mafbed")
     n <- Wlist$n
     m <- Wlist$m
     nbytes <- ceiling(n/4)
     if(is.null(cls)) cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,unlist(Wlist$rsids))
     nc <- length(cls)
     rws <- 1:n
     if(!is.null(ids)) rws <- match(ids,Wlist$ids)
     if(!is.null(Wlist$study_ids)) rws <- match(Wlist$study_ids,Wlist$ids)
     nr <- length(rws)
     af <-rep(0,nc)
     fnRAW <- Wlist$fnRAW
     qc <- .Fortran("mafbed", 
                    n = as.integer(n),
                    nr = as.integer(nr),
                    rws = as.integer(rws),
                    nc = as.integer(nc),
                    cls = as.integer(cls),
                    af = as.double(af),  
                    nbytes = as.integer(nbytes),  
                    fnRAW = as.character(fnRAW),                     
                    PACKAGE = 'qgg'
                    
     )
     #dyn.unload(dll)
     qc$maf <- qc$af
     qc$maf[qc$maf>0.5] <- 1-qc$maf[qc$maf>0.5]
     rsids <- unlist(Wlist$rsids)[cls]
     names(qc$maf) <- rsids
     return(qc$maf)
}


#######################################################################################

