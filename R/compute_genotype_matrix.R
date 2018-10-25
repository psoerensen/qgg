###############################################################################################
#  Prepare (processed genotypes) for target population
###############################################################################################

#' Prepare genotype data for all statistical analyses (initial step) 
#'
#' @description
#' The processed genotypes are stored in a matrix W. 
#' By default genotypes are allele counts of alternative allele which are centered and scaled. 
#' 
#' Currently the raw genotype format supported is plink bed/bim/fam files. 
#' To fascilitate a simple and common interface for using raw genotypes we prepare a Glist structure 
#' which contains information about the raw and processed genotypes. 
#' Glist is obtained using the prepG function based on plink bed/bim/fam files.
#' The genotype matrix, W, is too big to fit into memory and is therefore stored in a binary file.
#' By default the W matrix is stored in a binary file (n*m*4 bytes).
#' Information about the processed genotype matrix W is provided in a Glist structure. 
#' 
#' The output of the prepG functions is a Glist structure containing information 
#' about the genotype matrix W used in downstream analyses. This should be saved in Rdata file and used 
#' for downstream analyses. Most users do not use this directly. Furthermore this structure 
#' also allow us to change the underlying data structures without apparent changes in 
#' the Glist interface.  
#' 
#' @param study name of the study
#' @param fnRAW name of binary file used for storing genotypes on disk
#' @param bedfiles a vector names for the bed files
#' @param famfiles a vector names for the fam files
#' @param bimfiles a vector names for the bim files
#' @param ids individual IDs used in study
#' @param rsids marker rsids used in study
#' @param overwrite logical if TRUE overwite binary genotype file
#' @param ncores used to process genotypes
#' 
#' @return Returns a list structure with information about genotypes
#' 


#' @author Peter Soerensen

#' @examples
#

#' #Glist <- prepG( bedfiles, bimfiles, study, path, additional arguments...)
#' #W <- getW( Glist, ids, rsids, additional arguments...)


#' @export
#'

prepG <- function( study=NULL, fnRAW=NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, ids=NULL, rsids=NULL, overwrite=TRUE, ncores=1){
     
     nfiles <- length(bedfiles)
      
     if(nfiles>1) {
          nchr <- nfiles
          if(file.exists(fnRAW)) warning(paste("fnRAW allready exist"))
          
          Glist <- NULL
          if(!is.null(fnRAW)) Glist$fnRAW <- fnRAW
          Glist$rsids <- vector(mode="list",length=nchr)
          Glist$alleles <- vector(mode="list",length=nchr)
          Glist$position <- vector(mode="list",length=nchr)
          Glist$af <- vector(mode="list",length=nchr)
          Glist$maf <- vector(mode="list",length=nchr)
          Glist$nmiss <- vector(mode="list",length=nchr)
          Glist$het <- vector(mode="list",length=nchr)
          Glist$n0 <- vector(mode="list",length=nchr)
          Glist$n1 <- vector(mode="list",length=nchr)
          Glist$n2 <- vector(mode="list",length=nchr)
          Glist$chr <- 1:nchr
          Glist$chrnames <- 1:nchr
          Glist$nchr <- nchr
          
          Glist$study_ids <- NULL
          fam <- read.table(file=famfiles[1], header=FALSE)
          #fam <- fread(input=famfiles[1], header=FALSE)
          Glist$ids <- as.character(fam[,2])
          Glist$study_ids <- Glist$ids
          if(!is.null(ids)) {
               if(any(!ids%in%as.character(fam[,2]))) warning(paste("some ids not found in famfiles"))
               Glist$study_ids <- Glist$ids[Glist$ids%in%as.character(ids)]
          }
          if(any(duplicated(Glist$ids))) stop("Duplicated ids found in famfiles")
          for ( chr in 1:length(bedfiles) ) {
               #bim <- read.table(file=bimfiles[chr], header=FALSE)
               bim <- fread(input=bimfiles[chr], header=FALSE, data.table = FALSE)

               rsidsBIM <- as.character(bim[,2])
               if(!is.null(rsids)) bim <- droplevels(bim[rsidsBIM%in%rsids,])
               
               #fam <- read.table(file=famfiles[chr], header=FALSE)
               fam <- fread(input=famfiles[chr], header=FALSE, data.table = FALSE)

               if(any(!Glist$ids%in%as.character(fam[,2]))) stop(paste("some ids not found in famfiles"))
               Glist$alleles[[chr]] <- as.character(bim[,6])   
               Glist$position[[chr]] <- as.numeric(bim[,4])   
               Glist$rsids[[chr]] <- as.character(bim[,2])   
               print(paste("Finished processing bim file",bimfiles[chr]))
          }   
          Glist$study_rsids <- unlist(Glist$rsids)
          if(!is.null(rsids)) Glist$study_rsids <- as.character(rsids) 
          if(!is.null(rsids)) if( any(!rsids%in%unlist(Glist$rsids)) ) warning(paste("some rsids not found in bimfiles"))
          if (!is.null(rsids)) Glist$study_rsids <- unlist(Glist$rsids)[unlist(Glist$rsids)%in%rsids]
          

          Glist$mchr <- sapply(Glist$rsids,length)
          Glist$chr <- rep(1:nchr,times=Glist$mchr)
          Glist$rsids <- unlist(Glist$rsids)
          Glist$alleles <- unlist(Glist$alleles)
          Glist$position <- unlist(Glist$position)
          
          Glist$m <- sum(Glist$mchr)
          Glist$n <- length(Glist$ids)
          Glist$study <- study
          Glist$bedfiles <- bedfiles
          Glist$bimfiles <- bimfiles
          Glist$famfiles <- famfiles

          
          print("Preparing raw file")
          if(overwrite) computeRAW( Glist=Glist, ids=Glist$ids, rsids=Glist$rsids, overwrite=overwrite)  # write genotypes to .raw file 

          print("Computing allele frequencies, missingness")
          if(overwrite) Glist <- summaryRAW(Glist=Glist, ids=Glist$study_ids, rsids=Glist$rsids, ncores=ncores) # compute allele frequencies, missingness, ....    
    

     }

     if(nfiles==1) {

          #bim <- read.table(file=bimfiles[1], header=FALSE)
          #fam <- read.table(file=famfiles[1], header=FALSE)
          bim <- fread(input=bimfiles[1], header=FALSE, data.table = FALSE)
          fam <- fread(input=bimfiles[1], header=FALSE, data.table = FALSE)
          
          Glist <- NULL
          
          Glist$chr <- bim[!duplicated(bim[,1]),1]
          nchr <- length(Glist$chr)
          Glist$nchr <- length(Glist$chr)
          Glist$rsids <- split(as.character(bim[,2]),f=as.factor(bim[,1]))
          Glist$alleles <- split(as.character(bim[,6]),f=as.factor(bim[,1]))
          Glist$position <- split(bim[,4],f=as.factor(bim[,1]))
          
          if(file.exists(fnRAW)) warning(paste("fnRAW allready exist"))
          Glist$fnRAW <- fnRAW
          
          Glist$af <- vector(mode="list",length=nchr)
          Glist$maf <- vector(mode="list",length=nchr)
          Glist$nmiss <- vector(mode="list",length=nchr)
          Glist$nhet <- vector(mode="list",length=nchr)
          Glist$n0 <- vector(mode="list",length=nchr)
          Glist$n1 <- vector(mode="list",length=nchr)
          Glist$n2 <- vector(mode="list",length=nchr)
          Glist$chr <- 1:nchr
          Glist$chrnames <- 1:nchr
          
          Glist$study_ids <- NULL
          #fam <- read.table(file=famfiles[1], header=FALSE)
          fam <- fread(input=bimfiles[1], header=FALSE, data.table = FALSE)
          #Glist$ids <- as.character(fam[,2])
          #Glist$study_ids <- Glist$ids
          #if(!is.null(ids)) Glist$study_ids <- as.character(ids) 
          #if(!is.null(ids)) if(any(!ids%in%as.character(fam[,2]))) stop(paste("some ids not found in famfiles"))
          Glist$ids <- as.character(fam[,2])
          Glist$study_ids <- Glist$ids
          if(!is.null(ids)) {
               if(any(!ids%in%as.character(fam[,2]))) warning(paste("some ids not found in famfiles"))
               Glist$study_ids <- Glist$ids[Glist$ids%in%as.character(ids)]
          }
          
          if(any(duplicated(Glist$ids))) stop("Duplicated ids found in famfiles")

          Glist$study_rsids <- unlist(Glist$rsids)
          if(!is.null(rsids)) Glist$study_rsids <- as.character(rsids) 
          if(!is.null(rsids)) if( any(!rsids%in%unlist(Glist$rsids)) ) warning(paste("some rsids not found in bimfiles"))
          if (!is.null(rsids)) Glist$study_rsids <- unlist(Glist$rsids)[unlist(Glist$rsids)%in%rsids]
          
          Glist$mchr <- sapply(Glist$rsids,length)
          Glist$chr <- rep(1:nchr,times=Glist$mchr)
          Glist$rsids <- unlist(Glist$rsids)
          Glist$alleles <- unlist(Glist$alleles)
          Glist$position <- unlist(Glist$position)
          
          Glist$m <- sum(Glist$mchr)
          Glist$n <- length(Glist$ids)
          Glist$study <- study
          Glist$bedfiles <- bedfiles
          Glist$bimfiles <- bimfiles
          Glist$famfiles <- famfiles

          print("Preparing raw file")
          computeRAW( Glist=Glist, ids=ids, rsids=Glist$rsids, overwrite=overwrite)  # write genotypes to .raw file 
          
          print("Computing allele frequencies, missingness")
          Glist <- summaryRAW(Glist=Glist, ids=ids, rsids=Glist$rsids, ncores=ncores) # compute allele frequencies, missingness, ....    
          
     }
     return(Glist)
}

#' @export
#'

computeRAW <- function(Glist=NULL, ids=NULL, rsids=NULL, overwrite=FALSE) {
     bed2raw( fnRAW=Glist$fnRAW, bedfiles=Glist$bedfiles, bimfiles=Glist$bimfiles, famfiles=Glist$famfiles, ids=ids, rsids=rsids, overwrite=overwrite)
}



#' @export
#'

bed2raw <- function(fnRAW=NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, ids=NULL, rsids=NULL,overwrite=FALSE) {

     if(file.exists(fnRAW)) {
          warning(paste("fnRAW file allready exist"))
          if(!overwrite) stop(paste("fnRAW file allready exist"))
     }
     #bfRAW <- file(fnRAW,"wb")
     for ( chr in 1:length(bedfiles)) {
          print(paste("Processing bedfile:",bedfiles[chr]))
          #bim <- read.table(file=bimfiles[chr], header=FALSE)
          #fam <- read.table(file=famfiles[chr], header=FALSE)
          bim <- fread(input=bimfiles[chr], header=FALSE, data.table = FALSE)
          fam <- fread(input=famfiles[chr], header=FALSE, data.table = FALSE)

          n <- nrow(fam)
          m <- nrow(bim)
          rsidsBIM <- as.character(bim[,2])
          keep <- rep(TRUE,m)
          if(!is.null(rsids)) keep <- rsidsBIM%in%rsids
          #indx <- seq(1,n*2,2)
          nbytes <- ceiling(n/4)
          printmarker <- rep(F,m)
          printmarker[seq(1,m,10000)] <- T
          fnBED <- bedfiles[chr]  
          bfBED <- file(fnBED,"rb")
          magic <- readBin(bfBED, "raw", n=3)
          if (!all(magic[1] == "6c", magic[2] == "1b", magic[3] == "01"))
               stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
          #for ( j in 1:m) {
          #     raw <- readBin(bfBED, "raw", n=nbytes)
          #     if(keep[j]) writeBin( raw, bfRAW, size = 1, endian = "little")
          #     #raw <- as.logical(rawToBits(readBin(bfBED, "raw", bsize)))
          #     #raw1 <- raw[indx]
          #     #raw2 <- raw[indx+1]
          #     #isNA <- raw1==1 & raw2==0
          #     #g <- raw1 + raw2 + 1
          #     #g[isNA] <- 0
          #     #writeBin( as.raw(g[rws]), bfRAW, size = 1, endian = "little")
          #     if(printmarker[j]) print(paste("Finished marker",j))
          #}
          close(bfBED)
          append <- 1 
          if(chr==1) append <- 0
          res <- .Fortran("bed2raw", 
                          n = as.integer(n),
                          m = as.integer(m),
                          cls = as.integer(keep),
                          nbytes = as.integer(nbytes),
                          append = as.integer(append),
                          fnBED = as.character(fnBED),
                          fnRAW = as.character(fnRAW),
                          PACKAGE = 'qgg'
                          
          )
          print(paste("Finished processing bedfile:",bedfiles[chr]))
     }
     #close(bfRAW)
     
}

#' @export
#'
 
readbed <- function(Glist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL,scaled=TRUE, method="direct", ncores=1) { 
     n <- Glist$n
     m <- Glist$m
     nbytes <- ceiling(n/4)
     if(is.null(cls)) cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,Glist$rsids)
     nc <- length(cls)
     if (is.null(rws)) 1:n
     if(!is.null(ids)) rws <- match(ids,Glist$ids)
     #if (is.null(rws)) {
     #     rws <- 1:n
     #     if(!is.null(ids)) rws <- match(ids,Glist$ids)
     #     if(!is.null(Glist$study_ids)) rws <- match(Glist$study_ids,Glist$ids)
     #}
     nr <- length(rws)
     fnRAW <- Glist$fnRAW
     OS <- .Platform$OS.type
     #if(OS=="windows") fnRAW <- tolower(gsub("/","\\",fnRAW,fixed=T)) 
     if (method=="direct") {
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
     }
     if (method=="stream") {
          res <- .Fortran("readbedstream", 
                          n = as.integer(n),
                          nr = as.integer(nr),
                          rws = as.integer(rws),
                          nc = as.integer(nc),
                          cls = as.integer(cls),
                          scaled = as.integer(scaled),
                          W = matrix(as.double(0),nrow=nr,ncol=nc),
                          nbytes = as.integer(nbytes),
                          fnRAW = as.character(fnRAW),
                          ncores = as.integer(ncores),
                          PACKAGE = 'qgg'
                          
          )
     }
     
     #rownames(res$W) <- Glist$ids[rws]
     #colnames(res$W) <- Glist$ids[cls]
     return(res$W)
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

readbed.R <- function( bedfiles=NULL, bimfiles=NULL, famfiles=NULL, chr=NULL, rsids=NULL, alleles=NULL, ids=NULL){
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


getW <- function(Glist=NULL, ids=NULL, rsids=NULL, rws=NULL,cls=NULL, scaled=FALSE) {
     if(is.null(ids)) ids <- Glist$ids
     if(is.null(cls)) cls <- match(rsids,Glist$rsids)
     W <- readbed(Glist=Glist,ids=ids,rsids=rsids,rws=rws,cls=cls,scaled=scaled, method="direct") 
     rownames(W) <- ids
     colnames(W) <- Glist$rsids[cls]
     
     #if(is.null(rws)) rws <- match(ids,Glist$ids)
     #maf <- unlist(Glist$maf)[cls]
     #meanW <- 2*maf
     #sdW <- sqrt(2*maf*(1-maf))
     #n <- Glist$n
     #W <- matrix(logical(0),nrow=length(rws),ncol=length(cls))
     #W <- rep(0,length(rws)*length(cls))
     #dim(W) <- c(length(rws),length(cls))
     #bfW <- file(Glist$fnRAW,"rb")
     #current <- 0
     #for (i in 1:length(cls) ) {
     #     where <- (cls[i]-current-1)*n
     #     current <- cls[i]
     #     seek(bfW, where=where, origin="current", rw="read")
     #     w <- as.double(readBin( bfW, "raw", n=n, size = 1, endian = "little"))
     #     if(scaled) w[w>0] <- (w[w>0]-1-meanW[i])/sdW[i]
     #     W[,i] <- w[rws]
     #}
     #close(bfW)
     return(W)
}



#' @export
#'

summaryRAW <- function(Glist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL, ncores=1) {
     
     #qc <- qcraw(Glist=Glist,ids=ids, rsids=rsids, rws=rws, cls=cls)    
     qc <- qcbed(Glist=Glist,ids=ids, rsids=rsids, rws=rws, cls=cls, ncores=ncores)    
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

qcbed <- function(Glist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL,ncores=1) { 
     n <- Glist$n
     m <- Glist$m
     nbytes <- ceiling(n/4)
     if(is.null(cls)) cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,Glist$rsids)
     nc <- length(cls)
     rws <- 1:n
     if(!is.null(ids)) rws <- match(ids,Glist$ids)
     if(!is.null(Glist$study_ids)) rws <- match(Glist$study_ids,Glist$ids)
     nr <- length(rws)
     af <- nmiss <- n0 <- n1 <- n2 <- rep(0,nc)
     fnRAW <- Glist$fnRAW
     OS <- .Platform$OS.type
     #if(OS=="windows") fnRAW <- tolower(gsub("/","\\",fnRAW,fixed=T))    
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
                    ncores = as.integer(ncores),  
                    PACKAGE = 'qgg'
                    
     )
     qc$hom <- (qc$n0+qc$n2)/(qc$nr-qc$nmiss)
     qc$het <- qc$n1/(qc$nr-qc$nmiss)
     qc$maf <- qc$af
     qc$maf[qc$maf>0.5] <- 1-qc$maf[qc$maf>0.5]
     rsids <- Glist$rsids[cls]
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

mafbed <- function(Glist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL,ncores=1) { 
     n <- Glist$n
     m <- Glist$m
     nbytes <- ceiling(n/4)
     if(is.null(cls)) cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,Glist$rsids)
     nc <- length(cls)
     rws <- 1:n
     if(!is.null(ids)) rws <- match(ids,Glist$ids)
     if(!is.null(Glist$study_ids)) rws <- match(Glist$study_ids,Glist$ids)
     nr <- length(rws)
     af <- nmiss <- n0 <- n1 <- n2 <- rep(0,nc)
     fnRAW <- Glist$fnRAW
     OS <- .Platform$OS.type
     #if(OS=="windows") fnRAW <- tolower(gsub("/","\\",fnRAW,fixed=T))    
     qc <- .Fortran("mafbed", 
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
                    ncores = as.integer(ncores),  
                    PACKAGE = 'qgg'
                    
     )
     qc$hom <- (qc$n0+qc$n2)/(qc$nr-qc$nmiss)
     qc$het <- qc$n1/(qc$nr-qc$nmiss)
     qc$maf <- qc$af
     qc$maf[qc$maf>0.5] <- 1-qc$maf[qc$maf>0.5]
     rsids <- Glist$rsids[cls]
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




#######################################################################################
