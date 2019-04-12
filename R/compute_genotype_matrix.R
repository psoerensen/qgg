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

gprep <- function( study=NULL, fnRAW=NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, ids=NULL, rsids=NULL, overwrite=FALSE, ncores=1){
     
     nfiles <- length(bedfiles)

     Glist <- NULL
     Glist$study <- study
     
     Glist$fnRAW <- fnRAW
     if(file.exists(fnRAW)) warning(paste("fnRAW allready exist"))
     
     Glist$bedfiles <- bedfiles
     Glist$bimfiles <- bimfiles
     Glist$famfiles <- famfiles

     # Read fam information
     fam <- fread(input=famfiles[1], header=FALSE, data.table = FALSE, colClasses="character")
     Glist$ids <- as.character(fam[,2])
     Glist$study_ids <- Glist$ids
     if(!is.null(ids)) {
          if(any(!ids%in%as.character(fam[,2]))) warning(paste("some ids not found in famfiles"))
          Glist$study_ids <- Glist$ids[Glist$ids%in%as.character(ids)]
     }
     if(any(duplicated(Glist$ids))) stop("Duplicated ids found in famfiles")

     Glist$rsids <- vector(mode="list",length=nfiles)
     Glist$a1 <- vector(mode="list",length=nfiles)
     Glist$a2 <- vector(mode="list",length=nfiles)
     Glist$position <- vector(mode="list",length=nfiles)
     Glist$chr <- vector(mode="list",length=nfiles)
     
     for ( chr in 1:length(bedfiles) ) {
          bim <- fread(input=bimfiles[chr], header=FALSE, data.table = FALSE, colClasses="character")
          rsidsBIM <- as.character(bim[,2])
          if(!is.null(rsids)) bim <- droplevels(bim[rsidsBIM%in%rsids,])
          fam <- fread(input=famfiles[chr], header=FALSE, data.table = FALSE, colClasses="character")
          if(any(!Glist$ids%in%as.character(fam[,2]))) stop(paste("some ids not found in famfiles"))
          Glist$a1[[chr]] <- as.character(bim[,5])   
          Glist$a2[[chr]] <- as.character(bim[,6])   
          Glist$position[[chr]] <- as.numeric(bim[,4])   
          Glist$rsids[[chr]] <- as.character(bim[,2])   
          Glist$chr[[chr]] <- as.character(bim[,1])   
          print(paste("Finished processing bim file",bimfiles[chr]))
     }   
     Glist$study_rsids <- unlist(Glist$rsids)
     if(!is.null(rsids)) Glist$study_rsids <- as.character(rsids) 
     if(!is.null(rsids)) if( any(!rsids%in%unlist(Glist$rsids)) ) warning(paste("some rsids not found in bimfiles"))
     if (!is.null(rsids)) Glist$study_rsids <- unlist(Glist$rsids)[unlist(Glist$rsids)%in%rsids]
     
     Glist$mchr <- sapply(Glist$rsids,length)
     Glist$rsids <- unlist(Glist$rsids)
     Glist$a1 <- unlist(Glist$a1)
     Glist$a2 <- unlist(Glist$a2)
     Glist$position <- unlist(Glist$position)
     Glist$chr <- unlist(Glist$chr)
     Glist$nchr <- length(unique(Glist$chr))
     
     Glist$n <- length(Glist$ids)
     Glist$m <- length(Glist$rsids)
     
     print("Preparing raw file")
     writeRAW( Glist=Glist, ids=ids, rsids=Glist$rsids, overwrite=overwrite)  # write genotypes to .raw file 
     
     print("Computing allele frequencies, missingness")
     Glist <- summaryRAW(Glist=Glist, ids=ids, rsids=Glist$rsids, ncores=ncores) # compute allele frequencies, missingness, ....    
     Glist$af1 <- 1-Glist$af
     Glist$af2 <- Glist$af
     
     return(Glist)
}

#' @export
#'
writeRAW <- function(Glist=NULL, ids=NULL, rsids=NULL, overwrite=FALSE) {
     bed2raw( fnRAW=Glist$fnRAW, bedfiles=Glist$bedfiles, bimfiles=Glist$bimfiles, famfiles=Glist$famfiles, ids=ids, rsids=rsids, overwrite=overwrite)
}

#' @export
#'

bed2raw <- function(fnRAW=NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, ids=NULL, rsids=NULL,overwrite=FALSE) {
     if(file.exists(fnRAW)) {
          warning(paste("fnRAW file allready exist"))
          if(!overwrite) stop(paste("fnRAW file allready exist"))
     }
     #bim0 <- NULL
     for ( chr in 1:length(bedfiles)) {
          print(paste("Processing bedfile:",bedfiles[chr]))
          bim <- fread(input=bimfiles[chr], header=FALSE, data.table = FALSE, colClasses="character")
          fam <- fread(input=famfiles[chr], header=FALSE, data.table = FALSE)
          n <- nrow(fam)
          m <- nrow(bim)
          rsidsBIM <- as.character(bim[,2])
          keep <- rep(TRUE,m)
          if(!is.null(rsids)) keep <- rsidsBIM%in%rsids
          nbytes <- ceiling(n/4)
          printmarker <- rep(F,m)
          printmarker[seq(1,m,10000)] <- T
          fnBED <- bedfiles[chr]  
          bfBED <- file(fnBED,"rb")
          magic <- readBin(bfBED, "raw", n=3)
          if (!all(magic[1] == "6c", magic[2] == "1b", magic[3] == "01"))
               stop("Wrong magic number for bed file; should be -- 0x6c 0x1b 0x01 --.")
          close(bfBED)
          append <- 1 
          if(chr==1) append <- 0
          #if(chr==1) {
          #     append <- 0
          #     if(filetype=="bed") {
          #          bfRAW <- file(fnRAW,"wb")
          #          writeBin( as.raw(magic), bfRAW, size = 1, endian = "little")
          #          close(bfRAW)
          #          append <- 1
          #     }
          #}
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
          #bim0 <- rbind(bim0,bim) 
     }
     #fnBIM <- gsub(".raw",".bim",fnRAW)
     #fnFAM <- gsub(".raw",".fam",fnRAW)
     #write.table(bim0, file="fnBIM", col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
     #write.table(fam, file="fnFAM", col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
     
}


#' @export
#'
summaryRAW <- function(Glist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL, ncores=1) {
     
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
     if(OS=="windows") fnRAW <- tolower(gsub("/","\\",fnRAW,fixed=T))    
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
 
readbed <- function(Glist=NULL, bedfiles=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL,impute=TRUE, ncores=1) {
     if(!is.null(Glist)) {     
       n <- Glist$n
       m <- Glist$m
       nbytes <- ceiling(n/4)
       if(is.null(cls)) cls <- 1:m
       if(!is.null(rsids)) cls <- match(rsids,Glist$rsids)
       nc <- length(cls)
       if (is.null(rws)) rws <- 1:n
       if(!is.null(ids)) rws <- match(ids,Glist$ids)
       nr <- length(rws)
       fnRAW <- Glist$fnRAW
       OS <- .Platform$OS.type
       if(OS=="windows") fnRAW <- tolower(gsub("/","\\",fnRAW,fixed=T))
       ids <- Glist$ids[rws]
       rsids <- Glist$rsids[cls]
     }
     if(!is.null(bedfiles)) {
          fnRAW <- bedfiles
          bimfiles <- gsub(".bed",".bim",bedfiles)
          famfiles <- gsub(".bed",".fam",bedfiles)
          bim <- fread(input=bimfiles, header=FALSE, data.table = FALSE,showProgress=FALSE,colClasses="character")
          fam <- fread(input=famfiles, header=FALSE, data.table = FALSE,showProgress=FALSE,colClasses="character")
          n <- nrow(fam)
          m <- nrow(bim)
          nbytes <- ceiling(n/4)
          cls <- match(rsids,as.character(bim[,2]))
          if(any(is.na(cls))) {
               warning(paste("some rsids not found in bimfiles"))
               print(rsids[is.na(cls)])  
          }
          cls <- cls[!is.na(cls)]
          if(length(cls)==0) stop("No rsids found in bimfiles") 
          nc <- length(cls)
          if (is.null(rws)) rws <- 1:n
          if(!is.null(ids)) rws <- match(ids,as.character(fam[,2]))
          nr <- length(rws)
          OS <- .Platform$OS.type
          if(OS=="windows") fnRAW <- tolower(gsub("/","\\",fnRAW,fixed=T))
          ids <- as.character(fam[rws,2])
          rsids <- as.character(bim[cls,2])
     }
     W <- .Fortran("readbed", 
                     n = as.integer(n),
                     nr = as.integer(nr),
                     rws = as.integer(rws),
                     nc = as.integer(nc),
                     cls = as.integer(cls),
                     impute = as.integer(impute),
                     W = matrix(as.double(0),nrow=nr,ncol=nc),
                     nbytes = as.integer(nbytes),
                     fnRAW = as.character(fnRAW),
                     PACKAGE = 'qgg')$W
     rownames(W) <- ids
     colnames(W) <- rsids
     return(W)
}


#' @export
#'
getW <- function(Glist=NULL, ids=NULL, rsids=NULL, rws=NULL,cls=NULL, impute=FALSE) {
     if(is.null(ids)) ids <- Glist$ids
     if(is.null(cls)) cls <- match(rsids,Glist$rsids)
     W <- readbed(Glist=Glist,ids=ids,rsids=rsids,rws=rws,cls=cls,impute=impute) 
     rownames(W) <- ids
     colnames(W) <- Glist$rsids[cls]
     return(W)
}


