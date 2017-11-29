###############################################################################################
#  Compute LD for target population
###############################################################################################
#
# Description:
#
# Linkage disequilibrium (LD) information is often used in genomic analyses. LD information 
# may be obtained from the target population or from external ressources such as the 1000 
# human genomes project. Below are functions to compute LD based on on the genotype matrix W 
# and functions to access LD information as implemented in the qgg package
#
# The LD matrix is computed using the computeLD function based on the genotype matrix W
# (i.e. Wlist structure).
#
# The r2 metric used is the pairwise correlation between markers (allele count alternative allele) 
# in a specified region of the genome. The marker genotype is allele count of the alternative allele
# which is assumed to be centered and scaled. This is the default output format of 
# the computeW function.
#
# The LD matrix can be too big to fit into memory and is therefore stored in a binary file. 
# Information about the LD matrix is provided in a LDlist structure. The LD matrix can be access 
# using the getLD function. Sets of markers in LD with a specific marker can be obtained using 
# the getLDsets function. The markers sets can be mapped onto the markers in W using mapLDSets function. 
#
# Examples of simple usage:
#
# 	LDlist <- prepLD( Wlist, additional arguments...)
# 	LDlist <- computeLD( Wlist, msize, additional arguments...)
# 	LD <- getLD( LDlist, additional arguments...)
# 	ldSets <- getLDSets( LDlist, r2, additional arguments...)
#       ldSets <- mapLDSets( ldSets=ldSets, Wlist=Wlist)
#
# Required input (arguments):
#
# Wlist is a list structure providing information about the genotype matrix W
# LDlist is a list structure providing information about the LD matrix
# ldSets is a list of markers sets in LD at a given r2 threshold and genome region (i.e. msize)
# msize is the number of markers used up/down stream of the focus marker (i.e. LD region) 
# r2 is the threshold for pairwise correlation for including the marker in the LD sets  
#
# Optional input:
#
# ids is a vector of subject ids used in W 
# rsids is a vector of variant ids used in W 
#
# Output:
#
# The computeLD function computes the LD matrix and save it in a binary format. Furthermore the 
# computeLD and prepLD functions returns a LDlist structure containing information about the LD 
# matrix which can be used in downstream analyses. Most users do not use this directly. 
# Furthermore this structure also allow us to change the underlying data structures without noticeable 
# changes in the LDlist interface.  
# 
# The LDlist contains the following slots (which may change in future releases):
#
#   LDlist$fnLD vector of file names 
#   LDlist$rsids vector of rsids
#   LDlist$alleles vector of allele used as reference 
#   LDlist$chr vector chromosome
#   LDlist$nchr total number of chromosomes
#   LDlist$ids vector of subject ids used in W
#   LDlist$m total number of markers in W
#   LDlist$n total number if subjects
#   LDlist$msize number of markers in LD region
#   LDlist$study study name
#
###############################################################################################


#######################################################################################
# compute LD functions
#######################################################################################

#' @export
#'

prepLD <- function(Wlist=NULL,fnLD=NULL, msize=NULL) {
     
     if(is.null(msize)) msize <- 100
     
     LDlist <- NULL   
     LDlist$study <- Wlist$study
     LDlist$fnLD <- fnLD
     LDlist$ids <- Wlist$ids
     LDlist$msize <- msize   
     LDlist$n <- Wlist$n   
     LDlist$m <- Wlist$m 
     LDlist$chr <- 1:length(fnLD)
     LDlist$rsids <- Wlist$rsids
     LDlist$nchr <- length(LDlist$chr)
     LDlist$mchr <- sapply(LDlist$rsids,length)
     LDlist$fnRAW <- Wlist$fnRAW
     return(LDlist)

}

#' @export
#'

computeLD <- function(LDlist=NULL, chr=NULL, ids=NULL, scaled=TRUE) {
     
     if(is.null(chr)) chr <- LDlist$chr
     if(any(file.exists(fnLD[chr]))) stop("LD files allready exists - please specify other file names")
     
     n <- LDlist$n
     rws <- 1:n
     if(!is.null(ids)) rws <- match(ids,LDlist$ids)
     
     for ( i in chr) {

       rsids <- LDlist$rsids[[i]]
       cls <- match(rsids,unlist(LDlist$rsids))
       cls <- split(cls, ceiling(seq_along(cls)/msize))
       msets <- sapply(cls,length)
       nsets <- length(msets)
       m <- length(rsids)
       
       W1 <- matrix(logical(0),nrow=n,ncol=msize)
       W2 <- matrix(logical(0),nrow=n,ncol=msize)
       W3 <- matrix(logical(0),nrow=n,ncol=msize)
          
       bfLD <- file(LDlist$fnLD[i],"wb")

       for ( j in 1:nsets) {
           W1 <- W2
           W2 <- W3
           W3 <- readraw(Wlist=LDlist,cls=cls[[j]],scaled=scaled)
           #W3 <- W3[rws,]
           WW <- t(crossprod(cbind(W1,W2,W3),W2))
           N <- t(crossprod(!cbind(W1,W2,W3)==0,!W2==0))
           WW <- WW/(N-1)  # N-1 accounts for sample mean is estimated
           WW[is.na(WW)] <- 0
           if (j>1) {
            for ( k in 1:msize ) { 
              ld <- as.vector(WW[k,k:(k+2*msize)]) 
              writeBin( ld, bfLD, size = 8, endian = "little")
            }
           }
           if (j==nsets) {
            W1 <- W2
            W2 <- W3
            W3 <- matrix(0,nrow=n,ncol=msize)
            WW <- t(crossprod(cbind(W1,W2,W3),W2))
            N <- t(crossprod(!cbind(W1,W2,W3)==0,!W2==0))
            WW <- WW/(N-1)  # N-1 accounts for sample mean is estimated
            WW[is.na(WW)] <- 0
            for ( k in 1:msets[j] ) { 
               ld <- as.vector(WW[k,k:(k+2*msize)]) 
               writeBin( ld, bfLD, size = 8, endian = "little")
            }
           }
           print(paste("Finished block",j,"chromosome",i))
         }
         close(bfLD)
       }
}

#' @export
#'

getLDSets <- function(LDlist=NULL, chr=NULL, r2=0.5) {
     
     msize <- LDlist$msize
     
     ldSets <- NULL
     
     for ( chr in 1:LDlist$nchr ) {
          
          n <- LDlist$n
          mchr <- LDlist$mchr[chr]
          rsidsChr <- LDlist$rsids[[chr]]
          rsidsLD <- c(rep("start",msize),rsidsChr,rep("end",msize))
          
          fnLD <- LDlist$fnLD[chr]
          bfLD <- file(fnLD,"rb")
          nld <- as.integer(mchr*msize*2+1)
          ld <- readBin( bfLD, "double", n=nld, size = 8, endian = "little")
          ld <- matrix(ld,nrow=mchr,byrow=TRUE) 
          close(bfLD)
          ld[,msize+1] <- 1  
          
          ldSetsChr <- vector(length=mchr,mode="list")
          names(ldSetsChr) <- rsidsChr
          
          for ( i in 1:mchr) { 
               cls <- which(abs(ld[i,])>r2) + i - 1
               ldSetsChr[[i]] <- rsidsLD[cls]
          } 
          
          ldSets[[chr]] <- ldSetsChr
          
          print(paste("Finished chromosome",chr))
     }
     names(ldSets) <- 1:LDlist$nchr
     
     return(ldSets)
}

#' @export
#'

mapLDSets <- function( ldSets=NULL, rsids=NULL, Wlist=NULL, index=TRUE ) { 
     mpsets <- NULL
     if(!is.null(Wlist)) rsids <- unlist(Wlist$rsids)
     for ( chr in 1:length(ldSets) ) {
          rsSets <- ldSets[[chr]]
          nsets <- sapply(rsSets,length)
          rsChr <- rep(names(rsSets),times=nsets)
          rsSets <- unlist(rsSets,use.names=FALSE)
          rsSets <- match(rsSets,rsids)
          inW <- !is.na(rsSets)
          rsSets <- rsSets[inW]
          if(!index) rsSets <- rsids[rsSets]
          rsChr <-  rsChr[inW]
          rsSets <- split(rsSets,f=as.factor(rsChr))
          mpsets[[chr]] <- rsSets[Wlist$rsids[[chr]]]
     } 
     return(mpsets)
}

#' @export
#'

mapSets <- function( sets=NULL, rsids=NULL, Wlist=NULL, index=TRUE ) { 
     if(!is.null(Wlist)) rsids <- unlist(Wlist$rsids)
     nsets <- sapply(sets,length)
     rs <- rep(names(sets),times=nsets)
     rsSets <- unlist(sets,use.names=FALSE)
     rsSets <- match(rsSets,rsids)
     inW <- !is.na(rsSets)
     rsSets <- rsSets[inW]
     if(!index) rsSets <- rsids[rsSets]
     rs <-  rs[inW]
     rsSets <- split(rsSets,f=as.factor(rs))
     return(rsSets)
}

