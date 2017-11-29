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
     
     LDlist <- NULL   
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
     LDlist$nchr <- length(LDlist$rsids)
     LDlist$mchr <- sapply(LDlist$rsids,length)
     return(LDlist)
     
}

#' @export
#'

computeLD <- function(Wlist=NULL, chr=NULL, fnLD=NULL, path=NULL, msize=100) {
     
     if(is.null(chr)) chr <- unique(Wlist$chr)
     study <- Wlist$study
     if(is.null(study)) study <- "STUDY_UNKNOWN"
     if (is.null(fnLD)) fnLD <- paste(path,"/LD_CHR",chr,"_",study,sep="")
     if(any(file.exists(fnLD[chr]))) stop("LD files allready exists - please specify other file names")
     
     n <- Wlist$n
     
     for ( i in chr) {

       fnW <- Wlist$fnRAWCHR[i]
          
       rsids <- Wlist$rsids[[i]]
       W1 <- matrix(logical(0),nrow=n,ncol=msize)
       W2 <- matrix(logical(0),nrow=n,ncol=msize)
       W3 <- matrix(logical(0),nrow=n,ncol=msize)
          
       bfLD <- file(fnLD[i],"wb")
       
       m <- length(rsids)
       
       msets <- sapply(split(1:m, ceiling(seq_along(1:m)/msize)),length)
       
       nsets <- length(msets)
       bfW <- file(fnW,"rb")
         for ( j in 1:nsets) {
           W1 <- W2
           W2 <- W3
           W3 <- matrix(0,nrow=n,ncol=msize)
           for ( k in 1:msets[j] ) {
             w <- as.double(readBin( bfW, "raw", n=n, size = 1, endian = "little"))
             w[w==0] <- NA
             w <- scale(w)[,1]
             w[is.na(w)] <- 0
             W3[,k] <- w
           }
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
         close(bfW)
         close(bfLD)
       }
     }
     
#     LDlist <- NULL
#     LDlist$study <- Wlist$study
#     LDlist$fnLD <- fnLD
#     LDlist$rsids <- unlist(Wlist$rsids) 
#     LDlist$ids <- Wlist$ids
#     LDlist$msize <- msize   
#     LDlist$n <- Wlist$n   
#     LDlist$m <- Wlist$m 
#     LDlist$chr <- unique(Wlist$chr)
#     LDlist$nchr <- length(LDlist$chr)
#     LDlist$rsids <- NULL
#     for ( chr in LDlist$chr) {
#          LDlist$rsids[[chr]] <- unlist(Wlist$rsids[Wlist$chr==chr])
#     }
#     LDlist$mchr <- sapply(LDlist$rsids,length)
#     LDlist$fnLD <- paste(path,"/LD_CHR",LDlist$chr,"_",LDlist$study,sep="")
#     return(LDlist)
     
#}

# computeLD <- function(Wlist=NULL, chr=NULL, fnLD=NULL, path=NULL, msize=100) {
#      
#      if(is.null(chr)) chr <- unique(Wlist$chr)
#      study <- Wlist$study
#      if(is.null(study)) study <- "STUDY_UNKNOWN"
#      if (is.null(fnLD)) fnLD <- paste(path,"/LD_CHR",chr,"_",study,sep="")
#      if(any(file.exists(fnLD))) stop("LD files allready exists - please specify other file names")
#      
#      n <- Wlist$n
#      
#      for ( chrom in chr) {
#           rws <- Wlist$chr==chrom 
#           
#           rsids <- Wlist$rsids[rws]
#           fnWChr <- Wlist$fnW[rws]
#           
#           nbchr <- length(fnWChr) 
#           
#           W1 <- matrix(logical(0),nrow=n,ncol=msize)
#           W2 <- matrix(logical(0),nrow=n,ncol=msize)
#           W3 <- matrix(logical(0),nrow=n,ncol=msize)
#           
#           bfLD <- file(fnLD[chrom],"wb")
#           for ( i in 1:nbchr ) {
#                m <- length(rsids[[i]])
#                msets <- sapply(split(1:m, ceiling(seq_along(1:m)/msize)),length)
#                nb <- length(msets)
#                fnW <- fnWChr[i]
#                bfW <- gzfile(fnW,"rb")
#                for ( j in 1:nb) {
#                     W1 <- W2
#                     W2 <- W3
#                     W3 <- matrix(0,nrow=n,ncol=msize)
#                     for ( k in 1:msets[j] ) {
#                          W3[,k] <- readBin( bfW, "double", n=n, size = 8, endian = "little")
#                     }
#                     #print(paste("Finished loading block",j))
#                     WW <- t(crossprod(cbind(W1,W2,W3),W2))
#                     N <- t(crossprod(!cbind(W1,W2,W3)==0,!W2==0))
#                     WW <- WW/(N-1)  # N-1 accounts for sample mean is estimated
#                     WW[is.na(WW)] <- 0
#                     if (i>1 | j>1) {
#                          for ( k in 1:msize ) { 
#                               ld <- as.vector(WW[k,k:(k+2*msize)]) 
#                               writeBin( ld, bfLD, size = 8, endian = "little")
#                          }
#                     }
#                     if (i==nbchr & j==nb) {
#                          W1 <- W2
#                          W2 <- W3
#                          W3 <- matrix(0,nrow=n,ncol=msize)
#                          WW <- t(crossprod(cbind(W1,W2,W3),W2))
#                          N <- t(crossprod(!cbind(W1,W2,W3)==0,!W2==0))
#                          WW <- WW/(N-1)  # N-1 accounts for sample mean is estimated
#                          WW[is.na(WW)] <- 0
#                          for ( k in 1:msets[j] ) { 
#                               ld <- as.vector(WW[k,k:(k+2*msize)]) 
#                               writeBin( ld, bfLD, size = 8, endian = "little")
#                          }
#                     }
#                }
#                close(bfW)
#                print(paste("Finished block",i,"chromosome",chrom))
#           }
#           close(bfLD)
#      }
#      
#      LDlist <- NULL
#      LDlist$study <- Wlist$study
#      LDlist$fnLD <- fnLD
#      LDlist$rsids <- unlist(Wlist$rsids) 
#      LDlist$ids <- Wlist$ids
#      LDlist$msize <- msize   
#      LDlist$n <- Wlist$n   
#      LDlist$m <- Wlist$m 
#      LDlist$chr <- unique(Wlist$chr)
#      LDlist$nchr <- length(LDlist$chr)
#      LDlist$rsids <- NULL
#      for ( chr in LDlist$chr) {
#           LDlist$rsids[[chr]] <- unlist(Wlist$rsids[Wlist$chr==chr])
#      }
#      LDlist$mchr <- sapply(LDlist$rsids,length)
#      LDlist$fnLD <- paste(path,"/LD_CHR",LDlist$chr,"_",LDlist$study,sep="")
#      return(LDlist)
#      
# }



#' @export
#'

getLDSets <- function(LDlist=NULL, chr=NULL, r2=0.5) {
     
     msize <- LDlist$msize
     
     ldSets <- NULL
     
     for ( chr in 1:LDlist$nchr ) {
          
          rsidsChr <- LDlist$rsids[[chr]]
          
          n <- LDlist$n
          mchr <- LDlist$mchr[[chr]]
          rsidsChr <- LDlist$rsids[[chr]]
          rsidsLD <- c(rep("start",msize),rsidsChr,rep("end",msize))
          
          fnLD <- LDlist$fnLD[chr]
          bfLD <- file(fnLD,"rb")
          ld <- readBin( bfLD, "double", n=mchr*5000, size = 8, endian = "little")
          ld <- matrix(ld,nrow=mchr,byrow=TRUE) 
          close(bfLD)
          ld[,msize+1] <- 1  
          
          ldSetsChr <- vector(length=mchr,mode="list")
          names(ldSetsChr) <- rsidsChr
          
          for ( i in 1:mchr) { 
               cls <- which(abs(ld[i,])>r2) + i - 1
               ldSetsChr[[i]] <- rsidsLD[cls]
          } 
          
          #nsets <- sapply(ldSetsChr,length)
          #fsets <- as.factor(rep(names(ldSetsChr),times=nsets))
          #rsidsChrSets <- unlist(ldSetsChr)
          #indxSets <- match(rsidsChrSets,unlist(LDlist$rsids))
          #mydf <- data.frame(indxSets,fsets)
          #ldSetsChr <- unstack( mydf, indxSets ~ fsets )
          #ldSetsChr <- split(indxSets,f=fsets)  
          #ldSetsChr <- ldSetsChr[rsidsChr]
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


#' @export
#'

adjustStat <- function( stat=NULL, ldSets=NULL, threshold=1) {
     
     for ( i in 1:nrow(stat$P) ) {
          rsidsStat <- colnames(stat$P)
          mStat <- length(rsidsStat)
          indx1 <- rep(T,mStat)
          indx2 <- rep(F,mStat)
          for ( chr in 1:length(ldSets) ) {
               setsChr <- ldSets[[chr]]
               setsChr <- setsChr[names(setsChr)%in%rsidsStat]
               rsidsChr <- names(setsChr)
               clsChr <- match(rsidsChr,rsidsStat)
               p <- stat$P[i,clsChr]
               s <- stat$S[i,clsChr]
               o <- order(p, decreasing=FALSE)
               for ( j in o) {
                    if (p[j]<=threshold) { 
                         if (indx1[clsChr[j]]) {
                              rws <- setsChr[[j]]
                              indx1[rws] <- F
                              indx2[clsChr[j]] <- T
                         }
                    }
               }
          }
          stat$S[i,!indx2] <- 0
     }
     return(stat)
}


#' @export
#'

computePRS <- function(Wlist=NULL,S=NULL,msize=100, scaled=TRUE) {
     PRS <- matrix(0,nrow=Wlist$n,ncol=nrow(S))
     rownames(PRS) <- Wlist$ids
     colnames(PRS) <- rownames(S)
     rsidsS <- colnames(S)
     if(any(!unlist(Wlist$rsids)==rsidsS)) stop("rsids not found in Wlist")
     maf <- unlist(Wlist$maf)
     meanW <- 2*maf
     sdW <- sqrt(2*maf*(1-maf))
     n <- Wlist$n
     m <- Wlist$m
     W <- matrix(double(0),nrow=n,ncol=msize) 
     sets <- split(1:m, ceiling(seq_along(1:m)/msize))
     nsets <- length(sets)
     msets <- sapply(sets,length)
     #bfRAW <- file(Wlist$fnRAW,"rb")
     for ( i in 1:nsets ) {
          cls <- sets[[i]]
          W <- readraw(Wlist=Wlist,cls=cls)
     #     for (j in 1:msets[i]) {
     #          w <- as.double(readBin( bfRAW, "raw", n=n, size = 1, endian = "little"))
     #          if(scaled) w[w>0] <- (w[w>0]-1-meanW[cls[j]])/sdW[cls[j]]
     #          if(!scaled) w[w>0] <- w[w>0]-1
     #          W[,j] <- w
     #     }
          if(nrow(S)>1) PRS <- PRS + tcrossprod(W[,1:msets[i]],S[,cls])
          if(nrow(S)==1) PRS <- PRS + tcrossprod(W[,1:msets[i]],t(S[,cls]))
          print(paste("Finished block",i,"out of",nsets))
     }
     #close(bfRAW)
     PRS
}



#######################################################################################


