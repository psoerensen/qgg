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


adjustStat <- function( stat=NULL, ldSets=NULL, threshold=1) {
     
     for ( i in 1:ncol(stat$P) ) {
          rsidsStat <- rownames(stat$P)
          mStat <- length(rsidsStat)
          indx1 <- rep(T,mStat)
          indx2 <- rep(F,mStat)
          for ( chr in 1:length(ldSets) ) {
               setsChr <- ldSets[[chr]]
               setsChr <- setsChr[names(setsChr)%in%rsidsStat]
               rsidsChr <- names(setsChr)
               rwsChr <- match(rsidsChr,rsidsStat)
               p <- stat$P[rwsChr,i]
               s <- stat$S[rwsChr,i]
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
          stat$S[!indx2,i] <- 0
     }
     return(stat)
}




#' @export
#'

computePRS <- function(Wlist=NULL,S=NULL,msize=100, scaled=TRUE) {
     if (is.vector(S)) S <- as.matrix(S)
     rsids <- rownames(S)
     if(any( !rsids%in%unlist(Wlist$rsids) )) stop("rsids not found in Wlist")
     PRS <- matrix(0,nrow=Wlist$n,ncol=ncol(S))
       rownames(PRS) <- Wlist$ids
       colnames(PRS) <- colnames(S)
     cls <- match(rsids,unlist(Wlist$rsids))
     m <- length(cls)
     cls <- split(cls, ceiling(seq_along(cls)/msize))
     rws <- split(1:m, ceiling(seq_along(1:m)/msize))
     nsets <- length(rws)
     msets <- sapply(rws,length)
     for ( i in 1:nsets ) {
          W <- readraw(Wlist=Wlist,cls=cls[[i]], scaled=scaled)
          PRS <- PRS + tcrossprod(W,t(S[rws[[i]],]))
          print(paste("Finished block",i,"out of",nsets))
     }
     PRS
}


