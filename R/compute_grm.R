###############################################################################################
#  Compute GRM 
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


computeG <- function(Wlist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL,scaled=TRUE, msize=100, ncores=1, fnG=NULL, updateG=FALSE) {

 
     n <- Wlist$n
     m <- Wlist$m
     nbytes <- ceiling(n/4)
     if(is.null(cls)) cls <- 1:m
     if(!is.null(rsids)) cls <- match(rsids,unlist(Wlist$rsids))
     nc <- length(cls)
     if (is.null(rws)) {
          rws <- 1:n
          if(!is.null(ids)) rws <- match(ids,Wlist$ids)
     }
     nr <- length(rws)
     fnRAW <- Wlist$fnRAW

     # Initiate Glist
     idsG <- Wlist$ids[rws]
     rsidsG <- Wlist$rsids[cls] 
     nG <- length(idsG)
     mG <- length(rsidsG) 
     Glist <- list(fnG=fnG,ids=idsG,rsids=rsidsG,n=nG,m=mG)
     
     # Initiate G file
     if(!is.null(fnG)) {
          if (!updateG)  {
               if(file.exists(fnG)) stop("G file name allready exist")
               bfG <- file(fnG,"wb")
               seek(bfG,8*(nG**2))
               writeBin(raw(8),bfG)
               close(bfG)
          }
     }

     res <- .Fortran("grmbed", 
                     n = as.integer(n),
                     nr = as.integer(nr),
                     rws = as.integer(rws),
                     nc = as.integer(nc),
                     cls = as.integer(cls),
                     scaled = as.integer(scaled),
                     nbytes = as.integer(nbytes),
                     fnRAW = as.character(fnRAW),
                     msize = as.integer(msize),
                     ncores = as.integer(ncores),
                     G = matrix(as.double(0),nrow=nr,ncol=nr),
                     PACKAGE = 'qgg'
                     
     )
     rownames(res$G) <- colnames(res$G) <- idsG
          if(!is.null(fnG)) {
               paste("Writing G to file")
               idsG1 <- rownames(G)
               idsG2 <- colnames(G)
               n1 <- length(idsG1)     
               n2 <- length(idsG2)
               indx1 <- match(idsG1,Glist$idsG)		# no reorder needed
               indx2 <- match(idsG2,Glist$idsG)		# no reorder needed
               bfG <- file(Glist$fnG,"r+b")
               for (j in 1:n1) {
                    k <- indx1[j]
                    where <- (k-1)*Glist$n + indx2[1]-1
                    seek(bfG,where=where*8, rw="write")
                    writeBin( G[j,1:n2]/m, bfG, size = 8, endian = "little")
                    gc()
               }
               close(bfG)
          }

     if(is.null(fnG)) return(res$G)
     if(!is.null(fnG)) return(Glist)
}
