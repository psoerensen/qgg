####################################################################################################################
#    Module 6: GSCORE
####################################################################################################################
#'
#' Genomic prediction based on single marker summary statistics
#'
#'
#' @description
#' The gscore function is used for genomic predictions based on single marker summary statistics
#' (coefficients, log-odds ratios, z-scores) and observed genotypes.
#'

#' @param stat matrix of single marker effects
#' @param Glist list of information about genotype matrix
#' @param bedfiles name of the PLINK bed-files
#' @param famfiles name of the PLINK fam-files
#' @param bimfiles name of the PLINK bim-files
#' @param ids vector of individuals used in the analysis
#' @param scale logical if TRUE the genotype markers have been scale to mean zero and variance one
#' @param impute logical if TRUE missing genotypes are set to its expected value (2*af where af is allele frequency)
#' @param msize number of genotype markers used for batch processing
#' @param ncores number of cores used in the analysis

#' @author Peter Soerensen

#' @examples
#'

#' bedfiles <- system.file("extdata", "sample_22.bed", package = "qgg")
#' bimfiles <- system.file("extdata", "sample_22.bim", package = "qgg")
#' famfiles <- system.file("extdata", "sample_22.fam", package = "qgg")
#' 
#' Glist <- gprep(study="1000G", bedfiles=bedfiles, bimfiles=bimfiles,
#'                famfiles=famfiles, overwrite=TRUE)
#' 
#' rsids <- Glist$rsids
#' stat <- data.frame(rsids=Glist$rsids,alleles=Glist$a2, af=Glist$af, effect=rnorm(Glist$m))
#' 
#' W <- getW(Glist=Glist,rsids=Glist$rsids)
#' pgs1 <- W%*%stat[,4]
#' 
#' pgs2 <- gscore(Glist = Glist, stat = stat) 
#' 
#' pgs3 <- gscore(bedfiles=bedfiles, stat = stat) 
#' 
#' pgs4 <- gscore(bedfiles=bedfiles,bimfiles=bimfiles,famfiles=famfiles, stat = stat) 
#' 
#' 
#' cor(cbind(pgs1,pgs2,pgs3,pgs4))
#'

#' @export
#'

gscore <- function(Glist = NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, stat = NULL, ids = NULL, scale = TRUE, impute = TRUE, msize = 100, ncores = 1) {
     
     if ( !is.null(Glist))  {
          for (chr in 1:length(Glist$bedfiles)) {
               prschr <- run_gscore(bedfiles=Glist$bedfiles[chr], bimfiles=Glist$bimfiles[chr], famfiles=Glist$famfiles[chr], stat = stat, 
                                    ids = ids, scale = scale, impute = impute, msize = msize, ncores = ncores)
               if (chr==1) prs <- prschr
               if (chr>1) prs <- prs + prschr
          }
     }
     if ( !is.null(bedfiles))  {
          prs <- run_gscore(bedfiles=bedfiles, bimfiles=bimfiles, famfiles=famfiles, stat = stat, 
                            ids = ids, scale = scale, impute = impute, msize = msize, ncores = ncores)
     }   
     return(prs)
}


run_gscore <- function(Glist = NULL, bedfiles=NULL, bimfiles=NULL, famfiles=NULL, stat = NULL, ids = NULL, scale = TRUE, impute = TRUE, msize = 100, ncores = 1) {
     
     if(sum(is.na(stat))>0) stop("stat object contains NAs") 
     if(is.null(Glist) & is.null(bedfiles)) stop("Please provide Glist or bedfile")
     
     if (!is.null(bedfiles)) {
          if(!file.exists(bedfiles)) stop(paste("bedfiles does not exists:"),bedfiles) 
          Glist <- NULL
          Glist$bedfiles <- bedfiles[1]
          if (!is.null(bimfiles)) Glist$bimfiles <- bimfiles[1]
          if (!is.null(famfiles)) Glist$famfiles <- famfiles[1]
          if (is.null(bimfiles)) Glist$bimfiles <- gsub(".bed", ".bim", bedfiles[1])
          if (is.null(famfiles)) Glist$famfiles <- gsub(".bed", ".fam", bedfiles[1])
          if(!file.exists(Glist$bedfiles)) stop(paste("bedfiles does not exists:"),Glist$famfiles) 
          if(!file.exists(Glist$famfiles)) stop(paste("famfiles does not exists:"),Glist$famfiles) 
          if(!file.exists(Glist$bimfiles)) stop(paste("bimfiles does not exists:"),Glist$bimfiles) 
          
          # Read fam information
          fam <- fread(input = Glist$famfiles[1], header = FALSE, data.table = FALSE, colClasses = "character")
          Glist$ids <- as.character(fam[, 2])
          if (any(duplicated(Glist$ids))) stop("Duplicated ids found in famfiles")
          message(paste("Finished processing fam file", famfiles[1]))
          
          # Read bim information
          bim <- fread(input = Glist$bimfiles[1], header = FALSE, data.table = FALSE, colClasses = "character")
          Glist$rsids <- as.character(bim[, 2])
          Glist$a1 <- as.character(bim[, 5])
          Glist$a2 <- as.character(bim[, 6])
          Glist$position <- as.numeric(bim[, 4])
          Glist$chr <- as.character(bim[, 1])
          message(paste("Finished processing bim file", bimfiles[1]))
          
          Glist$n <- length(Glist$ids)
          Glist$m <- length(Glist$rsids)
          
     }
     
     # Prepase summary stat
     if (!sum(colnames(stat)[1:3] == c("rsids", "alleles", "af")) == 3) {
          stop("First three columns in data frame stat should be: rsids, alleles, af ")
     }
     rsidsOK <- stat$rsids %in% Glist$rsids
     if (any(!rsidsOK)) {
          warning("Some variants not found in genotype files")
          message(paste("Number of variants missing;", sum(!rsidsOK)))
          message(paste("Number of variants used;", sum(rsidsOK)))
          stat <- stat[rsidsOK, ]
          stat$rsids <- as.character(stat$rsids)
          stat$alleles <- as.character(stat$alleles)
     }
     S <- stat[, -c(1:3)]
     if (is.vector(S)) S <- as.matrix(S)
     S <- apply(S, 2, as.numeric)
     colnames(S) <- colnames(stat)[-c(1:3)]
     rsids <- as.character(stat$rsids)
     af <- stat$af
     
     # Prepare input data for mpgrs
     rws <- 1:Glist$n
     if (!is.null(ids)) rws <- match(ids, Glist$ids)
     nr <- length(rws)
     cls <- match(rsids, Glist$rsids)
     nc <- length(cls)
     direction <- as.integer(stat$alleles == Glist$a2[cls])
     
     if(!file.exists(Glist$bedfiles)) stop(paste("bed file does not exists:"),Glist$bedfiles) 
     

     # single core
     if(ncores==1) {
          Slist <- vector(ncol(S),mode="list")
          for (j in 1:ncol(S)) {
               Slist[[j]] <- S[,j]
          }
          #af <- 1-Glist$af[cls]
          grs <- .Call("_qgg_mtgrsbed", Glist$bedfiles, Glist$n, cls, af, Slist)
          grs <- as.matrix(as.data.frame(grs[[1]]))
          rownames(grs) <- Glist$ids
          colnames(grs) <- colnames(S)
     }

          
     # multicore     
     if(ncores>1) {
          size <- ceiling(length(cls)/ncores)
          af <- 1-Glist$af[cls]
          af <- split(af, ceiling(seq_along(af) / size))
          cls <- split(cls, ceiling(seq_along(cls) / size))
          rwsS <- nrow(S)
          rwsS <- split(rwsS, ceiling(seq_along(rwsS) / size))
          Slist <- vector(length(rwsS),mode="list")
          for (j in 1:length(rwsS)) {
               Slist[[j]] <- lapply(seq_len(ncol(S)), function(i) S[rwsS[[j]],i])
          }
          
          grslist <- mclapply(1:length(cls), function(set) { .Call("_qgg_mtgrsbed", Glist$bedfiles, Glist$n, cls[[set]], af[[set]], Slist[[set]]  ) }, mc.cores = ncores)
          grs <- as.matrix(as.data.frame(grslist[[1]]))
          for (j in 2:length(grslist)) {
               grs <- grs + as.matrix(as.data.frame(grslist[[j]]))      
          }
          rownames(grs) <- Glist$ids
          colnames(grs) <- colnames(S)
          # end multicore
     }
     
     return(grs)
}
