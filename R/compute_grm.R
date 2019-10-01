#######################################################################################
# compute GRM functions
#######################################################################################
#'
#' Computing the genomic relationship matrix (GRM)
#'
#' @description
#' The grm function is used to compute a genomic relationship matrix (GRM) based on all,
#' or a subset of marker genotypes. GRM for additive, and non-additive (dominance and
#' epistasis) genetic models can be constructed. The output of the grm function can either be a
#' within-memory GRM object (n x n matrix), or a GRM-list which is a list structure that
#' contains information about the GRM stored in a binary file on the disk.
#'
#' @param Glist list providing information about genotypes stored on disk
#' @param GRMlist list providing information about GRM matrix stored in binary files on disk
#' @param ids vector of individuals used for computing GRM
#' @param rsids vector marker rsids used for computing GRM
#' @param rws rows in genotype matrix used for computing GRM
#' @param cls columns in genotype matrix used for computing GRM
#' @param W matrix of centered and scaled genotypes
#' @param scale logical if TRUE the genotypes in Glist has been scaled to mean zero and variance one
#' @param method indicator of method used for computing GRM: additive (add, default), dominance (dom) or epistasis (epi-pairs or epi-hadamard (all genotype markers))
#' @param msize number of genotype markers used for batch processing
#' @param ncores number of cores used to compute the GRM
#' @param fnG name of the binary file used for storing the GRM on disk
#' @param overwrite logical if TRUE the binary file fnG will be overwritten
#' @param returnGRM logical if TRUE function returns the GRM matrix to the R environment
#' @param task either computation of GRM (task="grm"  which is default) or eigenvalue decomposition of GRM (task="eigen")
#' @param miss the missing code (miss=0 is default) used for missing values in the genotype data
#'
#'
#' @return Returns a genomic relationship matrix (GRM) if returnGRM=TRUE else a list structure (GRMlist) with information about the GRM  stored on disk

#' @author Peter Soerensen

#' @examples
#'
#' # Simulate data
#' W <- matrix(rnorm(1000000), ncol = 1000)
#' 	colnames(W) <- as.character(1:ncol(W))
#' 	rownames(W) <- as.character(1:nrow(W))
#'
#' # Compute GRM
#' GRM <- grm(W = W)
#'
#' \donttest{
#' 
#' # Eigen value decompostion GRM
#' eig <- grm(GRM=GRM, task="eigen")
#'
#' }

#' @export
#'

grm <- function(Glist = NULL, GRMlist = NULL, ids = NULL, rsids = NULL, rws = NULL, cls = NULL,
                W = NULL, method = "add", scale = TRUE, msize = 100, ncores = 1, fnG = NULL,
                overwrite = FALSE, returnGRM = FALSE, miss = 0, task = "grm") {
  if (task == "grm") {
    GRM <- computeGRM(
      Glist = Glist, ids = ids, rsids = rsids, rws = rws, cls = cls,
      W = W, method = method, scale = scale, msize = msize, ncores = ncores,
      fnG = fnG, overwrite = overwrite, returnGRM = returnGRM, miss = miss
    )
    return(GRM)
  }
  if (task == "eigen") {
    eig <- eigenGRM(GRM = GRM, GRMlist = GRMlist, method = "default", ncores = ncores)
    return(eig)
  }
}




computeGRM <- function(Glist = NULL, ids = NULL, rsids = NULL, rws = NULL, cls = NULL, W = NULL, method = "add", scale = TRUE, msize = 100, ncores = 1, fnG = NULL, overwrite = FALSE, returnGRM = FALSE, miss = 0) {
  if (method == "add") gmodel <- 1
  if (method == "dom") gmodel <- 2
  if (method == "epi-pairs") gmodel <- 3
  if (method == "epi-hadamard") gmodel <- 4

  if (!is.null(W)) {
    SS <- tcrossprod(W) # compute crossproduct, all SNPs
    N <- tcrossprod(!W == miss) # compute number of observations, all SNPs
    G <- SS / N
    return(G)
  }
  if (is.null(W)) {
    n <- Glist$n
    m <- Glist$m
    nbytes <- ceiling(n / 4)

    if (is.null(cls)) cls1 <- cls2 <- 1:m
    if (is.list(cls)) {
      gmodel <- 3
      cls1 <- cls[[1]]
      cls2 <- cls[[2]]
    }
    if (is.list(rsids)) {
      gmodel <- 3
      cls1 <- match(rsids[[1]], Glist$rsids)
      cls2 <- match(rsids[[2]], Glist$rsids)
    }
    if (is.vector(rsids) & !is.list(rsids)) {
      cls1 <- cls2 <- match(rsids, Glist$rsids)
    }
    nc <- length(cls1)

    if (is.null(rws)) {
      rws <- 1:n
      if (!is.null(ids)) rws <- match(ids, Glist$ids)
    }
    nr <- length(rws)
    fnRAW <- Glist$fnRAW

    # Initiate GRMlist
    idsG <- Glist$ids[rws]
    rsidsG <- Glist$rsids[unique(c(cls1, cls2))]
    nG <- length(idsG)
    mG <- length(rsidsG)
    GRMlist <- list(fnG = fnG, idsG = idsG, rsids = rsidsG, n = nG, m = mG, method = method)

    # Initiate G file
    if (file.exists(fnG)) {
      if (!overwrite) stop("G file name allready exist")
    }
    fnG <- GRMlist$fnG
    OS <- .Platform$OS.type
    if (OS == "windows") fnRAW <- tolower(gsub("/", "\\", fnRAW, fixed = T))
    if (OS == "windows") fnG <- gsub("/", "\\", fnG, fixed = T)

    #write.table(c(as.character(fnRAW),as.character(fnG)), file = "param.qgg", quote = FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
    
    res <- .Fortran("grmbed",
      n = as.integer(n),
      nr = as.integer(nr),
      rws = as.integer(rws),
      nc = as.integer(nc),
      cls1 = as.integer(cls1),
      cls2 = as.integer(cls2),
      scale = as.integer(scale),
      nbytes = as.integer(nbytes),
      fnRAWCHAR = as.integer(unlist(sapply(as.character(fnRAW),charToRaw),use.names=FALSE)),
      nchars = nchar(as.character(fnRAW)),
      msize = as.integer(msize),
      ncores = as.integer(ncores),
      fnGCHAR = as.integer(unlist(sapply(as.character(fnG),charToRaw),use.names=FALSE)),
      ncharsg = nchar(as.character(fnG)),
      gmodel = as.integer(gmodel),
      # G = matrix(as.double(0),nrow=nr,ncol=nr),
      PACKAGE = "qgg"
    )
    #file.remove("param.qgg")
    if (!returnGRM) return(GRMlist)
    if (returnGRM) {
      GRM <- getGRM(GRMlist = GRMlist, ids = GRMlist$idsG)
      return(GRM)
    }
  }
}



writeGRM <- function(GRM = NULL) {
  if (!is.null(GRM)) {
    for (i in 1:length(GRM)) {
      fileout <- file(paste("G", i, sep = ""), "wb")
      nr <- nrow(GRM[[i]])
      for (j in 1:nr) {
        writeBin(as.double(GRM[[i]][1:nr, j]), fileout, size = 8, endian = "little")
      }
      close(fileout)
    }
  }
}


#' Extract elements from genomic relationship matrix (GRM) stored on disk
#'
#' @description
#' Extract elements from genomic relationship matrix (GRM) (whole or subset) stored on disk.

#' @param GRMlist list providing information about GRM matrix stored in binary files on disk
#' @param ids vector of ids in GRM to be extracted
#' @param idsRWS vector of row ids in GRM to be extracted
#' @param idsCLS vector of column ids in GRM to be extracted
#' @param rws vector of rows in GRM to be extracted
#' @param cls vector of columns in GRM to be extracted


#' @export
#'

getGRM <- function(GRMlist = NULL, ids = NULL, idsCLS = NULL, idsRWS = NULL, cls = NULL, rws = NULL) {

  # GRMlist(fnG=fnG,idsG=idsG,rsids=rsidsG,n=nG,m=mG)

  if (!is.null(ids)) {
    idsRWS <- idsCLS <- ids
  }
  if (!is.null(rws)) idsRWS <- GRMlist$idsG[rws]
  if (!is.null(cls)) idsCLS <- GRMlist$idsG[cls]
  if (is.null(idsRWS)) stop("Please specify ids or idsRWS and idsCLS or rws and cls")
  if (is.null(idsCLS)) stop("Please specify ids or idsRWS and idsCLS or rws and cls")

  if (sum(!idsCLS %in% GRMlist$idsG) > 0) stop("Error some ids not found in idsG")
  if (sum(!idsRWS %in% GRMlist$idsG) > 0) stop("Error some ids not found in idsG")

  rws <- match(idsRWS, GRMlist$idsG) # no reorder needed
  cls <- match(idsCLS, GRMlist$idsG)

  nG <- GRMlist$n # nG <- GRMlist$nG
  nr <- length(rws)
  nc <- length(cls)

  # Sub matrix
  G <- matrix(0, nrow = nr, ncol = nc)
  rownames(G) <- GRMlist$idsG[rws]
  colnames(G) <- GRMlist$idsG[cls]

  # If full stored project study matrix
  fileG <- file(GRMlist$fnG, "rb")
  for (i in 1:nr) {
    k <- rws[i]
    where <- (k - 1) * nG
    seek(fileG, where = where * 8)
    grws <- readBin(fileG, "double", n = nG, size = 8, endian = "little")
    G[i, ] <- grws[cls]
  }
  close(fileG)
  return(G)
}


#' Merge multiple GRMlist objects
#' 
#' @description
#' Merge multiple GRMlist objects each with information about a 
#' genomic rfelationship matrix stored on disk

#' @param GRMlist list providing information about GRM matrix stored in binary files on disk

#' @export
#'


mergeGRM <- function(GRMlist = NULL) {
  GRMlist <- do.call(function(...) mapply(c, ..., SIMPLIFY = FALSE), args = GRMlist)
  GRMlist$idsG <- unique(GRMlist$idsG)
  GRMlist$n <- length(GRMlist$idsG)
  GRMlist$m <- length(GRMlist$rsids)
  GRMlist
}




eigenGRM <- function(GRM = NULL, GRMlist = NULL, method = "default", ncores = 1) {
  n <- ncol(GRM)
  evals <- rep(0, n)
  res <- .Fortran("eiggrm",
    n = as.integer(n),
    GRM = matrix(as.double(GRM), nrow = n, ncol = n),
    evals = as.double(evals),
    ncores = as.integer(ncores),
    PACKAGE = "qgg"
  )
  list(values = res$evals, U = res$GRM)
}
