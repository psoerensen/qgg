#######################################################################################
# compute GRM functions
#######################################################################################

#' @export
#'


computeGRM <- function(Glist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL, W=NULL, method="add", scaled=TRUE, msize=100, ncores=1, fnG=NULL, overwrite=FALSE, returnGRM=FALSE, miss=0) {

     if(method=="add") gmodel <- 1 
     if(method=="dom") gmodel <- 2
     if(method=="epi") gmodel <- 3

     if (!is.null(W)) { 
          SS <- tcrossprod(W)                              # compute crossproduct, all SNPs
          N <- tcrossprod(!W == miss)                      # compute number of observations, all SNPs
          G <- SS /N
          return(G)
     }     
     if (is.null(W)) { 
          
     n <- Glist$n
     m <- Glist$m
     nbytes <- ceiling(n/4)
     
     if(is.null(cls)) cls1 <- cls2 <- 1:m
     if(is.list(cls)) {
          gmodel <- 3
          cls1 <- cls[[1]]
          cls2 <- cls[[2]]
     }
     if(is.list(rsids)) {
          gmodel <- 3
          cls1 <- match(rsids[[1]],Glist$rsids)
          cls2 <- match(rsids[[2]],Glist$rsids)
     }
     if(is.vector(rsids)&!is.list(rsids)) {
          cls1 <- cls2 <- match(rsids,Glist$rsids)
     }
     nc <- length(cls1)
     
     if (is.null(rws)) {
       rws <- 1:n
       if(!is.null(ids)) rws <- match(ids,Glist$ids)
     }
     nr <- length(rws)
     fnRAW <- Glist$fnRAW

     # Initiate GRMlist
     idsG <- Glist$ids[rws]
     rsidsG <- Glist$rsids[unique(c(cls1,cls2))] 
     nG <- length(idsG)
     mG <- length(rsidsG) 
     GRMlist <- list(fnG=fnG,idsG=idsG,rsids=rsidsG,n=nG,m=mG,method=method)
     
     # Initiate G file
     if(file.exists(fnG)) {
       if(!overwrite) stop("G file name allready exist")
     }
     fnG <- GRMlist$fnG
     OS <- .Platform$OS.type
     if(OS=="windows") fnRAW <- tolower(gsub( "/", "\\" , fnRAW, fixed=T ))    
     if(OS=="windows") fnG <- gsub( "/", "\\" , fnG, fixed=T )    
     res <- .Fortran("grmbed", 
           n = as.integer(n),
           nr = as.integer(nr),
           rws = as.integer(rws),
           nc = as.integer(nc),
           cls1 = as.integer(cls1),
           cls2 = as.integer(cls2),
           scaled = as.integer(scaled),
           nbytes = as.integer(nbytes),
           fnRAW = as.character(fnRAW),
           msize = as.integer(msize),
           ncores = as.integer(ncores),
           fnG = as.character(fnG),
           gmodel=as.integer(gmodel),
           #G = matrix(as.double(0),nrow=nr,ncol=nr),
           PACKAGE = 'qgg'
     )
     if(!returnGRM) return(GRMlist)
     if(returnGRM) {
          GRM <- getGRM( GRMlist=GRMlist, ids=GRMlist$idsG)
          return(GRM)
     }
     }
     
}



#' @export
#'

getGRM <- function( GRMlist=NULL,ids=NULL, idsCLS=NULL, idsRWS=NULL, cls=NULL,rws=NULL) {
     
     #GRMlist(fnG=fnG,idsG=idsG,rsids=rsidsG,n=nG,m=mG)
     
     if (!is.null(ids)) {
       idsRWS <- idsCLS <- ids     
     }
     if (!is.null(rws)) idsRWS <- GRMlist$idsG[rws]
     if (!is.null(cls)) idsCLS <- GRMlist$idsG[cls]
     if (is.null(idsRWS)) stop("Please specify ids or idsRWS and idsCLS or rws and cls")
     if (is.null(idsCLS)) stop("Please specify ids or idsRWS and idsCLS or rws and cls")
     
     if (sum(!idsCLS%in%GRMlist$idsG)>0) stop("Error some ids not found in idsG")
     if (sum(!idsRWS%in%GRMlist$idsG)>0) stop("Error some ids not found in idsG")
     
     rws <- match(idsRWS,GRMlist$idsG)		# no reorder needed
     cls <- match(idsCLS,GRMlist$idsG)
     
     nG <- GRMlist$n  # nG <- GRMlist$nG
     nr <- length(rws)
     nc <- length(cls)
     
     # Sub matrix
     G <- matrix(0,nrow=nr,ncol=nc)
     rownames(G) <- GRMlist$idsG[rws]
     colnames(G) <- GRMlist$idsG[cls]
     
     # If full stored project study matrix
     fileG <- file(GRMlist$fnG,"rb")
     for (i in 1:nr) {
          k <- rws[i]
          where <- (k-1)*nG 
          seek(fileG,where=where*8)
          grws <- readBin( fileG, "double", n=nG, size = 8, endian = "little")
          G[i,] <- grws[cls]
     }
     close(fileG)
     return(G)
}



#' @export
#'

   mergeGRM <- function(GRMlist=NULL) {
     GRMlist <- do.call(function(...) mapply(c,...,SIMPLIFY = FALSE),args = GRMlist)
     GRMlist$idsG <- unique(GRMlist$idsG)
     GRMlist$n <- length(GRMlist$idsG)
     GRMlist$m <- length(GRMlist$rsids)
     GRMlist
   }
