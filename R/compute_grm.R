#######################################################################################
# compute GRM functions
#######################################################################################

#' @export
#'


computeG <- function(Wlist=NULL,ids=NULL,rsids=NULL,rws=NULL,cls=NULL,scaled=TRUE, msize=100, ncores=1, fnG=NULL, overwrite=FALSE) {

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
     rsidsG <- unlist(Wlist$rsids)[cls] 
     nG <- length(idsG)
     mG <- length(rsidsG) 
     Glist <- list(fnG=fnG,idsG=idsG,rsids=rsidsG,n=nG,m=mG)
     
     # Initiate G file
     if(!is.null(fnG)) {
       if(file.exists(fnG)) {
         if(!overwrite) stop("G file name allready exist")
       }
       if(file.exists(fnG)) {
       if (overwrite)  {
         bfG <- file(fnG,"wb")
         seek(bfG,8*(nG**2))
         writeBin(raw(8),bfG)
         close(bfG)
       }
       }
       if(!file.exists(fnG)) {
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
     #idsG1 <- rownames(G)
     #idsG2 <- colnames(G)
     #n1 <- length(idsG1)     
     #n2 <- length(idsG2)
     #indx1 <- match(idsG1,Glist$idsG)		# no reorder needed
     #indx2 <- match(idsG2,Glist$idsG)		# no reorder needed
     #bfG <- file(Glist$fnG,"r+b")
     #for (j in 1:n1) {
     #  k <- indx1[j]
     #  where <- (k-1)*Glist$n + indx2[1]-1
     #  seek(bfG,where=where*8, rw="write")
     #  writeBin( res$G[j,1:n2], bfG, size = 8, endian = "little")
     #  gc()
     #}
     bfG <- file(Glist$fnG,"wb")
     for (j in 1:nG) {
       writeBin( as.double(res$G[1:nG,j]), bfG, size = 8, endian = "little")
     }
     close(bfG)
     }
     if(is.null(fnG)) return(res$G)
     if(!is.null(fnG)) return(Glist)
}



#' @export
#'

getG <- function( Glist=NULL,ids=NULL, idsCLS=NULL, idsRWS=NULL, cls=NULL,rws=NULL) {
     
     #Glist(fnG=fnG,idsG=idsG,rsids=rsidsG,n=nG,m=mG)
     
     if (!is.null(ids)) {
       idsRWS <- idsCLS <- ids     
     }
     if (!is.null(rws)) idsRWS <- Glist$idsG[rws]
     if (!is.null(cls)) idsCLS <- Glist$idsG[cls]
     if (is.null(idsRWS)) stop("Please specify ids or idsRWS and idsCLS or rws and cls")
     if (is.null(idsCLS)) stop("Please specify ids or idsRWS and idsCLS or rws and cls")
     
     if (sum(!idsCLS%in%Glist$idsG)>0) stop("Error some ids not found in idsG")
     if (sum(!idsRWS%in%Glist$idsG)>0) stop("Error some ids not found in idsG")
     
     rws <- match(idsRWS,Glist$idsG)		# no reorder needed
     cls <- match(idsCLS,Glist$idsG)
     
     nG <- Glist$n  # nG <- Glist$nG
     nr <- length(rws)
     nc <- length(cls)
     
     # Sub matrix
     G <- matrix(0,nrow=nr,ncol=nc)
     rownames(G) <- Glist$idsG[rws]
     colnames(G) <- Glist$idsG[cls]
     
     # If full stored project study matrix
     fileG <- file(Glist$fnG,"rb")
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

   mergeG <- function(Glist=NULL) {
     Glist <- do.call(function(...) mapply(c,...,SIMPLIFY = FALSE),args = Glist)
     Glist$idsG <- unique(Glist$idsG)
     Glist$n <- length(Glist$idsG)
     Glist$m <- length(Glist$rsids)
     Glist
   }
