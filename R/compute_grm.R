#######################################################################################
# compute GRM functions
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
     rsidsG <- unlist(Wlist$rsids)[cls] 
     nG <- length(idsG)
     mG <- length(rsidsG) 
     Glist <- list(fnG=fnG,idsG=idsG,rsids=rsidsG,n=nG,m=mG)
     
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
       writeBin( res$G[1:nG,j], bfG, size = 8, endian = "little")
     }
     close(bfG)
     }
     if(is.null(fnG)) return(res$G)
     if(!is.null(fnG)) return(Glist)
}
