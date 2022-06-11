####################################################################################################################
#    Module 8: GSIM
####################################################################################################################
#'
#' Genomic simulation
#' 
#'
#' @description
#' The gsolve function is used for solving of linear mixed model equations. The algorithm used to solve the equation
#' system is based on a Gauss-Seidel (GS) method (matrix-free with residual updates) that handles large data sets.
#'
#' The linear mixed model fitted can account for multiple traits, multiple genetic factors (fixed or random genetic
#' marker effects), adjust for complex family relationships or population stratification, and adjust for other
#' non-genetic factors including lifestyle characteristics. Different genetic architectures (infinitesimal,
#' few large and many small effects) is accounted for by modeling genetic markers in different sets as fixed or
#' random effects and by specifying individual genetic marker weights.

#'
#' @param y vector or matrix of phenotypes
#' @param X design matrix of fixed effects
#' @param W matrix of centered and scaled genotypes
#' @param Glist list of information about genotype matrix stored on disk


#' @author Peter Soerensen

#' @examples
#'
#' # Simulate data
#' W <- matrix(rnorm(1000000), ncol = 1000)
#' 	colnames(W) <- as.character(1:ncol(W))
#' 	rownames(W) <- as.character(1:nrow(W))
#' m <- ncol(W)
#' causal <- sample(1:ncol(W),50)
#' y <- rowSums(W[,causal]) + rnorm(nrow(W),sd=sqrt(50))
#'
#' X <- model.matrix(y~1)
#'
#' Sg <- 50
#' Se <- 50
#' h2 <- Sg/(Sg+Se)
#' lambda <- Se/(Sg/m)
#' lambda <- m*(1-h2)/h2
#'
#' # BLUP of single marker effects and total genomic effects based on Gauss-Seidel procedure
#' fit <- gsolve( y=y, X=X, W=W, lambda=lambda)
#'





#'
#' @export
#'

qcstat <- function(Glist=NULL, stat=NULL, filename=NULL, 
                   excludeMAF=0.01, excludeMAFDIFF=0.05, excludeINFO=0.8, 
                   excludeCGAT=TRUE, excludeINDEL=TRUE, excludeDUPS=TRUE, excludeMHC=FALSE,
                   excludeMISS=0.05, excludeHWE=1e-12) {
  
  # we use cpra to link sumstats and Glist
  cpra <- unlist(Glist$cpra)
  rsids <- unlist(Glist$rsids)

  # stat is a data.frame
  if(!is.data.frame(stat)) stop("stat should be  a data frame")
  #if(!is.null(stat$marker)) rownames(stat) <- stat$marker
  if(!is.null(stat$rsids)) rownames(stat) <- stat$rsids
  
  # internal summary statistic column format
  # data.frame(rsids, chr, pos, a1, a2, af, b, seb, stat, p, n)     (single trait)
  # list(marker=(rsids, chr, pos, a1, a2, af), b, seb, stat, p, n)  (multiple trait)
  
  fm_internal <- c("rsids","chr","pos","a1","a2","af","b","seb","p","n")
  
  fm_external1 <- c("marker","chromosome", "position", "effect_allele", "non_effect_allele", 
                   "effect_allele_freq","effect", "effect_se", "effect_p", "effect_n")
  
  fm_external2 <- c("marker","chromosome", "position", "effect_allele", "non_effect_allele", 
                    "effect_allele_freq","effect", "effect_se", "effect_p")


  #fm_external4 <- c("marker","chromosome", "position", "effect_allele", "non_effect_allele", 
  #                  "effect_allele_freq","OR", "OR_se", "OR_p", "OR_n")
  
  format <- "unknown"

  if(all(fm_internal%in%colnames(stat))) format <- "internal"

  if(all(fm_external1%in%colnames(stat))) {
    format <- "external"
    fm_external <- fm_external1
  }
  
  if(all(fm_external2%in%colnames(stat))) {
    format <- "external"
    fm_external <- fm_external2
    fm_internal <- fm_internal[1:9]
  }
  
  if(format=="unknown") {
    message("Column headings for stat object not found")
    message("Column headings for stat object should be:")
    print(fm_external1)
    message("or:")
    print(fm_external2)
    message("or:")
    print(fm_internal)
    stop("please revised your stat object according to these ")
  }
  
  if(format=="external") {
    stat <- stat[,fm_external]
    cpra1 <- paste(stat[,"chromosome"],stat[,"position"],stat[,"effect_allele"],stat[,"non_effect_allele"],sep="_")
    cpra2 <- paste(stat[,"chromosome"],stat[,"position"],stat[,"non_effect_allele"],stat[,"effect_allele"],sep="_")

    mapped <- cpra1%in%cpra | cpra2%in%cpra
    message("Map markers based on cpra")
    message(paste("Number of markers in stat mapped to marker ids in Glist:",sum(mapped)))
    message(paste("Number of markers in stat not mapped to marker ids in Glist:",sum(!mapped)))
    
    stat <- stat[mapped,]
    cpra1 <- cpra1[mapped]
    cpra2 <- cpra2[mapped]
    rws1 <- match(cpra1,cpra)
    rws2 <- match(cpra2,cpra)
    
    stat$marker[!is.na(rws1)] <- rsids[rws1[!is.na(rws1)]]
    stat$marker[!is.na(rws2)] <- rsids[rws2[!is.na(rws2)]]
    
    isdup <- duplicated(stat$marker)
    if(any(isdup)) message("Removing markers with duplicated ids")
    if(any(isdup)) message(paste("Number of markers duplicated in stat:",sum(isdup)))
    stat <- stat[!isdup,] 
    rownames(stat) <- stat$marker    
    
  }
  

  # external summary statistic column format
  # optimal format:
  # marker, chromosome, position, effect_allele, non_effect_allele, 
  # effect_allele_freq, effect, effect_se, statistic, p, n
  # (which will allow best quality control)
  #
  # minimal option 1:
  # marker, effect_allele, effect, effect_se, n   (limited quality control)
  #
  # minimal option 2:
  # marker, effect_allele, sign, p, n             (limited quality control)
  
  marker <- data.frame(rsids=unlist(Glist$rsids),cpra=unlist(Glist$cpra),
                       chr=unlist(Glist$chr), pos=unlist(Glist$position), 
                       a1=unlist(Glist$a1), a2=unlist(Glist$a2),
                       af=unlist(Glist$af),stringsAsFactors = FALSE)
  
  rownames(marker) <- marker$rsids
  
  message("Filtering markers based on information in Glist:")
  message("")
  
  rsids <-  gfilter(Glist = Glist,
                    excludeMAF=excludeMAF, 
                    excludeMISS=excludeMISS, 
                    excludeCGAT=excludeCGAT, 
                    excludeINDEL=excludeINDEL, 
                    excludeDUPS=excludeDUPS, 
                    excludeHWE=excludeHWE, 
                    excludeMHC=excludeMHC)
  marker <- marker[marker$rsids%in%rsids,]
  message("")
  
  if(!is.null(Glist$rsidsLD)) {
    rsids <- unlist(Glist$rsidsLD)
    message(paste("Number of markers in sparse LD matrices:", sum(marker$rsids%in%rsids)))
    message("")
    marker <- marker[marker$rsids%in%rsids,]
  }
  
  message("Filtering markers based on information in stat:")
  message("")
  
  if(!is.null(stat$rsids)) marker_in_stat <- marker$rsids%in%stat$rsids
  if(!is.null(stat$marker)) marker_in_stat <- marker$rsids%in%stat$marker
  message(paste("Number of markers in stat also found in bimfiles:", sum(marker_in_stat)))
  message("")
  if(sum(marker_in_stat)==0) stop("No marker ids found in bimfiles")
  
  # align marker and stat object
  marker <- marker[marker_in_stat,]
  stat <- stat[marker$rsids,]
  
  if(!is.null(stat$effect_allele)) aligned <- stat$effect_allele==marker$a1
  if(!is.null(stat$a1)) aligned <- stat$a1==marker$a1
  message(paste("Number of effect alleles aligned with first allele in bimfiles:", sum(aligned)))
  message(paste("Number of effect alleles not aligned with first allele in bimfiles:", sum(!aligned)))
  message("")


  if(format=="external") {
    
    #original
    effect <- stat[,"effect"]
    effect_allele <- stat[,"effect_allele"]
    non_effect_allele <- stat[,"non_effect_allele"]
    effect_allele_freq <- stat[,"effect_allele_freq"]
    
    # aligned
    stat[!aligned,"effect"] <- -effect[!aligned]
    stat[!aligned,"effect_allele"] <- non_effect_allele[!aligned]
    stat[!aligned,"non_effect_allele"] <- effect_allele[!aligned] 
    stat[!aligned,"effect_allele_freq"] <- 1-effect_allele_freq[!aligned]
    excludeMAFDIFF <- abs(marker$af-stat$effect_allele_freq) > excludeMAFDIFF

    message(paste("Number of markers excluded by large difference between MAF difference:", sum(excludeMAFDIFF)))
    message("")
    
    stat <- stat[!excludeMAFDIFF,]
    marker <- marker[!excludeMAFDIFF,]
    colnames(stat) <- fm_internal
    if(is.null(stat$n)) stat$n <- neff(seb=stat$seb,af=stat$af)
  }  
  
  if(format=="internal") {
    
    #original
    effect <- stat[,"b"]
    effect_allele <- stat[,"a1"]
    non_effect_allele <- stat[,"a2"]
    effect_allele_freq <- stat[,"af"]
    
    # aligned
    stat[!aligned,"b"] <- -effect[!aligned]
    stat[!aligned,"a1"] <- non_effect_allele[!aligned]
    stat[!aligned,"a2"] <- effect_allele[!aligned] 
    stat[!aligned,"af"] <- 1-effect_allele_freq[!aligned]
    excludeMAFDIFF <- abs(marker$af-stat$af) > excludeMAFDIFF
    
    message(paste("Number of markers excluded by large difference between MAF difference:", sum(excludeMAFDIFF)))
    message("")
    
    stat <- stat[!excludeMAFDIFF,]
    marker <- marker[!excludeMAFDIFF,]
    if(is.null(stat$n)) stat$n <- neff(seb=stat$seb,af=stat$af)
    
  }  

  if(-!is.null(stat$info)) {
    lowINFO <- stat$info < excludeINFO
    message(paste("Number of markers excluded by low INFO score:", sum(lowINFO)))
    message("")
    stat <- stat[!lowINFO,]
  }  
  return(stat)
}


# checkStat <- function(Glist=NULL, stat=NULL, filename=NULL, excludeMAF=0.01, excludeMAFDIFF=0.05, excludeINFO=0.8, 
#                       excludeCGAT=TRUE, excludeINDEL=TRUE, excludeDUPS=TRUE, excludeMHC=FALSE) {
#   # effect, effect_se, effect_allele, alternative_allele, effect_allele_freq, nobs
#   cpra <- paste(unlist(Glist$chr),unlist(Glist$position),unlist(Glist$a1),unlist(Glist$a2), sep="_")
#   df <- data.frame(rsids=unlist(Glist$rsids),cpra,
#                    chr=unlist(Glist$chr), position=unlist(Glist$position), 
#                    a1=unlist(Glist$a1), a2=unlist(Glist$a2),
#                    af=unlist(Glist$af))
#   rsidsDUPS <- df$rsids[duplicated(df$rsids)]
#   df <- df[!df$rsids%in%rsidsDUPS,]
#   rsidsLD <- unlist(Glist$rsidsLD)
#   df <- df[df$rsids%in%rsidsLD,]
#   rownames(df) <- df$rsids
#   
#   inGlist <- stat$rsids%in%df$rsids
#   message(paste("Number of markers in stat also found in bedfiles:", sum(inGlist)))
#   
#   stat <- stat[inGlist,]
#   df <- df[rownames(stat),]
#   
#   aligned <- stat$effect_allele==df$a1
#   message(paste("Number of effect alleles aligned with first allele in bimfiles:", sum(aligned)))
#   message(paste("Number of effect alleles not aligned with first allele in bimfiles:", sum(!aligned)))
#   
#   if(!is.null(filename)) png(file=filename)
#   
#   layout(matrix(1:6,ncol=2,byrow=TRUE))
#   
#   try(lm(stat$effect_allele_freq[aligned]~ df$af[aligned]))
#   plot(stat$effect_allele_freq[aligned],df$af[aligned], ylab="AF in Glist (allele matching)",xlab="AF in stat (allele matching)")
#   
#   try(lm(stat$effect_allele_freq[!aligned]~ df$af[!aligned]))
#   plot(stat$effect_allele_freq[!aligned],df$af[!aligned], ylab="AF in Glist (allele not matching)",xlab="AF in stat (allele not matching)")
#   
#   stat[!aligned,"effect_allele_freq"] <- 1 - stat[!aligned,"effect_allele_freq"]
#   effect <- stat[!aligned,"b"]
#   effect_allele <- stat[!aligned,"effect_allele"]
#   alternative_allele <- stat[!aligned,"alternative_allele"]
#   stat[!aligned,"effect_allele"] <- alternative_allele 
#   stat[!aligned,"alternative_allele"] <- effect_allele 
#   stat[!aligned,"b"] <- -effect 
#   
#   try(lm(stat$effect_allele_freq~ df$af))
#   plot(stat$effect_allele_freq,df$af, ylab="AF in Glist",xlab="AF in stat (after allele flipped)")
#   
#   isDUPS <- duplicated(stat$rsids)
#   a1 <- df$a1
#   a2 <- df$a2
#   isAT <- a1=="A" & a2=="T"
#   isTA <- a1=="T" & a2=="A"
#   isCG <- a1=="C" & a2=="G"
#   isGC <- a1=="G" & a2=="C"
#   isCGAT <- isAT | isTA | isCG | isGC
#   CGTA <- c("C","G","T","A")
#   isINDEL <- !((a1%in%CGTA) & (a2%in%CGTA))
#   
#   largeMAFDIFF <- abs(df$af-stat$effect_allele_freq) > excludeMAFDIFF
#   maf <- stat$effect_allele_freq
#   maf[maf>0.5] <- 1-maf[maf>0.5]
#   lowMAF <- maf < excludeMAF
#   
#   message(paste("Number of markers excluded by low MAF:", sum(lowMAF)))
#   message(paste("Number of markers excluded by large difference between MAF difference:", sum(largeMAFDIFF)))
#   
#   rsidsQC <- lowMAF | largeMAFDIFF
#   
#   if(!is.null(stat$info)) {
#     lowINFO <- stat$info < excludeINFO
#     rsidsQC <- rsidsQC | lowINFO
#     message(paste("Number of markers excluded by low INFO score:", sum(lowINFO)))
#   }
#   
#   if(excludeCGAT) {
#     rsidsQC <- rsidsQC | isCGAT
#     message(paste("Number of markers excluded by ambiguity (CG or AT):", sum(isCGAT)))
#   }
#   if(excludeDUPS) {
#     rsidsQC <- rsidsQC | isDUPS
#     message(paste("Number of markers excluded by duplicated rsids", sum(isDUPS)))
#   }
#   if(excludeINDEL) {
#     rsidsQC <- rsidsQC | isINDEL
#     message(paste("Number of markers excluded by being INDEL:", sum(isINDEL)))
#   }
#   
#   rsidsQC <- !rsidsQC
#   plot(stat$effect_allele_freq[rsidsQC],df$af[rsidsQC], ylab="AF in Glist",xlab="AF in stat (after qc check)")
#   stat <- stat[rsidsQC,]
#   maf <- stat$effect_allele_freq
#   maf[maf>0.5] <- 1-maf[maf>0.5]
#   seb <- stat$seb
#   plot(y=seb,x=maf, ylab="SEB",xlab="MAF")
#   if(!is.null(filename)) dev.off()
#   stat$af <- stat$effect_allele_freq  
#   stat$alleles <- stat$effect_allele  
#   if(is.null(stat$n)) stat$n <- neff(seb=stat$seb,af=stat$af)
#   return(stat)
# }



# adjStat <- function(Glist=NULL,stat=NULL,filename=NULL, chr=NULL){
#   chromosomes <- chr
#   if(is.null(chromosomes)) chromosomes <- 1:Glist$nchr
#   badj <- NULL
#   for ( chr in chromosomes) {
#     
#     LD <- getSparseLD(Glist = Glist, chr = chr)
#     rsidsLD <- Glist$rsidsLD[[chr]]
#     
#     zobs <- zpred <- rep(0,length(rsidsLD))
#     names(zobs) <- names(zpred) <- rsidsLD
#     rsidsSTAT <- rownames(stat)[rownames(stat)%in%rsidsLD]
#     zobs[rsidsSTAT] <- stat[rsidsSTAT,"b"]
#     #zobs[rsidsSTAT] <- stat[rsidsSTAT,"b"]/stat[rsidsSTAT,"seb"]
#     for (i in 1:length(LD$indices)){
#       #zsum <- sum(zobs[LD$indices[[i]]]*LD$values[[i]])-zobs[i]
#       #nsum <- sum(abs(LD$values[[i]])) - 1
#       zsum <- sum(zobs[LD$indices[[i]]]*LD$values[[i]])
#       nsum <- sum(abs(LD$values[[i]]))
#       #zpred[i] <- zsum/nsum
#       if(!zobs[i]==0.0) zpred[i] <- zsum/nsum
#     }
#     # quantile normalisation (https://academic.oup.com/bioinformatics/article/19/2/185/372664)
#     zobs_rank <- rank(zobs, ties.method = "min")
#     zpred_rank <- rank(zpred, ties.method = "min")
#     zobs_sort <- sort(zobs)
#     zpred_sort <- sort(zpred)
#     zmean <- (zobs_sort+zpred_sort)/2
#     zobs_adj <- zmean[zobs_rank]
#     zpred_adj <- zmean[zpred_rank]
#     
#     if(!is.null(filename[chr])) {
#       png(filename[chr])
#       layout(matrix(1:4,ncol=2,byrow=TRUE))
#       plot(y=zobs,x=zpred, main=paste("Chr",chr), ylab="Observed Z",xlab="Predicted Z")
#       plot(y=zobs_adj,x=zpred_adj, main=paste("Chr",chr), ylab="Observed Z (normalized)",xlab="Predicted Z (normalized)")
#       plot(y=zobs,x=zobs_adj, main=paste("Chr",chr), ylab="Observed Z", xlab="Observed Z (normalized)")
#       plot(y=zobs,x=zpred_adj, main=paste("Chr",chr), ylab="Observed Z", xlab="Predicted Z (normalized)")
#       dev.off()
#     }
#     print(paste("Finished chr:",chr))
#     badj <- c(badj,zpred_adj)
#   }
#   badj <- badj[names(badj)%in%rownames(stat)]
#   stat <- cbind(stat[names(badj),],badj=badj)
#   return(stat)
# }


#'
#' @export
#'

adjStat <- function(stat = NULL, Glist = NULL, chr=NULL, statistics = "b", 
                    r2 = 0.9, ldSets = NULL, threshold = 1, header=NULL,
                    method = "pruning") {
  p <- stat$p
  if(is.null(stat$p)) p <- pnorm(abs(stat$b/stat$seb),lower.tail=FALSE)
  if(is.data.frame(stat)) names(p) <- rownames(stat)
  
  if(is.null(ldSets)) {p <- adjLD(Glist=Glist, stat=p, r2=r2, threshold=threshold)}
  if(!is.null(ldSets)) {p <- adjLD(Glist=Glist, stat=p, r2=r2, threshold=threshold, ldSets=ldSets)}
  
  p[p>0] <- 1
  if(is.null(header)) header <- c("rsids","chr","pos","allele","a1","a2","af")
  
  
  if(is.data.frame(stat)) {
    if(statistics=="b") {
      b <- stat[rownames(p),"b"]
      badj <- p*b
      #badj <- p*stat[rownames(p),"b"]
      colnames(badj) <- paste0("b_",threshold)
      if(any(colnames(stat)%in%header)) statadj <- data.frame(stat[rownames(badj),colnames(stat)%in%header],b,badj)
      if(!any(colnames(stat)%in%header)) statadj <- as.matrix(data.frame(b,badj))
      return(statadj)
    }
    if(statistics=="z") {
      z <- stat[rownames(p),"b"]/stat[rownames(p),"seb"]
      zadj <- p*z
      #zadj <- p*stat[rownames(p),"b"]/stat[rownames(p),"seb"]
      colnames(zadj) <- paste0("z_",threshold)
      if(any(colnames(stat)%in%header)) statadj <- data.frame(stat[rownames(zadj),colnames(stat)%in%header],z,zadj)
      if(!any(colnames(stat)%in%header)) statadj <- as.matrix(data.frame(z,zadj))
      return(statadj)
    }
  }
  if(is.list(stat)) {
    cls <- rep(1:ncol(stat$b),times=length(threshold))
    if(statistics=="b") {
      b <- stat$b[rownames(p),]
      badj <- p*b[,cls]
      colnames(b) <- paste0("b_",colnames(b))
      colnames(badj) <- paste0("b_",colnames(p))
      if(!is.null(stat$marker)) statadj <- data.frame(stat$marker[rownames(badj),colnames(stat$marker)%in%header],b,badj)
      if(is.null(stat$marker)) statadj <- as.matrix(data.frame(b,badj))
      return(statadj)
    }
    if(statistics=="z") {
      z <- stat$z[rownames(p),]
      if(is.null(z)) z <- stat$b[rownames(p),]/stat$seb[rownames(p),]
      zadj <- p*z[,cls]
      colnames(z) <- paste0("z_",colnames(z))
      colnames(zadj) <- paste0("z_",colnames(p))
      if(!is.null(stat$marker)) statadj <- data.frame(stat$marker[rownames(zadj),colnames(stat$marker)%in%header],z,zadj)
      if(is.null(stat$marker)) statadj <- as.matrix(data.frame(z,zadj))
      return(statadj)
    }
  }
  
}

#'
#' @export
#'

getStat <- function(stat=NULL, cls=NULL, rws=NULL) {
  if(is.null(rws)) rws <- 1:nrow(stat[[1]])
  for(i in 1:7) { stat[[i]] <- stat[[i]][rws,cls] }
  if(length(cls)==1)  stat <- as.data.frame(stat[1:7], stringsAsFactors=FALSE)
  return(stat)
}



#' LD pruning of summary statistics
#'
#' @description
#' Perform LD pruning of summary statistics before they are used in gene set enrichment analyses.
#' @param stat vector or matrix of single marker statistics (e.g. coefficients, t-statistics, p-values)
#' @param statistics is the type of statistics used in stat (e.g. statistics="p-value")
#' @param ldSets list of marker sets - names corresponds to row names in stat
#' @param r2 threshold for r2 used in LD pruning
#' @param threshold p-value threshold used in LD pruning
#' @param Glist list providing information about genotypes stored on disk
#' @param method used including method="pruning" which is default or "clumping"


#' @export

adjLD <- function(stat = NULL, Glist = NULL, chr=NULL, statistics = "p-value", r2 = 0.9, ldSets = NULL, threshold = 1,
                  method = "pruning") {
  
  if (any(is.na(stat))) stop(paste("NAs found in stat"))
  if(is.vector(stat)) stat <- as.matrix(stat) 
  rsidsStat <- rownames(stat)
  if (is.null(rsidsStat)) stop(paste("please provide names/rownames for stat"))
  #if (statistics == "p-value") pstat <- stat[, "p"]
  rsidsMapped <- rsidsStat%in%unlist(Glist$rsidsLD)
  stat <- stat[rsidsMapped,]
  stat <- as.matrix(stat)
  rsidsStat <- rownames(stat)
  if(any(!rsidsMapped)) {
    message(paste("Number of rsids found in LD matrices:", sum(rsidsMapped)))
    message(paste("Number of rsids not found in LD matrices:", sum(!rsidsMapped)))
    if(nrow(stat)<1) stop("No rsids remaining for pruning")
  }
  #rownames(pstat) <- rsidsStat
  if(is.null(colnames(stat))) colnames(stat) <- paste0("stat",1:ncol(stat))
  res <- NULL
  if (method %in% c("pruning", "clumping")) {
    if(!is.null(Glist$ldSets[[as.character(r2)]])) {
      message("Using ldSets in Glist")
      ldSets <- Glist$ldSets[[as.character(r2)]]
    }
    
    if (!is.null(ldSets)) nchr <- length(ldSets)
    if (!is.null(Glist)) nchr <- Glist$nchr
    if (is.null(chr)) chromosomes <- 1:nchr
    if (!is.null(chr)) chromosomes <- chr
    
    if(is.null(ldSets)) {
      ldSets <- vector(length=Glist$nchr, mode="list")
      for (chr in chromosomes) {
        message(paste("Extract LD information for chromosome:", chr))
        ldSets[[chr]] <- getLDsets(Glist = Glist, r2 = r2, chr = chr)
        #ldSets[[chr]] <- mapSets(sets = ldSets[[chr]], rsids = rsidsStat)
      }
    }
    for (chr in chromosomes) {
      ldSets[[chr]] <- mapSets(sets = ldSets[[chr]], rsids = rsidsStat)
    }
    
    for (thold in threshold) {
      pstat <- stat
      for (i in 1:ncol(stat)) {
        m <- length(rsidsStat)
        indx1 <- rep(T, m)
        indx2 <- rep(F, m)
        message(paste("Pruning stat column:", colnames(pstat)[i],"using threshold:",thold))
        for (chr in chromosomes) {
          #if (!is.null(Glist)) {
          #message(paste("Pruning chromosome:", chr, "for stat column:", colnames(pstat)[i]))
          #if(i==1) setsChr[[chr]] <- getLDsets(Glist = Glist, r2 = r2, chr = chr)
          #if(i==1) setsChr[[chr]] <- mapSets(sets = setsChr[[chr]], rsids = rsidsStat)
          #}
          if (!is.null(ldSets)) setsChr <- ldSets[[chr]]
          setsChr <- setsChr[names(setsChr)%in%rsidsStat]
          rsidsChr <- names(setsChr)
          rwsChr <- match(rsidsChr, rsidsStat)
          p <- pstat[rwsChr, i]
          o <- order(p, decreasing = FALSE)
          for (j in o) {
            if (p[j] <= thold) {
              if (indx1[rwsChr[j]]) {
                rws <- setsChr[[j]]
                indx1[rws] <- F
                indx2[rwsChr[j]] <- T
              }
            }
          }
          #message(paste("Finished pruning chromosome:", chr, "for stat column:", colnames(pstat)[i]))
        }
        if (method == "clumping") {
          pstat[indx1, i] <- 0
          p <- pstat[, i]
          pstat[p > thold, i] <- 0
        }
        if (method == "pruning") pstat[!indx2, i] <- 0
      }
      colnames(pstat) <- paste0(colnames(stat),"_",thold)
      res <- cbind(res,pstat)
    }
    
    
  }
  #res <- res[!rowSums(res == 0) == ncol(res), ]
  if (length(chromosomes)==1) res <- res[rownames(res)%in%rsidsChr,]
  return(res)
}



#'
#' @export
#'

adjLDStat <- function(stat=NULL, Glist = NULL, chr = NULL, region=NULL, msize=NULL, threshold=1e-5, overlap=NULL, niter=5) {
  
  if(is.null(stat)) stop("Please provide stat")
  if(is.null(Glist)) stop("Please provide Glist")
  if(is.null(chr)) stop("Please provide chr")
  if(is.null(region)) stop("Please specify region, e.g. 1000000:1200000")
  if(is.null(msize)) stop("Please specify window size: e.g. msize=200")
  if(is.null(overlap)) stop("Please specify window size overlap: e.g. overlap=100")
  
  chrStat <- stat[stat$chr==chr,]
  chrStat <- chrStat[chrStat$rsids%in%Glist$rsidsLD[[chr]],]
  
  region_rsids <- getMarkers(Glist=Glist, chr=chr, region=region)
  region_rsids <- region_rsids[region_rsids%in%chrStat$rsids]
  problem_rsids <- rep(FALSE,length(region_rsids))
  names(problem_rsids) <- region_rsids
  
  zobs <- chrStat[chrStat$rsids,"b"]/chrStat[chrStat$rsids,"seb"]
  names(zobs) <- chrStat$rsids
  
  sets <- splitWithOverlap(region_rsids,msize,overlap)
  
  for (iter in 1:niter) {
    if(iter>1) {
      if(length(remove_rsids)>0) chrStat <- chrStat[!chrStat$rsids%in%remove_rsids,]
    }
    region_rsids <- region_rsids[region_rsids%in%chrStat$rsids]
    sets <- splitWithOverlap(region_rsids,msize,overlap)
    
    for( i in 1:length(sets) ) {
      z <- vz <- rep(0,length(sets[[i]]))
      W <- getG(Glist=Glist,chr=chr,rsids=sets[[i]], scale=TRUE)
      WW <- cor(W)
      zo <- zobs[sets[[i]]]
      for (j in 1:length(sets[[i]])) {
        z[j] <- sum(WW[j,-j]*solve(WW[-j,-j])%*%zo[-j])
        vz[j] <- 1- t(WW[j,-j])%*%solve(WW[-j,-j])%*%WW[j,-j]
      }
      #plot(y=z,x=zo,ylab="Predicted", xlab="Observed",  frame.plot=FALSE)
      #abline(0,1, lwd=2, col=2, lty=2)
      tstat <- ((zo-z)**2)/vz
      p <- 1-pchisq(tstat,1)
      problem_rsids[sets[[i]][p<threshold]] <- TRUE 
      print(c(i,sum(p<threshold)))
    }
    remove_rsids <- names(problem_rsids)[problem_rsids]
    keep <- which.min(chrStat[remove_rsids,"p"])
    #remove_rsids <- remove_rsids[-keep] 
  }
  rsidsStat <- rownames(stat)
  rsidsChr <- rownames(stat[stat$chr==chr,])
  remove_rsids <- rsidsChr[!rsidsChr%in%rownames(chrStat)]
  rws <- match(remove_rsids,rsidsStat)
  stat <- stat[-rws,]
  return(stat)
}

# suggest to revise function
# read files outside function
# merge 2 or more data frames and create a list
# mergeStat <- function(Glist=NULL, fnDir=NULL, tnames=NULL){
#   # fnDir: path to summary stats from qcstat.function
#   # tnames: vector of trait names
#   nt <- length(fnDir)
#   b <- seb <- p <- n <- vector(nt, mode="list")
#   for(i in 1:nt){
#     stat <- fread(fnDir[i],data.table=F, showProgress=F)
#     rownames(stat) <- stat$rsids
#     b[[i]] <- stat[,"b",drop=F]
#     seb[[i]] <- stat[,"seb",drop=F]
#     p[[i]] <- stat[,"p",drop=F]
#     n[[i]] <- stat[,"n",drop=F]
#     message("Read summary stat file no ", i, "...\n")
#   }
#   b <- do.call(cbind,b)
#   seb <- do.call(cbind, seb)
#   p <- do.call(cbind, p)
#   n <- do.call(cbind, n)
#   colnames(b) <- colnames(seb) <- colnames(p) <- colnames(n) <- tnames
#   rsids <- rownames(b)
#   marker <- data.frame(rsids = rsids, chr = unlist(Glist$chr)[rsids],
#                        pos = unlist(Glist$pos)[rsids], a1 = unlist(Glist$a1)[rsids],
#                        a2 = unlist(Glist$a2)[rsids], af = unlist(Glist$a2)[rsids],
#                        stringsAsFactors=F)
#   
#   ma <- list(b = b, seb = seb, p = p, n = n, marker=marker)
#   return(ma)
# }

#'
#' @export
#'

mapStat <- function(Glist=NULL, stat=NULL) {
  
  # we use cpra to link sumstats and Glist
  cpra <- unlist(Glist$cpra)
  rsids <- unlist(Glist$rsids)
  
  # stat is a data.frame
  if(!is.data.frame(stat)) stop("stat should be  a data frame")
  if(!is.null(stat$rsids)) rownames(stat) <- stat$rsids
  
  
  fm_map1 <- c("rsids","chr","pos","a1","a2")
  
  format <- "unknown"
  
  if(all(fm_map1%in%colnames(stat))) format <- "fm_map1"
  
  if(format=="unknown") {
    message("Column headings for stat object not found")
    message("Column headings for stat object should be:")
    print(fm_map1)
    stop("please revised your stat object according to these ")
  }
  
  if(format=="fm_map1") {
    
    cpra1 <- paste(stat[,"chr"],stat[,"pos"],stat[,"a1"],stat[,"a2"],sep="_")
    cpra2 <- paste(stat[,"chr"],stat[,"pos"],stat[,"a2"],stat[,"a1"],sep="_")
    
    mapped <- cpra1%in%cpra | cpra2%in%cpra
    message("Map markers based on cpra")
    message(paste("Number of markers in stat mapped to marker ids in Glist:",sum(mapped)))
    message(paste("Number of markers in stat not mapped to marker ids in Glist:",sum(!mapped)))
    
    stat <- stat[mapped,]
    cpra1 <- cpra1[mapped]
    cpra2 <- cpra2[mapped]
    rws1 <- match(cpra1,cpra)
    rws2 <- match(cpra2,cpra)
    
    stat$rsids[!is.na(rws1)] <- rsids[rws1[!is.na(rws1)]]
    stat$rsids[!is.na(rws2)] <- rsids[rws2[!is.na(rws2)]]
    
    isdup <- duplicated(stat$rsids)
    if(any(isdup)) message("Removing markers with duplicated ids")
    if(any(isdup)) message(paste("Number of markers duplicated in stat:",sum(isdup)))
    stat <- stat[!isdup,] 
    rownames(stat) <- stat$rsids    
    
  }
  
  marker <- data.frame(rsids=unlist(Glist$rsids),cpra=unlist(Glist$cpra),
                       chr=unlist(Glist$chr), pos=unlist(Glist$position), 
                       a1=unlist(Glist$a1), a2=unlist(Glist$a2),
                       af=unlist(Glist$af),stringsAsFactors = FALSE)
  
  rownames(marker) <- marker$rsids
  
  message("Filtering markers based on information in Glist:")
  message("")
  
  rsids <-  gfilter(Glist = Glist,
                    excludeMAF=excludeMAF, 
                    excludeMISS=excludeMISS, 
                    excludeCGAT=excludeCGAT, 
                    excludeINDEL=excludeINDEL, 
                    excludeDUPS=excludeDUPS, 
                    excludeHWE=excludeHWE, 
                    excludeMHC=excludeMHC)
  marker <- marker[marker$rsids%in%rsids,]
  message("")
  
  message("Filtering markers based on information in stat:")
  message("")
  
  if(!is.null(stat$rsids)) marker_in_stat <- marker$rsids%in%stat$rsids
  if(!is.null(stat$marker)) marker_in_stat <- marker$rsids%in%stat$marker
  message(paste("Number of markers in stat also found in bimfiles:", sum(marker_in_stat)))
  message("")
  if(sum(marker_in_stat)==0) stop("No marker ids found in bimfiles")
  
  # align marker and stat object
  marker <- marker[marker_in_stat,]
  stat <- stat[marker$rsids,]
  
  if(!is.null(stat$effect_allele)) aligned <- stat$effect_allele==marker$a1
  if(!is.null(stat$a1)) aligned <- stat$a1==marker$a1
  message(paste("Number of effect alleles aligned with first allele in bimfiles:", sum(aligned)))
  message(paste("Number of effect alleles not aligned with first allele in bimfiles:", sum(!aligned)))
  message("")
  
  if(format=="fm_map1") {
    
    effect_columns <- !colnames(stat)%in%fm_map1
    
    #original
    effect <- stat[,effect_columns, drop = FALSE]
    effect_allele <- stat[,"a1"]
    non_effect_allele <- stat[,"a2"]
    effect_allele_freq <- stat[,"af"]
    
    # aligned
    stat[!aligned,effect_columns] <- -effect[!aligned, ,drop = FALSE]
    stat[!aligned,"a1"] <- non_effect_allele[!aligned]
    stat[!aligned,"a2"] <- effect_allele[!aligned] 
    stat[!aligned,"af"] <- 1-effect_allele_freq[!aligned]
    
  }  
  
  return(stat)
}

