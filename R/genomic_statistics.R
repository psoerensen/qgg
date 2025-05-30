####################################################################################################################
#    Module 8: Genomic statistics
####################################################################################################################
#' 
#' Quality Control of Marker Summary Statistics
#'
#' Performs quality control on GWAS summary statistics, which includes:
#' - Mapping marker ids to LD reference panel data.
#' - Checking effect allele, frequency, and build version.
#' - Excluding based on various criteria like MAF, HWE, INDELS, and more.
#'
#' @description
#' Quality control is a fundamental step in GWAS summary statistics analysis.
#' The function is equipped to handle various tasks including mapping marker ids,
#' checking the effect allele and its frequency, determining build versions, and 
#' excluding data based on multiple criteria.
#'
#' @param Glist List containing information about genotype matrix stored on disk.
#' @param stat Data frame of marker summary statistics. It should either follow 
#' the "internal" or "external" format.
#' @param excludeMAF Numeric. Exclusion threshold for minor allele frequency. Default is 0.01.
#' @param excludeMAFDIFF Numeric. Threshold for excluding markers based on allele frequency difference. Default is 0.05.
#' @param excludeINFO Numeric. Exclusion threshold for info score. Default is 0.8.
#' @param excludeMISS Numeric. Exclusion threshold for sample missingness. Default is 0.05.
#' @param excludeHWE Numeric. Exclusion threshold for Hardy Weinberg Equilibrium test p-value. Default is 1e-12.
#' @param excludeCGAT Logical. Exclude ambiguous alleles (CG or AT). Default is TRUE.
#' @param excludeMHC Logical. Exclude markers located in MHC region. Default is FALSE.
#' @param excludeINDEL Logical. Exclude insertion/deletion markers. Default is TRUE.
#' @param excludeDUPS Logical. Exclude markers with duplicated ids. Default is TRUE.
#'
#' @details
#' The function works with both "internal" and "external" formats of summary statistics. 
#' When the summary statistics format is "external", the function maps marker ids based 
#' on chr-pos-ref-alt information. It also aligns the effect allele with the LD reference 
#' panel and flips effect sizes if necessary. When allele frequencies are not provided, 
#' it uses the frequencies from the genotype data.
#' 
#' Required headers for external summary statistics:
#'  marker, chr, pos, ea, nea, eaf, b, seb, stat, p, n    
#' 
#' Required headers for internal summary statistics:
#'  rsids, chr, pos, ea, nea, eaf, b, seb, stat, p, n
#'
#' @return 
#' A data frame with processed and quality-controlled summary statistics.
#'
#' @author Peter Soerensen
#' 
#' @keywords internal
#' 
#' @export
#' 
checkStat <- function(Glist=NULL, stat=NULL, excludeMAF=0.01, excludeMAFDIFF=0.05, 
                      excludeINFO=0.8, excludeCGAT=TRUE, excludeINDEL=TRUE, 
                      excludeDUPS=TRUE, excludeMHC=FALSE, excludeMISS=0.05, 
                      excludeHWE=1e-12) {
  
  # Some summary statistics provide ea and nea in lower case
  stat$ea <- toupper(stat$ea)
  stat$nea <- toupper(stat$nea)
  
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
  
  fm_internal <- c("rsids","chr","pos","ea","nea","eaf","b","seb","p","n","info")
  fm_external <- c("marker","chr","pos","ea","nea","eaf","b","seb","p","n","info")
  fm_gbayes <- c("rsids","chr","pos","ea","nea","eaf","bm","dm","vm")
  
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
  
  format <- "unknown"
  
  if(all(fm_internal[1:5]%in%colnames(stat))) format <- "internal"
  if(all(fm_external[1:5]%in%colnames(stat))) format <- "external"
  if(all(fm_gbayes%in%colnames(stat))) {
    stat <- stat[,1:7]
    colnames(stat) <- c("rsids","chr","pos","ea","nea","eaf","b")
    format <- "internal"
  }
  if(format=="unknown") {
    message("Column headings for stat object not found")
    message("Column headings for stat object should minimal be:")
    print(fm_external[1:5])
    message("or:")
    print(fm_external)
    stop("please revised your stat object accordingly")
  }
  if(format=="internal") stat <- stat[,colnames(stat)%in%fm_internal]
  if(format=="external") stat <- stat[,colnames(stat)%in%fm_external]
  
  # map external summary statistics  
  if(format=="external") {
    cpra1 <- paste(stat[,"chr"],stat[,"pos"],stat[,"ea"],stat[,"nea"],sep="_")
    cpra2 <- paste(stat[,"chr"],stat[,"pos"],stat[,"nea"],stat[,"ea"],sep="_")
    
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
  
  # marker information in Glist
  marker <- data.frame(rsids=unlist(Glist$rsids),cpra=unlist(Glist$cpra),
                       chr=unlist(Glist$chr), pos=unlist(Glist$pos), 
                       ea=unlist(Glist$a1), nea=unlist(Glist$a2),
                       eaf=unlist(Glist$af),stringsAsFactors = FALSE)
  
  message("Filtering markers based on information in Glist:")
  message("")
  
  if(is.null(Glist$rsidsLD)) {
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
  }
  
  if(!is.null(Glist$rsidsLD)) {
    rsids <- unlist(Glist$rsidsLD)
    message(paste("Number of markers in sparse LD matrices:", sum(marker$rsids%in%rsids)))
    message("")
    marker <- marker[marker$rsids%in%rsids,]
  }

  rownames(marker) <- marker$rsids
  
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
  aligned <- stat$ea==marker$ea
  message(paste("Number of effect alleles aligned with first allele in bimfiles:", sum(aligned)))
  message(paste("Number of effect alleles not aligned with first allele in bimfiles:", sum(!aligned)))
  message("")
  
  
  # align stat if format external
  #if(format=="external") {
  
  #original
  effect <- stat[,"b"]
  effect_allele <- stat[,"ea"]
  non_effect_allele <- stat[,"nea"]
  if(!is.null(stat$eaf)) effect_allele_freq <- stat[,"eaf"]

  # aligned
  stat[!aligned,"b"] <- -effect[!aligned]
  stat[!aligned,"ea"] <- non_effect_allele[!aligned]
  stat[!aligned,"nea"] <- effect_allele[!aligned] 
  if(!is.null(stat$eaf)) stat[!aligned,"eaf"] <- 1-effect_allele_freq[!aligned]
  
  if(is.null(stat$eaf)) {
    message("No effect allele frequency (eaf) provided - using eaf in Glist")
    stat$eaf <- marker$eaf
  }
  
  # exclude based on maf  
  excludeMAFDIFF <- abs(marker$eaf-stat$eaf) > excludeMAFDIFF
  message(paste("Number of markers excluded by large difference between MAF difference:", sum(excludeMAFDIFF)))
  message("")
  stat <- stat[!excludeMAFDIFF,]
  marker <- marker[!excludeMAFDIFF,]
  
  if(colnames(stat)[1]=="marker")  colnames(stat)[1] <- "rsids"
    
  #colnames(stat) <- fm_internal
  #if(is.null(stat$n)) stat$n <- neff(seb=stat$seb,af=stat$eaf)
  if(is.null(stat$n)) stat$n <- 4/((2*stat$eaf*(1-stat$eaf))*stat$seb^2)
  
  if(!is.null(excludeMAF)) stat <- stat[stat$eaf>excludeMAF,]
  if(!is.null(excludeMAF)) stat <- stat[stat$eaf<(1-excludeMAF),]
  
  if(!is.null(stat$info)) {
    lowINFO <- stat$info < excludeINFO
    message(paste("Number of markers excluded by low INFO score:", sum(lowINFO)))
    message("")
    stat <- stat[!lowINFO,]
  }
  
  clsstat <- c("rsids","chr","pos","ea","nea","eaf","b","seb","p","n","info")%in%colnames(stat)
  clsstat <- c("rsids","chr","pos","ea","nea","eaf","b","seb","p","n","info")[clsstat]
  #colnames(stat) <- clsstat
  #return(stat)
  return(stat[,clsstat])
}


#' Adjustment of marker summary statistics using clumping and thresholding
#'
#' Adjust marker summary statistics using linkage disequilibrium information from Glist.
#'
#' Required input format for summary statistics:
#' 
#' stat can be a data.frame(rsids, chr, pos, ea, nea, eaf, b, seb, stat, p, n) (single trait)
#' 
#' stat can be a list(marker=(rsids, chr, pos, ea, nea, eaf), b, seb, stat, p, n) (multiple trait)
#'
#' @param Glist List of information about genotype matrix stored on disk.
#' @param stat A data frame with marker summary statistics (see required format above).
#' @param chr Chromosome(s) being processed.
#' @param statistics Specify what type of statistics ("b" or "z") is being processed. Default is "b".
#' @param r2 Threshold used in clumping/pruning procedure. Default is 0.9.
#' @param threshold P-value threshold used in clumping procedure. Default is 1.
#' @param method Method used in adjustment for linkage disequilibrium. Default is "clumping".
#' @param ldSets List of marker sets - names correspond to row names in `stat`.
#' @param header Character vector with column names to be excluded in the LD adjustment.
#'
#' @details
#' For details about the summary statistics format, see the main function description.
#'
#' @author Peter Soerensen
#' 
#' @export
#'
adjStat <- function(stat = NULL, Glist = NULL, chr=NULL, statistics = "b", 
                    r2 = 0.9, ldSets = NULL, threshold = 1, header=NULL,
                    method = "pruning") {
  p <- stat$p
  if(is.null(stat$p)) p <- pnorm(abs(stat$b/stat$seb),lower.tail=FALSE)
  if(is.data.frame(stat)) names(p) <- rownames(stat)
  
  if(is.null(ldSets)) {p <- adjLD(Glist=Glist, stat=p, r2=r2, chr=chr, threshold=threshold)}
  if(!is.null(ldSets)) {p <- adjLD(Glist=Glist, stat=p, r2=r2, chr=chr, threshold=threshold, ldSets=ldSets)}
  
  p[p>0] <- 1
  if(is.null(header)) header <- c("rsids","chr","pos","ea","nea","eaf")
  
  if(is.data.frame(stat)) {
    if(statistics=="b") {
      b <- stat[rownames(p),"b"]
      badj <- p*b
      colnames(badj) <- paste0("b_",threshold)
      if(any(colnames(stat)%in%header)) statadj <- data.frame(stat[rownames(badj),colnames(stat)%in%header],b,badj)
      if(!any(colnames(stat)%in%header)) statadj <- as.matrix(data.frame(b,badj))
      return(statadj)
    }
    if(statistics=="z") {
      z <- stat[rownames(p),"b"]/stat[rownames(p),"seb"]
      zadj <- p*z
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


#' LD pruning of summary statistics
#'
#' Perform LD pruning of summary statistics before they are used in gene set enrichment analyses.
#'
#' @param Glist List of information about the genotype matrix stored on disk.
#' @param stat A data frame with marker summary statistics. Ensure that it is in the correct format.
#' @param statistics Specify what type of statistics ("b", "z", or "p-value") is being processed. Default is "p-value".
#' @param chr Chromosome(s) being processed.
#' @param ldSets List of marker sets - names correspond to row names in `stat`.
#' @param r2 Threshold used in the clumping/pruning procedure. Default is 0.9.
#' @param threshold P-value threshold used in LD pruning.
#' @param method Method used in adjustment for linkage disequilibrium. Options are "pruning" or "clumping". Default is "pruning".
#' 
#' @keywords internal
#' 
#' @export
#' 
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
        message(paste("Pruning stat column:", i,"using threshold:",thold))
        for (chr in chromosomes) {
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
  if (length(chromosomes)==1) res <- res[rownames(res)%in%rsidsChr,]
  return(res)
}

#' Check concordance between marker effect and sparse LD matrix.
#' 
#' @description
#' Check concordance between predicted and observed marker effect. Marker effect is predicted based on sparse LD matrix in Glist.
#' 
#' @param Glist list of information about genotype matrix stored on disk
#' @param stat data frame with marker summary statistics (see required format above)
#' @param chr chromosome for which marker effect is checked
#' @param region genome region (in base pairs) for which marker effect is checked
#' @param threshold p-value threshold used for chisquare test for difference between observed and predicted marker effects
#' @param msize is the number of markers used in the prediction
#' @param overlap is the number of markers overlapping between adjacent genome region
#' @param niter is the number of iteration used for detecting outliers
#' 
#' @keywords internal
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
      tstat <- ((zo-z)**2)/vz
      p <- 1-pchisq(tstat,1)
      problem_rsids[sets[[i]][p<threshold]] <- TRUE 
      print(c(i,sum(p<threshold)))
    }
    remove_rsids <- names(problem_rsids)[problem_rsids]
    keep <- which.min(chrStat[remove_rsids,"p"])
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


#' Map marker summary statistics to Glist
#'
#' @description
#' Quality control is a critical step for working with summary statistics (in particular
#'                                                                         for external). 
#' Processing and quality control of GWAS summary statistics includes:                                                                      
#'
#' - map marker ids (rsids/cpra (chr, pos, ref, alt)) to LD reference panel data 
#' 
#' - check effect allele (flip EA, EAF, Effect)
#' 
#' - check effect allele frequency
#' 
#' - thresholds for MAF and HWE
#' 
#' - exclude INDELS, CG/AT and MHC region
#' 
#' - remove duplicated marker ids
#' 
#' - check which build version
#' 
#' - check for concordance between marker effect and LD data
#'
#' @param Glist list of information about genotype matrix stored on disk
#' @param stat dataframe with marker summary statistics
#' @param excludeMAF exclude marker if minor allele frequency (MAF) is below threshold (0.01 is default)
#' @param excludeMAFDIFF exclude marker if minor allele frequency difference (MAFDIFF) between Glist$af and stat$af is above threshold (0.05 is default)
#' @param excludeINFO exclude marker if info score (INFO) is below threshold (0.8 is default)
#' @param excludeMISS exclude marker if missingness (MISS) is above threshold (0.05 is default)
#' @param excludeHWE exclude marker if p-value for Hardy Weinberg Equilibrium test is below threshold (0.01 is default)
#' @param excludeCGAT exclude marker if alleles are ambigous (CG or AT)
#' @param excludeMHC exclude marker if located in MHC region 
#' @param excludeINDEL exclude marker if it an insertion/deletion  
#' @param excludeDUPS exclude marker id if duplicated
#' 
#' @keywords internal
#' 
#' @author Peter Soerensen
#'
#' @export
#'
mapStat <- function(Glist=NULL, stat=NULL, excludeMAF=0.01, excludeMAFDIFF=0.05, 
                    excludeINFO=0.8, excludeCGAT=TRUE, excludeINDEL=TRUE, 
                    excludeDUPS=TRUE, excludeMHC=FALSE,excludeMISS=0.05, 
                    excludeHWE=1e-12) {
  
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
                       chr=unlist(Glist$chr), pos=unlist(Glist$pos), 
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


mapStatID <- function(Glist=NULL, stat=NULL, cpra=NULL) {
  # we use cpra to link sumstats and Glist
  if(!is.null(cpra)) rsids <- cpra 
  if(is.null(cpra)) {
    cpra <- unlist(Glist$cpra)
    rsids <- unlist(Glist$rsids)
  }
  # stat is a data.frame
  if(!is.data.frame(stat)) stop("stat should be  a data frame")
  fm_internal <- c("marker","chr","pos","ea","nea")
  fm_internal1 <- c("rsids","chr","pos","ea","nea")
  fm_internal2 <- c("marker","chr","pos","ea","nea")
  format <- "unknown"
  if(all(fm_internal1%in%colnames(stat))) format <- "internal"
  if(all(fm_internal2%in%colnames(stat))) format <- "internal"
  if(format=="unknown") {
    message("Column headings for stat object not found")
    message("First five elements in column heading for stat object should be:")
    print(fm_internal)
    stop("please revised your stat object accordingly")
  }
  cpra1 <- paste(stat[,"chr"],stat[,"pos"],toupper(stat[,"ea"]),toupper(stat[,"nea"]),sep="_")
  cpra2 <- paste(stat[,"chr"],stat[,"pos"],toupper(stat[,"nea"]),toupper(stat[,"ea"]),sep="_")
  mapped <- cpra1%in%cpra | cpra2%in%cpra
  message("Map markers based on cpra")
  message(paste("Number of markers in stat mapped to marker ids in Glist:",sum(mapped)))
  message(paste("Number of markers in stat not mapped to marker ids in Glist:",sum(!mapped)))
  stat <- stat[mapped,]
  cpra1 <- cpra1[mapped]
  cpra2 <- cpra2[mapped]
  rws1 <- match(cpra1,cpra)
  rws2 <- match(cpra2,cpra)
  
  if(colnames(stat)[1]=="marker"){
    stat$marker[!is.na(rws1)] <- rsids[rws1[!is.na(rws1)]]
    stat$marker[!is.na(rws2)] <- rsids[rws2[!is.na(rws2)]]  
  }
  if(colnames(stat)[1]=="rsids"){
    stat$rsids[!is.na(rws1)] <- rsids[rws1[!is.na(rws1)]]
    stat$rsids[!is.na(rws2)] <- rsids[rws2[!is.na(rws2)]]  
  }
  
  colnames(stat)[1] <- "rsids"
  
  isdup <- duplicated(stat$rsids)
  if(any(isdup)) message("Removing markers with duplicated ids")
  if(any(isdup)) message(paste("Number of markers duplicated in stat:",sum(isdup)))
  
  stat <- stat[!isdup,]
  
  rownames(stat) <- stat$rsids
  
  return(stat)
}
