####################################################################################################################
#    Module 8: Genomic simulation
####################################################################################################################
#'
#' Genomic simulation
#' 
#'
#' Simulate Genotype and Phenotype Data
#'
#' This function simulates genotype and phenotype data based on the `Glist`, which is 
#' information about the genotype matrix. 
#'
#' @param Glist A list of information about the genotype matrix. Default is `NULL`.
#' @param chr The chromosome(s) being used in the simulation. Default is 1.
#' @param nt Number of traits. Default is 1.
#' @param W Matrix of centered and scaled genotypes. Default is `NULL`.
#' @param n Number of individuals. Default is 1000.
#' @param m Number of markers. Default is 1000.
#' @param rsids A character vector of rsids. Default is `NULL`.
#'
#' @return 
#' A list containing:
#' \itemize{
#'   \item \code{y}: Phenotypes.
#'   \item \code{W}: Matrix of centered and scaled genotypes.
#'   \item \code{e}: Errors.
#'   \item \code{g}: Genotype effect.
#'   \item \code{b0}, \code{b1}: Coefficients.
#'   \item \code{set0}, \code{set1}: Selected markers.
#'   \item \code{causal}: Causal markers.
#' }
#' 
#' @examples
#' ## Plink bed/bim/fam files
#' bedfiles <- system.file("extdata", paste0("sample_chr",1:2,".bed"), package = "qgg")
#' bimfiles <- system.file("extdata", paste0("sample_chr",1:2,".bim"), package = "qgg")
#' famfiles <- system.file("extdata", paste0("sample_chr",1:2,".fam"), package = "qgg")
#'
#' # Summarize bed/bim/fam files
#' Glist <- gprep(study="Example", bedfiles=bedfiles, bimfiles=bimfiles, famfiles=famfiles)
#'
#' # Simulate phenotype
#' sim <- gsim(Glist=Glist, chr=1, nt=1)
#' head(sim$y)
#' head(sim$e)
#' head(sim$causal)
#' 
#' @author Peter Soerensen
#' 
#' @export
#' 
gsim <- function(Glist=NULL, chr=1, nt=1,W=NULL, n=1000, m=1000, rsids=NULL) {
  if(!is.null(Glist)) {
    if(is.null(rsids)) rsids <- Glist$rsidsLD[[chr]]
    W <- getG(Glist=Glist, chr=chr, rsids=rsids)
  }
  if(is.null(W)) {
    W <- matrix(runif(n),ncol=1)
    for (i in 2:m) {
      W <- cbind(W,scale(W[,i-1]) + runif(n))  
    }
  }
  n <- nrow(W)
  m <- ncol(W)
  if(is.null(colnames(W))) colnames(W) <- paste0("m",1:m)
  if(is.null(rownames(W))) rownames(W) <- paste0("id",1:n)
  
  y <- e <- vector(length=nt,mode="list")
  names(y) <- paste0("D",1:nt)
  set0 <- sample(1:ncol(W),2)
  set1 <- b1 <- g1 <- vector(length=nt,mode="list")
  g <- NULL
  b0 <- sample(c(0.25,-0.25,0.5,-0.5),2)
  for (i in 1:nt){
    g0 <- W[,set0]%*%b0
    set1[[i]] <- sample(1:ncol(W),2)
    b1[[i]] <- sample(c(0.25,-0.25,0.5,-0.5),length(set1[[i]]))
    g1[[i]] <- W[,set1[[i]]]%*%b1[[i]]
    e[[i]] <- rnorm(nrow(W),mean=0,sd=1)
    y[[i]] <- as.vector(g0+g1[[i]]+e[[i]])
    names(y[[i]]) <- rownames(W)
    g <- cbind(g,g0+g1[[i]])
  }
  colnames(g) <- paste0("D",1:nt) 
  causal <- c(set0,unlist(set1))
  causal <- colnames(W)[causal]
  if(nt==1) return( list(y=y[[1]],W=W,e=e[[1]],g=g,b0=b0,b1=b1,set0=set0,set1=set1,causal=causal))
  if(nt>1) return( list(y=as.matrix(as.data.frame(y)),W=W,e=as.matrix(as.data.frame(e)),g=g,b0=b0,b1=b1,set0=set0,set1=set1,causal=causal))
}


#' Simulate Genetic Data Based on Given Parameters
#'
#' This function simulates phenotype data by random sampling of markers available on `Glist`.
#' Default parameters for the simulated phenotype reflect the genetic architecture assumed by BayesC prior(Habier et al., 2011).
#' This function is under active development.
#'
#' @param Glist A list containing genetic data. If NULL, the function will stop with an error.
#' @param h2 Heritability. If NULL, heritability of 0.5 is assumed.
#' @param pv Prevalence. A value between >0 and <1. If NULL then trait is continuous.
#' @param m Number of causal markers. The values for either `m` or `prp.cau` should be provided at any given time. 
#' If the list of quality controlled markers is not available then list of raw markers is used.  
#' If `m` is NULL and `pi` is also NULL, `pi` will default to 0.001.
#' @param pi Proportion of causal markers. The values for either `m` or `pi` should be provided at any given time. 
#' @param n Number of individuals randomly sampled from `Glist`. If NULL, all the individuals on `Glist` is used.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{y}: Vector of simulated phenotypes.
#'   \item \code{g}: Vector of simulated genetic values.
#'   \item \code{e}: Vector of simulated residual effects.
#'   \item \code{b}: Vector of effect sizes of the simulated causal markers.
#'   \item \code{causal}: Vector of ids for the simulated causal markers.
#'   \item \code{h2}: Estimated heritability of the simulated phenotype.
#' }
#'
#' @author Peter Soerensen 
#' @author Merina Shrestha
#' 
#' @keywords internal
#' 
#' @export
#' 
gsimC <- function(Glist = NULL, h2 = NULL, pi = NULL, pv = NULL, m = NULL, n = NULL) {
  if (is.null(Glist)) stop("Error: Glist is NULL")
  
  # Set SNP heritability
  hsnp <- if (is.null(h2)) 0.5 else h2
  message("Heritability: ", hsnp)
  
  # Determine number of causal SNPs
  if (is.null(m) && is.null(pi)) {
    pi <- 0.001
  }
  
  if (!is.null(m) && is.null(pi)) {
    snp_cau <- m
  } else {
    tot_snps <- if (!is.null(Glist$rsidsLD)) length(unlist(Glist$rsidsLD)) else length(unlist(Glist$rsids))
    snp_cau <- round(tot_snps * pi)
  }
  
  message("Number of causal SNPs: ", snp_cau)
  
  # Simulate effect sizes
  pop_var <- hsnp / snp_cau
  effect_sd <- sqrt(pop_var)
  b <- rnorm(snp_cau, mean = 0, sd = effect_sd)
  
  # Select causal SNPs
  rsids_all <- if (!is.null(Glist$rsidsLD)) {
    message("QC'd SNPs used")
    unlist(Glist$rsidsLD)
  } else {
    message("Raw SNPs used")
    unlist(Glist$rsids)
  }
  rsids_causal <- sample(rsids_all, snp_cau)
  
  # Determine sample individuals
  ids_used <- if (is.null(n)) Glist$ids else sample(Glist$ids, n)
  message("Number of samples: ", length(ids_used))
  
  # Compute genetic values
  g_total <- 0
  pb <- txtProgressBar(min = 0, max = snp_cau, style = 3)
  
  for (i in seq_along(rsids_causal)) {
    setTxtProgressBar(pb, i)
    snp_id <- rsids_causal[i]
    
    for (chr in seq_along(Glist$bedfiles)) {
      idx <- match(snp_id, Glist$rsids[[chr]])
      if (!is.na(idx)) {
        geno <- getG(Glist, chr = chr, cls = idx, ids = ids_used,
                     impute = TRUE, scale = TRUE)
        g_total <- g_total + geno * b[i]
        break
      }
    }
  }
  close(pb)
  
  g <- g_total
  
  # Simulate environmental noise
  e_sd <- sqrt(var(g) * ((1 / hsnp) - 1))
  e <- rnorm(length(g), mean = 0, sd = e_sd)
  
  # Create phenotype
  y <- as.vector(g + e)
  names(y) <- ids_used
  
  # Optional binary transformation
  if (!is.null(pv)) {
    liability <- g + e
    threshold <- quantile(liability, probs = 1 - pv)
    y <- ifelse(liability >= threshold, 1, 0)
  }
  
  # Output
  return(list(
    y = y,
    g = g,
    e = e,
    b = b,
    causal = rsids_causal,
    h2 = round(var(g) / (var(g) + var(e)), 2)
  ))
}


#' Simulate Genetic Data Based on Given Parameters
#'
#' This function simulates phenotype data by random sampling of markers available on `Glist`.
#' Default parameters for the simulated phenotype reflect the genetic architecture assumed by BayesR prior(Erbe et al., 2012).
#' Marker effect are sampled from mixture distributions leading to small, moderate and large effect sizes. 
#' This function is under active development.
#' 
#' @param Glist A list containing genetic data. If NULL, the function will stop with an error.
#' @param h2 Heritability. If NULL, heritability of 0.5 is assumed.
#' @param pv Prevalence. A value between >0 and <1. If NULL then trait is continuous.
#' @param m Number of causal markers. The values for either `m` or `prp.cau` should be provided at any given time. 
#' If the list of quality controlled markers is not available then list of raw markers is used.  
#' If `m` is NULL and `pi` is also NULL, `pi` will default to 0.001.
#' @param pi Proportion of causal markers. The values for either `m` or `pi` should be provided at any given time. 
#' @param n Number of individuals randomly sampled from `Glist`. If NULL, all the individuals on `Glist` is used.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{y}: Vector of simulated phenotypes.
#'   \item \code{g}: Vector of simulated genetic values.
#'   \item \code{e}: Vector of simulated residual effects.
#'   \item \code{b}: list of vectors for effect sizes of the simulated causal markers in three classes.
#'   \item \code{causal}: list of vectors for ids for the simulated causal markers in three classes.
#'   \item \code{h2}: Estimated heritability of the simulated phenotype.
#' }
#' 
#' @author Peter Soerensen 
#' @author Merina Shrestha
#' 
#' @keywords internal
#' @export
#' 
gsimR <- function(Glist = NULL, h2 = NULL, pi = NULL, pv = NULL, m = NULL, n = NULL) {
  
  if (is.null(Glist)) stop("Error: Glist is NULL")
  
  # Set SNP heritability
  hsnp <- if (is.null(h2)) 0.5 else h2
  message("Heritability: ", hsnp)
  
  # Determine number of causal SNPs
  if (is.null(m) && is.null(pi)) {
    pi <- 0.001
  }
  if (!is.null(m) && is.null(pi)) {
    snp_cau <- m
  } else {
    tot_snps <- if (!is.null(Glist$rsidsLD)) length(unlist(Glist$rsidsLD)) else length(unlist(Glist$rsids))
    snp_cau <- round(tot_snps * pi)
  }
  
  message("Number of causal SNPs: ", snp_cau)
  
  # Divide causal SNPs into 3 classes
  snp_cau_1 <- max(round(snp_cau * 0.94), 1)
  snp_cau_2 <- max(round(snp_cau * 0.05), 1)
  snp_cau_3 <- max(round(snp_cau * 0.01), 1)
  
  for (i in 1:3) {
    message("Number of causal SNPs in class ", i, ": ", get(paste0("snp_cau_", i)))
  }
  
  # Define class-specific heritabilities
  h_class1 <- (snp_cau_1 * 0.001) / hsnp
  h_class2 <- (snp_cau_2 * 0.01) / hsnp
  h_class3 <- (snp_cau_3 * 0.1) / hsnp
  
  # Flatten list of SNPs
  rsids_all <- if (!is.null(Glist$rsidsLD)) unlist(Glist$rsidsLD) else unlist(Glist$rsids)
  
  # Sample causal SNPs for each class
  rsids_cau_1 <- sample(rsids_all, snp_cau_1)
  rsids_all <- setdiff(rsids_all, rsids_cau_1)
  
  rsids_cau_2 <- sample(rsids_all, snp_cau_2)
  rsids_all <- setdiff(rsids_all, rsids_cau_2)
  
  rsids_cau_3 <- sample(rsids_all, snp_cau_3)
  
  # Sample individual IDs
  ids_used <- if (is.null(n)) Glist$ids else sample(Glist$ids, n)
  message("Number of samples: ", length(ids_used))
  
  # Function to generate genetic values
  g_eff <- function(h_class, n_causal, rsids_causal, Glist) {
    pop_var <- h_class / n_causal
    effect_sizes <- rnorm(n_causal, mean = 0, sd = sqrt(pop_var))
    g_total <- 0
    
    pb <- txtProgressBar(min = 0, max = n_causal, style = 3)
    for (i in seq_along(rsids_causal)) {
      setTxtProgressBar(pb, i)
      snp_id <- rsids_causal[i]
      
      for (chr in seq_along(Glist$bedfiles)) {
        match_idx <- match(snp_id, Glist$rsids[[chr]])
        if (!is.na(match_idx)) {
          geno <- getG(Glist, chr = chr, cls = match_idx, ids = ids_used,
                       impute = TRUE, scale = TRUE)
          g_total <- g_total + geno * effect_sizes[i]
          break
        }
      }
    }
    close(pb)
    list(g = g_total, b = effect_sizes)
  }
  
  # Compute genetic values per class
  message("Generating effects for class 1...")
  res1 <- g_eff(h_class1, snp_cau_1, rsids_cau_1, Glist)
  
  message("Generating effects for class 2...")
  res2 <- g_eff(h_class2, snp_cau_2, rsids_cau_2, Glist)
  
  message("Generating effects for class 3...")
  res3 <- g_eff(h_class3, snp_cau_3, rsids_cau_3, Glist)
  
  # Combine genetic values
  g <- res1$g + res2$g + res3$g
  
  # Simulate environmental noise
  e_sd <- sqrt(var(g) * ((1 / hsnp) - 1))
  e <- rnorm(length(g), mean = 0, sd = e_sd)
  y <- g + e
  names(y) <- ids_used
  
  # Optional binary transformation based on quantile
  if (!is.null(pv)) {
    liability <- g + e
    threshold <- quantile(liability, probs = 1 - pv)
    y <- ifelse(liability >= threshold, 1, 0)
  }
  
  # Output list
  return(list(
    y = y,
    e = e,
    g = g,
    b = list(class1 = res1$b, class2 = res2$b, class3 = res3$b),
    causal = list(class1 = rsids_cau_1, class2 = rsids_cau_2, class3 = rsids_cau_3),
    h2 = round(var(g) / (var(g) + var(e)), 2)
  ))
}



checkSimulatedData <- function(sim_result) {
  stopifnot(all(c("y", "g", "e") %in% names(sim_result)))
  
  y <- sim_result$y
  g <- sim_result$g
  e <- sim_result$e
  
  # Compute variance components
  var_g <- var(g)
  var_e <- var(e)
  var_y <- var(y)
  h2_est <- round(var_g / var_y, 3)
  
  # Set up 2x2 plot layout
  old_par <- par(mfrow = c(2, 2))
  
  # Plot 1: Histogram of phenotype
  hist(y, breaks = 50, col = "skyblue", main = "Histogram of Phenotype (y)", xlab = "Phenotype")
  
  # Plot 2: Histogram of genetic value
  hist(g, breaks = 50, col = "lightgreen", main = "Histogram of Genetic Value (g)", xlab = "Genetic Value")
  
  # Plot 3: Scatterplot g vs y
  plot(g, y, pch = 20, col = "darkblue", xlab = "Genetic Value (g)", ylab = "Phenotype (y)",
       main = "Genetic Value vs Phenotype")
  abline(lm(y ~ g), col = "red")
  
  # Plot 4: QQ plot of residuals
  residuals <- y - g
  qqnorm(residuals, main = "QQ Plot of Residuals (e = y - g)")
  qqline(residuals, col = "red")
  
  par(old_par)
  
  # Optional barplot of variance components
  barplot(c(var_g, var_e, var_y),
          names.arg = c("Genetic", "Environmental", "Total"),
          col = c("green", "red", "blue"),
          main = paste("Variance Components (h2 ≈", h2_est, ")"))
  
  # Binary outcome summary if applicable
  if (all(y %in% c(0, 1))) {
    cat("\nBinary outcome detected:\n")
    print(table(y))
    barplot(table(y), col = c("gray80", "gray40"), main = "Binary Phenotype Distribution",
            names.arg = c("0", "1"))
  }
}

