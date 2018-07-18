#############################################################################################
# utility function for using DMU 
#############################################################################################
#' 
#' Genomic Feature Model analyses implemented using Restriced Likelihood Methods in DMU
#'
#' @description
#' Genomic Feature Best Linear Unbiased Prediction models implemented using REML. 
#'
#' @details 
#' The models are implemented using restricted maximum likelihood methods. 
#' Variance components estimated using REML and predictions are based on MME. 
#' Predicted random effects and single marker effects and statistics can be obtained.
#' Cross validation procedures for assessment of prediction accuracy and model selection. 
#' This is an interface to be used for DMU. 
#' 
#' @param fm a formula with model statement for the linear mixed model 
#' @param data a data frame containing the phenotypic observations and fixed factors specified in the model statements
#' @param GRMlist a list of relationship / correlation matrices corresponding to random effects specified in vfm
#' @param validate a matrix or a list with the ids of validation sets corresponding to the rows in data
#' @param bin is the directory for DMU binaries (dmu1 and dmuai1)
#' @return Returns results in a list structure including 
#' \item{f}{list of predicted random effects} 
#' \item{sigma}{estimated variance components} 
#' \item{asd}{asymptotic standard deviation for the estimated variance components} 
#' \item{fitted}{fitted values from linear mixed model fit} 
#' \item{residuals}{residuals from linear mixed model fit} 
#' @author Peter Sørensen
#' @references Edwards, S. M., Sørensen, I. F., Sarup, P., Mackay, T. F., & Sørensen, P. (2016). Genomic prediction for quantitative traits is improved by mapping variants to gene ontology categories in Drosophila melanogaster. Genetics, 203(4), 1871-1883.
#’ @references Rohde, P. D., Demontis, D., Cuyabano, B. C. D., Børglum, A. D., & Sørensen, P. (2016). Covariance Association Test (CVAT) Identifies Genetic Markers Associated with Schizophrenia in Functionally Associated Biological Processes. Genetics, 203(4), 1901-1913.
#’ @references Edwards, S. M., Thomsen, B., Madsen, P., & Sørensen, P. (2015). Partitioning of genomic variance reveals biological pathways associated with udder health and milk production traits in dairy cattle. Genetics Selection Evolution, 47(1), 60.
#’ @references Sarup, P., Jensen, J., Ostersen, T., Henryon, M., & Sørensen, P. (2016). Increased prediction accuracy using a genomic feature model including prior information on quantitative trait locus regions in purebred Danish Duroc pigs. BMC genetics, 17(1), 11.
#' @examples
#' 
#' library(qgg)
#' 
#' setwd("C:/Users/sor/Dropbox/GFBLUP DGRP/scripts used in the analyses/work")
#' #bin <- "C:/Program Files (x86)/QGG-AU/DMUv6/R5.2-EM64T/bin"
#' 
#' # Simulate data
#' W <- matrix(rnorm(4000000), ncol = 10000)
#'   colnames(W) <- as.character(1:ncol(W))
#'   rownames(W) <- as.character(1:nrow(W))
#' 
#' G <- W %*% t(W) / ncol(W)
#'
#' y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + 10 * rnorm(nrow(W))
#' 
#' data <- data.frame(f = factor(sample(1:2, nrow(W), replace = TRUE)), g = factor(1:nrow(W)), y = y)
#' 
#' fm <- y ~ f + (1 | g~G) 
#' GRMlist <- list(G = G)
#' 
#' 
#' fit <- remlDMU(fm = fm, GRMlist = GRMlist, data = data)
#'   str(fit)
#' 
#' @export
#'

# Main function for reml analyses suing DMU
remlDMU <- function(fm = NULL, vfm = NULL, GRMlist = NULL, restrict = NULL, data = NULL, validate = NULL, bin = NULL) {
     
     tol <- 0.001
     fit <- cvfit <- NULL
     #model <- extractModel(fm = fm, data = data)
     model <- lapply(fm, function(x) {extractModel(fm = x, data = data)})
     model <- modelDMU(model = model, restrict = restrict)
     flevels <- writeDMU(model = model, data = data, GRMlist = GRMlist)
     executeDMU(bin = bin)
     fit <- readDMU(model = model, flevels = flevels)
     if (!is.null(validate)) {
          for (i in 1:ncol(validate)) {
               #set data missing
               writeDMU(model = model, data = data, GRMlist = GRMlist)
               executeDMU(bin = bin)
               cvfit[[i]] <- readDMU(model = model, flevels = flevels)
          }
          fit$cv <- cvfit
     }
     
     return(fit = fit)
     
}

# Extract model information from fm and vfm 
extractModel <- function(fm = NULL, data = NULL) {
     
     vtype <- sapply(data, class)
     
     ffvar <- frvar <- vvar <- yvar <- NULL 
     yvar <- as.character(fm)[2]
     
     fmsplit <- unlist(strsplit(as.character(fm)[3], split = "+", fixed = TRUE))
     rwsR <- grep(")", fmsplit, fixed = TRUE)
     rwsF <- (1:length(fmsplit))[-rwsR]
     
     fvar <- fmsplit[rwsF]
     fvar <- gsub(" ", "", fvar, fixed = TRUE)
     ffvar <- fvar[vtype[fvar] == "factor"]
     frvar <- fvar[vtype[fvar] == "numeric"]
     
     vvar <- fmsplit[rwsR]
     vvar <- gsub("1 | ", "", vvar, fixed = TRUE)
     vvar <- gsub("+", "", vvar, fixed = TRUE)
     vvar <- gsub("(", "", vvar, fixed = TRUE)
     vvar <- gsub(")", "", vvar, fixed = TRUE)
     vvar <- gsub(" ", "", vvar, fixed = TRUE)
     
     vvar <- lapply(vvar, function(x) {x <- unlist(strsplit(as.character(x), split = "~", fixed = TRUE))})
     
     cvvar <- sapply(vvar, function(x) {x[2]})
     vvar <- sapply(vvar, function(x) {x[1]})
     names(cvvar) <- vvar
     cvvar <- cvvar[!is.na(cvvar)]    
     vvartype <- rep("I", length(vvar))
     names(vvartype) <- vvar
     vvartype[names(cvvar)] <- "COR"
     
     nreg <- length(frvar)
     nrandom <- length(vvar)
     nfixed <- length(ffvar)
     nfactors <- nrandom + nfixed 
     
     variables <- list(fixed = ffvar, 
                       regression = frvar, 
                       random = vvar, 
                       response = yvar, 
                       factors = c(ffvar, vvar), 
                       variables = c(ffvar, vvar, yvar, frvar))
     
     n <- as.data.frame(t(sapply(variables, length)))
     
     return(list(fixed = ffvar, 
                 regression = frvar, 
                 random = vvar, 
                 response = yvar, 
                 factors = c(ffvar, vvar), 
                 variables = c(ffvar, vvar, yvar, frvar), 
                 n = n, 
                 covtype = vvartype, 
                 covmat = cvvar)
     ) 
     
}

modelDMU <- function(model = NULL, restrict = NULL) {
     
     fixed <- unique(unlist(lapply(model, function(x) {x$fixed})))
     random <- unique(unlist(lapply(model, function(x) {x$random})))
     regression <- unique(unlist(lapply(model, function(x) {x$regression})))
     response <- unique(unlist(lapply(model, function(x) {x$response})))
     covmat <- NULL
     for (i in 1:length(model)) {
          covmat <- c(covmat, model[[i]]$covmat) 
     }
     covmat <- covmat[!duplicated(names(covmat))]
     model$nt <- length(model)
     model$absorb <- rep(0, model$nt)
     
     model$data <- NULL
     model$data$missing <- -9999
     model$data$variables <- c(fixed, random, response, regression)
     model$data$nvariables <- length(model$data$variables)
     model$data$nintegers <- sum(length(fixed) + length(random))
     model$data$nreals <- sum(length(response) + length(regression))
     
     model$data$integers <- 1:model$data$nintegers
     names(model$data$integers) <- c(fixed, random)
     model$data$reals <- 1:model$data$nreals
     names(model$data$reals) <- c(response, regression)
     model$data$random <- 1:length(random)
     names(model$data$random) <- random
     model$data$covmat <- covmat
     model$restrict <- NULL
     model$restrict$nresiduals <- model$restrict$residuals <- 0
     if(!is.null(restrict)) {
          model$restrict$nresiduals <- nrow(restrict$residuals)
          model$restrict$residuals <- t(restrict$residuals)
     }
     
     return(model)
     
}

# Recode factors for DMU 
recodeDMU <- function(data = NULL) {
     
     flevels <- rlevels <- NULL
     for (i in 1:ncol(data)) {
          f <- data[, i]
          flevels[[i]] <- levels(f)
          names(flevels[[i]]) <- 1:nlevels(f)
          rlevels[[i]] <- 1:nlevels(f)
          names(rlevels[[i]]) <- levels(f)
          data[, i] <- names(flevels[[i]])[f]
     }
     names(flevels) <- names(rlevels) <- colnames(data)
     head(data)
     
     return(list(rlevels = rlevels, flevels = flevels, data = data))
     
}

# Write DIR file, data file, and cor files for DMU 
writeDMU <- function(model = NULL, data = NULL, GRMlist = NULL, tol = 0.001) {
     
     # Write DMU DIR file
     dir.file <- "gfm.DIR"
     write("$COMMENT", file = dir.file)
     write("DIR file DMU generated from R ", file = dir.file, append = TRUE)
     write(" ", file = dir.file, append = TRUE)
     write(paste("$ANALYSE", 1, 1, 0, 0), file = dir.file, append = TRUE)
     write(" ", file = dir.file, append = TRUE)
     write(c("$DATA ASCII (", model$data$nintegers, ",", model$data$nreals, ",", model$data$missing, ") data.txt"),
           file = dir.file, append = TRUE, ncolumns = 12, sep = "")
     write(" ", file = dir.file, append = TRUE)
     write("$VARIABLE", file = dir.file, append = TRUE)
     write(model$data$variables, file = dir.file, append = TRUE, ncolumns = model$data$nvariables)
     write(" ", file = dir.file, append = TRUE)
     
     write("$MODEL", file = dir.file, append = TRUE)
     write(model$nt, file = dir.file, append = TRUE)     # Number of traits
     write(model$absorb, file = dir.file, append = TRUE, ncolumns = 1)     # Weights - one line for each trait
     
     for (i in 1:model$nt) {
          write(c(model$data$reals[model[[i]]$response], 0, model[[i]]$n$factors, model$data$integers[model[[i]]$factors]), 
                file = dir.file, append = TRUE, ncolumns = 3 + model[[i]]$n$factors)
     }
     
     for (i in 1:model$nt) {
          write(c(model[[i]]$n$random,model$data$random[model[[i]]$random]), file = dir.file, append = TRUE,
                ncolumns = 1 + model[[i]]$n$random)
     }
     
     for (i in 1:model$nt) {
          write(c(model[[i]]$n$regression, model$data$reals[model[[i]]$regression]), file = dir.file, append = TRUE, 
                ncolumns = 1 + model[[i]]$n$regression)
     }
     
     write(model$restrict$nresiduals, file = dir.file, append = TRUE)     # Number of residual covariances that are assumed to be zero
     #if (model$restrict$nresiduals == 0)
     #write(model$restrict$residuals, file = dir.file, append = TRUE, ncolumns = 1, sep = " ")    # Trait number combination for zero residual covariance 
     #if (model$restrict$nresiduals > 0)
     write(model$restrict$residuals, file = dir.file, append = TRUE, ncolumns = 2, sep = " ")    # Trait number combination for zero residual covariance 
     write(" ", file = dir.file, append = TRUE)
     
     if (length(model$data$covmat) > 0) {
       for (i in 1:length(model$data$covmat)) {
         vvarname <- names(model$data$covmat)[i]
         vvarfile <- paste(vvarname, ".txt", sep = "")
         write(c("$VAR_STR", model$data$random[vvarname], "COR", "ASCII", vvarfile), file = dir.file, append = TRUE, 
               ncolumns = 5, sep = " ")
       }  
     }
     write(" ", file = dir.file, append = TRUE)
     
     write(c("$DMUAI", format(c(10, 0.0000001, 0.000001, 1, 0, 0), scientific = FALSE)), file = dir.file, append = TRUE)
     write(" ", file = dir.file, append = TRUE)
     
     write("$RESIDUALS ASCII", file = dir.file, append = TRUE)
     
     # Write DMU data file
     data.file <- "data.txt"
     # Recode data (factors only) (check recoding this again)
     rec <- recodeDMU(data[, names(model$data$integers)])
     data[, names(model$data$integers)] <- rec$data
     write.table(format(data[, model$data$variables], scientific = FALSE), file = data.file, 
                 quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
     
     # Write DMU cor data files
     if (!is.null(GRMlist)) {
          for (i in 1:length(model$data$covmat)) {
               vvarname <- names(model$data$covmat)[i]
               vvarfile <- paste(vvarname, ".txt", sep = "")
               iG <- qggginv(GRMlist[[model$data$covmat[i]]], tol = tol)
               colnames(iG$G) <- rownames(iG$G) <- rec$rlevels[[vvarname]][rownames(iG$G)]
               writeGDMU(G = iG$G, filename = vvarfile, ldet = iG$ldet)
          }
     }
     
     return(flevels = rec$flevels)
     
}

# Remove DMU output files 
rename.and.clean.windows <- function(jobname = NULL) {
     
     ll.ff <- list.files()
     "my.system" <- function(cmd) {return(system(paste(Sys.getenv("COMSPEC"), "/c", cmd), show.output.on.console = FALSE))}
     ll.name <- c("SOL", "PAROUT", "PAROUT_STD", "PEDDATA", "RESIDUAL", 
                  "LLIK", "SOL_STD")
     ll <- ll.name[ll.name %in% ll.ff]
     for (kk in ll) {my.system(paste("move ", kk, " ", jobname, ".", kk, sep = ""))}
     junk.files <- c("DMU1.log", "DMUAI.log", paste("COR", 1:20, sep = ""), "CODE_TABLE", "DMU1.dir", "DMUAI.dir", "DMU_LOG", 
                     "DUMMY", "FSPAKWK", "Latest_parm", "LEVAL", "MODINF", 
                     "PARIN", "RCDATA_I", "RCDATA_R", "INBREED", "AINV1", 
                     "AINV2", "PEDFILE1", "PEDFILE2", "fort.81", "fort.66", "fort.99")
     del.files <- ll.ff[ll.ff %in% junk.files]
     if (length(del.files) > 0) { 
          for (kk in 1:length(del.files)) {my.system(paste("del ", del.files[kk], sep = ""))}
     }
     
}

rename.and.clean <- function (jobname = NULL) {
     
     ll.ff <- list.files()
     ll.name <- c("SOL", "PAROUT", "PAROUT_STD", "PEDDATA", "RESIDUAL", "LLIK", "SOL_STD")
     ll <- ll.name[ll.name %in% ll.ff]
     for (kk in ll) {system(paste("mv ", kk, " ", jobname, ".", kk, sep = ""))}
     junk.files <- c("DMU1.log", "DMUAI.log", paste("COR", 1:20, sep = ""), "CODE_TABLE", "DMU1.dir", "DMUAI.dir", "DMU_LOG", 
                     "DUMMY", "FSPAKWK", "Latest_parm", "LEVAL", "MODINF", 
                     "PARIN", "RCDATA_I", "RCDATA_R", "INBREED", "AINV1", 
                     "AINV2", "PEDFILE1", "PEDFILE2", "fort.81", "fort.66", "fort.99")
     del.files <- ll.ff[ll.ff %in% junk.files]
     if (length(del.files)) {
          for (kk in 1:length(del.files)) {system(paste("rm ", del.files[kk], sep = ""))}
     }

}

# Execute DMU 
executeDMU <- function(bin = NULL) {
     
     jobname <- "gfm"
     
     dmu1 <- "dmu1"
     dmuai <- "dmuai"
     
     if (!is.null(bin)) dmu1 <- paste(bin, "dmu1", sep = "/")
     if (!is.null(bin)) dmuai <- paste(bin, "dmuai", sep = "/")
     
     out <- paste(jobname, ".dmuai.lst", sep = "")
     dir <- paste(jobname, ".DIR", sep = "")

     HW <- Sys.info()["machine"]
     OS <- .Platform$OS.type
     
     if (OS == "windows") {
          "my.system" <- function(cmd) {return(system(paste(Sys.getenv("COMSPEC"), "/c", cmd)))}
          my.system("set MKL_NUM_THREADS = 1")
          test <- my.system(paste(shQuote(dmu1), " < ", dir, " > ", out, sep = ""))
          if (test == 0 & "MODINF" %in% list.files()) {
               test <- my.system(paste(shQuote(dmuai), " < ", dir, " >> ", out, sep = ""))
          }
          rename.and.clean.windows(jobname)  
     }
     
     if (!OS == "windows") {
          ll.ff <- list.files()
          if (!("dmu1" %in% ll.ff)) 
          system(paste("cp ", dmu1, " dmu1", sep = ""))
          if (!("dmuai" %in% ll.ff))
          system(paste("cp ", dmuai, " dmuai", sep = ""))
          system("export MKL_NUM_THREADS = 1")
          test <- system(paste("time ./dmu1 < ", dir, "> ", out))
          if (test == 0 & "MODINF" %in% list.files()) {
               test <- system(paste("time ./dmuai >> ", out))
          }
          rename.and.clean(jobname)  
     }
     
}

# Read DMU output files
readDMU <- function(model = NULL, flevels = NULL) {
     
     jobname <- "gfm"
     
     fit <- NULL
     
     llik <- scan(paste(jobname, ".LLIK", sep = ""), what = character(0), quiet = TRUE)
     
     fit$llik <- as.numeric(llik[12])
     cls1 <- grep("Theta", llik)
     cls2 <- grep("ASD", llik)
     #fit$sigma <- as.numeric(llik[(cls1 + 1):(cls2 - 1)])
     #fit$asd <- as.numeric(llik[(cls2 + 1):length(llik)])
     #names(fit$sigma) <- names(fit$asd) <- c(model$random, "e")
     
     sol <- as.matrix(read.table(paste(jobname, ".SOL", sep = ""), as.is = TRUE)[-1,  , drop = FALSE])
     blue <- sol[sol[, 1] == 2, c(2, 4:6, 8:9)]  # "== 2" is estimates effects for fixed factors
     blup <- sol[sol[, 1] == 3, c(2, 4:6, 8:9)]  # "== 3" is estimates effects for random factors

     f <- vector("list", length = model$nt)
     for (i in 1:model$nt) {
        names(f)[i] <- model[[i]]$response
          for (j in 1:model[[i]]$n$random) {
          f[[i]][[j]] <- blup[blup[, 1] == i & blup[, 2] == j, 5:6]
          rownames(f[[i]][[j]]) <- flevels[[model[[i]]$random[j]]][blup[blup[, 1] == i & blup[, 2] == j, 3]]
          colnames(f[[i]][[j]]) <- c("Estimate", "SE") 
          names(f[[i]])[j] <- model[[i]]$random[j]
          } 
     }
     fit$f <- f
     
     sigma <- as.matrix(read.table(paste(jobname, ".PAROUT", sep = "")))
     asd <- as.matrix(read.table(paste(jobname, ".PAROUT_STD", sep = ""), skip = 1, nrow = nrow(sigma)))
     fit$sigma <- cbind(sigma, asd[, 3])
     colnames(fit$sigma) <- c("Random ID", "Row ID", "Col ID", "Estimate", "ASD")
     rownames(fit$sigma) <- 1:nrow(sigma)

     fit$covsigma <- as.matrix(read.table(paste(jobname, ".PAROUT_STD", sep = ""), skip = 1 + nrow(sigma)))
     colnames(fit$covsigma) <- c("Random ID", "Random ID", "Correlation", "ASE")
     rownames(fit$covsigma) <- 1:nrow(fit$covsigma)
     
     resi <- as.matrix(read.table(paste(jobname, ".RESIDUAL", sep = "")))
     nc <- ncol(resi)
     nn.per.tr <- 4 #if (object$glmm) 7 else 4
     
     n.trait <- (nc - 1) / nn.per.tr
     if (n.trait != round(n.trait)) {stop("something wrong")}
     fit$residuals <- resi[, 1 + (((nn.per.tr - 2) * n.trait + 1):((nn.per.tr - 1) * n.trait))]
     fit$fitted <- resi[, 1 + ((1 * n.trait + 1):(2 * n.trait))]
     fit$hat.matrix <- resi[, 1 + (((nn.per.tr - 1) * n.trait + 1):(nn.per.tr * n.trait))]
     
     return(fit)
     
}

vec2mat <- function(vec = NULL, n = NULL, rcnames = NULL) {
     
     X <- diag(n)
     X[lower.tri(X, diag = TRUE)] <- vec
     X <- X + t(X) - diag(diag(X))  
     if(!is.null(rcnames)) {rownames(X) <- colnames(X) <- rcnames}
     X
     
}

writeGDMU <- function(G = NULL, filename = NULL, clear = TRUE, ldet = NULL) {
     
     if (clear) {file.remove(filename)}
     nr <- nrow(G) 
     if (!is.null(ldet)) {write.table(t(c(0, 0, ldet)), filename, quote = F, sep = " ", row.names = F, 
                                      col.names = F, append = TRUE)}
     for (i in 1:nr) { 
          out <- data.frame(rownames(G)[i], rownames(G)[i:nr], G[i:nr, i])
          write.table(out, filename, quote = F, sep = " ", row.names = F, col.names = F, append = TRUE)
     }
     
}

