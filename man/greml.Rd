% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomic_reml.R
\name{greml}
\alias{greml}
\title{Genomic REML analysis}
\usage{
greml(y = NULL, X = NULL, GRMlist = NULL, GRM = NULL,
  theta = NULL, ids = NULL, validate = NULL, maxit = 100,
  tol = 1e-05, bin = NULL, ncores = 1, wkdir = getwd(),
  verbose = FALSE, makeplots = FALSE, interface = "R", fm = NULL,
  data = NULL)
}
\arguments{
\item{y}{vector or matrix of phenotypes}

\item{X}{design matrix for factors modeled as fixed effects}

\item{GRMlist}{list providing information about GRM matrix stored in binary files on disk}

\item{GRM}{list of one or more genomic relationship matrices}

\item{theta}{vector of initial values of co-variance for REML estimation}

\item{ids}{vector of individuals used in the analysis}

\item{validate}{dataframe or list of individuals used in cross-validation (one column/row for each validation set)}

\item{maxit}{maximum number of iterations used in REML analysis}

\item{tol}{tolerance, i.e. convergence criteria used in REML}

\item{bin}{directory for fortran binaries (e.g. DMU binaries dmu1 and dmuai)}

\item{ncores}{number of cores used for the analysis}

\item{wkdir}{is the working directory used for REML}

\item{verbose}{logical if TRUE print more details during optimization}

\item{makeplots}{logical if TRUE makes some plots or parameter estimates and prediction accuracy during cross validation}

\item{interface}{used for specifying whether to use R or Fortran implementations of REML}

\item{fm}{formula with model statement for the linear mixed model}

\item{data}{data frame containing the phenotypic observations and fixed factors specified in the model statements}
}
\value{
Returns a list structure including
\item{llik}{log-likelihood at convergence}
\item{theta}{covariance estimates from REML}
\item{asd}{asymptotic standard deviation}
\item{b}{vector of fixed effect estimates}
\item{varb}{vector of variances of fixed effect estimates}
\item{g}{vector or matrix of random effect estimates}
\item{e}{vector or matrix of residual effects}
\item{accuracy}{matrix of prediction accuracies (only returned if validate is provided)}
}
\description{
The greml function is used for estimation of genomic parameters (co-variance, heritability and correlation)
for linear mixed models using restricted maximum likelihood estimation (REML) and genomic prediction using
best linear unbiased prediction (BLUP).

The linear mixed model can account for multiple genetic factors (fixed and random genetic marker effects),
adjust for complex family relationships or population stratification, and adjust for other non-genetic factors
including lifestyle characteristics. Different genetic architectures (infinitesimal, few large and many
small effects) is accounted for by modeling genetic markers in different sets as fixed or random effects
and by specifying individual genetic marker weights. Different genetic models (e.g. additive and non-additive)
can be specified by providing additive and non-additive genomic relationship matrices (GRMs) (constructed using grm).
The GRMs can be accessed from the R environment or from binary files stored on disk facilitating analyses of
large-scale genetic data.

The output contains estimates of variance components, fixed and random effects, first and second derivatives of
log-likelihood, and the asymptotic standard deviation of parameter estimates.

Assessment of predictive accuracy (including correlation and R2, and AUC for binary phenotypes) can be obtained
by providing greml with a dataframe or list containing sample IDs used in the validation, see examples for details.

Genomic parameters can also be estimated with DMU (http://www.dmu.agrsci.dk/DMU/) if interface =”DMU”.
This option requires DMU to be installed locally, and the path to the DMU binary files has to be specified
(see examples below for details).
}
\examples{

\dontrun{

# Simulate data
W <- matrix(rnorm(20000000), ncol = 10000)
	colnames(W) <- as.character(1:ncol(W))
	rownames(W) <- as.character(1:nrow(W))
y <- rowSums(W[, 1:10]) + rowSums(W[, 1001:1010]) + rnorm(nrow(W))

# Create model
data <- data.frame(y = y, mu = 1)
fm <- y ~ 0 + mu
X <- model.matrix(fm, data = data)

# Compute GRM
GRM <- grm(W = W)

# REML analyses
fitG <- greml(y = y, X = X, GRM = list(GRM))

# REML analyses and cross validation

# Create marker sets
setsGB <- list(A = colnames(W)) # gblup model
setsGF <- list(C1 = colnames(W)[1:1000], C2 = colnames(W)[1001:2000],
               C3 = colnames(W)[2000:10000]) # gfblup model
setsGT <- list(C1 = colnames(W)[1:10], C2 = colnames(W)[1001:1010],
               C3 = colnames(W)[1:10000]) # true model

GB <- lapply(setsGB, function(x) {grm(W = W[, x])})
GF <- lapply(setsGF, function(x) {grm(W = W[, x])})
GT <- lapply(setsGT, function(x) {grm(W = W[, x])})

n <- length(y)
fold <- 10
nvalid <- 5

validate <- replicate(nvalid, sample(1:n, as.integer(n / fold)))
cvGB <- greml(y = y, X = X, GRM = GB, validate = validate)
cvGF <- greml(y = y, X = X, GRM = GF, validate = validate)
cvGT <- greml(y = y, X = X, GRM = GT, validate = validate)

cvGB$accuracy
cvGF$accuracy
cvGT$accuracy

}

}
\references{
Lee, S. H., & van Der Werf, J. H. (2006). An efficient variance component approach implementing an average information REML suitable for combined LD and linkage mapping with a general complex pedigree. Genetics Selection Evolution, 38(1), 25.
}
\author{
Peter Soerensen
}
