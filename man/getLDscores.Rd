% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomic_matrix.R
\name{getLDscores}
\alias{getLDscores}
\title{Extract Linkage Disequilibrium (LD) Scores}
\usage{
getLDscores(Glist = NULL, chr = NULL, rsids = NULL)
}
\arguments{
\item{Glist}{A list structure containing genotypic data, including rsids for LD calculation (`rsidsLD`), LD file locations (`ldfiles`), and `msize` which indicates the size for surrounding region to consider for LD.}

\item{chr}{A specific chromosome from which LD sets need to be extracted.}

\item{rsids}{A vector of rsids that need to be included in the sparse LD matrix. Default is NULL, implying all rsids in the chromosome will be used.}
}
\value{
A vector containing LD scores.
}
\description{
This function retrieves LD scores from a structured list (`Glist`) that includes
LD scores, chromosome information, and SNP identifiers. You can extract all LD scores
or filter them by a specific chromosome and/or SNP identifiers (rsids).
}
\keyword{internal}
