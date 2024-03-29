% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomic_matrix.R
\name{gfilter}
\alias{gfilter}
\title{Filter genetic marker data based on different quality measures}
\usage{
gfilter(
  Glist = NULL,
  excludeMAF = 0.01,
  excludeMISS = 0.05,
  excludeINFO = NULL,
  excludeCGAT = TRUE,
  excludeINDEL = TRUE,
  excludeDUPS = TRUE,
  excludeHWE = 1e-12,
  excludeMHC = FALSE,
  assembly = "GRCh37"
)
}
\arguments{
\item{Glist}{A list containing information about the genotype matrix stored on disk.}

\item{excludeMAF}{A scalar threshold. Exclude markers with a minor allele frequency (MAF) below this threshold. Default is 0.01.}

\item{excludeMISS}{A scalar threshold. Exclude markers with missingness (MISS) above this threshold. Default is 0.05.}

\item{excludeINFO}{A scalar threshold. Exclude markers with an info score (INFO) below this threshold. Default is 0.8.}

\item{excludeCGAT}{A logical value; if TRUE exclude markers if the alleles are ambiguous (i.e., either CG or AT combinations).}

\item{excludeINDEL}{A logical value; if TRUE exclude markers that are insertions or deletions (INDELs).}

\item{excludeDUPS}{A logical value; if TRUE exclude markers if their identifiers are duplicated.}

\item{excludeHWE}{A scalar threshold. Exclude markers where the p-value for the Hardy-Weinberg Equilibrium test is below this threshold. Default is 0.01.}

\item{excludeMHC}{A logical value; if TRUE exclude markers located within the MHC region.}

\item{assembly}{A character string indicating the name of the genome assembly (e.g., "GRCh38").}
}
\description{
Quality control is a critical step for working with summary statistics (in particular
                                                                        for external). 
Processing and quality control of GWAS summary statistics includes:                                                                      

- map marker ids (rsids/cpra (chr, pos, ref, alt)) to LD reference panel data 
- check effect allele (flip EA, EAF, Effect)
- check effect allele frequency
- thresholds for MAF and HWE
- exclude INDELS, CG/AT and MHC region
- remove duplicated marker ids
- check which build version
- check for concordance between marker effect and LD data

External summary statistics format:
 marker, chr, pos, effect_allele, non_effect_allele, effect_allele_freq, effect, effect_se, stat, p, n    

Internal summary statistics format:
 rsids, chr, pos, a1, a2, af, b, seb, stat, p, n
}
\author{
Peter Soerensen
}
