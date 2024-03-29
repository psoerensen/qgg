% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomic_matrix.R
\name{getMap}
\alias{getMap}
\title{Retrieve the map for specified rsids on a given chromosome.}
\usage{
getMap(Glist = NULL, chr = NULL, rsids = NULL)
}
\arguments{
\item{Glist}{A list structure with information about genotypes stored on disk.}

\item{chr}{A chromosome from which the map is retrieved.}

\item{rsids}{A vector of rsids for which the map is needed.}
}
\value{
A vector containing the map corresponding to the specified rsids on the given chromosome.
}
\description{
Fetch the map associated with provided rsids for a given chromosome from the list `Glist`.
}
\keyword{internal}
