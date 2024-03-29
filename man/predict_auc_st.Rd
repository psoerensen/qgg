% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qgg_utility_functions.R
\name{predict_auc_st}
\alias{predict_auc_st}
\title{Expected AUC for prediction of a binary trait}
\usage{
predict_auc_st(h2x, Nx, Kx, Px, Me, M)
}
\arguments{
\item{h2x}{Heritability of the target trait.}

\item{Nx}{Number of samples for the target trait.}

\item{Kx}{Prevalence of the target trait.}

\item{Px}{Case-control proportion of the target trait.}

\item{Me}{Number of independent chromosome segments.}

\item{M}{Number of markers.}
}
\value{
A numeric value representing the expected AUC.
}
\description{
Computes the expected Area Under the Curve (AUC) for predicting a binary trait.
}
\keyword{internal}
