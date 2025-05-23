% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allelic_imbalance.R
\name{correct_theta}
\alias{correct_theta}
\title{Fits a local regression model with dispersion as a function of total gene counts.
Predicted values are the expected dispersion for genes with similar expression levels
The initial dispersion estimates are shrunk towards the common dispersion}
\usage{
correct_theta(estimates, delta_set = 50, N_set = 30, thetaFilter = NULL)
}
\arguments{
\item{estimates}{Output of estim_bbparams}

\item{delta_set}{Delta parameter}

\item{N_set}{N parameter}

\item{thetaFilter}{Minimum dispersion value, genes with dispersion below the set threshold
are excluded from the shrinking procedure.}
}
\description{
Fits a local regression model with dispersion as a function of total gene counts.
Predicted values are the expected dispersion for genes with similar expression levels
The initial dispersion estimates are shrunk towards the common dispersion
}
\examples{
correct_theta()
}
