% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allelic_imbalance.R
\name{estim_delta}
\alias{estim_delta}
\title{Estimating appropriate tuning paramater delta and the number of degrees of freedom, N.
We assume that dispersion follows Gamma distribution. Appropriate shrinkage is selected based on
the MLE of the the difference between fitted dispersion and shrunk dispersion}
\usage{
estim_delta(estimates)
}
\arguments{
\item{estimates}{Output of estim_bbparams}
}
\description{
Estimating appropriate tuning paramater delta and the number of degrees of freedom, N.
We assume that dispersion follows Gamma distribution. Appropriate shrinkage is selected based on
the MLE of the the difference between fitted dispersion and shrunk dispersion
}
\examples{
estim_delta()
}
