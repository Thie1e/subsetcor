% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_functions.R
\name{se_z_finite}
\alias{se_z_finite}
\title{Functions for Fisher's z-transformation}
\usage{
se_z_finite(Npop, n)
}
\arguments{
\item{Npop}{Population size in case of finite populations, e.g. when a subset was
extracted from a sample}

\item{n}{Sample size}

\item{r}{Sample correlation coefficient}

\item{rho}{Population correlation coefficient}
}
\description{
Convert a sample correlation coefficient (r) to Fisher's z or the other way
around, compute the expected value of z (E_z), compute the standard error of z
(se_z), and the standard error of z with a correction for finite samples (se_z_finite) with
Npop being the size of the full sample.
}
\seealso{
\code{\link{rtoz}}, \code{\link{ztor}}, \code{\link{E_z}}, \code{\link{se_z}}, \code{\link{se_z_finite}}
}

