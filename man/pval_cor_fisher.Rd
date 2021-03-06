% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_functions.R
\name{pval_cor_fisher}
\alias{pval_cor_fisher}
\title{p-value of sample correlation coefficient (also in finite populations)}
\usage{
pval_cor_fisher(r, rho, n, Npop, alternative = "two.sided")
}
\arguments{
\item{r}{sample correlation coefficient}

\item{rho}{population correlation coefficient}

\item{n}{sample size}

\item{Npop}{population size}

\item{alternative}{a character string specifying the alternative hypothesis,
must be one of "two.sided" (default), "greater" or "less"}
}
\description{
Based on Fisher's z-transformation and a correction of the standard deviation for
finite samples: (Npop - n) / (Npop - 1)
Calculates probability of observing a sample correlation coefficient r given
a sample size n, a population size Npop and the population correlation coefficient rho.
}

