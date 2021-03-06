% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{zipnegloglik_nocov_cont}
\alias{zipnegloglik_nocov_cont}
\title{negative log likelihood function for zero-inflated Poisson hidden Markov model without covariates,
where zero-inflation only happens in state 1}
\usage{
zipnegloglik_nocov_cont(parm, M, y, ntimes, timeindex, udiff)
}
\arguments{
\item{parm}{working parameters}

\item{M}{number of hidden states}

\item{y}{observed series}

\item{ntimes}{length of the observed series}

\item{timeindex}{vector of observed time points}

\item{udiff}{unique time intervals}
}
\value{
negative log likelihood
}
\description{
negative log likelihood function for zero-inflated Poisson hidden Markov model without covariates,
where zero-inflation only happens in state 1
}
