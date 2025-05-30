% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{fit_mcmc}
\alias{fit_mcmc}
\title{Fit the base imputation model using a Bayesian approach}
\usage{
fit_mcmc(designmat, outcome, group, subjid, visit, method, quiet = FALSE)
}
\arguments{
\item{designmat}{The design matrix of the fixed effects.}

\item{outcome}{The response variable. Must be numeric.}

\item{group}{Character vector containing the group variable.}

\item{subjid}{Character vector containing the subjects IDs.}

\item{visit}{Character vector containing the visit variable.}

\item{method}{A \code{method} object as generated by \code{\link[=method_bayes]{method_bayes()}}.}

\item{quiet}{Specify whether the stan sampling log should be printed to the console.}
}
\value{
A named list composed by the following:
\itemize{
\item \code{samples}: a named list containing the draws for each parameter. It corresponds to the output of \code{\link[=extract_draws]{extract_draws()}}.
\item \code{fit}: a \code{stanfit} object.
}
}
\description{
\code{fit_mcmc()} fits the base imputation model using a Bayesian approach.
This is done through a MCMC method that is implemented in \code{stan}
and is run by using the function \code{rstan::sampling()}.
The function returns the draws from the posterior distribution of the model parameters
and the \code{stanfit} object. Additionally it performs multiple diagnostics checks of the chain
and returns warnings in case of any detected issues.
}
\details{
The Bayesian model assumes a multivariate normal likelihood function and weakly-informative
priors for the model parameters: in particular, uniform priors are assumed for the regression
coefficients and inverse-Wishart priors for the covariance matrices.
The chain is initialized using the REML parameter estimates from MMRM as starting values.

The function performs the following steps:
\enumerate{
\item Fit MMRM using a REML approach.
\item Prepare the input data for the MCMC fit as described in the \verb{data\{\}}
block of the Stan file. See \code{\link[=prepare_stan_data]{prepare_stan_data()}} for details.
\item Run the MCMC according the input arguments and using as starting values the REML parameter estimates
estimated at point 1.
\item Performs diagnostics checks of the MCMC. See \code{\link[=check_mcmc]{check_mcmc()}} for details.
\item Extract the draws from the model fit.
}

The chains perform \code{method$n_samples} draws by keeping one every \code{method$burn_between} iterations. Additionally
the first \code{method$burn_in} iterations are discarded. The total number of iterations will
then be \code{method$burn_in + method$burn_between*method$n_samples}.
The purpose of \code{method$burn_in} is to ensure that the samples are drawn from the stationary
distribution of the Markov Chain.
The \code{method$burn_between} aims to keep the draws uncorrelated each from other.
}
