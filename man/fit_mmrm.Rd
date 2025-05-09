% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmrm.R
\name{fit_mmrm}
\alias{fit_mmrm}
\title{Fit a MMRM model}
\usage{
fit_mmrm(
  designmat,
  outcome,
  subjid,
  visit,
  group,
  cov_struct = c("us", "ad", "adh", "ar1", "ar1h", "cs", "csh", "toep", "toeph"),
  REML = TRUE,
  same_cov = TRUE
)
}
\arguments{
\item{designmat}{a \code{data.frame} or \code{matrix} containing the covariates to use in the MMRM model.
Dummy variables must already be expanded out, i.e. via \code{\link[stats:model.matrix]{stats::model.matrix()}}. Cannot contain
any missing values}

\item{outcome}{a numeric vector. The outcome value to be regressed on in the MMRM model.}

\item{subjid}{a character / factor vector. The subject identifier used to link separate visits
that belong to the same subject.}

\item{visit}{a character / factor vector. Indicates which visit the outcome value occurred on.}

\item{group}{a character / factor vector. Indicates which treatment group the patient belongs to.}

\item{cov_struct}{a character value. Specifies which covariance structure to use. Must be one of \code{"us"} (default),
\code{"ad"}, \code{"adh"}, \code{"ar1"}, \code{"ar1h"}, \code{"cs"}, \code{"csh"}, \code{"toep"}, or \code{"toeph"})}

\item{REML}{logical. Specifies whether restricted maximum likelihood should be used}

\item{same_cov}{logical. Used to specify if a shared or individual covariance matrix should be
used per \code{group}}
}
\description{
Fits a MMRM model allowing for different covariance structures using \code{\link[mmrm:mmrm]{mmrm::mmrm()}}.
Returns a \code{list} of key model parameters \code{beta}, \code{sigma} and an additional element \code{failed}
indicating whether or not the fit failed to converge. If the fit did fail to converge
\code{beta} and \code{sigma} will not be present.
}
