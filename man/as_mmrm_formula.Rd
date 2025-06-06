% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmrm.R
\name{as_mmrm_formula}
\alias{as_mmrm_formula}
\title{Create MMRM formula}
\usage{
as_mmrm_formula(mmrm_df, cov_struct)
}
\arguments{
\item{mmrm_df}{an mmrm \code{data.frame} as created by \code{\link[=as_mmrm_df]{as_mmrm_df()}}}

\item{cov_struct}{Character - The covariance structure to be used, must be one of \code{"us"} (default),
\code{"ad"}, \code{"adh"}, \code{"ar1"}, \code{"ar1h"}, \code{"cs"}, \code{"csh"}, \code{"toep"}, or \code{"toeph"})}
}
\description{
Derives the MMRM model formula from the structure of mmrm_df.
returns a formula object of the form:
}
\details{
\if{html}{\out{<div class="sourceCode">}}\preformatted{outcome ~ 0 + V1 + V2 + V4 + ... + us(visit | group / subjid)
}\if{html}{\out{</div>}}
}
