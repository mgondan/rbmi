% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsmeans.R
\name{lsmeans}
\alias{lsmeans}
\title{Least Square Means}
\usage{
lsmeans(
  model,
  ...,
  .weights = c("counterfactual", "equal", "proportional_em", "proportional")
)
}
\arguments{
\item{model}{A model created by \code{lm}.}

\item{...}{Fixes specific variables to specific values i.e.
\code{trt = 1} or \code{age = 50}. The name of the argument must be the name
of the variable within the dataset.}

\item{.weights}{Character, either \code{"counterfactual"} (default), \code{"equal"},
\code{"proportional_em"} or \code{"proportional"}.
Specifies the weighting strategy to be used when calculating the lsmeans.
See the weighting section for more details.}
}
\description{
Estimates the least square means from a linear model. The exact implementation
/ interpretation depends on the weighting scheme; see the weighting section for more
information.
}
\section{Weighting}{

\subsection{Counterfactual}{

For \code{weights = "counterfactual"} (the default) the lsmeans are obtained by
taking the average of the predicted values for each patient after assigning all patients
to each arm in turn.
This approach is equivalent to standardization or g-computation.
In comparison to \code{emmeans} this approach is equivalent to:

\if{html}{\out{<div class="sourceCode">}}\preformatted{emmeans::emmeans(model, specs = "<treatment>", counterfactual = "<treatment>")
}\if{html}{\out{</div>}}

Note that to ensure backwards compatibility with previous versions of \code{rbmi}
\code{weights = "proportional"} is an alias for \code{weights = "counterfactual"}.
To get results consistent with \code{emmeans}'s \code{weights = "proportional"}
please use \code{weights = "proportional_em"}.
}

\subsection{Equal}{

For \code{weights = "equal"} the lsmeans are obtained by taking the model fitted
value of a hypothetical patient whose covariates are defined as follows:
\itemize{
\item Continuous covariates are set to \code{mean(X)}
\item Dummy categorical variables are set to \code{1/N} where \code{N} is the number of levels
\item Continuous * continuous interactions are set to \code{mean(X) * mean(Y)}
\item Continuous * categorical interactions are set to \code{mean(X) * 1/N}
\item Dummy categorical * categorical interactions are set to \code{1/N * 1/M}
}

In comparison to \code{emmeans} this approach is equivalent to:

\if{html}{\out{<div class="sourceCode">}}\preformatted{emmeans::emmeans(model, specs = "<treatment>", weights = "equal")
}\if{html}{\out{</div>}}
}

\subsection{Proportional}{

For \code{weights = "proportional_em"} the lsmeans are obtained as per \code{weights = "equal"}
except instead of weighting each observation equally they are weighted by the proportion
in which the given combination of categorical values occurred in the data.
In comparison to \code{emmeans} this approach is equivalent to:

\if{html}{\out{<div class="sourceCode">}}\preformatted{emmeans::emmeans(model, specs = "<treatment>", weights = "proportional")
}\if{html}{\out{</div>}}

Note that this is not to be confused with \code{weights = "proportional"} which is an alias
for \code{weights = "counterfactual"}.
}
}

\section{Fixing}{


Regardless of the weighting scheme any named arguments passed via \code{...} will
fix the value of the covariate to the specified value.
For example, \code{lsmeans(model, trt = "A")} will fix the dummy variable \code{trtA} to 1
for all patients (real or hypothetical) when calculating the lsmeans.

See the references for similar implementations as done in SAS and
in R via the \code{emmeans} package.
}

\examples{
\dontrun{
mod <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
lsmeans(mod)
lsmeans(mod, Species = "virginica")
lsmeans(mod, Species = "versicolor")
lsmeans(mod, Species = "versicolor", Petal.Length = 1)
}
}
\references{
\url{https://CRAN.R-project.org/package=emmeans}

\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_glm_details41.htm}
}
