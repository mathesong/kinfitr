% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_miscfuncs.R
\name{coef_kinfit}
\alias{coef_kinfit}
\title{Extract coefficients}
\usage{
coef_kinfit(x, ...)
}
\arguments{
\item{x}{The output object of the model fitting procedure.}

\item{...}{Additional optional arguments.}
}
\description{
Function to extract the coefficients from kinfit objects
}
\examples{
\dontrun{
loganout <- Loganplot(t_tac, tac, input, 10, weights)
coef(loganout)
}

}
\author{
Granville J Matheson, \email{mathesong@gmail.com}
}
