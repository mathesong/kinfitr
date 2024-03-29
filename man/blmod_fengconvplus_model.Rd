% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_bloodmodels.R
\name{blmod_fengconvplus_model}
\alias{blmod_fengconvplus_model}
\title{Convolved Feng model for fitting of arterial input functions with an additional rise term}
\usage{
blmod_fengconvplus_model(
  time,
  t0,
  A,
  alpha,
  B,
  beta,
  C,
  gamma,
  ti,
  asymptote,
  slope
)
}
\arguments{
\item{time}{Time of each sample.}

\item{t0}{The delay time. This is the point at which the linear rise begins.}

\item{A}{The multiplier of the first exponential.}

\item{alpha}{The rate of the first exponential.}

\item{B}{The multiplier of the second exponential.}

\item{beta}{The rate of the second exponential.}

\item{C}{The multiplier of the third exponential.}

\item{gamma}{The rate of the third exponential.}

\item{ti}{The infusion time of the bolus.}

\item{asymptote}{The value of the asymptote to which the rise rises.}

\item{slope}{The slope of the rise term, i.e. the rate at which it rises.}
}
\value{
Model predictions
}
\description{
This is the model itself for the Fengconvplus model of the AIF with an extended
injection and additional rising term.
}
\examples{
blmod_fengconvplus_model(1:1000, 30, 220, 0.4, 100, 0.05, 22, 0.003, 30, 40, 0.001)
}
