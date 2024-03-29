% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_pfmodels.R
\name{metab_gamma}
\alias{metab_gamma}
\title{Fit the Gamma Function for Modelling Parent Fraction.}
\usage{
metab_gamma(
  time,
  parentFraction,
  fit_ppf0 = FALSE,
  fit_delay = FALSE,
  lower = list(a = 0.8, b = 0.3, c = 0, d = 0, delay = -30),
  upper = list(a = 1.1, b = 3, c = 20, d = 1e+06, delay = 30),
  multstart_lower = NULL,
  multstart_upper = NULL,
  multstart_iter = 1000
)
}
\arguments{
\item{time}{Time in seconds.}

\item{parentFraction}{Measured values of parent fraction.}

\item{fit_ppf0}{Should the starting plasma parent fraction be fitted? Otherwise, it is set to 1. Defaults to FALSE.}

\item{fit_delay}{Should the delay of the plasma parent fraction be fitted? Otherwise, it is set to 0. Defaults to FALSE.}

\item{lower}{Named list of the lower limits.}

\item{upper}{Named list of the upper limits.}

\item{multstart_lower}{Named list of the lower starting limits.}

\item{multstart_upper}{Named list of the upper starting limits.}

\item{multstart_iter}{Number of fits to perform before deciding on an optimal.}
}
\value{
An nls fit object.
}
\description{
This function fits the integrated gamma function to parent fraction data.
}
\examples{
\dontrun{
pf <- bd_getdata(blooddata, output = "parentFraction")
metab_gamma(pf$time, pf$parentFraction)
}
}
\references{
Naganawa M, Jacobsen LK, Zheng MQ, Lin SF, Banerjee A, Byon W, Weinzimmer D, Tomasi G, Nabulsi N, Grimwood S, Badura LL. Evaluation of the agonist PET radioligand [11C] GR103545 to image kappa opioid receptor in humans: Kinetic model selection, test–retest reproducibility and receptor occupancy by the antagonist PF-04455242. Neuroimage. 2014 Oct 1;99:69-79.
}
