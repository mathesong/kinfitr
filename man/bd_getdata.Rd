% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_blooddata.R
\name{bd_getdata}
\alias{bd_getdata}
\title{Create an input object from a blooddata object.}
\usage{
bd_getdata(
  blooddata,
  startTime = 0,
  stopTime = NULL,
  interpPoints = 6000,
  output = c("input", "Blood", "BPR", "parentFraction", "AIF")
)
}
\arguments{
\item{blooddata}{A blooddata object, to which all the desired fits have been applied and added.}

\item{startTime}{The starting time for the interpolation. Defaults to zero. If, after application of the TimeShift value in the blooddata object, the startTime is still after zero, it will be set to zero.}

\item{stopTime}{The end time for the interpolation. Defaults to the maximum measured time.}

\item{interpPoints}{The number of points to interpolate over between the start and stop times.}

\item{output}{The output. This defaults to an "input" object, which can be used in a kinetic model fit. But if set to "Blood", "BPR", "parentFraction" or "AIF", it yields the appropriate input for the function which will be used to model these.}
}
\value{
A tibble containing the output specified.
}
\description{
DEPRECATION WARNING: this function will be slowly phased out of future
releases as it is based on the old PET BIDS standard. in favour of
bd_extract() and bd_create_input(). blooddata objects can be updated to the
new format using update_blooddata(), or (better) by creating it afresh from
the source.
}
\details{
This function extracts data from blooddata objects to create either an
input object for kinetic modelling, or sets of values for modelling of
blood-related curves.
}
\examples{
\dontrun{
bd_getdata(blooddata)
bd_getdata(blooddata, output = "parentFraction")
}
}
\author{
Granville J Matheson, \email{mathesong@gmail.com}
}
