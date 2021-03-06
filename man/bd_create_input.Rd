% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_blooddata.R
\name{bd_create_input}
\alias{bd_create_input}
\title{Create an input object from a blooddata object}
\usage{
bd_create_input(blooddata, startTime = 0, stopTime = NULL, interpPoints = 6000)
}
\arguments{
\item{blooddata}{The blooddata object.}

\item{startTime}{The starting time for the interpolation. Defaults to zero.
If, after application of the TimeShift value in the blooddata object, the
startTime is still after zero, it will be set to zero.}

\item{stopTime}{The end time for the interpolation. Defaults to the maximum
measured time.}

\item{interpPoints}{The number of points to interpolate over between the
start and stop times. Defaults to 6000.}
}
\value{
An input object
}
\description{
Get an input object for kinetic modelling
}
\examples{
\dontrun{
bd_create_input(blooddata)
}
}
\author{
Granville J Matheson, \email{mathesong@gmail.com}
}
