% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_blooddata.R
\name{create_blooddata_bids}
\alias{create_blooddata_bids}
\title{Create a blooddata object from BIDS data}
\usage{
create_blooddata_bids(bids_data, TimeShift = 0)
}
\arguments{
\item{bids_data}{The filename of a PET BIDS json sidecar, or a list
containing the information contained within a PET BIDS json sidecar.}

\item{TimeShift}{The extent to which all the times in the data should be
shifted (in seconds). Defaults to 0.}
}
\value{
A blooddata object
}
\description{
Deprecated. This function creates a blooddata object from JSON data structured according to
the old PET BIDS standard.
}
\examples{
\dontrun{
a <- create_blooddata_bids("bids_sidecar.json")
}
}
\author{
Granville J Matheson, \email{mathesong@gmail.com}
}
