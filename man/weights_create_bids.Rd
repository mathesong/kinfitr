% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_weights.R
\name{weights_create_bids}
\alias{weights_create_bids}
\title{Create weights from BIDS data}
\usage{
weights_create_bids(
  petinfo,
  tac,
  method = 2,
  minweight = 0.7,
  minweight_risetopeak = FALSE,
  weight_checkn = 5
)
}
\arguments{
\item{petinfo}{The parsed information from the PET JSON sidecar. This can be
extracted simply using jsonlite::fromJSON() on the JSON file, or else at a
study level using bids_parse_study()}

\item{tac}{The time activity curve: preferably for a region with high signal-to-noise ratio.}

\item{method}{Which method should be used? 1 represents duration^2 /
(tac_uncorrected). 2 represents sqrt(durations*tac_uncorrected). 3
represents duration / tac. 4 represents sqrt(durations). 5 represents
durations * exp((-ln(2)) / halflife ). 6 represents durations /
tac. 7 represents durations.}

\item{minweight}{The minimum weight. Weights will be calculated as a fraction
between this value and 1. A zero frame with duration=0 will be set to 0
though.}

\item{minweight_risetopeak}{Should there be a linear rise of the minimum weight to
the peak weight value? This means that values before the maximum weight can
be below the minweight specified. This is helpful for downweighting noisy
values towards the start of the TAC. Defaults to FALSE.}

\item{weight_checkn}{The number of values of the weights to check to make
sure that things haven't gone terribly wrong. It will check that this
number of weights are all above half of the maximum weight.}
}
\value{
A vector of weights.
}
\description{
This uses the BIDS PET JSON sidecar information and a TAC to create a
weights vector.
}
\examples{
\dontrun{
bids_weights_create(studydata$petinfo[[1]],tacdata$Neo.cx)
}
}
\author{
Granville J Matheson, \email{mathesong@gmail.com}
}
