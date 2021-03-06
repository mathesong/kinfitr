% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_bids.R
\name{bids_parse_blood}
\alias{bids_parse_blood}
\title{Extract blood data from BIDS study folder}
\usage{
bids_parse_blood(filedata)
}
\arguments{
\item{filedata}{A table of the file data, created using bids_parse_files()}
}
\value{
Data from these files: the Data section of the blooddata object
}
\description{
Extracts the Data section of the blooddata object
}
\examples{
\dontrun{
bd_dat <- bids_parse_blood(filedata)
}
}
