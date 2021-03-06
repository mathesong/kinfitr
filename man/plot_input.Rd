% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_bloodfuncs.R
\name{plot_input}
\alias{plot_input}
\title{Plot input}
\usage{
plot_input(input)
}
\arguments{
\item{input}{Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
using the \code{blood_interp} function.}
}
\value{
A ggplot2 object
}
\description{
Provides a plot of the input data to get an idea of how it looks. Blood, plasma and AIF data are scaled to a common value.
}
\examples{

data(pbr28)

input <- blood_interp(
  pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
  pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
  t_parentfrac = 1, parentfrac = 1
)

plot_input(input)
}
\author{
Granville J Matheson, \email{mathesong@gmail.com}
}
