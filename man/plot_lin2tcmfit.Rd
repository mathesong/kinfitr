% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_lin2tcm.R
\name{plot_lin2tcmfit}
\alias{plot_lin2tcmfit}
\title{Plot: Linear 2TCM}
\usage{
plot_lin2tcmfit(lin2tcmout, roiname = NULL)
}
\arguments{
\item{lin2tcmout}{The output object of the linearised 2TCM fitting procedure.}

\item{roiname}{Optional. The name of the Target Region to see it on the plot.}
}
\value{
A ggplot2 object of the plot.
}
\description{
Function to visualise the fit of the linearised 2TCM model to data.
}
\examples{
data(pbr28)

t_tac <- pbr28$tacs[[2]]$Times / 60
tac <- pbr28$tacs[[2]]$FC
weights <- pbr28$tacs[[2]]$Weights

input <- blood_interp(
  pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
  pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
  t_parentfrac = 1, parentfrac = 1
)

fit <- lin2tcm(t_tac, tac, input, weights, inpshift=0.1)
plot_lin2tcmfit(fit)
}
\author{
Granville J Matheson, \email{mathesong@gmail.com}
}
