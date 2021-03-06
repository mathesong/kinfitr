% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_refmllogan.R
\name{refmlLogan}
\alias{refmlLogan}
\title{Non-Invasive Multilinear Logan Plot}
\usage{
refmlLogan(
  t_tac,
  reftac,
  roitac,
  k2prime,
  tstarIncludedFrames,
  weights = NULL,
  dur = NULL,
  frameStartEnd = NULL
)
}
\arguments{
\item{t_tac}{Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
zero. If a time zero frame is not included, it will be added.}

\item{reftac}{Numeric vector of radioactivity concentrations in the reference tissue for each frame. We include zero at
time zero: if not included, it is added.}

\item{roitac}{Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
zero: if not included, it is added.}

\item{k2prime}{Value of k2prime to be used for the fitting, i.e. the average tissue-to-plasma clearance rate. This can be
obtained from another model, or set at a specified value. If using SRTM to estimate this value, it is equal to k2 / R1.}

\item{tstarIncludedFrames}{The number of frames to be used in the regression model, i.e. the number of frames for which
the function is linear after pseudo-equilibrium is reached. This is a count from the end of the measurement, so a value of
10 means that last 10 frames will be used. This value can be estimated using \code{refmlLogan_tstar}.}

\item{weights}{Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
if not included, it is added. If not specified, uniform weights will be used.}

\item{dur}{Optional. Numeric vector of the time durations of the frames. If
not included, the integrals will be calculated using trapezoidal integration.}

\item{frameStartEnd}{Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
This is to assess time stability.}
}
\value{
A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit}, a dataframe
containing the TACs of the data \code{out$tacs}, a dataframe containing the TACs of the fitted values \code{out$fitvals},
a vector of the weights \code{out$weights}, the specified k2prime value \code{out$k2prime}, and the specified
tstarIncludedFrames value \code{out$tstarIncludedFrames}
}
\description{
Function to fit the non-invasive modification of the multilinear Logan plot model of Turkheimer et al. (2003) to data.
}
\examples{
data(simref)

t_tac <- simref$tacs[[2]]$Times
reftac <- simref$tacs[[2]]$Reference
roitac <- simref$tacs[[2]]$ROI1
weights <- simref$tacs[[2]]$Weights

fit <- refmlLogan(t_tac, reftac, roitac, k2prime = 0.1, tstarIncludedFrames = 10, weights = weights)
}
\references{
Turkheimer FE, Aston JA, Banati RB, Riddell C, Cunningham VJ. A linear wavelet filter for parametric imaging with dynamic PET. IEEE transactions on medical imaging. 2003 Mar;22(3):289-301.
}
\author{
Granville J Matheson, \email{mathesong@gmail.com}
}
