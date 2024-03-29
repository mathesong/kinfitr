% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_bloodfuncs.R
\name{blood_interp}
\alias{blood_interp}
\title{Interpolate Blood Curves}
\usage{
blood_interp(
  t_blood,
  blood,
  t_plasma = NULL,
  plasma = NULL,
  t_parentfrac = NULL,
  parentfrac = 1,
  t_aif = NULL,
  aif = NULL,
  interpPoints = 6000
)
}
\arguments{
\item{t_blood}{Numeric vector of times for each blood measurement in minutes.}

\item{blood}{Numeric vector of the radioactivity concentration in each blood
measurement.}

\item{t_plasma}{Numeric vector of times for each plasma measurement in
minutes.}

\item{plasma}{Numeric vector of the radioactivity concentration in each
plasma measurement.}

\item{t_parentfrac}{Numeric vector of times for each parent fraction
measurement in minutes.}

\item{parentfrac}{Numeric vector of the radioactivity concentration in each
plasma measurement.}

\item{t_aif}{Optional. Numeric vector of times for each modelled AIF
value in minutes.}

\item{aif}{Optional. Numeric vector of the radioactivity concentration
for each modelled AIF value.}

\item{interpPoints}{The number of points to interpolate into.}
}
\value{
A dataframe containing the time, blood, plasma and parent fraction
  interpolated into the same times, with \code{interpPoints} number of
  points.
}
\description{
This function sorts out all the blood, plasma and parent
  fraction measurements into one convenient data frame for arterial models.
  It makes several 'editorial decisions' in the process. i) The data is
  interpolated into 6000 points by default (after some trial and error, the
  number of points can have dramatic implications for the speed of the
  function: 6000 and 1024 are very fast). ii) The different measurements are
  set to have the same times, thus if one measurement is taken for a shorter
  period than the others, the others are extended to that time point. This
  extension is performed by keeping the same value as the previous recorded
  value at that point. Again, get in touch if you want a better method
  introduced here. iii) This function sets the blood concentration, plasma
  concentration and parent fraction to 0, 0, and 1 respectively at time 0.
  Further, it removes any measurements at time <= 0. Best to add a bit to all
  measurements if you have time<=0 values, which can be fixed in the time
  shifting.
}
\details{
Function to interpolate the blood, plasma and parent fraction data in order
to create an \code{input} object for use with models requiring arterial input
functions. This function is a poor stand-in for the create_blooddata_*
functions.
}
\examples{
\dontrun{
input <- blood_interp(
  t_blood = blooddata$Time.sec. / 60,
  blood = blooddata$Cbl.nCi.cc.,
  t_plasma = plasmadata$Time.sec. / 60,
  plasma = plasmadata$Cpl.nCi.cc.,
  t_parentfrac = parentdata$Times / 60,
  parentfrac = parentdata$Fraction
)
}

}
\author{
Granville J Matheson, \email{mathesong@gmail.com}
}
