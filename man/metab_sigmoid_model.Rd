% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kinfitr_pfmodels.R
\name{metab_sigmoid_model}
\alias{metab_sigmoid_model}
\title{Sigmoid Model for Parent Fraction}
\usage{
metab_sigmoid_model(time, a, b, c, ppf0 = 1, delay = 0)
}
\arguments{
\item{time}{Time in seconds.}

\item{a}{Parameter A.}

\item{b}{Parameter B.}

\item{c}{Parameter C.}

\item{ppf0}{The starting point of the parent fraction curve.}

\item{delay}{The delay of the metabolism curve.}
}
\value{
Predicted values
}
\description{
This is the model function for fitting of the sigmoid function by Guo et al. (2013).
}
\examples{
metab_sigmoid_model(seq(0, 60 * 60, by = 120), 7, 0.6, 0.04, 1, 0)
}
\references{
Guo Q, Colasanti A, Owen DR, et al. Quantification of the specific translocator protein signal of 18F-PBR111 in healthy humans: a genetic polymorphism effect on in vivo binding. J Nucl Med 2013; 54: 1915–1923.
}
