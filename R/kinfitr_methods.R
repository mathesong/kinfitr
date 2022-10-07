#' Plot kinetic model fit
#'
#' Function to plot the output of a kinetic model.
#'
#' @param x The output object of the model fitting procedure.
#' @param ... Additional optional arguments.
#'
#' @details This function uses the \code{out$model} name to call the correct function to plot the model fit.
#'
#' @examples
#' \dontrun{
#' loganout <- Loganplot(t_tac, tac, input, 10, weights)
#' plot(loganout)
#'
#' srtmout <- srtm(t_tac, reftac, roitac)
#' plot(srtmout, roiname = "Putamen", refname = "Cerebellum")
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
plot.kinfit <- function(x, ...) {
  plot_kinfit(x, ...)
}

#' Predict Blood Model Values
#'
#' Method to predict values from blood models
#'
#' @param object The output object of the model fitting procedure.
#' @param ... Additional optional arguments.
#'
#' @details This function uses the \code{class(model)} name to call the correct function to plot the model fit.
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' blood <- bd_extract(blooddata, output = "Blood")
#' blood_fit <- blmod_splines(blood$time,
#'                            blood$activity,
#'                            Method = blood$Method)
#'
#' predict(blood_fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
predict.blmod <- function(object, ...) {
  predict_blmod(object, ...)
}


#' Plot interpolated blood
#'
#' Function to plot the interpolated blood input (\code{interpblood} object).
#'
#' @param x The interpolated blood \code{input} object.
#' @param ... Additional optional arguments.
#'
#' @details This function provides an idea of how the input object looks.
#'
#' @examples
#' \dontrun{
#' input <- bd_getdata(blooddata)
#' plot(input)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
plot.interpblood <- function(x, ...) {
  plot_input(x, ...)
}


#' Plot blood data
#'
#' Function to plot the \code{blooddata} object.
#'
#' @param x The \code{blooddata} object.
#' @param ... Additional optional arguments.
#'
#' @details This function provides an idea of how the blood data and their fitted data look.
#'
#' @examples
#' \dontrun{
#' plot(blooddata)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
plot.blooddata <- function(x, ...) {
  plot_blooddata(x, ...)
}
