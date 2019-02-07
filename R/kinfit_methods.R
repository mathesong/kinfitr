#' Plot kinetic model fit
#'
#' Function to plot the output of a kinetic model.
#'
#' @param x The output object of the model fitting procedure.
#' @param roiname The name of the ROI.
#' @param refname The name of the reference region.
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
#' plot(srtmout, roiname = 'Putamen', refname = 'Cerebellum')
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
#'
#' @examples
plot.kinfit <- function(x, ...) {
  plot_kinfit(x, ...)
}


