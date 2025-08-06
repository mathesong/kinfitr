#' Simplified Reference Tissue Model
#'
#' Function to fit the SRTM model of Lammertsma and Hume (1996) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param roitac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param weights Optional. Numeric vector of the weights assigned to each frame
#'   in the fitting. We include zero at time zero: if not included, it is added.
#'   If not specified, uniform weights will be used.
#' @param frameStartEnd Optional. This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes.
#' @param R1.start Optional. Starting parameter for fitting of R1. Default is 1.
#' @param R1.lower Optional. Lower bound for the fitting of R1. Default is 0.
#' @param R1.upper Optional. Upper bound for the fitting of R1. Default is 10.
#' @param k2.start Optional. Starting parameter for fitting of k2. Default is
#'   0.1.
#' @param k2.lower Optional. Lower bound for the fitting of k2. Default is 0.
#' @param k2.upper Optional. Upper bound for the fitting of k2. Default is 1.
#' @param bp.start Optional. Starting parameter for fitting of bp. Default is
#'   1.5.
#' @param bp.lower Optional. Lower bound for the fitting of bp. Default is -10.
#' @param bp.upper Optional. Upper bound for the fitting of bp. Default is 15.
#' @param multstart_iter Number of iterations for starting parameters. Default
#'   is 1. For more information, see \code{\link[nls.multstart]{nls_multstart}}.
#'   If specified as 1 for any parameters, the original starting value will be
#'   used, and the multstart_lower and multstart_upper values ignored.
#' @param multstart_lower Optional. Lower bounds for starting parameters.
#'   Defaults to the lower bounds.  Named list of whichever parameters' starting
#'   bounds should be altered.
#' @param multstart_upper Optional. Upper bounds for starting parameters.
#'   Defaults to the upper bounds.  Named list of whichever parameters' starting
#'   bounds should be altered.
#' @param printvals Optional. This displays the parameter values for each
#'   iteration of the model. This is useful for debugging and changing starting
#'   values and upper and lower bounds for parameters.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par},
#'   their percentage standard errors (scaled so that 1 represents 100\%)
#'   \code{out$par.se}, the model fit object \code{out$fit}, the model weights
#'   \code{out$weights}, and a dataframe containing the TACs both of the data
#'   and the fitted values \code{out$tacs}.
#'
#' @examples
#'
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' reftac <- simref$tacs[[2]]$Reference
#' roitac <- simref$tacs[[2]]$ROI1
#' weights <- simref$tacs[[2]]$Weights
#'
#' fit1 <- srtm(t_tac, reftac, roitac)
#' fit2 <- srtm(t_tac, reftac, roitac, weights, frameStartEnd = c(1, 33))
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Lammertsma AA, Hume SP. Simplified reference tissue model for PET
#'   receptor studies. Neuroimage. 1996 Dec 31;4(3):153-8.
#'
#' @export

srtm <- function(t_tac, reftac, roitac, weights = NULL, frameStartEnd = NULL, timeStartEnd = NULL,
                 R1.start = 1, R1.lower = 0, R1.upper = 10,
                 k2.start = 0.1, k2.lower = 0, k2.upper = 1,
                 bp.start = 1.5, bp.lower = 0, bp.upper = 15,
                 multstart_iter = 1, multstart_lower = NULL, multstart_upper = NULL,
                 printvals = F) {

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1], 
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  # Tidying

  tidyinput <- tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd)

  modeldata <- tidyinput


  # Parameters

  start <- c(R1 = R1.start, k2 = k2.start, bp = bp.start)
  lower <- c(R1 = R1.lower, k2 = k2.lower, bp = bp.lower)
  upper <- c(R1 = R1.upper, k2 = k2.upper, bp = bp.upper)

  multstart_pars <- fix_multstartpars(
    start, lower, upper, multstart_iter,
    multstart_lower, multstart_upper
  )
  multstart_upper <- multstart_pars$multstart_upper
  multstart_lower <- multstart_pars$multstart_lower


  # Solution

  if (prod(multstart_iter) == 1) {
    output <- minpack.lm::nlsLM(
      roitac ~ srtm_model(t_tac, reftac, R1, k2, bp),
      data = modeldata,
      start = start, lower = lower, upper = upper,
      weights = weights,
      control = minpack.lm::nls.lm.control(maxiter = 200),
      trace = printvals
    )
  } else {
    output <- nls.multstart::nls_multstart(
      roitac ~ srtm_model(t_tac, reftac, R1, k2, bp),
      data = modeldata,
      supp_errors = "Y",
      start_lower = multstart_lower,
      start_upper = multstart_upper,
      iter = multstart_iter, convergence_count = FALSE,
      lower = lower, upper = upper, modelweights = weights
    )
  }

  # Check for parameters hitting limits

  limcheck_u <- purrr::map2_lgl(round(upper,3), round(coef(output),3), identical)
  limcheck_l <- purrr::map2_lgl(round(lower,3), round(coef(output),3), identical)
  limcheck <- limcheck_u + limcheck_l
  limcheck <- limcheck==1

  if(
    any(limcheck)
  ) {
    warning(
      paste0(
        "\nFitted parameters are hitting upper or lower limit bounds. Consider \n",
        "either modifying the upper and lower limit boundaries, or else using \n",
        "multstart when fitting the model (see the function documentation).\n") )
  }

  # Output

  tacs <- data.frame(
    Time = modeldata$t_tac,
    Reference = modeldata$reftac,
    Target = modeldata$roitac,
    Target_fitted = as.numeric(fitted(output))
  )

  par <- as.data.frame(as.list(coef(output)))

  par.se <- par
  par.se[1,] <- purrr::map_dbl(names(coef(output)), ~ get_se(output, .x))
  names(par.se) <- paste0(names(par.se), ".se")

  out <- list(
    par = par, par.se = par.se,
    fit = output, weights = modeldata$weights, tacs = tacs,
    model = "srtm"
  )

  class(out) <- c("srtm", "kinfit")

  return(out)
}


#' Model: Simplified Reference Tissue Model
#'
#' This is the SRTM model itself by which predicted values are generated.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a zero.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param R1 Parameter value for R1
#' @param k2 Parameter value for k2
#' @param bp Parameter value for bp
#'
#' @return A numeric vector of the predicted values of the TAC in the target region.
#'
#' @examples
#' \dontrun{
#' srtm_model(t_tac, reftac, R1 = 0.9, k2 = 0.1, bp = 1.5)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Lammertsma AA, Hume SP. Simplified reference tissue model for PET receptor studies. Neuroimage. 1996 Dec 31;4(3):153-8.
#'
#' @export


srtm_model <- function(t_tac, reftac, R1, k2, bp) {
  interptime <- pracma::linspace(min(t_tac), max(t_tac), 1024)
  step <- interptime[2] - interptime[1]

  iref <- pracma::interp1(t_tac, reftac, interptime, method = "linear")

  a <- (k2 - (R1 * k2 / (1 + bp))) * iref
  b <- exp((-k2 / (1 + bp)) * interptime)

  ND <- R1 * iref
  BOUND <- kinfit_convolve(a, b, step)
  i_outtac <- ND + BOUND

  outtac <- pracma::interp1(interptime, i_outtac, t_tac)

  return(outtac)
}


#' Plot: Simplified Reference Tissue Model
#'
#' Function to visualise the fit of the SRTM model to data.
#'
#' @param srtmout The output object of the SRTM fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#' @param refname Optional. The name of the Reference Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' reftac <- simref$tacs[[2]]$Reference
#' roitac <- simref$tacs[[2]]$ROI1
#' weights <- simref$tacs[[2]]$Weights
#'
#' fit <- srtm(t_tac, reftac, roitac, weights, frameStartEnd = c(1, 33))
#'
#' plot_srtmfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_srtmfit <- function(srtmout, roiname = NULL, refname = NULL) {
  measured <- data.frame(
    Time = srtmout$tacs$Time,
    Reference = srtmout$tacs$Reference,
    ROI.measured = srtmout$tacs$Target,
    Weights = weights(srtmout$fit)
  )

  fitted <- data.frame(
    Time = srtmout$tacs$Time,
    ROI.fitted = srtmout$tacs$Target_fitted,
    Weights = weights(srtmout$fit)
  )

  if (is.null(roiname)) {
    roiname <- "ROI"
  }
  if (is.null(refname)) {
    refname <- "Reference"
  }

  measured <- plyr::rename(measured, c(
    "ROI.measured" = paste0(roiname, ".measured"),
    "Reference" = refname
  ))

  fitted <- plyr::rename(fitted, c("ROI.fitted" = paste0(roiname, ".fitted")))

  tidymeasured <- tidyr::gather(
    measured,
    key = Region, value = Radioactivity,
    -Time, -Weights, factor_key = F
  )

  tidyfitted <- tidyr::gather(
    fitted,
    key = Region, value = Radioactivity,
    -Time, -Weights, factor_key = F
  )


  Region <- forcats::fct_inorder(factor(c(tidymeasured$Region, tidyfitted$Region)))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  outplot <- ggplot(tidymeasured, aes(x = Time, y = Radioactivity, colour = Region)) +
    geom_point(data = tidymeasured, aes(shape = "a", size = Weights)) +
    geom_line(data = tidyfitted) +
    guides(shape = "none", color = guide_legend(order = 1)) + colScale +
    scale_size(range = c(1, 3)) +
    coord_cartesian(ylim = c(0, max(tidymeasured$Radioactivity)))

  return(outplot)
}
