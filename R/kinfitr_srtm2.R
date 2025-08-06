#' Simplified Reference Tissue Model 2
#'
#' Function to fit the SRTM2 model of Wu and Carson (2002) to data. Note that
#' without setting the k2prime, the model is effectively equivalent to a
#' conventional SRTM model. This configuration is allowed for estimating an
#' appropriate k2prime value, but without re-fitting with a specified k2prime
#' value, the model is not really SRTM2.
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
#' @param k2prime Optional. If empty, then the model will fit a value of
#'   k2prime. If specified, the model will be fitted with this parameter set
#'   (i.e. as a 2 parameter model).
#' @param weights Optional. Numeric vector of the weights assigned to each frame
#'   in the fitting. We include zero at time zero: if not included, it is added.
#'   If not specified, uniform weights will be used.
#' @param frameStartEnd Optional. This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes.
#' @param R1.start Optional. Starting parameter for fitting of R1. Default is 1.
#' @param R1.lower Optional. Lower bound for the fitting of R1. Default is 0.
#' @param R1.upper Optional. Upper bound for the fitting of R1. Default is 10.
#' @param k2prime.start Optional. Starting parameter for fitting of k2prime.
#'   Default is 0.1.
#' @param k2prime.lower Optional. Lower bound for the fitting of k2prime.
#'   Default is 0.001.
#' @param k2prime.upper Optional. Upper bound for the fitting of k2prime.
#'   Default is 1.
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
#' fit_setk2prime <- srtm2(t_tac, reftac, roitac, k2prime=0.1)
#'
#' # Note: this is not really SRTM2 because the k2prime is not specified
#' fitk2prime <- srtm2(t_tac, reftac, roitac)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Wu Y, Carson RE. Noise reduction in the simplified reference
#'   tissue model for neuroreceptor functional imaging. J Cereb Blood Flow
#'   Metab. 2002;22:1440-1452.
#'
#'
#'
#' @export

srtm2 <- function(t_tac, reftac, roitac, k2prime=NULL, weights = NULL, frameStartEnd = NULL, timeStartEnd = NULL,
                 R1.start = 1, R1.lower = 0, R1.upper = 10,
                 k2prime.start = 0.1, k2prime.lower = 0.001, k2prime.upper = 1,
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

  if(is.null(k2prime)) {
    start <- c(R1 = R1.start, k2prime = k2prime.start, bp = bp.start)
    lower <- c(R1 = R1.lower, k2prime = k2prime.lower, bp = bp.lower)
    upper <- c(R1 = R1.upper, k2prime = k2prime.upper, bp = bp.upper)

    rlang::inform("Note: Without specifying a k2prime value for SRTM2, it is effectively
            equivalent to the conventional SRTM model. This can be useful for
            selecting an appropriate k2prime value, but without re-fitting the
            model with a specified k2prime value, the model is not really SRTM2.",
            .frequency = "once", .frequency_id = "srtm2_message")

  } else {

    if(length(k2prime) > 1) {
      stop("k2prime must be specified by a single value.")
    }

    start <- c(R1 = R1.start, bp = bp.start)
    lower <- c(R1 = R1.lower, bp = bp.lower)
    upper <- c(R1 = R1.upper, bp = bp.upper)
  }

  multstart_pars <- fix_multstartpars(
    start, lower, upper, multstart_iter,
    multstart_lower, multstart_upper
  )
  multstart_upper <- multstart_pars$multstart_upper
  multstart_lower <- multstart_pars$multstart_lower


  # Solution

  formula <- paste0("roitac ~ srtm2_model(t_tac, reftac, R1, k2prime",
                    ifelse(is.null(k2prime),
                           yes = "",
                           no = paste0("=", k2prime)),
                    ", bp)")

  if (prod(multstart_iter) == 1) {
    output <- minpack.lm::nlsLM(
      formula = as.formula(formula),
      data = modeldata,
      start = start, lower = lower, upper = upper,
      weights = weights,
      control = minpack.lm::nls.lm.control(maxiter = 200),
      trace = printvals
    )
  } else {
    output <- nls.multstart::nls_multstart(
      formula = as.formula(formula),
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


  # Coefficients
  par <- as.data.frame(as.list(coef(output)))

  par.se <- par
  par.se[1,] <- purrr::map_dbl(names(coef(output)), ~ get_se(output, .x))
  names(par.se) <- paste0(names(par.se), ".se")

  if(!is.null(k2prime)) {
    par$k2prime = k2prime
    par.se$k2prime=0
  }

  # Derived
  par$k2a = with(par, (R1 * k2prime) / (bp + 1)  )

  par.se$k2a.se <- get_se(output, "(R1 * k2prime) / (bp + 1)")

  out <- list(
    par = par, par.se = par.se,
    fit = output, weights = modeldata$weights, tacs = tacs,
    model = "srtm2"
  )

  class(out) <- c("srtm2", "kinfit")

  return(out)
}


#' Model: Simplified Reference Tissue Model 2
#'
#' This is the SRTM2 model itself by which predicted values are generated.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a zero.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param R1 Parameter value for R1
#' @param k2prime Parameter value for k2prime
#' @param bp Parameter value for bp
#'
#' @return A numeric vector of the predicted values of the TAC in the target region.
#'
#' @examples
#' \dontrun{
#' srtm2_model(t_tac, reftac, R1 = 0.9, k2prime = 0.1, bp = 0.1)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Wu Y, Carson RE. Noise reduction in the simplified reference tissue
#'   model for neuroreceptor functional imaging. J Cereb Blood Flow Metab.
#'   2002;22:1440-1452.
#'
#' @export


srtm2_model <- function(t_tac, reftac, R1, k2prime, bp) {

  interptime <- pracma::linspace(min(t_tac), max(t_tac), 1024)
  step <- interptime[2] - interptime[1]

  iref <- pracma::interp1(t_tac, reftac, interptime, method = "linear")

  k2a <- (R1 * k2prime) / (bp + 1)

  a <- R1 * (k2prime - k2a) * iref
  b <- exp(-k2a * interptime)

  ND <- R1 * iref
  BOUND <- kinfit_convolve(a, b, step)
  i_outtac <- ND + BOUND

  outtac <- pracma::interp1(interptime, i_outtac, t_tac)

  return(outtac)
}


#' Plot: Simplified Reference Tissue Model 2
#'
#' Function to visualise the fit of the SRTM2 model to data.
#'
#' @param srtm2out The output object of the SRTM2 fitting procedure.
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
#' fit <- srtm2(t_tac, reftac, roitac, weights=weights)
#'
#' plot_srtm2fit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_srtm2fit <- function(srtm2out, roiname = NULL, refname = NULL) {

  measured <- data.frame(
    Time = srtm2out$tacs$Time,
    Reference = srtm2out$tacs$Reference,
    ROI.measured = srtm2out$tacs$Target,
    Weights = weights(srtm2out$fit)
  )

  fitted <- data.frame(
    Time = srtm2out$tacs$Time,
    ROI.fitted = srtm2out$tacs$Target_fitted,
    Weights = weights(srtm2out$fit)
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
