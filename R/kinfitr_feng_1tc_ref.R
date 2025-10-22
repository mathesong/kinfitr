#' Feng input model convolved with a 1TC IRF
#'
#' This function is intended for describing a TAC, such as the reference region
#' TAC. The parameters are not meant to be interpreted as true values, but
#' merely as a guide for interpolating the measured curve through a parametric
#' description. For this reason, there are no upper and lower limits, because
#' completely incorrect values can be helpful for interpolating our measured
#' TACs more closely.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the tissue for
#'   each frame. We include zero at time zero: if not included, it is added.
#' @param weights Optional. Numeric vector of the weights assigned to each frame
#'   in the fitting. We include zero at time zero: if not included, it is added.
#'   If not specified, uniform weights will be used.
#' @param fit_t0 Should a time zero point be fitted? If TRUE, the model can
#'   accommodate zero values in the TAC before rising. If FALSE, the TAC must
#'   rise at time = 0.
#' @param frameStartEnd Optional. This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This can be used to assess
#'   time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and
#'   end time point instead of defining the frame numbers using frameStartEnd.
#'   This function will restrict the model to all time frames whose t_tac is
#'   between the values, i.e. c(0,5) will select all frames with midtimes
#'   during the first 5 minutes.
#' @param multstart_iter Number of iterations for starting parameters. Default
#'   is 500. For more information, see
#'   \code{\link[nls.multstart]{nls_multstart}}.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par},
#'   their percentage standard errors \code{out$par.se}, the model fit object
#'   \code{out$fit}, the model weights \code{out$weights}, and a dataframe
#'   containing the TACs both of the data and the fitted values \code{out$tacs}.
#'
#' @examples
#'
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' tac <- simref$tacs[[2]]$Reference
#' weights <- simref$tacs[[2]]$Weights
#'
#' fit <- feng_1tc_tac(t_tac, tac, weights)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Jiao, J. et al, 2023. NiftyPAD-Novel Python Package for
#'   Quantitative Analysis of Dynamic PET Data. Neuroinformatics, pp.1-12.
#'   Matheson, G.J & Ogden, R.T., in preparation. SiMBA for Reference Tissue
#'   Models.
#'
#' @export
feng_1tc_tac <- function(t_tac, tac, weights = NULL,
                         fit_t0 = TRUE,
                         frameStartEnd = NULL,
                         timeStartEnd = NULL,
                         multstart_iter = 500) {

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  # Tidying

  tidyinput <- tidyinput_ref(t_tac, tac, tac, weights, frameStartEnd)

  modeldata <- tidyinput

  maxtac <- max(modeldata$reftac)
  maxtactime <- modeldata$t_tac[which.max(modeldata$TAC)]
  meantac <- mean(modeldata$reftac)


  # Parameters

  start <- c(t0 = 0.1,
             A = maxtac * 100,
             B = maxtac * 10,
             C = maxtac * 1,
             alpha = 1,
             beta=0.1,
             gamma = 0.01,
             Ph1 = 0.1,
             Th1 = 0.1)

  upper <- start * 10

  lower <- c(t0 = 0,
             A = maxtac,
             B = 0,
             C = 0,
             alpha = 0,
             beta=0,
             gamma = 0,
             Ph1 = 0.001,
             Th1 = 0.001)

  if( !fit_t0 ) {
    start <- start[-1]
    upper <- upper[-1]
    lower <- lower[-1]
  }

  multstart_pars <- fix_multstartpars(
    start, lower, upper, multstart_iter,
    NULL, NULL  )

  multstart_upper <- multstart_pars$multstart_upper
  multstart_lower <- multstart_pars$multstart_lower


  # Solution

  formula <- paste("reftac ~ feng_1tc_tac_model(t_tac, ",
                   "t0",ifelse(fit_t0, yes = "", no = "=0"),
                   ", A, B, C, alpha, beta, gamma, Ph1, Th1)",
                   sep="")

  output <- nls.multstart::nls_multstart(
    formula = as.formula(formula),
    data = modeldata,
    supp_errors = "Y",
    start_lower = multstart_lower,
    start_upper = multstart_upper,
    iter = multstart_iter,
    convergence_count = FALSE,
    #lower = lower, upper = upper,  # These are removed to allow "wrong"
                                    # parameter values with better fits,
                                    # since we don't intend to actually
                                    # interpret the parameter values.
    modelweights = weights
  )


  # Output

  tacs <- data.frame(
    Time = modeldata$t_tac,
    Reference = modeldata$reftac,
    Reference_fitted = as.numeric(fitted(output))
  )

  par <- as.data.frame(as.list(coef(output)))

  out <- list(
    par = par,
    fit = output,
    weights = modeldata$weights,
    tacs = tacs,
    model = "feng_1tc_tac"
  )

  class(out) <- c("feng_1tc_tac", "kinfit")

  return(out)
}



#' Model: Feng input model convolved with a 1TC IRF
#'
#' This is the model definition for the Feng input model convolved with a 1TC
#' IRF.
#'
#' @param t_tac Numeric vector of time values in minutes.
#' @param t0 The time point at which the curve begins to increase
#' @param A Feng AIF model A parameter
#' @param B Feng AIF model B parameter
#' @param C Feng AIF model C parameter
#' @param alpha Feng AIF model alpha parameter
#' @param beta Feng AIF model beta parameter
#' @param gamma Feng AIF model gamma parameter
#' @param Ph1 1TC IRF model K1 parameter
#' @param Th1 1TC IRF model k2 parameter
#'
#' @return The predicted values
#' @export
#'
#' @examples
#' data(simref)
#' t_tac <- simref$tacs[[2]]$Times
#'
#' feng_1tc_tac_model(t_tac, 0, 3, 1, 0.2, 0.6, 0.2, 0.01, 0.2, 0.1)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Jiao, J. et al, 2023. NiftyPAD-Novel Python Package for
#'   Quantitative Analysis of Dynamic PET Data. Neuroinformatics, pp.1-12.
#'   Matheson, G.J & Ogden, R.T., in preparation. SiMBA for Reference Tissue
#'   Models.
feng_1tc_tac_model <- function(t_tac, t0,
                                  A, B, C,
                                  alpha, beta, gamma,
                                  Ph1, Th1) {

  tcorr <- t_tac - t0

  (tcorr > 0) * (
    (Ph1*((B*(-1 + exp(tcorr*(-alpha + Th1))))/(alpha - Th1) +
            (C*(-1 + exp(tcorr*(-alpha + Th1))))/(alpha - Th1) +
            (B*(-1 + exp(tcorr*(-beta + Th1))))/(-beta + Th1) +
            (C*(-1 + exp(tcorr*(-gamma + Th1))))/(-gamma + Th1) +
            (A*(1 + exp(tcorr*(-alpha + Th1))*(-1 - alpha*tcorr + tcorr*Th1)))/
            (alpha - Th1)^(2)))/exp(tcorr*Th1)
  )
}



#' Plot Feng input model convolved with a 1TC IRF fit
#'
#' This function fits the results of the feng_1tc_tac model. For this
#' particular plotting function, I've added little crosses on the TAC
#' points in case the interpolation misses between points. The specific
#' values are important for SRTM and FRTM which use the "raw" reference
#' TAC, so extreme values at the exact times of the data points can
#' be problematic.
#'
#' @param feng_1tc_tac_fit_out The model output fit object.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' tac <- simref$tacs[[2]]$Reference
#' weights <- simref$tacs[[2]]$Weights
#'
#' fit <- feng_1tc_tac(t_tac, tac, weights)
#'
#' plot_feng_1tc_tacfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export
plot_feng_1tc_tacfit <- function(feng_1tc_tac_fit_out) {

  measured <- data.frame(
    Time = feng_1tc_tac_fit_out$tacs$Time,
    Radioactivity = feng_1tc_tac_fit_out$tacs$Reference,
    Weights = weights(feng_1tc_tac_fit_out$fit)
  )

  measured$Pred = predict(feng_1tc_tac_fit_out$fit)

  ifitted <- data.frame(
    Time = seq(0, max(feng_1tc_tac_fit_out$tacs$Time), length.out=1000)
  )

  ifitted$Radioactivity = predict(feng_1tc_tac_fit_out$fit,
                                  newdata=list(t_tac = ifitted$Time))



  outplot <- ggplot(measured, aes(x = Time, y = Radioactivity)) +
    geom_point(data = measured, aes(shape = "a", size = Weights)) +
    geom_line(data = ifitted, colour="red") +
    guides(shape = "none") +
    scale_size(range = c(1, 3)) +
    coord_cartesian(ylim = c(0, max(c(measured$Radioactivity,
                                      measured$Pred)*1.2))) +
    geom_point(data = measured, aes(y=Pred), shape = 3, size=1, colour="orange")

  return(outplot)
}
