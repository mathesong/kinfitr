#' Spline-based Reference TAC Fitting
#'
#' This function fits a reference TAC using a piecewise approach: NLS to optimize
#' t0 (the time when the curve starts rising from zero), and GAM to fit the curve
#' after t0. The parameters are not meant to be interpreted as true physiological
#' values, but merely as a guide for denoising.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the tissue for
#'   each frame. We include zero at time zero: if not included, it is added.
#' @param weights Optional. Numeric vector of the weights assigned to each frame
#'   in the GAM fitting (not used in NLS optimisation). We include zero at time
#'   zero: if not included, it is added. If not specified, uniform weights will
#'   be used.
#' @param frameStartEnd Optional. This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This can be used to assess
#'   time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and
#'   end time point instead of defining the frame numbers using frameStartEnd.
#'   This function will restrict the model to all time frames whose t_tac is
#'   between the values, i.e. c(0,5) will select all frames with midtimes
#'   during the first 5 minutes.
#' @param t0.lower Lower bound for the fitting of t0. Default is 0.
#' @param t0.upper Upper bound for the fitting of t0. Default is 5.
#' @param k Basis dimension for the smooth term, controlling wiggliness. Higher
#'   values allow more flexible (wiggly) fits. Default is -1, which allows mgcv
#'   to choose the dimension automatically. Typical values range from 5 (smooth)
#'   to 15 (flexible).
#'
#' @return A list with a data frame of the fitted parameters \code{out$par},
#'   the optimization fit object \code{out$fit}, the final GAM fit object
#'   \code{out$gam_fit}, the model weights \code{out$weights}, and a dataframe
#'   containing the TACs both of the data and the fitted values \code{out$tacs}.
#'
#' @details The model fits the curve in two pieces:
#'   \itemize{
#'     \item Before t0: TAC = 0
#'     \item After t0: TAC ~ s(sqrt(t_tac - t0), k=k) using mgcv::gam with REML
#'   }
#'   Any predicted values below zero are set to 0. The k parameter controls the
#'   basis dimension of the smooth, which affects how wiggly the fitted curve can
#'   be. Lower values produce smoother fits, higher values allow more flexibility.
#'
#' @examples
#'
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' tac <- simref$tacs[[2]]$Reference
#' weights <- simref$tacs[[2]]$Weights
#'
#' fit <- spline_tac(t_tac, tac, weights)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
spline_tac <- function(t_tac, tac, weights = NULL,
                       frameStartEnd = NULL,
                       timeStartEnd = NULL,
                       t0.lower = 0,
                       t0.upper = 5,
                       k = -1) {

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  # Tidying - use tac for both reftac and roitac
  tidyinput <- tidyinput_ref(t_tac, tac, tac, weights, frameStartEnd)

  modeldata <- tidyinput

  # Store original data for closure
  full_t_tac <- modeldata$t_tac
  full_tac <- modeldata$roitac
  full_weights <- modeldata$weights

  # Define objective function for optimize()
  objective_fn <- function(t0_val) {

    predictions <- rep(0, length(full_t_tac))
    after_t0_mask <- full_t_tac >= t0_val

    # Need at least 3 points to fit a spline
    if(sum(after_t0_mask) < 3) {
      return(1e10)
    }

    t_adj <- full_t_tac[after_t0_mask] - t0_val
    tac_subset <- full_tac[after_t0_mask]
    weights_subset <- full_weights[after_t0_mask]

    gam_data <- data.frame(
      tac = tac_subset,
      t_adj = t_adj
    )

    # Fit GAM with weights
    gam_fit <- try(mgcv::gam(
      tac ~ s(sqrt(t_adj), k = k),
      data = gam_data,
      method = "REML",
      weights = weights_subset
    ), silent = TRUE)

    if(inherits(gam_fit, "try-error")) {
      return(1e10)
    }

    # Predict for all time points
    predictions[after_t0_mask] <- predict(gam_fit, newdata = gam_data, type = "response")

    # Set negative predictions to 0
    predictions[predictions < 0] <- 0

    # Calculate residual sum of squares
    residuals <- full_tac - predictions
    rss <- sum(residuals^2)

    return(rss)
  }

  # Optimize t0 using base R optimize()
  opt_result <- optimize(
    f = objective_fn,
    interval = c(t0.lower, t0.upper),
    tol = 0.0001
  )

  # Extract fitted t0
  fitted_t0 <- opt_result$minimum

  # Refit final GAM with optimal t0 for output
  after_t0_mask <- modeldata$t_tac >= fitted_t0

  if(sum(after_t0_mask) >= 3) {
    t_adj <- modeldata$t_tac[after_t0_mask] - fitted_t0
    tac_subset <- modeldata$roitac[after_t0_mask]
    weights_subset <- modeldata$weights[after_t0_mask]

    gam_data <- data.frame(
      tac = tac_subset,
      t_adj = t_adj
    )

    gam_fit <- mgcv::gam(
      tac ~ s(sqrt(t_adj), k = k),
      data = gam_data,
      method = "REML",
      weights = weights_subset
    )
  } else {
    gam_fit <- NULL
  }

  # Calculate fitted values using the optimal t0
  fitted_vals <- rep(0, nrow(modeldata))
  if(!is.null(gam_fit) && sum(after_t0_mask) > 0) {
    gam_data_final <- data.frame(
      t_adj = modeldata$t_tac[after_t0_mask] - fitted_t0
    )
    fitted_vals[after_t0_mask] <- predict(gam_fit, newdata = gam_data_final, type = "response")
  }
  fitted_vals[fitted_vals < 0] <- 0

  # Create a pseudo-fit object for compatibility with plotting functions
  fit_obj <- list(
    par = c(t0 = fitted_t0),
    value = opt_result$objective,
    convergence = 0,
    fitted_values = fitted_vals
  )
  class(fit_obj) <- "spline_tac_fit"

  # Output
  tacs <- data.frame(
    Time = modeldata$t_tac,
    TAC = modeldata$roitac,
    TAC_fitted = fitted_vals
  )

  par <- data.frame(t0 = fitted_t0)

  out <- list(
    par = par,
    fit = fit_obj,
    gam_fit = gam_fit,
    weights = modeldata$weights,
    tacs = tacs,
    model = "spline_tac"
  )

  class(out) <- c("spline_tac", "kinfit")

  return(out)
}



#' Predict method for spline_tac
#'
#' Generate predictions from a fitted spline_tac model at specified time points.
#'
#' @param object A spline_tac fit object.
#' @param newdata Optional. A list or data frame with variable \code{t_tac}
#'   containing new time points at which to predict. If NULL, predictions
#'   are generated at the original time points.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of predicted TAC values.
#'
#' @examples
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' tac <- simref$tacs[[2]]$Reference
#' weights <- simref$tacs[[2]]$Weights
#'
#' fit <- spline_tac(t_tac, tac, weights)
#'
#' # Predict at original times
#' pred1 <- predict(fit)
#'
#' # Predict at new times
#' new_times <- seq(0, max(t_tac), length.out = 100)
#' pred2 <- predict(fit, newdata = list(t_tac = new_times))
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
predict.spline_tac <- function(object, newdata = NULL, ...) {

  fitted_t0 <- object$par$t0

  # Default: predict at original time points
  if(is.null(newdata)) {
    t_tac <- object$tacs$Time
  } else {
    t_tac <- newdata$t_tac
  }

  # Initialize predictions to zero
  predictions <- rep(0, length(t_tac))
  after_t0_mask <- t_tac >= fitted_t0

  # Predict from GAM for times after t0
  if(sum(after_t0_mask) > 0 && !is.null(object$gam_fit)) {
    t_adj <- t_tac[after_t0_mask] - fitted_t0
    newdata_gam <- data.frame(t_adj = t_adj)

    predictions[after_t0_mask] <- predict(
      object$gam_fit,
      newdata = newdata_gam,
      type = "response"
    )
  }

  # Set negative predictions to 0
  predictions[predictions < 0] <- 0

  return(predictions)
}



#' Plot Spline TAC fit
#'
#' This function plots the results of the spline_tac model with interpolated
#' fitted curve. Small crosses are shown on the measured TAC points to highlight
#' the exact measured values, which is important since the interpolated curve
#' may miss between-point extrema.
#'
#' @param spline_tac_fit_out The model output fit object.
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
#' fit <- spline_tac(t_tac, tac, weights)
#'
#' plot_spline_tacfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export
plot_spline_tacfit <- function(spline_tac_fit_out) {

  measured <- data.frame(
    Time = spline_tac_fit_out$tacs$Time,
    Radioactivity = spline_tac_fit_out$tacs$TAC,
    Weights = spline_tac_fit_out$weights
  )

  measured$Pred = spline_tac_fit_out$tacs$TAC_fitted

  # Create interpolated fitted curve
  ifitted <- data.frame(
    Time = seq(0, max(spline_tac_fit_out$tacs$Time), length.out=1000)
  )

  ifitted$Radioactivity = predict(spline_tac_fit_out,
                                  newdata = list(t_tac = ifitted$Time))

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
