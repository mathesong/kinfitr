#' Likelihood Estimation in Graphical Analysis (LEGA)
#'
#' Function to fit the LEGA approach of Ogden (2003) to
#' data. LEGA models the original (untransformed) tissue time activity curve
#' directly under a Gaussian likelihood, and as a result is (asymptotically)
#' free of the measurement bias of the Logan plot.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it  will be added.
#' @param dur Numeric vector of the time durations of the frames in minutes. The
#' LEGA estimator is explicitly frame-aware (it integrates over each acquisition
#' interval), so frame durations must be provided.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param tstar The t* specification for fitting. If tstar_type="frames",
#' this is the number of frames from the end to include (e.g., 10 means last 10 frames).
#' If tstar_type="time", this is the time point (in minutes) after which all frames
#' with midpoints later than this time are included. This value can be estimated using \code{LEGA_tstar}.
#' @param tstar_type Either "frames" (default) or "time", specifying how to interpret tstar.
#' @param tstarIncludedFrames Deprecated. Use 'tstar' with 'tstar_type="frames"' instead.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param inpshift Optional. The number of minutes by which to shift the timing of the input data frame forwards or backwards.
#' If not specified, this will be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will be ignored and assumed to be 0%. If specified, it
#' will be corrected for prior to parameter estimation using the following equation:
#' \deqn{C_{T}(t) = \frac{C_{Measured}(t) - vB\times C_{B}(t)}{1-vB}}
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes.
#' @param Vt.start Optional. Reference value for \eqn{V_T}. If not supplied, it is
#' taken from the reference fit (see \code{fallback}). This value is used as the
#' reference against which the stability fallback is judged, and as the value
#' returned if the fallback is triggered.
#' @param gamma.start Optional. Starting value for the intercept (\eqn{\gamma})
#' around which the one-dimensional optimisation is centred. If not supplied, it
#' is taken from the reference fit (see \code{fallback}).
#' @param fallback The LEGA estimator can occasionally be numerically unstable in
#' high-noise data. Following the recommendation of Ogden (2003), if the LEGA
#' estimate of \eqn{V_T} is non-finite, non-positive, or more than twice the size
#' of the reference estimate (\code{Vt.start}, or the reference-method estimate),
#' a more numerically stable estimate is returned instead (with a message). One
#' of: \code{"MA1"} (default; revert to the MA1 estimate, which targets
#' essentially the same quantity as LEGA), \code{"Logan"} (revert to the ordinary
#' Logan estimate, as in Ogden 2003), or \code{"none"} (never revert; the raw
#' LEGA estimate is returned, with a warning if it looks unreliable).
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the
#' standard error of \eqn{V_T} as a fraction of the estimate \code{out$par.se},
#' the model fit object \code{out$fit}, a dataframe containing the TACs of the
#' data \code{out$tacs}, a dataframe containing the fitted values
#' \code{out$fitvals}, the blood input data frame after time shifting
#' \code{input}, a vector of the weights \code{out$weights}, the inpshift value
#' used \code{inpshift}, the specified vB value \code{out$vB} and the specified
#' tstarIncludedFrames value \code{out$tstarIncludedFrames}.
#'
#' @details
#' LEGA is fitted using the efficient profiled algorithm of Ogden (2003,
#' equations 11-13): for a fixed value of the intercept (\eqn{\gamma}) the
#' noise-free tissue values are linear in the slope (\eqn{\beta = V_T}), so the
#' slope is available in closed form, and only the one-dimensional intercept is
#' optimised numerically. At the optimum the estimate coincides with the
#' weighted least-squares fit of the original tissue values onto the
#' model-derived terms, which is the \code{lm} object returned in \code{out$fit}.
#'
#' Note that the standard error returned here is conditional on the estimated
#' intercept. This is a reasonable preliminary approximation; the full
#' (unconditional) non-linear least squares standard errors described by Ogden
#' (2003) could be added in future.
#'
#' @examples
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' tac <- pbr28$tacs[[2]]$FC
#' weights <- pbr28$tacs[[2]]$Weights
#' dur <- pbr28$tacs[[2]]$Duration/60
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1
#' )
#'
#' fit1 <- LEGA(t_tac, dur, tac, input, 10, weights)
#' fit2 <- LEGA(t_tac, dur, tac, input, 10, weights, inpshift = 0.1, vB = 0.05)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Ogden RT. Estimation of kinetic parameters in graphical analysis of PET imaging data. Statistics in Medicine. 2003 Nov 30;22(22):3557-68.
#'
#' @export

LEGA <- function(t_tac, dur, tac, input, tstar, weights = NULL,
                 inpshift = 0, vB = 0, frameStartEnd = NULL, timeStartEnd = NULL,
                 tstar_type = "frames", tstarIncludedFrames = NULL,
                 Vt.start = NULL, gamma.start = NULL,
                 fallback = c("MA1", "Logan", "none")) {

  if (missing(dur) || is.null(dur)) {
    stop("'dur' (frame durations) is a required input for LEGA.")
  }

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  # Handle deprecated parameter
  if (!is.null(tstarIncludedFrames)) {
    warning("tstarIncludedFrames is deprecated and will be removed in a future version. Use 'tstar' with 'tstar_type=\"frames\"' instead", call. = FALSE)
    if (!missing(tstar)) {
      stop("Cannot specify both 'tstar' and 'tstarIncludedFrames'")
    }
    tstar <- tstarIncludedFrames
    tstar_type <- "frames"
  }

  # Validate tstar_type and fallback method
  tstar_type <- match.arg(tstar_type, c("frames", "time"))
  fallback <- match.arg(fallback, c("MA1", "Logan", "none"))

  # Reference fit. This serves two purposes: it provides default starting values
  # for the intercept (gamma) and the reference slope (Vt), and it supplies the
  # numerically stable estimate that the stability fallback reverts to. When the
  # fallback is disabled a reference is only computed if it is still needed to
  # seed the optimisation (i.e. gamma.start was not supplied), and is never used
  # to override the result. If the user supplies everything that would otherwise
  # be derived, no reference fit is computed at all.
  ref_method <- if (fallback != "none") fallback else "Logan"
  need_ref <- is.null(gamma.start) ||
    (fallback != "none" && is.null(Vt.start))

  if (need_ref) {
    if (ref_method == "MA1") {
      reffit <- ma1(t_tac, tac, input, tstar = tstar, weights = weights,
                    inpshift = inpshift, vB = vB, dur = dur,
                    frameStartEnd = frameStartEnd, tstar_type = tstar_type)
      Vt_ref_default <- reffit$par$Vt
      gamma_ref_default <- as.numeric(1 / stats::coef(reffit$fit)[2])
      ref_se <- reffit$par.se$Vt.se
    } else {
      reffit <- Loganplot(t_tac, tac, input, tstar = tstar, weights = weights,
                          inpshift = inpshift, vB = vB, dur = dur,
                          frameStartEnd = frameStartEnd, tstar_type = tstar_type)
      Vt_ref_default <- reffit$par$Vt
      gamma_ref_default <- as.numeric(stats::coef(reffit$fit)[1])
      ref_se <- reffit$par.se$Vt.se
    }
  } else {
    Vt_ref_default <- NA_real_
    gamma_ref_default <- NA_real_
    ref_se <- NA_real_
  }

  gamma_start <- if (!is.null(gamma.start)) gamma.start else gamma_ref_default
  Vt_ref <- if (!is.null(Vt.start)) Vt.start else Vt_ref_default

  # Tidying

  tidyinput <- tidyinput_art(t_tac, tac, weights, frameStartEnd)

  tidyinput_dur <- tidyinput_art(dur, tac, weights, frameStartEnd)
  dur <- tidyinput_dur$t_tac

  t_tac <- tidyinput$t_tac
  tac <- tidyinput$tac
  weights <- tidyinput$weights


  # Convert tstar based on type
  if (tstar_type == "time") {
    frames_after_tstar <- which(t_tac >= tstar)
    tstarIncludedFrames <- length(frames_after_tstar)
  } else {
    tstarIncludedFrames <- tstar
  }


  newvals <- shift_timings(
    t_tac = t_tac,
    tac = tac,
    input = input,
    inpshift = inpshift
  )

  t_tac <- newvals$t_tac
  tac <- newvals$tac

  t_inp <- newvals$input$Time
  blood <- newvals$input$Blood
  aif <- newvals$input$AIF

  # Parameters

  interptime <- newvals$input$Time

  i_tac <- pracma::interp1(t_tac, tac, interptime, method = "linear")

  # Blood Volume Correction (nothing happens if vB = 0)
  i_tac <- (i_tac - vB * blood) / (1 - vB)
  tac_uncor <- tac
  tac <- pracma::interp1(interptime, i_tac, t_tac, method = "linear")

  # Integrated plasma up to each frame midtime (x_i in Ogden 2003)
  x <- as.numeric(pracma::cumtrapz(interptime, aif))
  x <- pracma::interp1(interptime, x, t_tac, method = "linear")

  # Z_i is the (blood-volume-corrected) measured tissue concentration
  Z <- tac
  n <- length(Z)

  # First equilibrium frame. The added time-zero frame (dur = 0) guarantees that
  # at least one pre-equilibrium frame exists for conditioning.
  k <- n - tstarIncludedFrames + 1
  if (k < 2) k <- 2
  eq_idx <- k:n

  d_eq <- dur[eq_idx]
  x_eq <- x[eq_idx]
  Z_eq <- Z[eq_idx]
  w_eq <- weights[eq_idx]
  m <- length(eq_idx)

  # W = sum of Z_j * d_j over pre-equilibrium frames (observed values, j < k)
  W <- sum(Z[1:(k - 1)] * dur[1:(k - 1)])

  # For a fixed gamma, compute the A_i and B_i coefficients of the recursion
  # Z*_i = A_i + B_i * beta  (Ogden 2003, equations 11 and 12)
  lega_AB <- function(gamma) {
    A <- numeric(m)
    B <- numeric(m)
    denom <- gamma - (3 / 8) * d_eq

    A[1] <- (W + (1 / 8) * d_eq[1] * Z[k - 1]) / denom[1]
    B[1] <- (-x_eq[1]) / denom[1]

    SA <- d_eq[1] * A[1] # running sum_{j=k}^{i-1} d_j A_j
    SB <- d_eq[1] * B[1]

    if (m >= 2) {
      for (p in 2:m) {
        A[p] <- (W + SA + (1 / 8) * d_eq[p] * A[p - 1]) / denom[p]
        B[p] <- (SB + (1 / 8) * d_eq[p] * B[p - 1] - x_eq[p]) / denom[p]
        SA <- SA + d_eq[p] * A[p]
        SB <- SB + d_eq[p] * B[p]
      }
    }

    list(A = A, B = B)
  }

  # Closed-form (weighted) slope for a given gamma (Ogden 2003, equation 13)
  lega_beta <- function(AB) {
    sum(w_eq * AB$B * (Z_eq - AB$A)) / sum(w_eq * AB$B^2)
  }

  # One-dimensional profiled objective in gamma
  lega_sse <- function(gamma) {
    AB <- lega_AB(gamma)
    beta <- lega_beta(AB)
    fitted <- AB$A + AB$B * beta
    sum(w_eq * (Z_eq - fitted)^2)
  }

  # Search interval for gamma, centred on the starting intercept (supplied via
  # gamma.start, or the OLS Logan intercept by default). The upper bound is kept
  # strictly below (3/8) * min(dur) so the denominator (gamma - 3/8 * d_i) never
  # changes sign within the search.
  gamma_span <- max(abs(gamma_start), 1) * 5
  g_lower <- gamma_start - gamma_span
  g_upper <- min(gamma_start + gamma_span, (3 / 8) * min(d_eq) - 1e-6)
  if (g_upper <= g_lower) {
    g_lower <- g_upper - gamma_span
  }

  opt <- stats::optimize(lega_sse, interval = c(g_lower, g_upper))
  gamma_hat <- opt$minimum

  # Solution: at the optimal gamma the estimate equals the WLS fit of Z onto B
  # with offset A, which gives an lm object compatible with the rest of kinfitr
  AB <- lega_AB(gamma_hat)
  A_eq <- AB$A
  B_eq <- AB$B

  lega_model <- lm(Z_eq ~ B_eq - 1, offset = A_eq, weights = w_eq)

  Vt <- as.numeric(stats::coef(lega_model)["B_eq"])

  par.se <- data.frame(Vt.se = get_se(lega_model, "B_eq"))

  # Numerical-stability check (Ogden 2003, p.3564): the LEGA estimate is
  # considered unreliable if it is non-finite, non-positive, or more than twice
  # the reference estimate (Vt.start, or the reference-method estimate).
  unreliable <- !is.finite(Vt) || Vt <= 0 ||
    (is.finite(Vt_ref) && Vt_ref > 0 && Vt > 2 * Vt_ref)

  fallback_used <- FALSE
  fallback_method <- NA_character_
  if (unreliable) {
    if (fallback != "none" && is.finite(Vt_ref)) {
      message("LEGA estimate of Vt was non-positive or more than twice the reference (",
              ref_method, ") estimate: returning the more numerically stable ",
              ref_method, " estimate.")
      Vt <- Vt_ref
      par.se <- data.frame(Vt.se = if (!is.null(Vt.start)) NA_real_ else ref_se)
      fallback_used <- TRUE
      fallback_method <- ref_method
    } else {
      # fallback == "none" (or no usable reference): keep the raw LEGA estimate
      # but warn so a pathological value is not returned silently.
      warning("LEGA estimate of Vt looks unreliable (non-positive, non-finite, ",
              "or implausibly large), but fallback is disabled: returning the raw ",
              "LEGA estimate.", call. = FALSE)
    }
  }

  # Output

  par <- data.frame(Vt = Vt)

  tacs <- data.frame(Time = t_tac, Target = tac, Target_uncor = tac_uncor, # uncorrected for blood volume
                     Duration = dur)

  fitvals <- data.frame(
    Time = t_tac[eq_idx], Target = Z_eq, A = A_eq, B = B_eq,
    Target_fitted = as.numeric(fitted(lega_model)), Weights = w_eq,
    gamma = gamma_hat
  )

  input <- newvals$input

  out <- list(
    par = par, par.se = par.se, fit = lega_model, tacs = tacs, fitvals = fitvals,
    input = input, weights = weights, inpshift = inpshift, vB = vB,
    tstarIncludedFrames = tstarIncludedFrames, gamma = gamma_hat,
    fallback_used = fallback_used, fallback_method = fallback_method,
    model = "LEGA"
  )

  class(out) <- c("LEGA", "kinfit")

  return(out)
}

#' Plot: Likelihood Estimation in Graphical Analysis (LEGA)
#'
#' Function to visualise the fit of the LEGA model to data. Because LEGA models
#' the original (untransformed) tissue time activity curve, the fit is shown in
#' the original TAC domain.
#'
#' @param LEGAout The output object of the LEGA fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' tac <- pbr28$tacs[[2]]$FC
#' weights <- pbr28$tacs[[2]]$Weights
#' dur <- pbr28$tacs[[2]]$Duration/60
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1
#' )
#'
#' fit <- LEGA(t_tac, dur, tac, input, 10, weights)
#' plot_LEGAfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_LEGAfit <- function(LEGAout, roiname = NULL) {
  measured <- data.frame(
    Time = LEGAout$tacs$Time,
    Target.measured = LEGAout$tacs$Target,
    Weights = LEGAout$weights
  )

  fitted <- data.frame(
    Time = LEGAout$fitvals$Time,
    Target.fitted = LEGAout$fitvals$Target_fitted,
    Weights = LEGAout$fitvals$Weights
  )

  if (is.null(roiname)) {
    roiname <- "ROI"
  }

  measured <- dplyr::rename(measured, !!paste0(roiname, ".measured") := Target.measured)

  fitted <- dplyr::rename(fitted, !!paste0(roiname, ".fitted") := Target.fitted)

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
    geom_line(data = tidyfitted) + colScale +
    guides(shape = "none", color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

  return(outplot)
}


#' Tstar Finder: Likelihood Estimation in Graphical Analysis (LEGA)
#'
#' Function to identify where t* is for LEGA.
#'
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param dur Numeric vector of the time durations of the frames in minutes. The
#' LEGA estimator is explicitly frame-aware, so frame durations must be provided.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param filename The name of the output image: filename_LEGA.jpeg
#' @param inpshift Optional. The number of minutes by which to shift the timing of the input data frame forwards or backwards.
#' If not specified, this will be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will
#'   be ignored and assumed to be 0%. If specified, it will be corrected for
#'   prior to parameter estimation using the following equation: \deqn{C_{T}(t)
#'   = \frac{C_{Measured}(t) - vB\times C_{B}(t)}{1-vB}}
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_LEGA.jpeg
#'
#' @examples
#' \dontrun{
#' LEGA_tstar(t_tac, dur, lowroi, medroi, highroi, input,
#'   filename = "demonstration",
#'   inpshift = onetcmout$par$inpshift, frameStartEnd, gridbreaks = 4
#' )
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export


LEGA_tstar <- function(t_tac, dur, lowroi, medroi, highroi, input, filename = NULL, inpshift = 0, vB = 0, frameStartEnd = NULL, timeStartEnd = NULL, gridbreaks = 2) {
  frameStartEnd <- tstar_frameStartEnd(t_tac, frameStartEnd, timeStartEnd)

  frames <- length(t_tac)
  lowroi_fit <- LEGA(t_tac, dur, lowroi, input, tstar = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  medroi_fit <- LEGA(t_tac, dur, medroi, input, tstar = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  highroi_fit <- LEGA(t_tac, dur, highroi, input, tstar = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)

  low_linplot <- plot_LEGAfit(lowroi_fit) + ggtitle("Low") + ylim(0, max(lowroi_fit$tacs$Target * 1.1)) + theme(legend.position = "none")
  med_linplot <- plot_LEGAfit(medroi_fit) + ggtitle("Medium") + ylim(0, max(medroi_fit$tacs$Target * 1.1)) + theme(legend.position = "none")
  high_linplot <- plot_LEGAfit(highroi_fit) + ggtitle("High") + ylim(0, max(highroi_fit$tacs$Target * 1.1)) + theme(legend.position = "none")

  tstarInclFrames <- 3:frames

  fitfunc <- function(roitac, tstar) {
    LEGA(t_tac, dur, roitac, input, tstar = tstar, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  }

  comp <- tstar_compute(t_tac, lowroi, medroi, highroi, fitfunc,
                        function(f) f$par$Vt, tstarInclFrames)

  outcome_vals <- c(comp$outcome_df$Low, comp$outcome_df$Medium, comp$outcome_df$High)
  outcome_ylim <- c(min(outcome_vals), max(outcome_vals))
  if (outcome_ylim[2] > 20 || outcome_ylim[1] < 0) {
    outcome_ylim <- c(0, 20)
  }

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow = 1)

  totalplot <- tstar_finalise_plot(
    linrow, lowroi_fit, medroi_fit, highroi_fit,
    comp$r2_df, comp$maxperc_df, comp$outcome_df,
    outcome_ylab = expression(V[T]), t_tac = t_tac,
    tstarInclFrames = tstarInclFrames, gridbreaks = gridbreaks,
    outcome_ylim = outcome_ylim, filename = filename, modelname = "LEGA"
  )

  return(totalplot)
}
