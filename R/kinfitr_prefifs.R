#' Pseudo-Reference-based Input Function Shape (pRef-IFS) 1TC model
#'
#' Implements the pRef-IFS method of Volpi et al. (2026) for the 1-tissue
#' compartment (1TC) case. The shape of the metabolite-corrected plasma input
#' function is recovered from a pseudo-reference (pRef) region TAC by smoothing
#' with the Feng-AIF / 1TC IRF model, computing the first derivative of the
#' smoothed curve, and combining it with a supplied \code{k2prime}:
#' \deqn{pRef\_IFS(t) = dC_T'(t)/dt + k2' \cdot C_T'(t)}
#' The shape is then scaled by matching its early AUC to an early
#' image-derived blood time-activity curve (ID-BTAC), giving an estimated
#' input function. This is used to fit a 1TC model to the target ROI TAC.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame. If a time-zero frame is not included it
#'   will be added.
#' @param pref_par The fitted parameters of a \code{\link{feng_1tc_tac}}
#'   smoother applied to the pseudo-reference TAC. Either the 1-row data frame
#'   returned in \code{$par} or the full fit object (in which case \code{$par}
#'   is extracted automatically). Must contain the columns \code{A, B, C,
#'   alpha, beta, gamma, Ph1, Th1}; \code{t0} is optional and defaults to 0 if
#'   absent. Fitting the smoother once externally and re-using \code{pref_par}
#'   across many target ROIs avoids re-fitting the (expensive) smoother for
#'   every ROI in a brain.
#' @param roitac Numeric vector of radioactivity concentrations in the target
#'   region for each frame.
#' @param t_blood Numeric vector of times for each early blood sample in
#'   minutes. Used for the scaling step.
#' @param blood Numeric vector of early blood radioactivity values (e.g. an
#'   image-derived blood TAC).
#' @param k2prime Required. Numeric scalar. The pRef clearance rate constant,
#'   either estimated externally via SRTM/SRTMC or supplied as a
#'   population-average value.
#' @param weights Optional. Numeric vector of weights for each frame in the
#'   target 1TC fit. Defaults to uniform.
#' @param scale_time Numeric. Upper limit (minutes) of the early window over
#'   which the scale factor is computed (paper default is 5 minutes).
#' @param frameStartEnd Optional. First and last frame indices to use in the
#'   target 1TC fit.
#' @param timeStartEnd Optional. First and last times (in minutes) to use.
#' @param K1.start,K1.lower,K1.upper Starting value and bounds for K1.
#' @param k2.start,k2.lower,k2.upper Starting value and bounds for k2.
#' @param multstart_iter Number of starting parameter iterations for the 1TC
#'   fit. See \code{\link[nls.multstart]{nls_multstart}}.
#' @param multstart_lower,multstart_upper Optional named lists for multstart
#'   starting bounds.
#' @param derivative Method for computing \eqn{dC_T'(t)/dt}. Either
#'   "analytical" (default; uses the 1TC ODE directly:
#'   \eqn{dC_T'/dt = Ph_1 \cdot Feng\_AIF(t) - Th_1 \cdot C_T'(t)}) or
#'   "symbolic" (uses R's \code{D()} on the smoother expression). Both should
#'   give numerically identical results; the option is exposed mainly for
#'   verification.
#' @param printvals If \code{TRUE}, the fitter prints iteration values.
#'
#' @return A list of class \code{c("prefifs_1tcm", "kinfit")} containing:
#'   \itemize{
#'     \item \code{par} — a data frame with K1, k2, Vt
#'     \item \code{par.se} — percentage standard errors
#'     \item \code{fit} — the underlying \code{nls} fit
#'     \item \code{tacs} — Time, Target, Target_fitted
#'     \item \code{input} — the constructed input tibble (Time, Blood, Plasma,
#'       ParentFraction, AIF)
#'     \item \code{prefifs} — Time, unscaled and scaled pRef-IFS curves
#'     \item \code{pref_par} — the smoother parameters used (echoed back)
#'     \item \code{k2prime}, \code{scale_factor}, \code{scale_time},
#'       \code{derivative_method}, \code{weights}, \code{model}
#'   }
#'
#' @examples
#' \dontrun{
#' data(pbr28)
#' tac_pref <- pbr28$tacs[[1]]$CBL
#' tac_tgt  <- pbr28$tacs[[1]]$FC
#' t_tac    <- pbr28$tacs[[1]]$Times / 60
#' weights  <- pbr28$tacs[[1]]$Weights
#' bd       <- pbr28$procblood[[1]]
#'
#' # Fit the smoother on the pseudo-reference TAC once
#' pref_par <- feng_1tc_tac(t_tac, tac_pref, weights)$par
#'
#' # ... and reuse for any number of target ROIs
#' fit <- prefifs_1tcm(
#'   t_tac, pref_par, tac_tgt,
#'   t_blood = bd$Time / 60, blood = bd$Cbl_dispcorr,
#'   k2prime = 0.05, weights = weights
#' )
#' plot(fit)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Volpi T, Naganawa M, Carson RE, Hillmer AT. "Pseudo-reference"-based
#'   Input Function Shape (pRef-IFS): Towards a True image-derived input
#'   function for PET kinetic modeling. J Cereb Blood Flow Metab. 2026;46(5):1238-1252.
#'
#' @export
prefifs_1tcm <- function(t_tac, pref_par, roitac,
                         t_blood, blood,
                         k2prime,
                         weights = NULL,
                         scale_time = 5,
                         frameStartEnd = NULL, timeStartEnd = NULL,
                         K1.start = 0.1, K1.lower = 0.0001, K1.upper = 1,
                         k2.start = 0.1, k2.lower = 0.0001, k2.upper = 0.5,
                         multstart_iter = 1,
                         multstart_lower = NULL, multstart_upper = NULL,
                         derivative = c("analytical", "symbolic"),
                         printvals = FALSE) {

  if (missing(k2prime) || is.null(k2prime)) {
    stop("k2prime must be supplied (e.g. from SRTM/1TC or a population value).")
  }
  if (length(k2prime) != 1 || !is.finite(k2prime) || k2prime <= 0) {
    stop("k2prime must be a single positive numeric value.")
  }
  derivative <- match.arg(derivative)

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    if (timeStartEnd[1] > max(t_tac) || timeStartEnd[2] < min(t_tac)) {
      stop("timeStartEnd = c(", timeStartEnd[1], ", ", timeStartEnd[2],
           ") lies outside the t_tac range [", min(t_tac), ", ", max(t_tac),
           "]. Did you mix up seconds and minutes?")
    }
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  # 1. Build the scaled pRef-IFS input function from the pre-fit smoother and blood
  shape <- prefifs_shape(
    t_tac = t_tac, pref_par = pref_par,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime,
    scale_time = scale_time,
    derivative = derivative
  )

  input <- shape$input

  # 2. Fit 1TCM to the target TAC using the constructed input
  modeldata <- tidyinput_art(t_tac, roitac, weights, frameStartEnd)
  weights_full <- modeldata$weights
  modeldata <- as.list(modeldata)
  modeldata$input <- input
  names(modeldata)[names(modeldata) == "tac"] <- "tac"

  start <- c(K1 = K1.start, k2 = k2.start)
  lower <- c(K1 = K1.lower, k2 = k2.lower)
  upper <- c(K1 = K1.upper, k2 = k2.upper)

  multstart_pars <- fix_multstartpars(
    start, lower, upper, multstart_iter,
    multstart_lower, multstart_upper
  )
  multstart_upper <- multstart_pars$multstart_upper
  multstart_lower <- multstart_pars$multstart_lower

  # vB is fixed at 0 in this implementation; bake it into the formula
  formula <- as.formula(
    "tac ~ onetcm_model(t_tac, input, K1, k2, vB = 0)"
  )

  if (prod(multstart_iter) == 1) {
    output <- minpack.lm::nlsLM(
      formula, data = modeldata,
      start = start, lower = lower, upper = upper,
      weights = weights_full,
      control = minpack.lm::nls.lm.control(maxiter = 200),
      trace = printvals
    )
  } else {
    # nls.multstart evaluates `modelweights` as a symbol inside `data`, so we
    # must hand it the name of a column in `modeldata` (here, `weights`, which
    # tidyinput_art() guarantees is present and length-matched).
    output <- nls.multstart::nls_multstart(
      formula, data = modeldata,
      supp_errors = "Y",
      start_lower = multstart_lower, start_upper = multstart_upper,
      iter = multstart_iter, convergence_count = FALSE,
      lower = lower, upper = upper, modelweights = weights
    )
  }

  # Limit-hit warning, matching onetcm()
  limcheck_u <- purrr::map2_lgl(round(upper, 3), round(coef(output), 3), identical)
  limcheck_l <- purrr::map2_lgl(round(lower, 3), round(coef(output), 3), identical)
  limcheck <- (limcheck_u + limcheck_l) == 1
  if (any(limcheck)) {
    warning(
      "\nFitted parameters are hitting upper or lower limit bounds. Consider \n",
      "either modifying the upper and lower limit boundaries, or else using \n",
      "multstart when fitting the model (see the function documentation).\n"
    )
  }

  tacs <- data.frame(
    Time = modeldata$t_tac,
    Target = modeldata$tac,
    Target_fitted = as.numeric(fitted(output))
  )

  par <- as.data.frame(as.list(coef(output)))
  par.se <- par
  par.se[1, ] <- purrr::map_dbl(names(par), ~ get_se(output, .x))
  names(par.se) <- paste0(names(par.se), ".se")

  par$Vt <- par$K1 / par$k2
  par.se$Vt.se <- get_se(output, "K1/k2")

  # Diagnostic: AUC of the target TAC over the same early window used for
  # scaling. Lets the user eyeball whether the scaled AIF is the right
  # magnitude relative to the tissue uptake (useful when SF looks suspect).
  early_tac_idx <- tacs$Time <= shape$scale_time
  if (sum(early_tac_idx) >= 2) {
    auc_target_early <- pracma::trapz(tacs$Time[early_tac_idx],
                                      tacs$Target[early_tac_idx])
  } else {
    auc_target_early <- NA_real_
  }

  out <- list(
    par = par, par.se = par.se,
    fit = output,
    tacs = tacs,
    input = input,
    prefifs = shape$prefifs,
    pref_par = shape$pref_par,
    k2prime = k2prime,
    scale_factor = shape$scale_factor,
    scale_time = shape$scale_time,
    auc_blood_early  = shape$auc_blood_early,
    auc_pref_early   = shape$auc_pref_early,
    auc_target_early = auc_target_early,
    frac_clamped     = shape$frac_clamped,
    derivative_method = derivative,
    weights = weights_full,
    model = "prefifs_1tcm"
  )
  class(out) <- c("prefifs_1tcm", "kinfit")
  out
}


#' Recover the scaled pRef-IFS input function
#'
#' Standalone helper that runs the pRef-IFS recovery (smooth the pRef TAC,
#' compute its derivative, form
#' \eqn{dC_T'(t)/dt + k_2' \cdot C_T'(t)}, and scale by matching the early
#' AUC of an image-derived blood TAC). Returns the constructed input object
#' alongside the unscaled and scaled curves so the shape can be inspected
#' separately from the 1TC fit.
#'
#' @inheritParams prefifs_1tcm
#' @param interpPoints Integer. Number of points in the uniform fine time grid
#'   used for the derivative, scaling and convolution. Defaults to 6000 to
#'   match \code{\link{blood_interp}}.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{input} — tibble with Time, Blood, Plasma, ParentFraction, AIF
#'     \item \code{prefifs} — data frame with Time, prefifs_unscaled, prefifs_scaled
#'     \item \code{pref_par} — the smoother parameters used (echoed back)
#'     \item \code{scale_factor}, \code{scale_time}, \code{k2prime},
#'       \code{derivative_method}
#'   }
#' @export
prefifs_shape <- function(t_tac, pref_par,
                          t_blood, blood,
                          k2prime,
                          scale_time = 5,
                          derivative = c("analytical", "symbolic"),
                          interpPoints = 6000) {

  if (missing(k2prime) || is.null(k2prime)) {
    stop("k2prime must be supplied (e.g. from SRTM/1TC or a population value).")
  }
  if (length(k2prime) != 1 || !is.finite(k2prime) || k2prime <= 0) {
    stop("k2prime must be a single positive numeric value.")
  }
  derivative <- match.arg(derivative)

  # Strip NA pairs from the blood data (cryptic 'missing value where TRUE/FALSE
  # needed' downstream otherwise).
  lengths_match <- length(t_blood) == length(blood)
  if (!lengths_match) {
    stop("t_blood and blood must have the same length (",
         length(t_blood), " vs ", length(blood), ").")
  }
  na_pairs <- !is.finite(t_blood) | !is.finite(blood)
  if (any(na_pairs)) {
    warning(sum(na_pairs), " NA/non-finite sample(s) in t_blood / blood were dropped.")
    t_blood <- t_blood[!na_pairs]
    blood   <- blood[!na_pairs]
  }
  if (length(t_blood) < 2) {
    stop("After dropping NA samples, t_blood / blood has fewer than 2 points.")
  }

  coefs <- coerce_pref_par(pref_par)

  # Build a uniform fine time grid spanning the TAC
  t_max <- max(t_tac)
  grid <- seq(0, t_max, length.out = interpPoints)

  # Evaluate the smoothed C_T'(t) and dC_T'/dt on the grid directly from coefs
  Ct_grid <- with(coefs, feng_1tc_tac_model(
    grid, t0, A, B, C, alpha, beta, gamma, Ph1, Th1
  ))
  dCt_grid <- feng_1tc_tac_deriv(grid, coefs, method = derivative)

  # 4. Unscaled pRef-IFS
  prefifs_unscaled <- dCt_grid + k2prime * Ct_grid
  prefifs_unscaled[!is.finite(prefifs_unscaled)] <- 0

  # 5. Scale by matching the early AUC to the early blood AUC
  early_end <- min(scale_time, max(t_blood))
  if (early_end < scale_time) {
    warning(sprintf(
      "max(t_blood) = %.3f is below scale_time = %.3f; using %.3f as the early window upper limit.",
      max(t_blood), scale_time, early_end
    ))
  }

  early_idx <- grid <= early_end
  if (sum(early_idx) < 2) {
    stop("The early scaling window contains fewer than 2 grid points. ",
         "Check t_blood and scale_time.")
  }

  # Pad with (t=0, blood=0) before approx() so that grid times below
  # min(t_blood) interpolate from zero rather than flat-extrapolating the
  # first sample (which would inflate auc_blood when the BTAC starts post-zero
  # on the rising edge). Keep rule=2 so the LATE end still uses the last
  # observed value (extrapolating to zero on the tail would be wrong).
  if (min(t_blood) > 0) {
    t_blood_pad <- c(0, t_blood)
    blood_pad   <- c(0, blood)
  } else {
    t_blood_pad <- t_blood
    blood_pad   <- blood
  }
  blood_on_grid_early <- stats::approx(
    x = t_blood_pad, y = blood_pad,
    xout = grid[early_idx],
    rule = 2
  )$y

  auc_blood <- pracma::trapz(grid[early_idx], blood_on_grid_early)
  auc_pref  <- pracma::trapz(grid[early_idx], prefifs_unscaled[early_idx])

  if (!is.finite(auc_blood) || auc_blood <= 0) {
    stop("Computed early AUC of the supplied blood is non-positive (",
         signif(auc_blood, 3), "); cannot derive a scale factor. ",
         "Inspect the blood TAC over [0, scale_time].")
  }
  if (!is.finite(auc_pref) || auc_pref <= 0) {
    stop("Computed early AUC of the unscaled pRef-IFS is non-positive (",
         signif(auc_pref, 3), "); cannot derive a scale factor. ",
         "Inspect the pRef smoothing fit.")
  }
  SF <- auc_blood / auc_pref

  prefifs_scaled <- SF * prefifs_unscaled

  # Clamp negatives. Under the 1T assumption (Volpi et al. eq. 1) the bracketed
  # term dC_T'/dt + k2' * C_T' is K1' * C_p(t) and so must be non-negative;
  # negatives only appear when the smoother's Th1 differs materially from the
  # user-supplied k2prime. Inform the user how much was clamped so they can
  # judge whether k2prime is compatible with the smoother fit.
  neg_idx <- prefifs_scaled < 0
  if (any(neg_idx)) {
    frac_neg <- mean(neg_idx)
    if (frac_neg > 0.01) {
      message(sprintf(
        "prefifs_shape: %.1f%% of grid points had negative pRef-IFS values and were clamped to 0. This usually means k2prime (%.4g) differs from the smoother's Th1 (%.4g) more than the 1T assumption tolerates.",
        100 * frac_neg, k2prime, coefs$Th1
      ))
    }
    prefifs_scaled[neg_idx] <- 0
  }

  # 6. Build the input object. vB = 0 in this implementation, but we still
  # populate the Blood column with the scaled pRef-IFS so that a future
  # variant supporting vB fitting has the closest blood-like signal available.
  input <- tibble::tibble(
    Time           = grid,
    Blood          = prefifs_scaled,
    Plasma         = prefifs_scaled,
    ParentFraction = rep(1, length(grid)),
    AIF            = prefifs_scaled
  )
  class(input) <- c("interpblood", class(input))

  list(
    input = input,
    prefifs = data.frame(
      Time              = grid,
      prefifs_unscaled  = prefifs_unscaled,
      prefifs_scaled    = prefifs_scaled
    ),
    pref_par = as.data.frame(coefs),
    scale_factor = SF,
    scale_time = early_end,
    auc_blood_early = auc_blood,
    auc_pref_early  = auc_pref,
    frac_clamped    = if (any(neg_idx)) mean(neg_idx) else 0,
    k2prime = k2prime,
    derivative_method = derivative
  )
}


# Internal: coerce a user-supplied `pref_par` into a named list of coefficients
# for the feng_1tc_tac smoother. Accepts either a `feng_1tc_tac` fit object
# (extracts $par) or a data.frame / named list with the required columns.
# Defaults `t0` to 0 if absent (the fit_t0 = FALSE case in feng_1tc_tac).
coerce_pref_par <- function(pref_par) {
  if (missing(pref_par) || is.null(pref_par)) {
    stop("pref_par must be supplied (the $par data frame from feng_1tc_tac).")
  }
  if (inherits(pref_par, "kinfit") || inherits(pref_par, "feng_1tc_tac")) {
    if (is.null(pref_par$par)) {
      stop("pref_par looks like a kinfit object but has no $par element.")
    }
    pref_par <- pref_par$par
  }
  if (is.data.frame(pref_par)) {
    if (nrow(pref_par) != 1) {
      stop("pref_par must be a 1-row data frame (one fitted smoother).")
    }
    coefs <- as.list(pref_par)
  } else if (is.list(pref_par)) {
    coefs <- pref_par
  } else {
    stop("pref_par must be a data frame, a named list, or a feng_1tc_tac fit object.")
  }

  required <- c("A", "B", "C", "alpha", "beta", "gamma", "Ph1", "Th1")
  missing_cols <- setdiff(required, names(coefs))
  if (length(missing_cols) > 0) {
    stop("pref_par is missing required column(s): ",
         paste(missing_cols, collapse = ", "), ".")
  }
  if (is.null(coefs$t0)) {
    message("coerce_pref_par: no `t0` column in pref_par; defaulting to t0 = 0.")
    coefs$t0 <- 0
  }
  # Flatten any 1-length vectors to scalars (data.frame cols come through as length-1)
  coefs <- lapply(coefs, function(x) if (length(x) == 1) as.numeric(x) else x)
  coefs
}


# Internal: derivative of feng_1tc_tac_model at a given grid using fitted coefs.
# method = "analytical" uses the 1TC ODE
#   dC_T'/dt = Ph1 * Feng_AIF(t - t0) - Th1 * C_T'(t)
# method = "symbolic" uses R's D() on the smoother body.
feng_1tc_tac_deriv <- function(t_tac, coefs, method = c("analytical", "symbolic")) {
  method <- match.arg(method)
  t0    <- coefs$t0
  A     <- coefs$A
  B     <- coefs$B
  C     <- coefs$C
  alpha <- coefs$alpha
  beta  <- coefs$beta
  gamma <- coefs$gamma
  Ph1   <- coefs$Ph1
  Th1   <- coefs$Th1

  if (method == "analytical") {
    tcorr <- t_tac - t0
    pos <- tcorr > 0
    # Feng AIF on the corrected time axis
    feng_aif <- (A * tcorr - B - C) * exp(-alpha * tcorr) +
                B * exp(-beta * tcorr) +
                C * exp(-gamma * tcorr)
    # C_T' from the smoother (re-use the model expression)
    Ct <- feng_1tc_tac_model(t_tac, t0, A, B, C, alpha, beta, gamma, Ph1, Th1)
    out <- Ph1 * feng_aif - Th1 * Ct
    out[!pos] <- 0
    return(out)
  }

  # symbolic
  dfun <- feng_1tc_tac_deriv_symbolic_expr()
  env <- list(
    t_tac = t_tac,
    t0 = t0, A = A, B = B, C = C,
    alpha = alpha, beta = beta, gamma = gamma,
    Ph1 = Ph1, Th1 = Th1
  )
  out <- eval(dfun, envir = env)
  out * as.numeric(t_tac > t0)
}


# Builds the symbolic derivative expression of the smoother body (without the
# t > t0 step function, which is non-differentiable). Cached on first call.
feng_1tc_tac_deriv_symbolic_expr <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    inner <- quote(
      Ph1 * (
        (B * (-1 + exp((t_tac - t0) * (-alpha + Th1)))) / (alpha - Th1) +
        (C * (-1 + exp((t_tac - t0) * (-alpha + Th1)))) / (alpha - Th1) +
        (B * (-1 + exp((t_tac - t0) * (-beta + Th1))))  / (-beta + Th1) +
        (C * (-1 + exp((t_tac - t0) * (-gamma + Th1)))) / (-gamma + Th1) +
        (A * (1 + exp((t_tac - t0) * (-alpha + Th1)) *
              (-1 - alpha * (t_tac - t0) + (t_tac - t0) * Th1))) /
          (alpha - Th1)^2
      ) / exp((t_tac - t0) * Th1)
    )
    cached <<- stats::D(inner, "t_tac")
    cached
  }
})


#' Plot the fit of a pRef-IFS 1TC model
#'
#' Single-panel plot showing the measured target TAC, the 1TC-fitted target
#' curve, and the recovered scaled pRef-IFS used as the input function. The
#' y-axis is capped a little above the target peak so the target fit is easy
#' to inspect; the scaled pRef-IFS typically peaks well above this cap by
#' construction (it is matched to early blood AUC, not to tissue uptake).
#'
#' @param prefifs_1tcmout The output object from \code{\link{prefifs_1tcm}}.
#' @param roiname Optional. Display name for the target region.
#'
#' @return A ggplot2 object.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#' @export
plot_prefifs_1tcmfit <- function(prefifs_1tcmout,
                                 roiname = NULL) {
  if (is.null(roiname)) roiname <- "ROI"

  measured <- data.frame(
    Time = prefifs_1tcmout$tacs$Time,
    Radioactivity = prefifs_1tcmout$tacs$Target,
    Weights = weights(prefifs_1tcmout$fit),
    Region = paste0(roiname, ".Measured")
  )

  # Restrict the fine-grid curves (AIF and fitted) to the fit window so that
  # frameStartEnd / timeStartEnd subsets plot correctly (the fine grid in
  # input$Time spans the full t_tac, but tacs$Time only spans the fit subset).
  fit_window <- range(measured$Time)
  fine_idx   <- prefifs_1tcmout$input$Time >= fit_window[1] &
                prefifs_1tcmout$input$Time <= fit_window[2]
  fine_time  <- prefifs_1tcmout$input$Time[fine_idx]
  fine_aif   <- prefifs_1tcmout$input$AIF[fine_idx]

  aifdf <- data.frame(
    Time = fine_time,
    Radioactivity = fine_aif,
    Weights = 1,
    Region = "pRefIFS (scaled)"
  )

  # Re-predict the fitted curve on the (clipped) fine grid. The formula's
  # response variable `tac` is not actually used by predict() — only RHS
  # variables matter — but the captured-data nls predict still requires the
  # symbol to resolve, so we pass a same-length dummy.
  i_fit <- predict(prefifs_1tcmout$fit, newdata = list(
    t_tac = fine_time,
    tac   = rep(0, length(fine_time))
  ))

  fitdf <- data.frame(
    Time = fine_time,
    Radioactivity = i_fit,
    Weights = 1,
    Region = paste0(roiname, ".Fitted")
  )

  plotdf <- rbind(aifdf, measured, fitdf)
  plotdf$Region <- forcats::fct_inorder(factor(plotdf$Region))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Region)
  colScale <- scale_colour_manual(name = "Curve", values = myColors)

  ggplot(plotdf, aes(x = Time, y = Radioactivity, colour = Region)) +
    colScale +
    geom_point(
      data = subset(plotdf, plotdf$Region == paste0(roiname, ".Measured")),
      aes(shape = "a", size = Weights)
    ) +
    geom_line(
      data = subset(plotdf, plotdf$Region != paste0(roiname, ".Measured"))
    ) +
    guides(shape = "none", color = guide_legend(order = 1)) +
    scale_size(range = c(1, 3)) +
    coord_cartesian(ylim = c(0, max(measured$Radioactivity) * 1.3))
}
