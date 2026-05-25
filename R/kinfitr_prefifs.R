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
#' @param reftac Numeric vector of radioactivity concentrations in the
#'   pseudo-reference region for each frame.
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
#'   fitting. Defaults to uniform.
#' @param scale_time Numeric. Upper limit (minutes) of the early window over
#'   which the scale factor is computed (paper default is 5 minutes).
#' @param frameStartEnd Optional. First and last frame indices to use.
#' @param timeStartEnd Optional. First and last times (in minutes) to use.
#' @param K1.start,K1.lower,K1.upper Starting value and bounds for K1.
#' @param k2.start,k2.lower,k2.upper Starting value and bounds for k2.
#' @param multstart_iter Number of starting parameter iterations for the 1TC
#'   fit. See \code{\link[nls.multstart]{nls_multstart}}.
#' @param multstart_lower,multstart_upper Optional named lists for multstart
#'   starting bounds.
#' @param fit_t0 Should a t0 parameter be fitted in the pRef smoother? Passed
#'   through to \code{\link{feng_1tc_tac}}.
#' @param pref_multstart_iter Number of starting parameter iterations for the
#'   pRef smoothing step.
#' @param derivative Method for computing \eqn{dC_T'(t)/dt}. Either
#'   "analytical" (default; uses the 1TC ODE directly:
#'   \eqn{dC_T'/dt = Ph_1 \cdot Feng\_AIF(t) - Th_1 \cdot C_T'(t)}) or
#'   "symbolic" (uses R's \code{D()} on the smoother expression). Both should
#'   give numerically identical results; the option is exposed mainly for
#'   verification.
#' @param printvals If \code{TRUE}, the fitter prints iteration values.
#'
#' @return A list of class \code{c("pRefIFS_1tcm", "kinfit")} containing:
#'   \itemize{
#'     \item \code{par} — a data frame with K1, k2, Vt
#'     \item \code{par.se} — percentage standard errors
#'     \item \code{fit} — the underlying \code{nls} fit
#'     \item \code{tacs} — Time, Target, Target_fitted
#'     \item \code{input} — the constructed input tibble (Time, Blood, Plasma,
#'       ParentFraction, AIF)
#'     \item \code{pRefIFS} — Time, unscaled and scaled pRef-IFS curves
#'     \item \code{pref_fit} — the \code{feng_1tc_tac} result for the pRef TAC
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
#' fit <- pRefIFS_1tcm(
#'   t_tac, tac_pref, tac_tgt,
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
pRefIFS_1tcm <- function(t_tac, reftac, roitac,
                         t_blood, blood,
                         k2prime,
                         weights = NULL,
                         scale_time = 5,
                         frameStartEnd = NULL, timeStartEnd = NULL,
                         K1.start = 0.1, K1.lower = 0.0001, K1.upper = 1,
                         k2.start = 0.1, k2.lower = 0.0001, k2.upper = 0.5,
                         multstart_iter = 1,
                         multstart_lower = NULL, multstart_upper = NULL,
                         fit_t0 = TRUE,
                         pref_multstart_iter = 500,
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
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  # 1. Build the scaled pRef-IFS input function from the pRef TAC and early blood
  shape <- pRefIFS_shape(
    t_tac = t_tac, reftac = reftac,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime,
    weights = weights,
    scale_time = scale_time,
    frameStartEnd = frameStartEnd,
    fit_t0 = fit_t0,
    pref_multstart_iter = pref_multstart_iter,
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

  out <- list(
    par = par, par.se = par.se,
    fit = output,
    tacs = tacs,
    input = input,
    pRefIFS = shape$pRefIFS,
    pref_fit = shape$pref_fit,
    k2prime = k2prime,
    scale_factor = shape$scale_factor,
    scale_time = shape$scale_time,
    derivative_method = derivative,
    weights = weights_full,
    model = "pRefIFS_1tcm"
  )
  class(out) <- c("pRefIFS_1tcm", "kinfit")
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
#' @inheritParams pRefIFS_1tcm
#' @param interpPoints Integer. Number of points in the uniform fine time grid
#'   used for the derivative, scaling and convolution. Defaults to 6000 to
#'   match \code{\link{blood_interp}}.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{input} — tibble with Time, Blood, Plasma, ParentFraction, AIF
#'     \item \code{pRefIFS} — data frame with Time, pRefIFS_unscaled, pRefIFS_scaled
#'     \item \code{pref_fit} — the underlying \code{feng_1tc_tac} fit
#'     \item \code{scale_factor}, \code{scale_time}, \code{k2prime},
#'       \code{derivative_method}
#'   }
#' @export
pRefIFS_shape <- function(t_tac, reftac,
                          t_blood, blood,
                          k2prime,
                          weights = NULL,
                          scale_time = 5,
                          frameStartEnd = NULL,
                          fit_t0 = TRUE,
                          pref_multstart_iter = 500,
                          derivative = c("analytical", "symbolic"),
                          interpPoints = 6000) {

  if (missing(k2prime) || is.null(k2prime)) {
    stop("k2prime must be supplied (e.g. from SRTM/SRTMC or a population value).")
  }
  if (length(k2prime) != 1 || !is.finite(k2prime) || k2prime <= 0) {
    stop("k2prime must be a single positive numeric value.")
  }
  derivative <- match.arg(derivative)

  # 1. Smooth the pRef TAC
  pref_fit <- feng_1tc_tac(
    t_tac = t_tac, tac = reftac, weights = weights,
    fit_t0 = fit_t0,
    frameStartEnd = frameStartEnd,
    multstart_iter = pref_multstart_iter
  )

  # 2. Build a uniform fine time grid spanning the TAC
  t_max <- max(t_tac)
  grid <- seq(0, t_max, length.out = interpPoints)

  # 3. Predict smoothed C_T'(t) and dC_T'/dt on the grid
  coefs <- as.list(coef(pref_fit$fit))
  if (is.null(coefs$t0)) coefs$t0 <- 0    # fit_t0 = FALSE case

  Ct_grid  <- predict(pref_fit$fit, newdata = list(t_tac = grid))
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

  blood_on_grid_early <- stats::approx(
    x = t_blood, y = blood,
    xout = grid[early_idx],
    rule = 2  # extrapolate flat at the ends
  )$y

  auc_blood <- pracma::trapz(grid[early_idx], blood_on_grid_early)
  auc_pref  <- pracma::trapz(grid[early_idx], prefifs_unscaled[early_idx])

  if (!is.finite(auc_pref) || auc_pref <= 0) {
    stop("Computed early AUC of the unscaled pRef-IFS is non-positive (",
         signif(auc_pref, 3), "); cannot derive a scale factor. ",
         "Inspect the pRef smoothing fit.")
  }
  SF <- auc_blood / auc_pref

  prefifs_scaled <- SF * prefifs_unscaled

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
    pRefIFS = data.frame(
      Time              = grid,
      pRefIFS_unscaled  = prefifs_unscaled,
      pRefIFS_scaled    = prefifs_scaled
    ),
    pref_fit = pref_fit,
    scale_factor = SF,
    scale_time = early_end,
    k2prime = k2prime,
    derivative_method = derivative
  )
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
#' Two panels stacked: (top) the pRef TAC with its Feng+1TC smoother; (bottom)
#' the target TAC with its 1TC fit, alongside the recovered scaled pRef-IFS.
#'
#' @param pRefIFS_1tcmout The output object from \code{\link{pRefIFS_1tcm}}.
#' @param roiname Optional. Display name for the target region.
#' @param refname Optional. Display name for the pseudo-reference region.
#'
#' @return A ggplot2 object.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#' @export
plot_pRefIFS_1tcmfit <- function(pRefIFS_1tcmout,
                                 roiname = NULL,
                                 refname = NULL) {
  if (is.null(roiname)) roiname <- "ROI"
  if (is.null(refname)) refname <- "pRef"

  # Bottom panel data (target + AIF)
  measured <- data.frame(
    Time = pRefIFS_1tcmout$tacs$Time,
    Radioactivity = pRefIFS_1tcmout$tacs$Target,
    Weights = weights(pRefIFS_1tcmout$fit),
    Region = paste0(roiname, ".Measured")
  )

  aifdf <- data.frame(
    Time = pRefIFS_1tcmout$input$Time,
    Radioactivity = pRefIFS_1tcmout$input$AIF,
    Weights = 1,
    Region = "pRefIFS (scaled)"
  )

  i_fit <- predict(pRefIFS_1tcmout$fit, newdata = list(
    t_tac = pRefIFS_1tcmout$input$Time,
    tac = pracma::interp1(
      pRefIFS_1tcmout$tacs$Time,
      pRefIFS_1tcmout$tacs$Target,
      pRefIFS_1tcmout$input$Time,
      method = "linear"
    )
  ))

  fitdf <- data.frame(
    Time = pRefIFS_1tcmout$input$Time,
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
    coord_cartesian(ylim = c(0, max(measured$Radioactivity) * 1.5))
}
