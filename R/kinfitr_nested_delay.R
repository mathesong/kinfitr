#' Nested 1TCM with Delay Estimation
#'
#' Fits the One Tissue Compartment Model to multiple regions simultaneously,
#' estimating a shared input delay (inpshift) across all regions using optim()
#' in the outer loop, while fitting K1 and k2 per region using nlsLM in the
#' inner loop.
#'
#' @param t_tac Numeric vector of times for each frame in minutes (repeated for
#'   each region). If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations, stacked across
#'   regions.
#' @param region Character vector identifying the region for each observation.
#' @param input Data frame containing the blood, plasma, and parent fraction
#'   concentrations over time. Generated using \code{blood_interp}.
#' @param vB The blood volume fraction. Default is 0.05.
#' @param weights Optional. Numeric vector of weights.
#' @param roiweights Optional. Numeric vector of ROI weights. Can be a named
#'   vector (one per region), an unnamed vector of the same length as the
#'   number of regions, or a long vector the same length as \code{tac}
#'   (one per observation; the first value per region is used). ROI size or
#'   relative size is a useful quantity to supply here, as larger regions
#'   provide more precise mean TAC estimates. Normalised internally so the
#'   maximum is 1.
#' @param inpshift.start Starting value for the delay parameter. Default is 0.
#' @param inpshift.lower Lower bound for the delay parameter. Default is -0.5.
#' @param inpshift.upper Upper bound for the delay parameter. Default is 0.5.
#' @param K1.start Starting parameter for K1. Default is 0.1.
#' @param K1.lower Lower bound for K1. Default is 0.0001.
#' @param K1.upper Upper bound for K1. Default is 1.
#' @param k2.start Starting parameter for k2. Default is 0.1.
#' @param k2.lower Lower bound for k2. Default is 0.0001.
#' @param k2.upper Upper bound for k2. Default is 0.5.
#' @param optim_method Optimization method for optim(). Default is "L-BFGS-B".
#'   Only "L-BFGS-B" supports bounds on the shared (outer) parameters; if
#'   another method is selected, the lower and upper limits for shared
#'   parameters are ignored.
#' @param optim_control List of control parameters for optim().
#' @param multstart_iter Number of multistart iterations for the inner
#'   (per-region) fits. Default is 1. A single number, or one value per fitted
#'   parameter. Where more than one iteration is used, starting values are
#'   chosen by improved Latin hypercube sampling over the multstart bounds
#'   rather than at random.
#' @param multstart_lower Optional. Lower bounds for multistart starting params.
#' @param multstart_upper Optional. Upper bounds for multistart starting params.
#' @param frameStartEnd Optional. Frame range c(start, end). Applied per region.
#' @param timeStartEnd Optional. Time range for frame selection.
#' @param keep_inner_fits Logical. If TRUE, store per-region nlsLM fit objects
#'   in \code{$fits_inner}. Default is FALSE.
#'
#' @details The recommended usage for delay estimation is to fit the first few
#'   minutes of the acquisition (5 minutes is usually a good default) and not
#'   to use frame weights, so as to prioritise the low radioactivity values
#'   where the delay signal is most apparent.
#'
#'   The standard errors in \code{par.se} follow the usual kinfitr
#'   convention of being expressed as a fraction of the estimate. Those for the
#'   shared (outer) parameter are derived from the curvature of the profile
#'   objective at its optimum, obtained by finite differences, and are therefore
#'   more approximate than those of the per-region parameters. The per-region
#'   standard errors are conditional on the shared parameter: they treat
#'   the fitted delay as known rather than estimated, and so are somewhat
#'   optimistic.
#'
#' @return A list with class c("nested_1tcm_delay", "kinfit").
#'   \code{par} and \code{par.se} have one row per region, with the shared
#'   parameter values repeated. \code{weights} holds the full stacked vector
#'   of frame weights, aligned with the rows of \code{tacs}.
#'
#' @examples
#' \dontrun{
#' data("pbr28")
#' tac_wide <- pbr28$tacs[[2]]
#' t_tac <- tac_wide$Times / 60
#' input <- pbr28$input[[2]]
#'
#' regions <- c("FC", "TC", "THA")
#' long_data <- do.call(rbind, lapply(regions, function(r) {
#'   data.frame(t_tac = t_tac, tac = tac_wide[[r]], region = r)
#' }))
#'
#' out <- nested_1tcm_delay(
#'   long_data$t_tac, long_data$tac, long_data$region, input,
#'   timeStartEnd = c(0, 5)
#' )
#' plot(out)
#' }
#'
#' @export
nested_1tcm_delay <- function(
    t_tac, tac, region, input,
    vB = 0.05,
    weights = NULL, roiweights = NULL,
    inpshift.start = 0, inpshift.lower = -0.5, inpshift.upper = 0.5,
    K1.start = 0.1, K1.lower = 0.0001, K1.upper = 1,
    k2.start = 0.1, k2.lower = 0.0001, k2.upper = 0.5,
    optim_method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN"),
    optim_control = list(),
    multstart_iter = 1,
    multstart_lower = NULL, multstart_upper = NULL,
    frameStartEnd = NULL, timeStartEnd = NULL,
    keep_inner_fits = FALSE) {

  optim_method <- match.arg(optim_method)

  model_name <- "nested_1tcm_delay"
  unnested_name <- "onetcm"

  region <- as.character(region)
  regions_raw <- .nested_regions(region, model_name, unnested_name)

  # Resolve roiweights against the original (pre-tidy) region vector, since a
  # per-observation roiweights vector is sized to the input TACs. tidyinput_long
  # may prepend a zero frame per region, which would otherwise make it too short.
  roiweights <- .nested_roiweights(roiweights, region, regions_raw)

  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    t_first <- t_tac[region == regions_raw[1]]
    frameStartEnd <- c(which(t_first >= timeStartEnd[1])[1],
                       tail(which(t_first <= timeStartEnd[2]), 1))
  }

  tidied <- tidyinput_long(t_tac, tac, region, weights, frameStartEnd)
  t_tac   <- tidied$t_tac
  tac     <- tidied$tac
  region  <- tidied$region
  weights <- tidied$weights
  regions <- unique(region)

  weights_by_region <- .nested_weights_by_region(weights, region, regions)

  start <- c(K1 = K1.start, k2 = k2.start)
  lower <- c(K1 = K1.lower, k2 = k2.lower)
  upper <- c(K1 = K1.upper, k2 = k2.upper)

  multstart_bounds <- .nested_multstart_bounds(start, lower, upper,
                                               multstart_iter,
                                               multstart_lower, multstart_upper)
  multstart_lower <- multstart_bounds$lower
  multstart_upper <- multstart_bounds$upper

  formula_str <- paste0("tac ~ onetcm_model(t_tac, input, K1, k2, vB=", vB, ")")

  outer_objective <- function(outer_vals) {
    inpshift_val <- outer_vals[["inpshift"]]
    shifted <- shift_timings_long(t_tac, tac, region, input, inpshift_val)

    total_rss <- 0
    for (r in regions) {
      region_mask <- shifted$region == r
      modeldata <- list(
        tac = shifted$tac[region_mask],
        t_tac = shifted$t_tac[region_mask],
        weights = weights_by_region[[r]],
        input = shifted$input
      )
      fit <- .nested_fit_region(formula_str, modeldata, start, lower, upper,
                                multstart_iter, multstart_lower, multstart_upper)
      if (is.null(fit)) {
        total_rss <- total_rss + 1e10 * unname(roiweights[r])
      } else {
        total_rss <- total_rss +
          sum(weights(fit) * residuals(fit)^2) * unname(roiweights[r])
      }
    }
    return(total_rss)
  }

  optim_args <- list(
    par = c(inpshift = inpshift.start),
    fn = outer_objective, method = optim_method, control = optim_control
  )
  if (optim_method == "L-BFGS-B") {
    optim_args$lower <- c(inpshift = inpshift.lower)
    optim_args$upper <- c(inpshift = inpshift.upper)
  }
  optim_result <- do.call(optim, optim_args)

  # Final pass
  optimal_inpshift <- optim_result$par[["inpshift"]]
  shifted <- shift_timings_long(t_tac, tac, region, input, optimal_inpshift)

  par_list <- list()
  par_se_list <- list()
  fitted_tacs <- list()
  inner_fits <- list()

  for (r in regions) {
    region_mask <- shifted$region == r
    region_t <- shifted$t_tac[region_mask]
    region_tac <- shifted$tac[region_mask]

    modeldata <- list(
      tac = region_tac, t_tac = region_t,
      weights = weights_by_region[[r]], input = shifted$input
    )
    fit <- .nested_fit_region(formula_str, modeldata, start, lower, upper,
                              multstart_iter, multstart_lower, multstart_upper)
    .nested_require_fit(fit, r, model_name)

    coefs <- as.list(coef(fit))
    par_list[[r]] <- data.frame(
      region = r,
      inpshift = optimal_inpshift,
      K1 = coefs$K1, k2 = coefs$k2,
      vB = vB,
      Vt = coefs$K1 / coefs$k2,
      stringsAsFactors = FALSE
    )

    par_se_list[[r]] <- data.frame(
      region = r,
      inpshift.se = NA_real_,
      K1.se = get_se(fit, "K1"),
      k2.se = get_se(fit, "k2"),
      Vt.se = get_se(fit, "K1/k2"),
      stringsAsFactors = FALSE
    )

    fitted_tacs[[r]] <- data.frame(
      Time = region_t, Region = r,
      Radioactivity = region_tac,
      Fitted = as.numeric(fitted(fit)),
      Weights = weights_by_region[[r]],
      stringsAsFactors = FALSE
    )

    if (keep_inner_fits) inner_fits[[r]] <- fit
  }

  par <- do.call(rbind, par_list)
  par.se <- do.call(rbind, par_se_list)
  tacs <- do.call(rbind, fitted_tacs)

  # Approximate profile-based SE for the shared delay
  par.se$inpshift.se <- .nested_outer_se(
    optim_result, outer_objective,
    n_obs = sum(weights > 0),
    n_par = 1 + length(start) * length(regions)
  )[["inpshift"]]

  rownames(par) <- NULL
  rownames(par.se) <- NULL
  rownames(tacs) <- NULL

  out <- list(
    par = par,
    par.se = par.se,
    outer_pars = "inpshift",
    fit = optim_result,
    tacs = tacs,
    input = shifted$input,
    weights = weights,
    roiweights = roiweights,
    vB = vB,
    model = "nested_1tcm_delay"
  )

  if (keep_inner_fits) {
    out$fits_inner <- tibble::tibble(
      region = regions,
      fit = inner_fits[regions]
    )
  }

  class(out) <- c("nested_1tcm_delay", "kinfit")
  return(out)
}


#' Nested 2TCM with Delay Estimation
#'
#' Fits the Two Tissue Compartment Model to multiple regions simultaneously,
#' estimating a shared input delay (inpshift) across all regions using optim()
#' in the outer loop, while fitting K1, k2, k3, k4 per region using nlsLM in
#' the inner loop.
#'
#' @inheritParams nested_1tcm_delay
#' @param k3.start Starting parameter for k3. Default is 0.1.
#' @param k3.lower Lower bound for k3. Default is 0.0001.
#' @param k3.upper Upper bound for k3. Default is 0.5.
#' @param k4.start Starting parameter for k4. Default is 0.1.
#' @param k4.lower Lower bound for k4. Default is 0.0001.
#' @param k4.upper Upper bound for k4. Default is 0.5.
#'
#' @details The recommended usage for delay estimation is to fit the first few
#'   minutes of the acquisition (5 minutes is usually a good default) and not
#'   to use frame weights, so as to prioritise the low radioactivity values
#'   where the delay signal is most apparent.
#'
#'   The standard errors in \code{par.se} follow the usual kinfitr
#'   convention of being expressed as a fraction of the estimate. Those for the
#'   shared (outer) parameter are derived from the curvature of the profile
#'   objective at its optimum, obtained by finite differences, and are therefore
#'   more approximate than those of the per-region parameters. The per-region
#'   standard errors are conditional on the shared parameter: they treat
#'   the fitted delay as known rather than estimated, and so are somewhat
#'   optimistic.
#'
#' @return A list with class c("nested_2tcm_delay", "kinfit").
#'   \code{par} and \code{par.se} have one row per region, with the shared
#'   parameter values repeated. \code{weights} holds the full stacked vector
#'   of frame weights, aligned with the rows of \code{tacs}.
#'
#' @examples
#' \dontrun{
#' data("pbr28")
#' tac_wide <- pbr28$tacs[[2]]
#' t_tac <- tac_wide$Times / 60
#' input <- pbr28$input[[2]]
#'
#' regions <- c("FC", "TC", "THA")
#' long_data <- do.call(rbind, lapply(regions, function(r) {
#'   data.frame(t_tac = t_tac, tac = tac_wide[[r]], region = r)
#' }))
#'
#' out <- nested_2tcm_delay(
#'   long_data$t_tac, long_data$tac, long_data$region, input,
#'   timeStartEnd = c(0, 5)
#' )
#' plot(out)
#' }
#'
#' @export
nested_2tcm_delay <- function(
    t_tac, tac, region, input,
    vB = 0.05,
    weights = NULL, roiweights = NULL,
    inpshift.start = 0, inpshift.lower = -0.5, inpshift.upper = 0.5,
    K1.start = 0.1, K1.lower = 0.0001, K1.upper = 1,
    k2.start = 0.1, k2.lower = 0.0001, k2.upper = 0.5,
    k3.start = 0.1, k3.lower = 0.0001, k3.upper = 0.5,
    k4.start = 0.1, k4.lower = 0.0001, k4.upper = 0.5,
    optim_method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN"),
    optim_control = list(),
    multstart_iter = 1,
    multstart_lower = NULL, multstart_upper = NULL,
    frameStartEnd = NULL, timeStartEnd = NULL,
    keep_inner_fits = FALSE) {

  optim_method <- match.arg(optim_method)

  model_name <- "nested_2tcm_delay"
  unnested_name <- "twotcm"

  region <- as.character(region)
  regions_raw <- .nested_regions(region, model_name, unnested_name)

  # Resolve roiweights against the original (pre-tidy) region vector, since a
  # per-observation roiweights vector is sized to the input TACs. tidyinput_long
  # may prepend a zero frame per region, which would otherwise make it too short.
  roiweights <- .nested_roiweights(roiweights, region, regions_raw)

  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    t_first <- t_tac[region == regions_raw[1]]
    frameStartEnd <- c(which(t_first >= timeStartEnd[1])[1],
                       tail(which(t_first <= timeStartEnd[2]), 1))
  }

  tidied <- tidyinput_long(t_tac, tac, region, weights, frameStartEnd)
  t_tac   <- tidied$t_tac
  tac     <- tidied$tac
  region  <- tidied$region
  weights <- tidied$weights
  regions <- unique(region)

  weights_by_region <- .nested_weights_by_region(weights, region, regions)

  start <- c(K1 = K1.start, k2 = k2.start, k3 = k3.start, k4 = k4.start)
  lower <- c(K1 = K1.lower, k2 = k2.lower, k3 = k3.lower, k4 = k4.lower)
  upper <- c(K1 = K1.upper, k2 = k2.upper, k3 = k3.upper, k4 = k4.upper)

  multstart_bounds <- .nested_multstart_bounds(start, lower, upper,
                                               multstart_iter,
                                               multstart_lower, multstart_upper)
  multstart_lower <- multstart_bounds$lower
  multstart_upper <- multstart_bounds$upper

  formula_str <- paste0("tac ~ twotcm_model(t_tac, input, K1, k2, k3, k4, vB=", vB, ")")

  outer_objective <- function(outer_vals) {
    inpshift_val <- outer_vals[["inpshift"]]
    shifted <- shift_timings_long(t_tac, tac, region, input, inpshift_val)

    total_rss <- 0
    for (r in regions) {
      region_mask <- shifted$region == r
      modeldata <- list(
        tac = shifted$tac[region_mask],
        t_tac = shifted$t_tac[region_mask],
        weights = weights_by_region[[r]],
        input = shifted$input
      )
      fit <- .nested_fit_region(formula_str, modeldata, start, lower, upper,
                                multstart_iter, multstart_lower, multstart_upper)
      if (is.null(fit)) {
        total_rss <- total_rss + 1e10 * unname(roiweights[r])
      } else {
        total_rss <- total_rss +
          sum(weights(fit) * residuals(fit)^2) * unname(roiweights[r])
      }
    }
    return(total_rss)
  }

  optim_args <- list(
    par = c(inpshift = inpshift.start),
    fn = outer_objective, method = optim_method, control = optim_control
  )
  if (optim_method == "L-BFGS-B") {
    optim_args$lower <- c(inpshift = inpshift.lower)
    optim_args$upper <- c(inpshift = inpshift.upper)
  }
  optim_result <- do.call(optim, optim_args)

  # Final pass
  optimal_inpshift <- optim_result$par[["inpshift"]]
  shifted <- shift_timings_long(t_tac, tac, region, input, optimal_inpshift)

  par_list <- list()
  par_se_list <- list()
  fitted_tacs <- list()
  inner_fits <- list()

  for (r in regions) {
    region_mask <- shifted$region == r
    region_t <- shifted$t_tac[region_mask]
    region_tac <- shifted$tac[region_mask]

    modeldata <- list(
      tac = region_tac, t_tac = region_t,
      weights = weights_by_region[[r]], input = shifted$input
    )
    fit <- .nested_fit_region(formula_str, modeldata, start, lower, upper,
                              multstart_iter, multstart_lower, multstart_upper)
    .nested_require_fit(fit, r, model_name)

    coefs <- as.list(coef(fit))
    par_list[[r]] <- data.frame(
      region = r,
      inpshift = optimal_inpshift,
      K1 = coefs$K1, k2 = coefs$k2, k3 = coefs$k3, k4 = coefs$k4,
      vB = vB,
      Vt = (coefs$K1 / coefs$k2) * (1 + coefs$k3 / coefs$k4),
      Vnd = coefs$K1 / coefs$k2,
      BPnd = coefs$k3 / coefs$k4,
      BPp = (coefs$K1 / coefs$k2) * (coefs$k3 / coefs$k4),
      stringsAsFactors = FALSE
    )

    par_se_list[[r]] <- data.frame(
      region = r,
      inpshift.se = NA_real_,
      K1.se = get_se(fit, "K1"),
      k2.se = get_se(fit, "k2"),
      k3.se = get_se(fit, "k3"),
      k4.se = get_se(fit, "k4"),
      Vt.se = get_se(fit, "(K1/k2) * (1 + k3/k4)"),
      Vnd.se = get_se(fit, "K1/k2"),
      BPnd.se = get_se(fit, "k3/k4"),
      BPp.se = get_se(fit, "(K1/k2) * (k3/k4)"),
      stringsAsFactors = FALSE
    )

    fitted_tacs[[r]] <- data.frame(
      Time = region_t, Region = r,
      Radioactivity = region_tac,
      Fitted = as.numeric(fitted(fit)),
      Weights = weights_by_region[[r]],
      stringsAsFactors = FALSE
    )

    if (keep_inner_fits) inner_fits[[r]] <- fit
  }

  par <- do.call(rbind, par_list)
  par.se <- do.call(rbind, par_se_list)
  tacs <- do.call(rbind, fitted_tacs)

  # Approximate profile-based SE for the shared delay
  par.se$inpshift.se <- .nested_outer_se(
    optim_result, outer_objective,
    n_obs = sum(weights > 0),
    n_par = 1 + length(start) * length(regions)
  )[["inpshift"]]

  rownames(par) <- NULL
  rownames(par.se) <- NULL
  rownames(tacs) <- NULL

  out <- list(
    par = par,
    par.se = par.se,
    outer_pars = "inpshift",
    fit = optim_result,
    tacs = tacs,
    input = shifted$input,
    weights = weights,
    roiweights = roiweights,
    vB = vB,
    model = "nested_2tcm_delay"
  )

  if (keep_inner_fits) {
    out$fits_inner <- tibble::tibble(
      region = regions,
      fit = inner_fits[regions]
    )
  }

  class(out) <- c("nested_2tcm_delay", "kinfit")
  return(out)
}


# --- Plot functions ---

#' Plot nested 1TCM delay fit
#'
#' @param nested_out Output from \code{nested_1tcm_delay}.
#' @param roiname Optional. Not used (for API consistency).
#'
#' @return A ggplot2 object.
#'
#' @export
plot_nested_1tcm_delayfit <- function(nested_out, roiname = NULL) {
  tacs <- nested_out$tacs
  regions <- unique(tacs$Region)
  fine_t <- nested_out$input$Time

  fitted_pieces <- lapply(regions, function(r) {
    rpar <- nested_out$par[nested_out$par$region == r, ]
    pred <- onetcm_model(fine_t, nested_out$input,
                         rpar$K1, rpar$k2, rpar$vB)
    data.frame(Time = fine_t, Radioactivity = pred, Region = r,
               stringsAsFactors = FALSE)
  })
  fitteddf <- do.call(rbind, fitted_pieces)

  .plot_nested_delayfit(nested_out, fitteddf)
}


#' Plot nested 2TCM delay fit
#'
#' @param nested_out Output from \code{nested_2tcm_delay}.
#' @param roiname Optional. Not used (for API consistency).
#'
#' @return A ggplot2 object.
#'
#' @export
plot_nested_2tcm_delayfit <- function(nested_out, roiname = NULL) {
  tacs <- nested_out$tacs
  regions <- unique(tacs$Region)
  fine_t <- nested_out$input$Time

  fitted_pieces <- lapply(regions, function(r) {
    rpar <- nested_out$par[nested_out$par$region == r, ]
    pred <- twotcm_model(fine_t, nested_out$input,
                         rpar$K1, rpar$k2, rpar$k3, rpar$k4, rpar$vB)
    data.frame(Time = fine_t, Radioactivity = pred, Region = r,
               stringsAsFactors = FALSE)
  })
  fitteddf <- do.call(rbind, fitted_pieces)

  .plot_nested_delayfit(nested_out, fitteddf)
}


#' Internal: shared delay plot
#'
#' @param nested_out The nested delay fit output.
#' @param fitteddf Data frame with fine-grid fitted values (Time, Radioactivity,
#'   Region).
#'
#' @return A ggplot2 object.
#' @keywords internal
.plot_nested_delayfit <- function(nested_out, fitteddf) {
  tacs <- nested_out$tacs
  regions <- unique(tacs$Region)

  measureddf <- data.frame(
    Time = tacs$Time,
    Radioactivity = tacs$Radioactivity,
    Weights = tacs$Weights,
    Region = tacs$Region,
    stringsAsFactors = FALSE
  )

  inputdf <- data.frame(
    Time = nested_out$input$Time,
    Radioactivity = nested_out$input$AIF,
    stringsAsFactors = FALSE
  )

  # Region color palette
  n_regions <- length(regions)
  if (n_regions <= 8) {
    region_colors <- RColorBrewer::brewer.pal(max(3, n_regions), "Set2")[seq_len(n_regions)]
  } else {
    region_colors <- grDevices::hcl.colors(n_regions, palette = "Set 2")
  }
  names(region_colors) <- regions
  colScale <- ggplot2::scale_colour_manual(name = "Region", values = region_colors)

  outplot <- ggplot2::ggplot() +
    # AIF: fixed red dashed line, outside legend
    ggplot2::geom_line(data = inputdf,
                       ggplot2::aes(x = Time, y = Radioactivity),
                       colour = "#E41A1C", linetype = "dashed") +
    # Measured points: colored by region
    ggplot2::geom_point(data = measureddf,
                        ggplot2::aes(x = Time, y = Radioactivity,
                                     colour = Region, shape = "a",
                                     size = Weights)) +
    # Fitted lines: colored by region, interpolated
    ggplot2::geom_line(data = fitteddf,
                       ggplot2::aes(x = Time, y = Radioactivity,
                                    colour = Region)) +
    colScale +
    ggplot2::guides(shape = "none",
                    color = ggplot2::guide_legend(order = 1)) +
    ggplot2::scale_size(range = c(1, 3)) +
    ggplot2::coord_cartesian(ylim = c(0, max(measureddf$Radioactivity) * 1.5))

  return(outplot)
}
