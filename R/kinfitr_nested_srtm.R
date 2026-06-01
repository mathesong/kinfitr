#' Nested SRTM with Shared k2prime Estimation
#'
#' Fits the Simplified Reference Tissue Model 2 (SRTM2) to multiple regions
#' simultaneously, estimating a shared k2prime across all regions using optim()
#' in the outer loop, while fitting R1 and bp per region using nlsLM in the
#' inner loop. This is the proper way to estimate k2prime for SRTM2 — by
#' leveraging information across multiple regions.
#'
#' @param t_tac Numeric vector of times (minutes), repeated per region.
#' @param roitac Numeric vector of target ROI radioactivity, stacked per region.
#' @param reftac Numeric vector of reference tissue radioactivity. Should be the
#'   same reference region values repeated for each target region (i.e., same
#'   length as roitac). Alternatively, can be a single vector of per-frame
#'   values that will be repeated for each region.
#' @param region Character vector identifying the region for each observation.
#' @param weights Optional. Numeric vector of weights.
#' @param roiweights Optional. Numeric vector of ROI weights. Can be a named
#'   vector (one per region), an unnamed vector of the same length as the
#'   number of regions, or a long vector the same length as \code{roitac}
#'   (one per observation; the first value per region is used). ROI size or
#'   relative size is a useful quantity to supply here, as larger regions
#'   provide more precise mean TAC estimates. Normalised internally so the
#'   maximum is 1.
#' @param k2prime.start Starting value for k2prime. Default is 0.1.
#' @param k2prime.lower Lower bound for k2prime. Default is 0.001.
#' @param k2prime.upper Upper bound for k2prime. Default is 1.
#' @param R1.start Starting parameter for R1. Default is 1.
#' @param R1.lower Lower bound for R1. Default is 0.
#' @param R1.upper Upper bound for R1. Default is 10.
#' @param bp.start Starting parameter for bp. Default is 1.5.
#' @param bp.lower Lower bound for bp. Default is 0.
#' @param bp.upper Upper bound for bp. Default is 15.
#' @param optim_method Optimization method. Default is "L-BFGS-B".
#'   Only "L-BFGS-B" supports bounds on the shared (outer) parameters; if
#'   another method is selected, the lower and upper limits for shared
#'   parameters are ignored.
#' @param optim_control Control parameters for optim().
#' @param multstart_iter Multistart iterations for inner fits. Default is 1.
#' @param multstart_lower Optional. Lower bounds for multistart.
#' @param multstart_upper Optional. Upper bounds for multistart.
#' @param frameStartEnd Optional. Frame range c(start, end).
#' @param timeStartEnd Optional. Time range for frame selection.
#' @param keep_inner_fits Logical. If TRUE, store per-region nlsLM fit objects
#'   in \code{$fits_inner}. Default is FALSE.
#'
#' @return A list with class c("nested_srtm", "kinfit").
#'
#' @export
nested_srtm <- function(
    t_tac, roitac, reftac, region,
    weights = NULL, roiweights = NULL,
    k2prime.start = 0.1, k2prime.lower = 0.001, k2prime.upper = 1,
    R1.start = 1, R1.lower = 0, R1.upper = 10,
    bp.start = 1.5, bp.lower = 0, bp.upper = 15,
    optim_method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN"),
    optim_control = list(),
    multstart_iter = 1,
    multstart_lower = NULL, multstart_upper = NULL,
    frameStartEnd = NULL, timeStartEnd = NULL,
    keep_inner_fits = FALSE) {

  optim_method <- match.arg(optim_method)

  # Convert timeStartEnd to frameStartEnd
  region <- as.character(region)
  regions_raw <- unique(region)

  # Resolve roiweights against the original (pre-tidy) region vector, since a
  # per-observation roiweights vector is sized to the input TACs. tidyinput_long
  # may prepend a zero frame per region, which would otherwise make it too short.
  roiweights <- .nested_roiweights(roiweights, region, regions_raw)

  # Handle reftac: if single vector, repeat for each region
  n_per_region <- length(roitac) / length(regions_raw)
  if (length(reftac) == n_per_region && length(reftac) != length(roitac)) {
    reftac <- rep(reftac, times = length(regions_raw))
  }

  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    t_first <- t_tac[region == regions_raw[1]]
    frameStartEnd <- c(which(t_first >= timeStartEnd[1])[1],
                       tail(which(t_first <= timeStartEnd[2]), 1))
  }

  # Tidy roitac and reftac together using the same transformations
  tidied_roi <- tidyinput_long(t_tac, roitac, region, weights, frameStartEnd)
  tidied_ref <- tidyinput_long(t_tac, reftac, region, weights = NULL, frameStartEnd)
  t_tac   <- tidied_roi$t_tac
  roitac  <- tidied_roi$tac
  reftac  <- tidied_ref$tac
  region  <- tidied_roi$region
  weights <- tidied_roi$weights
  regions <- unique(region)

  weights_per_region <- weights[region == regions[1]]

  # Inner fit parameters
  start <- c(R1 = R1.start, bp = bp.start)
  lower <- c(R1 = R1.lower, bp = bp.lower)
  upper <- c(R1 = R1.upper, bp = bp.upper)

  if (multstart_iter > 1) {
    if (is.null(multstart_lower)) multstart_lower <- lower
    if (is.null(multstart_upper)) multstart_upper <- upper
  }

  # Outer objective function
  outer_objective <- function(outer_vals) {
    k2prime_val <- outer_vals[["k2prime"]]

    formula_str <- paste0("roitac ~ srtm2_model(t_tac, reftac, R1, k2prime=",
                          k2prime_val, ", bp)")

    total_rss <- 0
    for (r in regions) {
      region_mask <- region == r
      modeldata <- list(
        roitac = roitac[region_mask],
        reftac = reftac[region_mask],
        t_tac = t_tac[region_mask],
        weights = weights_per_region
      )

      fit <- .nested_fit_region(formula_str, modeldata, start, lower, upper,
                                multstart_iter, multstart_lower, multstart_upper)

      if (is.null(fit)) {
        total_rss <- total_rss + 1e10 * roiweights[r]
      } else {
        total_rss <- total_rss + sum(weights(fit) * residuals(fit)^2) * roiweights[r]
      }
    }
    return(total_rss)
  }

  # Run optim
  optim_args <- list(
    par = c(k2prime = k2prime.start),
    fn = outer_objective,
    method = optim_method,
    control = optim_control
  )
  if (optim_method == "L-BFGS-B") {
    optim_args$lower <- c(k2prime = k2prime.lower)
    optim_args$upper <- c(k2prime = k2prime.upper)
  }

  optim_result <- do.call(optim, optim_args)

  # Final pass at optimal parameters
  optimal_k2prime <- optim_result$par[["k2prime"]]
  formula_str <- paste0("roitac ~ srtm2_model(t_tac, reftac, R1, k2prime=",
                        optimal_k2prime, ", bp)")

  par_list <- list()
  par_se_list <- list()
  fitted_tacs <- list()
  inner_fits <- list()

  for (r in regions) {
    region_mask <- region == r
    region_t <- t_tac[region_mask]
    region_roitac <- roitac[region_mask]
    region_reftac <- reftac[region_mask]

    modeldata <- list(
      roitac = region_roitac,
      reftac = region_reftac,
      t_tac = region_t,
      weights = weights_per_region
    )

    fit <- .nested_fit_region(formula_str, modeldata, start, lower, upper,
                              multstart_iter, multstart_lower, multstart_upper)

    coefs <- as.list(coef(fit))
    k2a_val <- (coefs$R1 * optimal_k2prime) / (coefs$bp + 1)

    par_list[[r]] <- data.frame(
      region = r,
      k2prime = optimal_k2prime,
      R1 = coefs$R1,
      bp = coefs$bp,
      k2a = k2a_val,
      stringsAsFactors = FALSE
    )

    k2prime_str <- format(optimal_k2prime, digits = 10)
    par_se_list[[r]] <- data.frame(
      region = r,
      R1.se = get_se(fit, "R1"),
      bp.se = get_se(fit, "bp"),
      k2a.se = get_se(fit, paste0("(R1 * ", k2prime_str, ") / (bp + 1)")),
      stringsAsFactors = FALSE
    )

    fitted_tacs[[r]] <- data.frame(
      Time = region_t,
      Region = r,
      Radioactivity = region_roitac,
      Reference = region_reftac,
      Fitted = as.numeric(fitted(fit)),
      Weights = weights_per_region,
      stringsAsFactors = FALSE
    )

    if (keep_inner_fits) inner_fits[[r]] <- fit
  }

  par <- do.call(rbind, par_list)
  par.se <- do.call(rbind, par_se_list)
  tacs <- do.call(rbind, fitted_tacs)
  rownames(par) <- NULL
  rownames(par.se) <- NULL
  rownames(tacs) <- NULL

  out <- list(
    par = par,
    par.se = par.se,
    outer_pars = "k2prime",
    fit = optim_result,
    tacs = tacs,
    reftac = reftac[region == regions[1]],
    weights = weights_per_region,
    roiweights = roiweights,
    model = "nested_srtm"
  )

  if (keep_inner_fits) {
    out$fits_inner <- tibble::tibble(
      region = regions,
      fit = inner_fits[regions]
    )
  }

  class(out) <- c("nested_srtm", "kinfit")
  return(out)
}


#' Plot nested SRTM fit (faceted)
#'
#' Faceted plot with one panel per region. Each panel shows measured target ROI
#' points (sized by weight), fitted line, and the shared reference TAC,
#' following the conventions of the unnested SRTM2 plot.
#'
#' @param nested_out Output from \code{nested_srtm}.
#' @param roiname Optional. Not used (for API consistency).
#'
#' @return A single ggplot2 object if <= 3 regions, otherwise a list of
#'   ggplot objects with up to 3 facets each.
#'
#' @export
plot_nested_srtmfit <- function(nested_out, roiname = NULL) {
  tacs <- nested_out$tacs

  make_plot <- function(tacs_subset) {
    regions_sub <- unique(tacs_subset$Region)

    plot_pieces <- lapply(regions_sub, function(r) {
      rmask <- tacs_subset$Region == r

      # Reference first so fct_inorder gives it Set1[1] = red
      ref_r <- data.frame(
        Time = tacs_subset$Time[rmask],
        Radioactivity = tacs_subset$Reference[rmask],
        Weights = 1,
        Type = "Reference",
        Facet = r,
        stringsAsFactors = FALSE
      )

      measured <- data.frame(
        Time = tacs_subset$Time[rmask],
        Radioactivity = tacs_subset$Radioactivity[rmask],
        Weights = tacs_subset$Weights[rmask],
        Type = "ROI.measured",
        Facet = r,
        stringsAsFactors = FALSE
      )

      fitted_r <- data.frame(
        Time = tacs_subset$Time[rmask],
        Radioactivity = tacs_subset$Fitted[rmask],
        Weights = 1,
        Type = "ROI.fitted",
        Facet = r,
        stringsAsFactors = FALSE
      )

      rbind(ref_r, measured, fitted_r)
    })

    plotdf <- do.call(rbind, plot_pieces)
    plotdf$Type <- forcats::fct_inorder(factor(plotdf$Type))

    myColors <- RColorBrewer::brewer.pal(3, "Set1")
    names(myColors) <- levels(plotdf$Type)
    colScale <- ggplot2::scale_colour_manual(name = "Region", values = myColors)

    max_measured <- max(tacs_subset$Radioactivity)

    outplot <- ggplot2::ggplot(plotdf, ggplot2::aes(
      x = Time, y = Radioactivity, colour = Type
    )) +
      colScale +
      ggplot2::geom_point(
        data = plotdf[plotdf$Type == "ROI.measured", ],
        ggplot2::aes(shape = "a", size = Weights)
      ) +
      ggplot2::geom_line(
        data = plotdf[plotdf$Type != "ROI.measured", ]
      ) +
      ggplot2::facet_wrap(~Facet) +
      ggplot2::guides(shape = "none",
                      color = ggplot2::guide_legend(order = 1)) +
      ggplot2::scale_size(range = c(1, 3)) +
      ggplot2::coord_cartesian(ylim = c(0, max_measured)) +
      ggplot2::labs(title = bquote("Nested SRTM (Shared" ~ k[2]*"'" * ")"))

    return(outplot)
  }

  .nested_paginate_facets(tacs, make_plot, max_facets = 3)
}
