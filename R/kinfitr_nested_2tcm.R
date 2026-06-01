#' Nested 2TCM with Shared Macroparameters
#'
#' Fits the Two Tissue Compartment Model (macro parameterisation) to multiple
#' regions simultaneously. Optimizes shared parameters (Vnd and/or k4) across
#' all regions using optim() in the outer loop, while fitting per-region
#' parameters using nlsLM in the inner loop. Uses
#' \code{twotcm_macro_model(t_tac, input, K1, Vnd, BPp, k4, vB)}.
#'
#' The \code{shared} argument controls which parameters are estimated jointly
#' across regions in the outer loop:
#' \itemize{
#'   \item \code{"Vnd"}: Vnd shared across regions. Inner fits K1, BPp, k4
#'     per region.
#'   \item \code{"k4"}: k4 shared across regions. Inner fits K1, Vnd, BPp
#'     per region.
#'   \item \code{"Vnd_k4"}: Both Vnd and k4 shared. Inner fits K1, BPp
#'     per region.
#' }
#'
#' @param t_tac Numeric vector of times (minutes), repeated per region.
#' @param tac Numeric vector of radioactivity, stacked per region.
#' @param region Character vector identifying regions.
#' @param input Data frame from blood_interp().
#' @param shared Character. Which parameters to share across regions in the
#'   outer loop. One of \code{"Vnd"}, \code{"k4"}, or \code{"Vnd_k4"}.
#' @param inpshift The delay in minutes. Must be pre-fitted (not estimated
#'   here). Default is 0.
#' @param vB Blood volume fraction. Default is 0.05.
#' @param weights Optional. Numeric vector of weights.
#' @param roiweights Optional. Numeric vector of ROI weights. Can be a named
#'   vector (one per region), an unnamed vector of the same length as the
#'   number of regions, or a long vector the same length as \code{tac}
#'   (one per observation; the first value per region is used). ROI size or
#'   relative size is a useful quantity to supply here, as larger regions
#'   provide more precise mean TAC estimates. Normalised internally so the
#'   maximum is 1.
#' @param Vnd.start Starting value for Vnd. Default is 1.
#' @param Vnd.lower Lower bound for Vnd. Default is 0.0001.
#' @param Vnd.upper Upper bound for Vnd. Default is 10.
#' @param k4.start Starting value for k4. Default is 0.1.
#' @param k4.lower Lower bound for k4. Default is 0.0001.
#' @param k4.upper Upper bound for k4. Default is 0.5.
#' @param K1.start Starting parameter for K1. Default is 0.1.
#' @param K1.lower Lower bound for K1. Default is 0.0001.
#' @param K1.upper Upper bound for K1. Default is 1.
#' @param BPp.start Starting parameter for BPp. Default is 1.
#' @param BPp.lower Lower bound for BPp. Default is 0.0001.
#' @param BPp.upper Upper bound for BPp. Default is 50.
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
#' @return A list with class c("nested_2tcm", "kinfit").
#'
#' @export
nested_2tcm <- function(
    t_tac, tac, region, input,
    shared = c("Vnd", "k4", "Vnd_k4"),
    inpshift = 0, vB = 0.05,
    weights = NULL, roiweights = NULL,
    Vnd.start = 1, Vnd.lower = 0.0001, Vnd.upper = 10,
    k4.start = 0.1, k4.lower = 0.0001, k4.upper = 0.5,
    K1.start = 0.1, K1.lower = 0.0001, K1.upper = 1,
    BPp.start = 1, BPp.lower = 0.0001, BPp.upper = 50,
    optim_method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN"),
    optim_control = list(),
    multstart_iter = 1,
    multstart_lower = NULL, multstart_upper = NULL,
    frameStartEnd = NULL, timeStartEnd = NULL,
    keep_inner_fits = FALSE) {

  shared <- match.arg(shared)
  optim_method <- match.arg(optim_method)

  # Convert timeStartEnd to frameStartEnd
  region <- as.character(region)
  regions_raw <- unique(region)

  # Resolve roiweights against the original (pre-tidy) region vector, since a
  # per-observation roiweights vector is sized to the input TACs. tidyinput_long
  # may prepend a zero frame per region, which would otherwise make it too short.
  roiweights <- .nested_roiweights(roiweights, region, regions_raw)

  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    t_first <- t_tac[region == regions_raw[1]]
    frameStartEnd <- c(which(t_first >= timeStartEnd[1])[1],
                       tail(which(t_first <= timeStartEnd[2]), 1))
  }

  # Tidy input
  tidied <- tidyinput_long(t_tac, tac, region, weights, frameStartEnd)
  t_tac   <- tidied$t_tac
  tac     <- tidied$tac
  region  <- tidied$region
  weights <- tidied$weights
  regions <- unique(region)

  weights_per_region <- weights[region == regions[1]]

  # Shift timings (inpshift is pre-fitted)
  shifted <- shift_timings_long(t_tac, tac, region, input, inpshift)

  # --- Determine outer and inner parameters based on shared mode ---

  # Outer parameters (optimized by optim)
  outer_start <- c()
  outer_lower <- c()
  outer_upper <- c()

  if (shared %in% c("Vnd", "Vnd_k4")) {
    outer_start <- c(outer_start, Vnd = Vnd.start)
    outer_lower <- c(outer_lower, Vnd = Vnd.lower)
    outer_upper <- c(outer_upper, Vnd = Vnd.upper)
  }
  if (shared %in% c("k4", "Vnd_k4")) {
    outer_start <- c(outer_start, k4 = k4.start)
    outer_lower <- c(outer_lower, k4 = k4.lower)
    outer_upper <- c(outer_upper, k4 = k4.upper)
  }

  # Inner parameters depend on shared mode
  build_inner_config <- function() {
    inner_start <- c(K1 = K1.start, BPp = BPp.start)
    inner_lower <- c(K1 = K1.lower, BPp = BPp.lower)
    inner_upper <- c(K1 = K1.upper, BPp = BPp.upper)

    if (shared == "Vnd") {
      inner_start <- c(inner_start, k4 = k4.start)
      inner_lower <- c(inner_lower, k4 = k4.lower)
      inner_upper <- c(inner_upper, k4 = k4.upper)
    } else if (shared == "k4") {
      inner_start <- c(K1 = K1.start, Vnd = Vnd.start, BPp = BPp.start)
      inner_lower <- c(K1 = K1.lower, Vnd = Vnd.lower, BPp = BPp.lower)
      inner_upper <- c(K1 = K1.upper, Vnd = Vnd.upper, BPp = BPp.upper)
    }

    list(start = inner_start, lower = inner_lower, upper = inner_upper)
  }

  inner_config <- build_inner_config()

  if (multstart_iter > 1) {
    if (is.null(multstart_lower)) multstart_lower <- inner_config$lower
    if (is.null(multstart_upper)) multstart_upper <- inner_config$upper
  }

  # Build formula string with outer params fixed from current outer values
  build_formula <- function(Vnd_val, k4_val) {
    vnd_part <- if (shared %in% c("Vnd", "Vnd_k4")) {
      paste0("Vnd=", Vnd_val)
    } else {
      "Vnd"
    }
    k4_part <- if (shared %in% c("k4", "Vnd_k4")) {
      paste0("k4=", k4_val)
    } else {
      "k4"
    }

    paste0("tac ~ twotcm_macro_model(t_tac, input, K1, ",
           vnd_part, ", BPp, ", k4_part, ", vB=", vB, ")")
  }

  # Outer objective function
  outer_objective <- function(outer_vals) {
    Vnd_val <- if ("Vnd" %in% names(outer_vals)) outer_vals[["Vnd"]] else NULL
    k4_val  <- if ("k4"  %in% names(outer_vals)) outer_vals[["k4"]]  else NULL

    formula_str <- build_formula(Vnd_val, k4_val)

    total_rss <- 0
    for (r in regions) {
      region_mask <- shifted$region == r
      modeldata <- list(
        tac = shifted$tac[region_mask],
        t_tac = shifted$t_tac[region_mask],
        weights = weights_per_region,
        input = shifted$input
      )

      fit <- .nested_fit_region(formula_str, modeldata,
                                inner_config$start, inner_config$lower,
                                inner_config$upper,
                                multstart_iter, multstart_lower,
                                multstart_upper)

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
    par = outer_start,
    fn = outer_objective,
    method = optim_method,
    control = optim_control
  )
  if (optim_method == "L-BFGS-B") {
    optim_args$lower <- outer_lower
    optim_args$upper <- outer_upper
  }

  optim_result <- do.call(optim, optim_args)

  # Extract optimal outer values
  opt_Vnd <- if ("Vnd" %in% names(optim_result$par)) optim_result$par[["Vnd"]] else NULL
  opt_k4  <- if ("k4"  %in% names(optim_result$par)) optim_result$par[["k4"]]  else NULL

  formula_str <- build_formula(opt_Vnd, opt_k4)

  # Determine outer_pars character vector
  outer_pars <- names(optim_result$par)

  # Final pass
  par_list <- list()
  par_se_list <- list()
  fitted_tacs <- list()
  inner_fits <- list()

  for (r in regions) {
    region_mask <- shifted$region == r
    region_t <- shifted$t_tac[region_mask]
    region_tac <- shifted$tac[region_mask]

    modeldata <- list(
      tac = region_tac,
      t_tac = region_t,
      weights = weights_per_region,
      input = shifted$input
    )

    fit <- .nested_fit_region(formula_str, modeldata,
                              inner_config$start, inner_config$lower,
                              inner_config$upper,
                              multstart_iter, multstart_lower,
                              multstart_upper)

    coefs <- as.list(coef(fit))

    # Resolve all macro params (some from outer, some from inner)
    this_K1  <- coefs$K1
    this_BPp <- coefs$BPp
    this_Vnd <- if (!is.null(opt_Vnd) && shared %in% c("Vnd", "Vnd_k4")) opt_Vnd else coefs$Vnd
    this_k4  <- if (!is.null(opt_k4) && shared %in% c("k4", "Vnd_k4"))  opt_k4  else coefs$k4

    # Derived params
    this_Vt   <- this_Vnd + this_BPp
    this_BPnd <- this_BPp / this_Vnd
    this_k2   <- this_K1 / this_Vnd
    this_k3   <- this_BPnd * this_k4

    # Build par row: region, outer params (repeated), inner params, fixed, derived
    region_par <- data.frame(region = r, stringsAsFactors = FALSE)
    if (shared == "Vnd") {
      region_par$Vnd <- this_Vnd
      region_par$K1  <- this_K1
      region_par$BPp <- this_BPp
      region_par$k4  <- this_k4
    } else if (shared == "k4") {
      region_par$k4  <- this_k4
      region_par$K1  <- this_K1
      region_par$Vnd <- this_Vnd
      region_par$BPp <- this_BPp
    } else {
      region_par$Vnd <- this_Vnd
      region_par$k4  <- this_k4
      region_par$K1  <- this_K1
      region_par$BPp <- this_BPp
    }
    region_par$vB       <- vB
    region_par$inpshift <- inpshift
    region_par$Vt       <- this_Vt
    region_par$BPnd     <- this_BPnd
    region_par$k2       <- this_k2
    region_par$k3       <- this_k3

    # SEs
    inner_names <- names(coef(fit))
    se_vals <- list(region = r)

    for (pname in inner_names) {
      se_vals[[paste0(pname, ".se")]] <- get_se(fit, pname)
    }

    if (shared == "Vnd_k4") {
      se_vals$Vt.se   <- NA
      se_vals$BPnd.se <- NA
      se_vals$k2.se   <- NA
      se_vals$k3.se   <- NA
    } else if (shared == "Vnd") {
      se_vals$Vt.se   <- NA
      se_vals$BPnd.se <- NA
      se_vals$k2.se   <- NA
      se_vals$k3.se   <- NA
    } else if (shared == "k4") {
      se_vals$Vt.se   <- get_se(fit, "Vnd + BPp")
      se_vals$BPnd.se <- get_se(fit, "BPp / Vnd")
      se_vals$k2.se   <- get_se(fit, "K1 / Vnd")
      se_vals$k3.se   <- NA
    }

    region_se <- as.data.frame(se_vals, stringsAsFactors = FALSE)

    fitted_tacs[[r]] <- data.frame(
      Time = region_t,
      Region = r,
      Radioactivity = region_tac,
      Fitted = as.numeric(fitted(fit)),
      Weights = weights_per_region,
      stringsAsFactors = FALSE
    )

    par_list[[r]] <- region_par
    par_se_list[[r]] <- region_se
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
    outer_pars = outer_pars,
    fit = optim_result,
    tacs = tacs,
    input = shifted$input,
    weights = weights_per_region,
    roiweights = roiweights,
    vB = vB,
    inpshift = inpshift,
    shared = shared,
    model = "nested_2tcm"
  )

  if (keep_inner_fits) {
    out$fits_inner <- tibble::tibble(
      region = regions,
      fit = inner_fits[regions]
    )
  }

  class(out) <- c("nested_2tcm", "kinfit")
  return(out)
}


#' Plot nested 2TCM fit (faceted)
#'
#' Faceted plot with one panel per region. Each panel shows measured points
#' (sized by weight), fitted line (interpolated at fine time resolution), and
#' AIF, following the conventions of the unnested 2TCM plot.
#'
#' @param nested_out Output from \code{nested_2tcm}.
#' @param roiname Optional. Not used (for API consistency).
#'
#' @return A single ggplot2 object if <= 3 regions, otherwise a list of
#'   ggplot objects with up to 3 facets each.
#'
#' @export
plot_nested_2tcmfit <- function(nested_out, roiname = NULL) {
  tacs <- nested_out$tacs

  make_plot <- function(tacs_subset) {
    regions_sub <- unique(tacs_subset$Region)
    fine_t <- nested_out$input$Time

    plot_pieces <- lapply(regions_sub, function(r) {
      rmask <- tacs_subset$Region == r

      # AIF (rbind first so fct_inorder gives it Set1[1] = red)
      aif_r <- data.frame(
        Time = fine_t,
        Radioactivity = nested_out$input$AIF,
        Weights = 1,
        Type = "AIF",
        Facet = r,
        stringsAsFactors = FALSE
      )

      measured <- data.frame(
        Time = tacs_subset$Time[rmask],
        Radioactivity = tacs_subset$Radioactivity[rmask],
        Weights = tacs_subset$Weights[rmask],
        Type = "Measured",
        Facet = r,
        stringsAsFactors = FALSE
      )

      # Fine-grid interpolation via model function
      rpar <- nested_out$par[nested_out$par$region == r, ]
      pred <- twotcm_macro_model(fine_t, nested_out$input,
                                  rpar$K1, rpar$Vnd, rpar$BPp,
                                  rpar$k4, rpar$vB)

      fitted_r <- data.frame(
        Time = fine_t,
        Radioactivity = pred,
        Weights = 1,
        Type = "Fitted",
        Facet = r,
        stringsAsFactors = FALSE
      )

      rbind(aif_r, measured, fitted_r)
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
        data = plotdf[plotdf$Type == "Measured", ],
        ggplot2::aes(shape = "a", size = Weights)
      ) +
      ggplot2::geom_line(
        data = plotdf[plotdf$Type != "Measured", ]
      ) +
      ggplot2::facet_wrap(~Facet) +
      ggplot2::guides(shape = "none",
                      color = ggplot2::guide_legend(order = 1)) +
      ggplot2::scale_size(range = c(1, 3)) +
      ggplot2::coord_cartesian(ylim = c(0, max_measured * 1.5)) +
      ggplot2::labs(title = switch(nested_out$shared,
        "Vnd"    = bquote("Nested 2TCM (Shared" ~ V[ND] * ")"),
        "k4"     = bquote("Nested 2TCM (Shared" ~ k[4] * ")"),
        "Vnd_k4" = bquote("Nested 2TCM (Shared" ~ V[ND] ~ "and" ~ k[4] * ")")
      ))

    return(outplot)
  }

  .nested_paginate_facets(tacs, make_plot, max_facets = 3)
}
