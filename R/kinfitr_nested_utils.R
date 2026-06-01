#' Fit a single region using nlsLM or nls_multstart
#'
#' Internal helper used by all nested model functions. Fits a nonlinear model
#' to a single region's data using either \code{nlsLM} (single start) or
#' \code{nls_multstart} (multiple random starts).
#'
#' @param formula_str Character string of the nlsLM formula.
#' @param modeldata Named list with t_tac, tac (or roitac), weights, and
#'   input (or reftac) as appropriate for the model.
#' @param start Named numeric vector of starting values.
#' @param lower Named numeric vector of lower bounds.
#' @param upper Named numeric vector of upper bounds.
#' @param multstart_iter Numeric. If > 1, use nls_multstart.
#' @param multstart_lower Named numeric vector of multstart lower bounds.
#' @param multstart_upper Named numeric vector of multstart upper bounds.
#'
#' @return The nls fit object, or NULL if fitting fails.
#' @keywords internal
.nested_fit_region <- function(formula_str, modeldata, start, lower, upper,
                               multstart_iter = 1,
                               multstart_lower = NULL,
                               multstart_upper = NULL) {
  tryCatch({
    if (prod(multstart_iter) == 1) {
      minpack.lm::nlsLM(
        as.formula(formula_str),
        data = modeldata,
        start = start, lower = lower, upper = upper,
        weights = weights,
        control = minpack.lm::nls.lm.control(maxiter = 200)
      )
    } else {
      nls.multstart::nls_multstart(
        as.formula(formula_str),
        data = modeldata,
        supp_errors = "Y",
        start_lower = multstart_lower,
        start_upper = multstart_upper,
        iter = multstart_iter, convergence_count = FALSE,
        lower = lower, upper = upper,
        modelweights = weights
      )
    }
  }, error = function(e) NULL)
}


#' Process ROI weights for nested models
#'
#' Handles roiweights provided as a short vector (one per region), a named
#' vector, or a long vector (same length as the stacked tac/roitac). Returns
#' a named, normalised vector with one weight per region.
#'
#' @param roiweights NULL, or a numeric vector of ROI weights.
#' @param region Character vector of region labels (long format).
#' @param regions Character vector of unique region names.
#'
#' @return Named numeric vector of length \code{length(regions)}, normalised
#'   so that the maximum value is 1.
#' @keywords internal
.nested_roiweights <- function(roiweights, region, regions) {
  if (is.null(roiweights)) {
    roiweights <- rep(1, length(regions))
    names(roiweights) <- regions
  } else if (length(roiweights) == length(region)) {
    roiweights <- vapply(regions, function(r) roiweights[region == r][1],
                         numeric(1))
  } else if (!is.null(names(roiweights))) {
    roiweights <- roiweights[regions]
  } else if (length(roiweights) == length(regions)) {
    names(roiweights) <- regions
  } else {
    stop("roiweights must be NULL, a named vector, a vector of length ",
         length(regions), " (one per region), or length ", length(region),
         " (one per observation)")
  }
  roiweights / max(roiweights)
}


#' Paginate faceted plots for nested models
#'
#' If more than \code{max_facets} regions, splits into multiple ggplot objects.
#'
#' @param tacs Data frame with Time, Region, Radioactivity, Fitted (and
#'   optionally Reference) columns.
#' @param plot_fn A function that takes a subset of tacs and returns a ggplot.
#' @param max_facets Maximum facets per plot. Default is 3 (single row).
#'
#' @return A single ggplot if <= max_facets regions, otherwise a list of
#'   ggplot objects.
#' @keywords internal
.nested_paginate_facets <- function(tacs, plot_fn, max_facets = 3) {
  regions <- unique(tacs$Region)

  if (length(regions) <= max_facets) {
    return(plot_fn(tacs))
  }

  region_groups <- split(regions, ceiling(seq_along(regions) / max_facets))

  plots <- lapply(region_groups, function(grp) {
    plot_fn(tacs[tacs$Region %in% grp, ])
  })

  return(plots)
}


# --- predict S3 methods for nested models ---

#' Predict from nested 1TCM delay fit
#'
#' @param object Output from \code{nested_1tcm_delay}.
#' @param newdata Optional data.frame with \code{t_tac} and \code{region}
#'   columns. If NULL, returns fitted values from the original data.
#' @param ... Not used.
#'
#' @return Numeric vector of predicted values.
#'
#' @export
predict.nested_1tcm_delay <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$tacs$Fitted)

  preds <- numeric(nrow(newdata))
  for (r in unique(newdata$region)) {
    mask <- newdata$region == r
    rpar <- object$par[object$par$region == r, ]
    preds[mask] <- onetcm_model(
      newdata$t_tac[mask], object$input,
      rpar$K1, rpar$k2, rpar$vB
    )
  }
  preds
}


#' Predict from nested 2TCM delay fit
#'
#' @param object Output from \code{nested_2tcm_delay}.
#' @param newdata Optional data.frame with \code{t_tac} and \code{region}
#'   columns. If NULL, returns fitted values from the original data.
#' @param ... Not used.
#'
#' @return Numeric vector of predicted values.
#'
#' @export
predict.nested_2tcm_delay <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$tacs$Fitted)

  preds <- numeric(nrow(newdata))
  for (r in unique(newdata$region)) {
    mask <- newdata$region == r
    rpar <- object$par[object$par$region == r, ]
    preds[mask] <- twotcm_model(
      newdata$t_tac[mask], object$input,
      rpar$K1, rpar$k2, rpar$k3, rpar$k4, rpar$vB
    )
  }
  preds
}


#' Predict from nested 2TCM fit
#'
#' @param object Output from \code{nested_2tcm}.
#' @param newdata Optional data.frame with \code{t_tac} and \code{region}
#'   columns. If NULL, returns fitted values from the original data.
#' @param ... Not used.
#'
#' @return Numeric vector of predicted values.
#'
#' @export
predict.nested_2tcm <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$tacs$Fitted)

  preds <- numeric(nrow(newdata))
  for (r in unique(newdata$region)) {
    mask <- newdata$region == r
    rpar <- object$par[object$par$region == r, ]
    preds[mask] <- twotcm_macro_model(
      newdata$t_tac[mask], object$input,
      rpar$K1, rpar$Vnd, rpar$BPp, rpar$k4, rpar$vB
    )
  }
  preds
}


#' Predict from nested SRTM fit
#'
#' @param object Output from \code{nested_srtm}.
#' @param newdata Optional data.frame with \code{t_tac}, \code{region}, and
#'   \code{reftac} columns. If NULL, returns fitted values from the original
#'   data.
#' @param ... Not used.
#'
#' @return Numeric vector of predicted values.
#'
#' @export
predict.nested_srtm <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$tacs$Fitted)

  preds <- numeric(nrow(newdata))
  for (r in unique(newdata$region)) {
    mask <- newdata$region == r
    rpar <- object$par[object$par$region == r, ]
    preds[mask] <- srtm2_model(
      newdata$t_tac[mask], newdata$reftac[mask],
      rpar$R1, rpar$k2prime, rpar$bp
    )
  }
  preds
}
