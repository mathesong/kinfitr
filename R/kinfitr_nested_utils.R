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
#' @param multstart_iter Numeric. A single number, or one value per fitted
#'   parameter. If the product exceeds 1, \code{nls_multstart} is used.
#' @param multstart_lower Named numeric vector of multstart lower bounds.
#' @param multstart_upper Named numeric vector of multstart upper bounds.
#'
#' @return The nls fit object, or NULL if fitting fails.
#' @keywords internal
.nested_fit_region <- function(formula_str, modeldata, start, lower, upper,
                               multstart_iter = 1,
                               multstart_lower = NULL,
                               multstart_upper = NULL) {
  fit <- tryCatch({
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
        lhstype = "improved",
        lower = lower, upper = upper,
        modelweights = weights
      )
    }
  }, error = function(e) NULL)

  # Record the failure so that callers which cannot see it directly -- notably
  # .nested_outer_se(), whose objective values pass through optim() -- can tell
  # that a fit failed rather than only seeing the penalty it produced.
  if (is.null(fit)) {
    .nested_fit_failures$n <- .nested_fit_failures$n + 1L
  }

  fit
}


#' Tally of failed inner fits
#'
#' A running count of the inner fits that have failed, incremented by
#' \code{.nested_fit_region}. Only differences in the count are meaningful:
#' callers record it before a block of fitting and compare afterwards.
#'
#' @keywords internal
.nested_fit_failures <- new.env(parent = emptyenv())
.nested_fit_failures$n <- 0L


#' Stop informatively when an inner fit fails
#'
#' The outer optimisation tolerates failed inner fits by penalising them, but
#' a failure in the final pass (i.e. at the optimal shared parameters) leaves
#' nothing to report for that region. This gives a legible error rather than
#' allowing a downstream \code{coef(NULL)} failure.
#'
#' @param fit The inner fit object, or NULL.
#' @param region The region label being fitted.
#' @param model The name of the calling nested model.
#'
#' @return \code{fit}, invisibly, if it is not NULL.
#' @keywords internal
.nested_require_fit <- function(fit, region, model) {
  if (is.null(fit)) {
    stop(model, ": the inner fit for region '", region, "' did not converge ",
         "at the optimised shared parameter values. Consider adjusting the ",
         "starting values or bounds for the per-region parameters, or using ",
         "multstart_iter > 1.", call. = FALSE)
  }
  invisible(fit)
}


#' Validate region labels for nested models
#'
#' @param region Character vector of region labels.
#' @param model The name of the calling nested model.
#' @param alternative The unnested model to suggest for a single region.
#'
#' @return The unique region labels, in order of appearance. Errors if fewer
#'   than two regions are supplied.
#' @keywords internal
.nested_regions <- function(region, model, alternative) {
  if (anyNA(region) || any(!nzchar(region))) {
    stop(model, ": region labels must all be non-missing and non-empty.",
         call. = FALSE)
  }

  regions <- unique(region)

  if (length(regions) < 2) {
    stop(model, " requires at least two regions, but was given only one ('",
         regions[1], "'). Nesting works by borrowing information across ",
         "regions: with a single region there is nothing to share, and the ",
         "result would be recorded as a nested fit while being no such thing. ",
         "Use ", alternative, "() instead.", call. = FALSE)
  }

  regions
}


#' Split frame weights by region
#'
#' Nested models fit each region separately, so each region must be given its
#' own frame weights rather than those of the first region.
#'
#' @param weights Numeric vector of frame weights (long format).
#' @param region Character vector of region labels (long format).
#' @param regions Character vector of unique region names.
#'
#' @return Named list of numeric vectors, one per region, in \code{regions}
#'   order.
#'
#' @details These weights stay aligned with the TACs after the timings are
#'   shifted because the frame count per region does not change:
#'   \code{shift_timings_df} only prepends a row when there is no zero frame,
#'   and \code{tidyinput_long} has already added one by the time the nested
#'   models call it. The alignment therefore depends on that ordering, and on
#'   each region's frames arriving in ascending time order.
#' @keywords internal
.nested_weights_by_region <- function(weights, region, regions) {
  out <- lapply(regions, function(r) weights[region == r])
  names(out) <- regions
  out
}


#' Align multstart bounds with the fitted parameters
#'
#' @param bounds Numeric vector of multstart bounds, named or unnamed.
#' @param start Named numeric vector of starting values.
#' @param argname Name of the argument, for error messages.
#'
#' @return Named numeric vector in the same order as \code{start}.
#' @keywords internal
.nested_align_bounds <- function(bounds, start, argname) {
  if (!is.null(names(bounds)) && all(names(start) %in% names(bounds))) {
    return(bounds[names(start)])
  }
  if (length(bounds) != length(start)) {
    stop(argname, " must have one value per fitted parameter (",
         paste(names(start), collapse = ", "), ")", call. = FALSE)
  }
  stats::setNames(as.numeric(bounds), names(start))
}


#' Resolve multstart bounds for nested models
#'
#' Validates \code{multstart_iter} against the fitted parameters and fills in
#' the multstart starting bounds from the fitting bounds where they have not
#' been supplied.
#'
#' @param start Named numeric vector of starting values.
#' @param lower Named numeric vector of lower bounds.
#' @param upper Named numeric vector of upper bounds.
#' @param multstart_iter Number of multstart iterations. A single number, or
#'   one value per fitted parameter.
#' @param multstart_lower Optional. Lower bounds for multstart.
#' @param multstart_upper Optional. Upper bounds for multstart.
#'
#' @return A list with \code{lower} and \code{upper} elements.
#' @keywords internal
.nested_multstart_bounds <- function(start, lower, upper, multstart_iter,
                                     multstart_lower, multstart_upper) {
  if (!is.numeric(multstart_iter) || anyNA(multstart_iter) ||
      any(multstart_iter < 1)) {
    stop("multstart_iter must be one or more numbers, each at least 1",
         call. = FALSE)
  }
  if (!(length(multstart_iter) %in% c(1L, length(start)))) {
    stop("multstart_iter must be either a single number, or one value per ",
         "fitted parameter (", length(start), ": ",
         paste(names(start), collapse = ", "), ")", call. = FALSE)
  }

  if (prod(multstart_iter) > 1) {
    if (is.null(multstart_lower)) multstart_lower <- lower
    if (is.null(multstart_upper)) multstart_upper <- upper
    multstart_lower <- .nested_align_bounds(multstart_lower, start,
                                            "multstart_lower")
    multstart_upper <- .nested_align_bounds(multstart_upper, start,
                                            "multstart_upper")
  }

  list(lower = multstart_lower, upper = multstart_upper)
}


#' Format a fixed parameter value for use in a delta method expression
#'
#' @param x Numeric value.
#'
#' @return A character string of the value, parenthesised.
#' @keywords internal
.nested_num <- function(x) {
  paste0("(", format(x, digits = 15), ")")
}


#' Standard errors for the shared (outer) parameters of a nested model
#'
#' The outer objective is a profile weighted residual sum of squares: the
#' per-region parameters are optimised out at every outer step. Its Hessian at
#' the optimum therefore carries the profile curvature of the shared
#' parameters, from which approximate standard errors follow in the usual
#' nonlinear least squares way, \code{2 * sigma^2 * H^-1}.
#'
#' Because the objective is only available through repeated inner fits, the
#' Hessian is obtained by finite differences and is correspondingly
#' approximate. Steps are scaled to each parameter rather than using optim()'s
#' fixed default, which is far too small for this kind of objective.
#'
#' @param optim_result The object returned by \code{optim()}.
#' @param objective The outer objective function.
#' @param n_obs Number of contributing observations across all regions.
#' @param n_par Total number of estimated parameters (shared plus per-region).
#'
#' @return Named numeric vector of standard errors as a fraction of the
#'   estimate (matching the convention of \code{get_se}), or NA values where
#'   they could not be derived -- including where any inner fit failed while
#'   the Hessian was being evaluated, since the penalty such a failure
#'   contributes would otherwise masquerade as very sharp curvature.
#'
#'   Where \code{optim()} reported a problem, a standard error is still given
#'   if the optimum survives a direct check: no neighbouring point along any
#'   axis is lower. This admits the common case of L-BFGS-B failing its line
#'   search on a slightly kinked objective at the optimum, while still refusing
#'   an optimisation that stopped part-way down a slope.
#' @keywords internal
.nested_outer_se <- function(optim_result, objective, n_obs, n_par) {
  pars <- optim_result$par
  failed <- stats::setNames(rep(NA_real_, length(pars)), names(pars))

  rss <- optim_result$value
  df <- n_obs - n_par
  if (!is.finite(rss) || df <= 0) return(failed)

  # Finite differences require a deterministic objective. With a scalar
  # multstart_iter the inner fits draw a fresh random start design on every
  # evaluation, so the objective is stochastic and its numerical derivatives are
  # noise. Rather than reason about which code path is stochastic, test the
  # property directly: evaluate twice at the optimum and compare.
  repeated <- tryCatch(c(objective(pars), objective(pars)),
                       error = function(e) NULL)
  if (is.null(repeated) || !all(is.finite(repeated)) ||
      !isTRUE(all.equal(repeated[1], repeated[2], tolerance = 1e-8))) {
    return(failed)
  }

  ndeps <- pmax(1e-3, abs(pars) * 0.01)

  # Curvature only describes uncertainty at an optimum, so a reported problem
  # from optim() is a reason for suspicion. It is not proof of one, though:
  # L-BFGS-B returns code 52 when its line search fails, which a profile
  # objective built from inner fits can provoke at the optimum itself, its
  # small numerical kinks defeating the line search while the estimate is
  # perfectly good. So rather than trust the code, check what it stands in for
  # -- that the search really did stop at the bottom, with no neighbouring
  # point along any axis lower than this one.
  if (!isTRUE(optim_result$convergence == 0)) {
    failures_before <- .nested_fit_failures$n

    at_minimum <- vapply(seq_along(pars), function(j) {
      step <- rep(0, length(pars))
      step[j] <- ndeps[j]
      neighbours <- tryCatch(c(objective(pars + step), objective(pars - step)),
                             error = function(e) NULL)
      !is.null(neighbours) && all(is.finite(neighbours)) &&
        all(neighbours >= repeated[1])
    }, logical(1))

    # A neighbour that only looks higher because an inner fit failed there
    # proves nothing, so treat that as a failure of the check too.
    if (.nested_fit_failures$n > failures_before || !all(at_minimum)) {
      return(failed)
    }
  }

  # A failed inner fit contributes a large but finite penalty to the objective.
  # If one occurs at a perturbed point, the resulting curvature is enormous and
  # the standard error comes out minute -- overconfident, and passing every
  # finiteness check on the way. Count the failures across the Hessian
  # evaluations and refuse to report a standard error if any occurred.
  failures_before <- .nested_fit_failures$n

  hess <- tryCatch(
    stats::optimHess(pars, objective, control = list(ndeps = ndeps)),
    error = function(e) NULL
  )
  if (is.null(hess) || !all(is.finite(hess))) return(failed)
  if (.nested_fit_failures$n > failures_before) return(failed)

  vcov <- tryCatch(2 * (rss / df) * solve(hess), error = function(e) NULL)
  if (is.null(vcov)) return(failed)

  variances <- diag(vcov)
  variances[!is.finite(variances) | variances < 0] <- NA_real_

  stats::setNames(abs(sqrt(variances) / pars), names(pars))
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
  } else if (!is.numeric(roiweights)) {
    stop("roiweights must be numeric", call. = FALSE)
  } else if (length(roiweights) == length(region)) {
    roiweights <- vapply(regions,
                         function(r) unname(roiweights[region == r][1]),
                         numeric(1))
  } else if (!is.null(names(roiweights))) {
    missing_regions <- setdiff(regions, names(roiweights))
    if (length(missing_regions) > 0) {
      stop("named roiweights is missing a value for the region",
           if (length(missing_regions) > 1) "s: " else ": ",
           paste(missing_regions, collapse = ", "), call. = FALSE)
    }
    roiweights <- roiweights[regions]
  } else if (length(roiweights) == length(regions)) {
    names(roiweights) <- regions
  } else {
    stop("roiweights must be NULL, a named vector, a vector of length ",
         length(regions), " (one per region), or length ", length(region),
         " (one per observation)", call. = FALSE)
  }

  if (anyNA(roiweights) || any(!is.finite(roiweights)) ||
      any(roiweights <= 0)) {
    stop("roiweights must all be finite and greater than zero", call. = FALSE)
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
#' @return A single ggplot if <= max_facets regions, otherwise a
#'   \code{nested_kinfit_plots} list of ggplot objects, which prints each of
#'   its plots in turn.
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
  names(plots) <- NULL

  class(plots) <- "nested_kinfit_plots"
  return(plots)
}


#' Print a paginated set of nested model plots
#'
#' @param x A \code{nested_kinfit_plots} object.
#' @param ... Passed on to the print method of each plot.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.nested_kinfit_plots <- function(x, ...) {
  for (p in x) print(p, ...)
  invisible(x)
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
