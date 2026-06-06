# Internal: run an mgcv fit while muffling the benign, expected warning that
# arises from the default global + factor-smooth structure
# (s(time) + s(time, group, bs = "fs")), where mgcv notes "repeated 1-d smooths
# of same variable". gam.side() resolves the identifiability via sum-to-zero
# constraints, so the message is informational; because metab_hgam refits many
# times during the IRLS/ramp it is otherwise needlessly noisy. Only this exact
# message is muffled -- all other warnings pass through untouched.
quiet_repeated_smooth <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("repeated 1-d smooths of same variable",
                conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Hierarchical GAM for Parent Fraction with Optional Monotonicity
#'
#' Fits a hierarchical generalised additive model (GAM) to time-varying
#' quantities such as metabolite parent fraction, pooling information across many
#' PET measurements at once, with optional enforcement of a monotonic decrease in
#' the predicted curves. This wraps \code{\link[mgcv]{gam}} and, when a constraint
#' is requested, adds an iteratively reweighted penalty on positive derivatives
#' evaluated on a (group, time) grid.
#'
#' Unlike the parametric parent fraction models in \pkg{kinfitr}
#' (\code{\link{metab_hill}}, \code{\link{metab_sigmoid}}, etc.), which fit a
#' single curve at a time, \code{metab_hgam} fits all curves jointly. The data
#' should therefore be in long format: one row per observation, with columns for
#' time, the measured parent fraction, and a grouping variable that uniquely
#' identifies each curve (typically the PET measurement).
#'
#' The monotonicity penalty acts on the derivative of the predicted curve for each
#' level of \code{group_var}. With \code{monotone = "soft"} the smoothing
#' parameter for the penalty is estimated by REML alongside the wiggliness
#' penalties, letting the data override the constraint where it strongly
#' disagrees. With \code{monotone = "hard"} the smoothing parameter is fixed and
#' ramped upward until the maximum positive slope on the grid falls below
#' \code{hard_tol}. With \code{monotone = "none"} (the default) a standard
#' unconstrained GAM is returned.
#'
#' If \code{formula} is not supplied, a default hierarchical formula is
#' constructed from the column-name arguments:
#' \code{parentFraction ~ s(time, k = 8) + s(time, pet, bs = "fs", k = 5)}.
#' Power users can pass any \pkg{mgcv}-style \code{formula} instead;
#' \code{time_var} and \code{group_var} are still required because they define the
#' grid on which the derivative (and hence the monotonicity penalty) is evaluated.
#'
#' @param data A data frame in long format containing (at least) the columns
#'   named by \code{time_var}, \code{parentFraction_var}, and \code{group_var}.
#' @param time_var Name of the time column. Time is conventionally in seconds in
#'   \pkg{kinfitr}. Default \code{"time"}.
#' @param parentFraction_var Name of the measured parent fraction column. Default
#'   \code{"parentFraction"}.
#' @param group_var Name of the finest-grain grouping column, which uniquely
#'   identifies each predicted curve. Default \code{"pet"}.
#' @param formula Optional \pkg{mgcv}-style formula. If \code{NULL} (default), a
#'   default hierarchical formula is built from the column-name arguments (see
#'   Details).
#' @param monotone One of \code{"none"} (default), \code{"soft"}, or
#'   \code{"hard"}. See Details.
#' @param family An \pkg{mgcv} family. Defaults to
#'   \code{\link[mgcv]{betar}(link = "logit")} for responses bounded in (0, 1).
#' @param n_grid Number of time points at which to evaluate the slope. Default 60.
#' @param max_iter Maximum IRLS iterations per penalty pass. Default 6.
#' @param tol Coefficient convergence tolerance for IRLS. Default 1e-4.
#' @param hard_tol Maximum allowed positive slope under \code{monotone = "hard"}.
#'   Default 1e-6.
#' @param hard_lambda_init Starting smoothing parameter for \code{"hard"} mode.
#'   Default 1e4.
#' @param hard_lambda_max Maximum smoothing parameter before giving up in
#'   \code{"hard"} mode. Default 1e12.
#' @param verbose If \code{TRUE}, print iteration diagnostics. Default
#'   \code{FALSE}.
#'
#' @return An object of class \code{"gam"}, identical in structure to the output
#'   of \code{\link[mgcv]{gam}}. Downstream tools such as \code{predict},
#'   \code{summary}, and \code{\link[mgcv]{k.check}} work unchanged.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @seealso \code{\link[mgcv]{gam}}, \code{\link[mgcv]{betar}},
#'   \code{\link{metab_hill}}
#'
#' @examples
#' \dontrun{
#' # Long-format parent fraction data with one row per (pet, time) observation
#' set.seed(42)
#' d <- expand.grid(
#'   subj_i = 1:4, pet_i = 1:2,
#'   time = seq(1, 60, length.out = 12)
#' )
#' d$pet <- factor(sprintf("sub-%02d_ses-%02d", d$subj_i, d$pet_i))
#' d$parentFraction <- plogis(2 - 0.07 * d$time + rnorm(nrow(d), 0, 0.1))
#'
#' # Unconstrained
#' f0 <- metab_hgam(d)
#'
#' # Soft and hard monotonicity
#' f1 <- metab_hgam(d, monotone = "soft")
#' f2 <- metab_hgam(d, monotone = "hard")
#' }
#'
#' @export
metab_hgam <- function(data,
                       time_var           = "time",
                       parentFraction_var = "parentFraction",
                       group_var          = "pet",
                       formula            = NULL,
                       monotone           = c("none", "soft", "hard"),
                       family             = mgcv::betar(link = "logit"),
                       n_grid = 60, max_iter = 6, tol = 1e-4,
                       hard_tol = 1e-6,
                       hard_lambda_init = 1e4,
                       hard_lambda_max  = 1e12,
                       verbose = FALSE) {

  monotone <- match.arg(monotone)

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame in long format.")
  }
  if (!time_var %in% names(data)) {
    stop(sprintf("Column '%s' not found in data. Set time_var explicitly.",
                 time_var))
  }
  if (!parentFraction_var %in% names(data)) {
    stop(sprintf(
      "Column '%s' not found in data. Set parentFraction_var explicitly.",
      parentFraction_var))
  }
  if (!group_var %in% names(data)) {
    stop(sprintf("Column '%s' not found in data. Set group_var explicitly.",
                 group_var))
  }

  if (is.null(formula)) {
    formula <- stats::as.formula(sprintf(
      "%s ~ s(%s, k = 8) + s(%s, %s, bs = 'fs', k = 5)",
      parentFraction_var, time_var, time_var, group_var))
  }

  if (!is.factor(data[[group_var]])) {
    data[[group_var]] <- factor(data[[group_var]])
  }

  fit <- quiet_repeated_smooth(
    mgcv::gam(formula, data = data, family = family, method = "REML"))
  if (monotone == "none") return(fit)

  # ---- Build derivative matrix D over (group, time) grid ----
  groups   <- levels(droplevels(data[[group_var]]))
  template <- data[!duplicated(data[[group_var]]), , drop = FALSE]
  template <- template[match(groups, template[[group_var]]), , drop = FALSE]
  time_seq <- seq(min(data[[time_var]]), max(data[[time_var]]),
                  length.out = n_grid)

  newdat <- do.call(rbind, lapply(seq_along(groups), function(i) {
    out <- template[rep(i, n_grid), , drop = FALSE]
    out[[time_var]] <- time_seq
    out
  }))
  newdat_hi <- newdat
  eps <- diff(range(time_seq)) * 1e-5
  newdat_hi[[time_var]] <- newdat_hi[[time_var]] + eps

  X_lo <- stats::predict(fit, newdata = newdat,    type = "lpmatrix")
  X_hi <- stats::predict(fit, newdata = newdat_hi, type = "lpmatrix")
  D <- (X_hi - X_lo) / eps

  # ---- IRLS pass at a given smoothing-parameter value ----
  # sp_value = -1 lets REML estimate it (soft); sp_value > 0 fixes it (hard).
  irls_fit <- function(beta_init, sp_value) {
    beta <- beta_init
    fit_local <- fit
    for (iter in seq_len(max_iter)) {
      derivs <- as.numeric(D %*% beta)
      w <- as.numeric(derivs > 0)
      n_pos <- sum(w)
      if (verbose) {
        max_slope <- if (n_pos > 0) max(derivs[w == 1]) else 0
        message(sprintf("  iter %d: %d/%d positive slopes, max %.4g",
                        iter, n_pos, length(derivs), max_slope))
      }
      if (n_pos == 0) {
        return(list(fit = fit_local, beta = beta, max_slope = 0))
      }

      Dw     <- D[w == 1, , drop = FALSE]
      P_mono <- crossprod(Dw) + 1e-10 * diag(ncol(D))
      rk     <- qr(Dw)$rank

      # sp_value >= 0 fixes the penalty's smoothing parameter (hard mode);
      # sp_value < 0 (i.e. -1) lets REML estimate it (soft mode).
      fixed_sp <- sp_value >= 0

      G <- quiet_repeated_smooth(
        mgcv::gam(formula, data = data, family = family,
                  fit = FALSE, method = "REML"))
      m      <- length(G$S)
      G$S    <- c(G$S,    list(P_mono))
      G$off  <- c(G$off,  1L)
      G$rank <- c(G$rank, rk)

      # mgcv maps the (log) smoothing parameters of all penalties from a set of
      # free working parameters via log(sp) = L %*% lsp + lsp0. We always supply
      # an explicit (L, lsp0): leaving them NULL lets mgcv auto-build a mapping
      # that, for extended families such as betar (which carry extra scale/shape
      # parameters), is inconsistent with our appended penalty and triggers
      # recycling warnings inside the REML optimiser.
      #   - Hard mode: *fix* the penalty by adding a zero row to L (no free
      #     working parameter) and carrying log(sp_value) in lsp0; G$sp (the free
      #     working parameters) is left unchanged.
      #   - Soft mode: *estimate* the penalty by REML by giving it its own free
      #     working parameter (an extra row and column in L). Left to REML, the
      #     monotonicity penalty is weighted by the data rather than driven to
      #     zero.
      baseL    <- if (is.null(G$L))    diag(m)   else G$L
      baseLsp0 <- if (is.null(G$lsp0)) rep(0, m) else G$lsp0
      if (fixed_sp) {
        G$L    <- rbind(baseL, rep(0, ncol(baseL)))
        G$lsp0 <- c(baseLsp0, log(sp_value))
      } else {
        G$L    <- rbind(cbind(baseL, rep(0, nrow(baseL))),
                        c(rep(0, ncol(baseL)), 1))
        G$lsp0 <- c(baseLsp0, 0)
        G$sp   <- c(G$sp, sp_value)
      }
      fit_local <- quiet_repeated_smooth(mgcv::gam(G = G))
      new_beta  <- stats::coef(fit_local)
      converged <- max(abs(new_beta - beta)) < tol
      beta      <- new_beta
      if (converged) break
    }
    derivs <- as.numeric(D %*% beta)
    list(fit = fit_local, beta = beta, max_slope = max(c(0, derivs)))
  }

  if (monotone == "soft") {
    res <- tryCatch(irls_fit(stats::coef(fit), sp_value = -1),
                    error = function(e) NULL)
    if (is.null(res)) {
      warning("Soft monotonicity fit failed numerically; returning the unconstrained fit.")
      return(fit)
    }
    return(res$fit)
  }

  # monotone == "hard". The penalised fit can become numerically unstable once
  # the fixed smoothing parameter is very large (particularly for extended
  # families such as betar, where the scale/shape parameter estimation diverges
  # as the linear predictor flattens). We therefore ramp lambda upward but stop
  # gracefully at the last stable fit if a step fails.
  lambda       <- hard_lambda_init
  beta_current <- stats::coef(fit)
  fit_current  <- fit
  res          <- NULL
  while (lambda <= hard_lambda_max) {
    if (verbose) message(sprintf("Hard: trying lambda = %.2g", lambda))
    res_try <- tryCatch(irls_fit(beta_current, sp_value = lambda),
                        error = function(e) e)
    if (inherits(res_try, "error")) {
      warning(sprintf(
        "Hard fit became numerically unstable at lambda = %.2g; returning the closest stable fit (max positive slope = %.2g).",
        lambda, if (is.null(res)) NA_real_ else res$max_slope))
      return(fit_current)
    }
    res          <- res_try
    fit_current  <- res$fit
    beta_current <- res$beta
    if (res$max_slope < hard_tol) {
      if (verbose) {
        message(sprintf(
          "Hard constraint satisfied (lambda = %.2g, max slope = %.2g)",
          lambda, res$max_slope))
      }
      return(fit_current)
    }
    lambda <- lambda * 10
  }
  warning(sprintf(
    "Hard monotonicity not achieved within lambda_max = %.2g (max positive slope = %.2g). Returning closest fit.",
    hard_lambda_max,
    if (is.null(res)) NA_real_ else res$max_slope))
  fit_current
}
