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

# Internal: beta regression whose log precision changes linearly in a supplied
# covariate, so that log(theta_i) = th[1] + th[2] * t_std[i].
#
# mgcv has no beta location-scale family (betar carries a single scalar theta),
# but betar stores and differentiates theta on the log scale and every slot is
# vectorised. Each betar theta-derivative therefore maps to the two-parameter
# case by multiplying by powers of t_std, so no new derivatives are needed.
# Verified against finite differences of the deviance to ~1e-8.
betar_lintheta <- function(t_std, ini = c(0, 0), link = "logit",
                           eps = .Machine$double.eps * 100) {

  stopifnot(is.numeric(t_std), all(is.finite(t_std)))
  base <- do.call(mgcv::betar, list(theta = NULL, link = link, eps = eps))

  env <- new.env(parent = emptyenv())
  assign(".Theta", ini, envir = env)
  assign(".t_std", t_std, envir = env)
  assign(".warned_recycle", FALSE, envir = env)

  # The replacement saturated.ll below was checked against optimize() up to a
  # log precision of 20; beyond about 25 its bounded Newton iteration is no
  # longer reliable at boundary responses. Cap well inside the verified range,
  # which is far above any plausible parent-fraction precision.
  .lt_min <- -20   # theta ~ 2e-9
  .lt_max <-  20   # theta ~ 4.9e8

  base_Dd  <- base$Dd
  base_dev <- base$dev.resids
  base_aic <- base$aic

  # Per-observation log theta, with a guard against silent row misalignment.
  logth <- function(theta, n) {
    u <- get(".t_std", envir = env)
    if (length(u) != n) {
      stop(sprintf(paste0("betar_lintheta: the precision covariate has length ",
                          "%d but the model has %d observations. Drop ",
                          "incomplete rows before building the family."),
                   length(u), n))
    }
    # Clamp the per-observation log precision. A scalar theta stays in range on
    # its own, but theta[1] + theta[2] * u can reach extremes during the line
    # search, and exp() then overflows so dbeta() returns NaN and the REML
    # optimiser fails with "missing value where TRUE/FALSE needed".
    pmax(pmin(theta[1] + theta[2] * u, .lt_max), .lt_min)
  }

  fam <- base
  fam$family   <- "Beta regression (linear log-precision)"
  fam$n.theta  <- 2
  fam$ini.theta <- ini

  fam$getTheta <- function(trans = FALSE) {
    th <- get(".Theta", envir = env)
    if (trans) th else th          # theta is a 2-vector; no scalar transform
  }
  fam$putTheta <- function(theta) assign(".Theta", theta, envir = env)

  fam$dev.resids <- function(y, mu, wt, theta = NULL) {
    if (is.null(theta)) theta <- get(".Theta", envir = env)
    base_dev(y, mu, wt, logth(theta, length(y)))
  }

  fam$aic <- function(y, mu, theta = NULL, wt, dev) {
    if (is.null(theta)) theta <- get(".Theta", envir = env)
    base_aic(y, mu, logth(theta, length(y)), wt, dev)
  }

  # Saturated log-likelihood is zero for betar; only the dimensions change.
  fam$ls <- function(y, w, theta, scale) {
    list(ls = 0, lsth1 = c(0, 0),
         LSTH1 = matrix(0, length(y), 2), lsth2 = matrix(0, 2, 2))
  }

  fam$Dd <- function(y, mu, theta, wt, level = 0) {
    n <- length(y)
    u <- get(".t_std", envir = env)
    r <- base_Dd(y, mu, logth(theta, n), wt, level)

    # d logtheta_i / d(th1, th2) = (1, u_i), so a first-order theta slot gains
    # a u-scaled column and a second-order slot gains u and u^2 columns
    # (upper-triangular order: (1,1), (1,2), (2,2)). betar populates these
    # slots only at the levels that need them, so expand whatever is present
    # and leave the rest absent.
    # Where logth() clamped the log precision, theta_i no longer varies with
    # the parameters, so every theta-derivative for that observation is zero.
    # Without this the analytic derivatives disagree with the likelihood the
    # clamp actually produces.
    lt_raw <- theta[1] + theta[2] * u
    act    <- as.numeric(lt_raw > .lt_min & lt_raw < .lt_max)

    ex1 <- function(v) {
      if (is.null(v)) NULL else cbind(v * act, v * u * act)
    }
    ex2 <- function(v) {
      if (is.null(v)) NULL else cbind(v * act, v * u * act, v * u^2 * act)
    }

    for (nm in c("Dth", "Dmuth", "Dmu2th", "EDmu2th", "Dmu3th")) {
      if (!is.null(r[[nm]])) r[[nm]] <- ex1(r[[nm]])
    }
    for (nm in c("Dth2", "Dmuth2", "Dmu2th2")) {
      if (!is.null(r[[nm]])) r[[nm]] <- ex2(r[[nm]])
    }
    r
  }

  fam$variance <- function(mu) {
    th <- get(".Theta", envir = env)
    mu * (1 - mu) / (1 + exp(logth_recycle(th, length(mu))))
  }

  # betar's own saturated.ll subsets y and eta but passes theta through
  # unsubsetted, which is correct for a scalar and produces NaN for a
  # per-observation vector. The beta log-density is concave in mu, so maximise
  # it per observation with a bounded Newton iteration instead.
  fam$saturated.ll <- function(y, wt, theta = NULL) {
    if (is.null(theta)) theta <- get(".Theta", envir = env)
    phi <- exp(logth(theta, length(y)))
    lo  <- 1e-10
    mu  <- pmin(pmax(y, lo), 1 - lo)
    lly <- log(y)
    ll1 <- log1p(-y)
    for (i in 1:60) {
      a  <- phi * mu
      b  <- phi * (1 - mu)
      g  <- phi * (-digamma(a) + digamma(b) + lly - ll1)
      h  <- -phi^2 * (trigamma(a) + trigamma(b))     # always negative
      st <- g / h
      st[!is.finite(st)] <- 0
      st <- pmax(pmin(st, 0.1), -0.1)                # limit the step
      mu <- pmin(pmax(mu - st, lo), 1 - lo)
      if (max(abs(st)) < 1e-12) break
    }
    LS <- stats::dbeta(y, phi * mu, phi * (1 - mu), log = TRUE)
    LS[!is.finite(LS)] <- 0
    list(f = sum(wt * LS), term = LS, mu = mu)
  }

  # betar's rd, qf and postproc read a scalar .Theta from betar's own
  # environment, which stays at its starting value here, so they would report a
  # precision of 1 regardless of the fit. postproc also pastes theta into the
  # family label, which turns family$family into a length-2 vector.
  # Unlike the fitting slots, qf and rd may be called with a vector shorter
  # than the data (for example a single mu), so recycle the covariate instead
  # of insisting on an exact length. They assume mu is in the data's row order.
  logth_recycle <- function(theta, n_out) {
    u_full <- get(".t_std", envir = env)
    if (n_out != length(u_full) && !get(".warned_recycle", envir = env)) {
      assign(".warned_recycle", TRUE, envir = env)
      warning("betar_lintheta: asked for ", n_out, " values while the ",
              "precision covariate has length ", length(u_full),
              ". Recycling from the first row, so the precision used may not ",
              "correspond to the observations you passed. Pass values in the ",
              "data's row order and at full length for exact results.",
              call. = FALSE)
    }
    u <- rep_len(u_full, n_out)
    pmax(pmin(theta[1] + theta[2] * u, .lt_max), .lt_min)
  }

  fam$rd <- function(mu, wt, scale) {
    Theta <- exp(logth_recycle(get(".Theta", envir = env), length(mu)))
    r <- stats::rbeta(length(mu), shape1 = Theta * mu,
                      shape2 = Theta * (1 - mu))
    r[r >= 1 - eps] <- 1 - eps
    r[r < eps] <- eps
    r
  }

  fam$qf <- function(p, mu, wt, scale) {
    Theta <- exp(logth_recycle(get(".Theta", envir = env), length(mu)))
    q <- stats::qbeta(p, shape1 = Theta * mu, shape2 = Theta * (1 - mu))
    q[q >= 1 - eps] <- 1 - eps
    q[q < eps] <- eps
    q
  }

  fam$postproc <- function(family, y, prior.weights, fitted,
                           linear.predictors, offset, intercept) {
    theta <- family$getTheta()
    lf <- family$saturated.ll(y, prior.weights, theta)
    l2 <- family$dev.resids(y, fitted, prior.weights)
    posr <- list()
    posr$deviance <- 2 * lf$f + sum(l2)
    wtdmu <- if (intercept) {
      sum(prior.weights * y) / sum(prior.weights)
    } else {
      family$linkinv(offset)
    }
    posr$null.deviance <- 2 * lf$f +
      sum(family$dev.resids(y, wtdmu, prior.weights))
    # A single string, so family$family stays length 1.
    posr$family <- paste0("Beta regression(log precision ",
                          round(theta[1], 3), " + ",
                          round(theta[2], 3), " * t_std)")
    posr
  }

  fam$lt_range <- c(.lt_min, .lt_max)
  fam$env <- env
  fam
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
#' A \code{formula} that does not depend on \code{time_var} is rejected in the
#' constrained modes, because the derivative grid would then be identically zero
#' and an increasing curve would be reported as monotone.
#'
#' \strong{What the constraint does and does not guarantee.} Three limits are
#' worth knowing:
#' \itemize{
#'   \item The slope is constrained at \code{n_grid} points, not everywhere, so
#'     a dense evaluation can show a small excess between grid points. Raise
#'     \code{n_grid} if that matters.
#'   \item The constraint acts on the linear predictor. With the default logit
#'     link the response is a monotone increasing function of it, so the two
#'     agree. With a decreasing link (for example \code{Gamma(link =
#'     "inverse")}) they do not, and the response can increase while the linear
#'     predictor decreases.
#'   \item Monotonicity holds over the observed range of \code{time_var} only.
#'     \code{predict} beyond that range extrapolates and can increase.
#'   \item The derivative grid takes every variable other than
#'     \code{time_var} from each group's first row. If a custom
#'     \code{formula} uses a predictor that varies \emph{within} a level of
#'     \code{group_var}, the constraint is applied only at that first row's
#'     value, and reordering the rows can change the result. Keep such
#'     predictors constant within \code{group_var}, or add them to
#'     \code{group_var}.
#' }
#'
#' The grid limit can bite hard with an oscillating custom \code{formula}: if
#' the oscillation happens to be non-increasing at all \code{n_grid} points, the
#' fit is reported as satisfying the constraint while rising steeply between
#' them. The default penalised-spline formula cannot alias this way, but an
#' arbitrary parametric one can.
#'
#' \strong{Time-varying precision.} Parent fraction is usually measured less
#' precisely late in a scan than early, but \code{\link[mgcv]{betar}} carries a
#' single scalar precision that has to compromise between the two. Setting
#' \code{theta_time = TRUE} replaces it with a family whose log precision changes
#' linearly across the scan, \code{log(theta_i) = a + b * t_std}, where
#' \code{t_std} runs from 0 at the first observed time to 1 at the last. A
#' negative \code{b} means the precision falls over the scan. Both parameters are
#' estimated by REML alongside everything else, and \code{b = 0} reproduces the
#' constant-precision fit, so the option only adds flexibility. The monotonicity
#' machinery is unaffected, because the precision does not enter the mean's linear
#' predictor.
#'
#' Two conditions apply. Every variable in \code{formula} must be a column of
#' \code{data}, because the precision covariate is built before fitting and has
#' to line up with the rows mgcv keeps. And the log precision is capped at plus
#' or minus 20; if the fitted precision reaches that cap for every observation
#' then the slope cannot affect the likelihood, so \code{metab_hgam} warns and
#' sets \code{$metab_hgam$theta_identified} to \code{FALSE}. Check that flag
#' before interpreting \code{theta}.
#'
#' \strong{What slope to expect.} Because \code{log(theta)} is linear in time,
#' \code{theta} itself changes by a constant proportion per unit time, which is
#' how radioactive decay behaves. If every sample has the same volume and
#' counting time, counts track activity, activity decays as
#' \code{exp(-lambda * t)}, and a counted proportion out of \code{n} has
#' precision of roughly \code{n}. The slope should therefore be about the
#' negative decay constant:
#'
#' \deqn{b \approx -\log(2) \times T / t_{1/2}}{b = -log(2) * T / half_life}
#'
#' where \code{T} is the time span of the data, \code{diff(range(time_var))},
#' in the same units as the half-life. For carbon-11 (20.4 min) over an hour
#' that is about -2, and for fluorine-18 (109.8 min) about -0.37. Comparing the
#' fitted slope with this value gives an external check that no model comparison
#' provides.
#'
#' Two things loosen the prediction. Simulations recovered it without bias at
#' moderate sample sizes, but showed a steepening of roughly 10 per cent with few
#' measurements or low late-scan counts, so do not read a small deviation as
#' meaningful. And HPLC noise from peak integration and baseline does not decay,
#' so once it dominates, the real \code{log(theta)} flattens rather than
#' continuing to fall. Unequal sample volumes or counting times break the
#' derivation outright.
#'
#' \strong{Standard errors.} \code{predict(se.fit = TRUE)}, and anything built
#' on it such as \code{gratia::fitted_values}, reports intervals from the
#' penalised covariance matrix \code{Vp}. Under \code{monotone = "hard"} those
#' intervals are conditional on a smoothing parameter that was \emph{chosen by
#' a search} in order to satisfy \code{hard_tol}, and they do not account for
#' that search. They are therefore too narrow wherever the constraint binds,
#' which is exactly where the data disagreed with monotonicity. Use
#' \code{monotone = "none"} for uncertainty, or resample, and keep the
#' constrained fit for the point estimate. The precision parameters' own
#' uncertainty is not propagated either, as with a scalar
#' \code{\link[mgcv]{betar}} theta.
#'
#' \strong{Parent fractions of exactly 0 or 1.} These are accepted, because
#' values of exactly 1 are common at the first sample, but a beta likelihood has
#' no density at the boundary. \code{\link[mgcv]{betar}} nudges them inwards by
#' its \code{eps}, which can leave the reported deviance negative and produce
#' repeated \dQuote{saturated likelihood may be inaccurate} warnings. If that
#' matters, nudge them yourself before fitting. Values outside [0, 1] and
#' non-finite values are rejected outright.
#'
#' \strong{Relationship to \code{bd_addfit}.} The returned object is not a
#' drop-in replacement for the single-curve parent fraction fits accepted by
#' \code{\link{bd_addfit}}. \code{\link{bd_getdata}} calls \code{predict} with
#' a time vector alone, which fails here because the model also needs
#' \code{group_var}. Predict per measurement and pass the fitted values on
#' instead.
#'
#' @param data A data frame in long format containing (at least) the columns
#'   named by \code{time_var}, \code{parentFraction_var}, and \code{group_var}.
#' @param time_var Name of the time column. Either seconds or minutes works,
#'   because the penalty is scaled to be dimensionless, so the
#'   \code{hard_lambda_*} defaults suit both. The ramp is not invariant to the
#'   unit, though: the same data may converge at a different lambda in seconds
#'   than in minutes. \code{hard_tol} is unit-dependent, being a slope. Default
#'   \code{"time"}.
#' @param parentFraction_var Name of the measured parent fraction column. Default
#'   \code{"parentFraction"}.
#' @param group_var Name of the finest-grain grouping column, which uniquely
#'   identifies each predicted curve. Default \code{"pet"}.
#' @param formula Optional \pkg{mgcv}-style formula. If \code{NULL} (default), a
#'   default hierarchical formula is built from the column-name arguments (see
#'   Details).
#' @param monotone One of \code{"none"} (default), \code{"soft"}, or
#'   \code{"hard"}. See Details. Note that \code{"soft"} guarantees nothing:
#'   REML weighs the monotonicity penalty against the data and routinely shrinks
#'   it, so on data that genuinely increases the effect is often under 1 per cent
#'   of the unconstrained slope. Use \code{"hard"} when you need the constraint
#'   to hold.
#' @param family An \pkg{mgcv} family. Defaults to
#'   \code{\link[mgcv]{betar}(link = "logit")} for responses bounded in (0, 1).
#' @param theta_time If \code{TRUE}, let the log of the beta precision change
#'   linearly across the scan instead of holding it constant. Requires a beta
#'   \code{family}. See Details. Default \code{FALSE}.
#' @param n_grid Number of time points at which to evaluate the slope, and hence
#'   the number of points at which monotonicity is enforced. Default 60.
#' @param max_iter Maximum IRLS iterations per penalty pass. Default 6.
#' @param tol Coefficient convergence tolerance for IRLS. Default 1e-4.
#' @param hard_tol Maximum allowed positive slope of the linear predictor under
#'   \code{monotone = "hard"}, in units of per unit of \code{time_var}. Scale it
#'   with your time unit: a tolerance of 1e-6 per minute is 1.67e-8 per second.
#'   Default 1e-6.
#' @param hard_lambda_init Starting smoothing parameter for \code{"hard"} mode.
#'   Default 1e4.
#' @param hard_lambda_max Maximum smoothing parameter before giving up in
#'   \code{"hard"} mode. Default 1e12.
#' @param verbose If \code{TRUE}, print iteration diagnostics. Default
#'   \code{FALSE}.
#'
#' @return An object of class \code{"gam"}, as returned by
#'   \code{\link[mgcv]{gam}}, so \code{predict}, \code{summary} and
#'   \code{\link[mgcv]{k.check}} work unchanged. One extra element,
#'   \code{$metab_hgam}, records what the constraint actually achieved, so code
#'   looping over many measurements can check the outcome without catching
#'   warnings:
#'   \itemize{
#'     \item \code{monotone}: the mode requested.
#'     \item \code{converged}: \code{TRUE} if \code{"hard"} met
#'       \code{hard_tol}; \code{FALSE} if it fell back to the closest fit it
#'       obtained; and \code{NA} for \code{"none"}, and for \code{"soft"} when
#'       it succeeded, since neither has a tolerance to meet. A \code{"soft"}
#'       fit that failed numerically reports \code{FALSE}.
#'     \item \code{max_slope}: the largest positive slope on the grid, for the
#'       fit that is actually returned.
#'     \item \code{lambda}: the final penalty smoothing parameter tried.
#'     \item \code{iterations}: the index of the IRLS iterate returned. Because
#'       the best iterate is kept rather than the last, this is not necessarily
#'       the number of iterations run.
#'     \item \code{hard_tol}, \code{n_grid}: the settings used.
#'     \item \code{theta_time}, \code{theta}: whether a linear precision was
#'       fitted, and its \code{intercept} and \code{slope} on the log scale.
#'     \item \code{theta_identified}: \code{FALSE} if the fitted precision sat
#'       at its cap for every observation, in which case \code{theta} is
#'       meaningless. \code{NA} when \code{theta_time = FALSE}.
#'     \item \code{error}: the underlying error message if a fit failed.
#'   }
#'   Always check \code{converged} before trusting a \code{"hard"} fit.
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
#'
#' # Always check the outcome before trusting a hard fit
#' f2$metab_hgam$converged
#' f2$metab_hgam$max_slope
#'
#' # Let the beta precision fall across the scan
#' f3 <- metab_hgam(d, monotone = "hard", theta_time = TRUE)
#' f3$metab_hgam$theta
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
                       theta_time = FALSE,
                       n_grid = 60, max_iter = 6, tol = 1e-4,
                       hard_tol = 1e-6,
                       hard_lambda_init = 1e4,
                       hard_lambda_max  = 1e12,
                       verbose = FALSE) {

  monotone <- match.arg(monotone)
  t_std_keep <- NULL

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

  # ---- Validate the control parameters ----
  # These all feed the derivative grid or the lambda ramp, where a bad value
  # otherwise surfaces as an opaque failure from inside mgcv (or, worse, as a
  # silent fallback to the unconstrained fit).
  chk_scalar <- function(x, nm, min_val, integer = FALSE) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < min_val ||
        (integer && x != round(x))) {
      stop(sprintf("'%s' must be a single %s value >= %g.", nm,
                   if (integer) "whole" else "finite", min_val))
    }
  }
  chk_scalar(n_grid,   "n_grid",   2, integer = TRUE)
  chk_scalar(max_iter, "max_iter", 1, integer = TRUE)
  # Tolerances and lambdas only need to be positive. Requiring at least
  # machine epsilon rejected legitimately tiny tolerances that the pre-change
  # API accepted.
  chk_scalar(tol,              "tol",              .Machine$double.xmin)
  chk_scalar(hard_tol,         "hard_tol",         .Machine$double.xmin)
  chk_scalar(hard_lambda_init, "hard_lambda_init", .Machine$double.xmin)
  chk_scalar(hard_lambda_max,  "hard_lambda_max",  .Machine$double.xmin)
  if (hard_lambda_max < hard_lambda_init) {
    stop("'hard_lambda_max' must be >= 'hard_lambda_init'.")
  }

  # ---- Validate and clean the data ----
  # Drop incomplete rows up front rather than letting mgcv's na.action drop
  # them later: the derivative grid is built from min()/max() of the time
  # column, so a retained NA would otherwise propagate into seq() and fail
  # with an error that does not name the cause.
  # Which variables does the model actually use, and are they all columns of
  # `data`?
  model_vars <- tryCatch(all.vars(mgcv::interpret.gam(formula)$fake.formula),
                         error = function(e) all.vars(formula))
  model_cols <- intersect(
    unique(c(time_var, parentFraction_var, group_var, model_vars)),
    names(data))
  external_vars <- setdiff(model_vars, names(data))

  # Drop incomplete rows only when every variable the model uses is a column of
  # `data`. A variable taken from the formula's environment keeps its original
  # length, so removing rows from `data` would leave the two disagreeing and
  # mgcv would stop with "variable lengths differ". In that case leave the data
  # alone and let mgcv's own na.action drop rows.
  if (!length(external_vars)) {
    keep <- stats::complete.cases(data[, model_cols, drop = FALSE])
    if (!all(keep)) {
      n_drop <- sum(!keep)
      warning(sprintf(
        "Dropped %d row%s with missing values in: %s.",
        n_drop, if (n_drop == 1L) "" else "s",
        paste(model_cols, collapse = ", ")))
      data <- data[keep, , drop = FALSE]
    }
  }
  if (!nrow(data)) {
    stop("No complete rows remain in 'data' after dropping missing values.")
  }

  tv <- data[[time_var]]
  if (!is.numeric(tv) || !all(is.finite(tv))) {
    stop(sprintf("Column '%s' must be numeric and finite.", time_var))
  }
  if (monotone != "none" && diff(range(tv)) <= 0) {
    stop(sprintf("Column '%s' takes only one distinct value, so no ",
                 time_var),
         "derivative grid can be built. Use monotone = \"none\", which needs ",
         "no grid.")
  }

  pf <- data[[parentFraction_var]]
  if (!is.numeric(pf) || !all(is.finite(pf))) {
    stop(sprintf(
      "Column '%s' must be numeric and finite. Non-finite values are ",
      parentFraction_var),
      "otherwise silently treated as the nearest boundary value.")
  }
  if (any(pf < 0 | pf > 1)) {
    rng <- range(pf)
    stop(sprintf(
      "Column '%s' has values outside [0, 1] (range %.4g to %.4g). A beta ",
      parentFraction_var, rng[1], rng[2]),
      "likelihood cannot represent these, and they are otherwise clamped ",
      "to the boundary without warning.")
  }


  if (!is.factor(data[[group_var]])) {
    data[[group_var]] <- factor(data[[group_var]])
  }

  # ---- Optional linear change in the beta precision ----
  # Parent fraction is measured less precisely late in a scan than early, so a
  # single scalar theta has to compromise between the two. theta_time = TRUE
  # lets log(theta) change linearly across the scan instead.
  if (!identical(theta_time, FALSE)) {
    if (!isTRUE(theta_time)) {
      stop("'theta_time' must be TRUE or FALSE.")
    }
    if (length(external_vars)) {
      stop("'theta_time = TRUE' needs every variable in 'formula' to be a ",
           "column of 'data', because the precision covariate is built before ",
           "fitting and must line up with the rows mgcv keeps. These come ",
           "from the formula's environment instead: ",
           paste(external_vars, collapse = ", "),
           ". Add them to 'data', or set 'theta_time = FALSE'.")
    }
    if (identical(as.numeric(family$n.theta), 0)) {
      stop("'theta_time = TRUE' cannot be combined with a fixed precision. ",
           "'family' was built with 'theta' supplied, which holds the ",
           "precision constant. Either drop 'theta' from the family or set ",
           "'theta_time = FALSE'.")
    }
    if (!grepl("^Beta regression", family$family)) {
      stop("'theta_time = TRUE' needs a beta family; 'family' is currently \"",
           family$family, "\". Either leave 'family' at its default or set ",
           "'theta_time = FALSE'.")
    }
    # Standardise time to [0, 1] over the observed range so that theta[2] is
    # unit-free and reads as the change in log precision across the scan.
    t_rng <- range(data[[time_var]])
    t_std <- (data[[time_var]] - t_rng[1]) / diff(t_rng)
    t_std_keep <- t_std
    family <- betar_lintheta(t_std = t_std, link = family$link)
  }

  fit <- quiet_repeated_smooth(
    mgcv::gam(formula, data = data, family = family, method = "REML"))
  # ---- Keep a stored fit's precision consistent with its coefficients ----
  # Every refit is handed the same `family` object, and an extended family
  # keeps its theta in a mutable environment. A fit stored from an earlier
  # iterate therefore reports a later iterate's theta, which makes its
  # deviance, residuals and family label disagree with its coefficients.
  # Snapshot theta whenever a fit is kept, and restore it before returning.
  snap_theta <- function(object) {
    if (is.function(object$family$getTheta)) object$family$getTheta() else NULL
  }
  restore_theta <- function(object, theta) {
    if (!is.null(theta) && is.function(object$family$putTheta)) {
      object$family$putTheta(theta)
    }
    object
  }
  best_out <- function(best, error) {
    c(list(fit = restore_theta(best$fit, best$theta)),
      best[c("beta", "max_slope", "iters")], list(error = error))
  }

  # A linear precision is only identified while some observations sit inside
  # the clamp: once every one is pinned to a bound, the slope cannot change the
  # likelihood and the reported value is meaningless.
  theta_is_identified <- function(object) {
    rng <- object$family$lt_range
    if (is.null(rng)) return(NA)
    th <- object$family$getTheta()
    lt <- th[1] + th[2] * t_std_keep
    ok <- any(lt > rng[1] & lt < rng[2])
    if (!ok) {
      warning("The fitted beta precision is at its limit for every ",
              "observation, so the linear precision is not identified: its ",
              "slope cannot affect the likelihood. Treat $metab_hgam$theta ",
              "as unreliable and consider theta_time = FALSE.")
    }
    ok
  }

  # ---- Record what the constraint machinery actually achieved ----
  # Every return path attaches this list, so a caller looping over many
  # measurements can tell a satisfied constraint from a fallback without having
  # to catch warnings.
  attach_status <- function(object, converged, lambda = NA_real_,
                            max_slope = NA_real_, iters = NA_integer_,
                            error = NULL) {
    object$metab_hgam <- list(
      monotone   = monotone,
      converged  = converged,
      lambda     = lambda,
      max_slope  = max_slope,
      iterations = iters,
      hard_tol   = if (monotone == "hard") hard_tol else NA_real_,
      n_grid     = n_grid,
      theta_time = !identical(theta_time, FALSE),
      theta      = if (identical(theta_time, FALSE)) NULL else
                     stats::setNames(object$family$getTheta(),
                                     c("intercept", "slope")),
      theta_identified = if (identical(theta_time, FALSE)) NA else
                           theta_is_identified(object),
      error      = if (is.null(error)) NULL else conditionMessage(error)
    )
    object
  }

  if (monotone == "none") {
    return(attach_status(fit, converged = NA, max_slope = NA_real_))
  }

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

  # If perturbing time_var does not move the linear predictor, the penalty is
  # identically zero and every fit trivially reports zero positive slopes. That
  # happens when a custom formula depends on some other time column, in which
  # case an increasing curve would otherwise be reported as monotone.
  # Test for exactly zero, not a small threshold: if the model does not use
  # time_var then X_hi and X_lo are identical and D is exactly 0, whereas a
  # genuinely time-dependent but small-scaled column (say I(1e-9 * time)) has
  # tiny but non-zero derivatives and must be allowed.
  if (max(abs(D)) <= 0) {
    stop(sprintf(
      "The model does not depend on '%s', so no monotonicity constraint can ",
      time_var),
      "be applied. Check that 'formula' uses the column named by 'time_var'.")
  }
  # ---- IRLS pass at a given smoothing-parameter value ----
  # sp_value = -1 lets REML estimate it (soft); sp_value > 0 fixes it (hard).
  # fit_init must be the fit whose coefficients are beta_init. Returning `fit`
  # here instead would hand back the globally unconstrained fit whenever a pass
  # achieved nothing, while reporting the previous pass's slope.
  irls_fit <- function(fit_init, beta_init, sp_value) {
    theta_init <- snap_theta(fit_init)
    beta <- beta_init
    fit_local <- fit_init
    best <- NULL
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
        return(list(fit = restore_theta(fit_init, theta_init), beta = beta,
                    max_slope = 0, iters = iter, error = NULL))
      }

      Dw <- D[w == 1, , drop = FALSE]

      # Scale the penalty by its largest eigenvalue so that lambda is
      # dimensionless. crossprod(Dw) carries the units of time_var (it shrinks
      # by 1/3600 when time is given in seconds rather than minutes), so
      # without this an equivalent lambda would have to be ~3600x larger and
      # the hard_lambda_* defaults would only ever suit one choice of unit.
      #
      # The rank declared to mgcv must describe the matrix actually handed to
      # it, so both come from the same eigendecomposition at the same
      # tolerance. An earlier version added a 1e-10 ridge here while declaring
      # rank(Dw): because the whole penalty is multiplied by lambda, that ridge
      # grew into a large ridge on every coefficient (the intercept included)
      # at the top of the ramp, and the full rank it implied disagreed with the
      # deficient rank declared, which drove mgcv's REML log-determinant to
      # NaN. mgcv supports rank-deficient penalties natively through G$rank, so
      # the ridge is not needed to keep the fit well posed.
      C_mono  <- crossprod(Dw)
      ev      <- eigen(C_mono, symmetric = TRUE, only.values = TRUE)$values
      scale_C <- max(ev)
      if (!is.finite(scale_C) || scale_C <= 0) {
        stop("The monotonicity penalty is degenerate: the fitted model has no ",
             "derivative with respect to '", time_var, "' on the evaluation ",
             "grid. Check that `formula` actually depends on `time_var`.")
      }
      P_mono <- C_mono / scale_C
      # Use the same threshold as mgcv's mroot(), which keeps eigenvalues
      # above max * .Machine$double.eps. sqrt(eps) drops one genuine direction
      # here, which measurably weakens the constraint at high lambda.
      rk     <- sum(ev / scale_C > .Machine$double.eps)

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

      # Catch a failed step per *iteration* rather than per pass: an earlier
      # version wrapped the whole pass, so a single bad step discarded every
      # iterate before it -- including iterates that had already met the
      # requested tolerance.
      step <- tryCatch(quiet_repeated_smooth(mgcv::gam(G = G)),
                       error = function(e) e)
      if (inherits(step, "error")) {
        if (is.null(best)) {
          return(list(fit = restore_theta(fit_init, theta_init),
                      beta = beta_init,
                      max_slope = max(c(0, as.numeric(D %*% beta_init))),
                      iters = iter, error = step))
        }
        return(best_out(best, step))
      }

      fit_local <- step
      new_beta  <- stats::coef(fit_local)
      new_slope <- max(c(0, as.numeric(D %*% new_beta)))

      # The active set can oscillate between iterations, so the last iterate is
      # not necessarily the best one. Keep the iterate with the smallest
      # positive slope and return that.
      if (is.null(best) || new_slope < best$max_slope) {
        best <- list(fit = fit_local, beta = new_beta,
                     max_slope = new_slope, iters = iter,
                     theta = snap_theta(fit_local))
      }

      converged <- max(abs(new_beta - beta)) < tol
      beta      <- new_beta
      if (converged) break
    }

    if (is.null(best)) {
      return(list(fit = restore_theta(fit_init, theta_init), beta = beta_init,
                  max_slope = max(c(0, as.numeric(D %*% beta_init))),
                  iters = max_iter, error = NULL))
    }
    best_out(best, NULL)
  }

  if (monotone == "soft") {
    res <- tryCatch(irls_fit(fit, stats::coef(fit), sp_value = -1),
                    error = function(e) e)
    if (inherits(res, "error")) {
      warning("Soft monotonicity fit failed numerically (",
              conditionMessage(res), "); returning the unconstrained fit.")
      return(attach_status(fit, converged = FALSE, error = res))
    }
    # A failed mgcv step is reported inside the returned list, not raised, so
    # it has to be surfaced here or it would pass silently.
    if (!is.null(res$error)) {
      warning("Soft monotonicity fit failed numerically (",
              conditionMessage(res$error),
              "); returning the closest fit obtained.")
      return(attach_status(res$fit, converged = FALSE,
                           max_slope = res$max_slope, iters = res$iters,
                           error = res$error))
    }
    # Soft mode has no tolerance to meet: REML weighs the monotonicity penalty
    # against the data and routinely shrinks it, so the constraint may be
    # barely active. Report the achieved slope and leave `converged` as NA
    # rather than implying the curve is monotone.
    return(attach_status(res$fit, converged = NA,
                         max_slope = res$max_slope, iters = res$iters,
                         error = res$error))
  }

  # monotone == "hard". Ramp the fixed smoothing parameter upward until the
  # largest positive slope on the grid falls below hard_tol. If a step fails
  # numerically we stop at the best fit obtained so far rather than continuing.
  lambda       <- hard_lambda_init
  beta_current <- stats::coef(fit)
  fit_current  <- fit
  res          <- NULL
  repeat {
    if (verbose) message(sprintf("Hard: trying lambda = %.2g", lambda))
    res_try <- tryCatch(irls_fit(fit_current, beta_current, sp_value = lambda),
                        error = function(e) e)
    if (inherits(res_try, "error")) {
      achieved <- max(c(0, as.numeric(D %*% stats::coef(fit_current))))
      warning(sprintf(
        "Hard fit became numerically unstable at lambda = %.2g (%s); returning the closest stable fit (max positive slope = %.2g).",
        lambda, conditionMessage(res_try), achieved))
      return(attach_status(fit_current, converged = FALSE, lambda = lambda,
                           max_slope = achieved, error = res_try))
    }

    # Judge the fit that is actually returned before reacting to a failed step:
    # a pass can meet hard_tol and then fail on a later iteration, and that is
    # still a success.
    if (!is.null(res_try$error) && res_try$max_slope >= hard_tol) {
      warning(sprintf(
        "Hard fit became numerically unstable at lambda = %.2g (%s); returning the closest stable fit (max positive slope = %.2g).",
        lambda, conditionMessage(res_try$error), res_try$max_slope))
      return(attach_status(res_try$fit, converged = FALSE, lambda = lambda,
                           max_slope = res_try$max_slope,
                           iters = res_try$iters, error = res_try$error))
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
      return(attach_status(fit_current, converged = TRUE, lambda = lambda,
                           max_slope = res$max_slope, iters = res$iters,
                           error = res$error))
    }

    # Always attempt hard_lambda_max itself: stepping by decades would
    # otherwise skip it whenever it does not sit on a power-of-ten rung.
    if (lambda >= hard_lambda_max) break
    lambda <- min(lambda * 10, hard_lambda_max)
  }

  warning(sprintf(
    "Hard monotonicity not achieved within lambda_max = %.2g (max positive slope = %.2g). Returning closest fit.",
    hard_lambda_max,
    if (is.null(res)) NA_real_ else res$max_slope))
  attach_status(fit_current, converged = FALSE, lambda = lambda,
                max_slope = if (is.null(res)) NA_real_ else res$max_slope,
                iters = if (is.null(res)) NA_integer_ else res$iters,
                error = if (is.null(res)) NULL else res$error)
}

