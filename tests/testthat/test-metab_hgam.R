context("test-metab_hgam.R")

# Helper: generate parent-fraction-like data with controllable monotonicity
# violation. `violate` adds an upward quadratic bump that pushes the curve
# non-monotone.
gen_pf_data <- function(n_subj = 4, n_pet_per_subj = 2, n_time = 12,
                        violate = 0, seed = 1) {
  set.seed(seed)
  d <- expand.grid(
    subj_i = seq_len(n_subj),
    pet_i  = seq_len(n_pet_per_subj),
    time   = seq(1, 60, length.out = n_time)
  )
  subj_eff <- stats::rnorm(n_subj, 0, 0.4)
  pet_eff  <- stats::rnorm(n_subj * n_pet_per_subj, 0, 0.2)
  d$subject <- factor(sprintf("sub-%02d", d$subj_i))
  d$pet     <- factor(sprintf("sub-%02d_ses-%02d", d$subj_i, d$pet_i))
  pet_idx   <- (d$subj_i - 1) * n_pet_per_subj + d$pet_i
  d$logit_pf <- 2 + subj_eff[d$subj_i] + pet_eff[pet_idx] -
    0.07 * d$time +
    violate * (d$time / 30)^2 +
    stats::rnorm(nrow(d), 0, 0.1)
  d$parentFraction <- pmin(pmax(stats::plogis(d$logit_pf), 0.001), 0.999)
  d[, c("subject", "pet", "time", "parentFraction")]
}

# Helper: maximum slope of fitted curves on a grid
max_slope <- function(fit, data, time_var = "time", group_var = "pet",
                      n_grid = 100) {
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
  eta_lo <- stats::predict(fit, newdata = newdat,    type = "link")
  eta_hi <- stats::predict(fit, newdata = newdat_hi, type = "link")
  max((eta_hi - eta_lo) / eps)
}


test_that("metab_hgam returns a gam object", {
  d <- gen_pf_data()
  fit <- metab_hgam(d)
  expect_s3_class(fit, "gam")
})

test_that("monotone = 'none' is identical to direct mgcv::gam call", {
  d <- gen_pf_data()
  fit1 <- metab_hgam(d, monotone = "none")
  fit2 <- mgcv::gam(
    parentFraction ~ s(time, k = 8) + s(time, pet, bs = "fs", k = 5),
    data = d, family = mgcv::betar(link = "logit"), method = "REML"
  )
  expect_equal(stats::coef(fit1), stats::coef(fit2), tolerance = 1e-6)
})

test_that("monotone = 'hard' produces a non-increasing fit when data violates", {
  d <- gen_pf_data(violate = 1)
  # Sanity: violating data really does produce a positive slope unconstrained
  fit_none <- metab_hgam(d, monotone = "none")
  expect_gt(max_slope(fit_none, d), 0)

  # A successful hard fit must not fall back, so it must not warn.
  expect_no_warning(
    fit_hard <- metab_hgam(d, monotone = "hard", hard_tol = 1e-5)
  )
  expect_true(fit_hard$metab_hgam$converged)

  # The constraint is enforced on the n_grid evaluation points, so the slope
  # recorded there must be below hard_tol exactly.
  expect_lt(fit_hard$metab_hgam$max_slope, 1e-5)

  # Between those points a small excess is possible. Hold it to a tight
  # multiple of hard_tol rather than an arbitrary loose threshold.
  expect_lt(max_slope(fit_hard, d, n_grid = 2000), 2e-5)
})

test_that("monotone = 'hard' works when time is in seconds", {
  # crossprod(D) carries the units of time_var, so before the penalty was
  # normalised this failed at the first lambda and silently returned the
  # unconstrained fit.
  d <- gen_pf_data(violate = 1)
  d$time <- d$time * 60

  expect_no_warning(
    fit_hard <- metab_hgam(d, monotone = "hard", hard_tol = 1e-5 / 60)
  )
  expect_true(fit_hard$metab_hgam$converged)
  expect_lt(fit_hard$metab_hgam$max_slope, 1e-5 / 60)

  # The result must differ from the unconstrained fit.
  fit_none <- metab_hgam(d, monotone = "none")
  expect_false(isTRUE(all.equal(stats::coef(fit_hard), stats::coef(fit_none))))
})

test_that("monotone = 'hard' constrains a single group", {
  d <- gen_pf_data(n_subj = 1, n_pet_per_subj = 1, n_time = 20, violate = 1)
  expect_gt(max_slope(metab_hgam(d, monotone = "none"), d), 0)

  expect_no_warning(
    fit_hard <- metab_hgam(d, monotone = "hard", hard_tol = 1e-5)
  )
  expect_true(fit_hard$metab_hgam$converged)
  expect_lt(fit_hard$metab_hgam$max_slope, 1e-5)
})

test_that("the returned object records what the constraint achieved", {
  d <- gen_pf_data(violate = 1)

  st_none <- metab_hgam(d, monotone = "none")$metab_hgam
  expect_identical(st_none$monotone, "none")
  expect_true(is.na(st_none$converged))

  # Soft mode has no tolerance to meet, so it does not claim convergence.
  st_soft <- metab_hgam(d, monotone = "soft")$metab_hgam
  expect_identical(st_soft$monotone, "soft")
  expect_true(is.na(st_soft$converged))
  expect_gt(st_soft$max_slope, 0)

  st_hard <- metab_hgam(d, monotone = "hard", hard_tol = 1e-5)$metab_hgam
  expect_identical(st_hard$monotone, "hard")
  expect_true(st_hard$converged)
  expect_gt(st_hard$lambda, 0)
  expect_identical(st_hard$hard_tol, 1e-5)
  expect_identical(st_hard$n_grid, 60)
  expect_null(st_hard$error)
})

test_that("a failed hard fit reports converged = FALSE rather than warning only", {
  d <- gen_pf_data(violate = 1)
  # An unreachable tolerance with a low lambda ceiling must not be reported as
  # success.
  expect_warning(
    fit <- metab_hgam(d, monotone = "hard", hard_tol = 1e-12,
                      hard_lambda_init = 1e2, hard_lambda_max = 1e3),
    "not achieved"
  )
  expect_false(fit$metab_hgam$converged)
  expect_gt(fit$metab_hgam$max_slope, 0)
})

test_that("hard_lambda_max is attempted even when off the decade ladder", {
  set.seed(1)
  tt <- rep(seq(0, 1, length.out = 30), 4)
  d <- data.frame(
    time = tt, pet = factor("p1"),
    parentFraction = stats::plogis(0.2 + 0.6 * tt + stats::rnorm(120, 0, 0.01))
  )
  expect_warning(
    fit <- metab_hgam(d, formula = parentFraction ~ time,
                      family = stats::gaussian(), monotone = "hard",
                      hard_tol = 1e-12, hard_lambda_init = 1,
                      hard_lambda_max = 5),
    "not achieved"
  )
  expect_identical(fit$metab_hgam$lambda, 5)
})

test_that("invalid control parameters are rejected", {
  d <- gen_pf_data()
  expect_error(metab_hgam(d, monotone = "hard", n_grid = 1), "n_grid")
  expect_error(metab_hgam(d, monotone = "hard", n_grid = 0), "n_grid")
  expect_error(metab_hgam(d, monotone = "hard", max_iter = 0), "max_iter")
  expect_error(metab_hgam(d, monotone = "hard", hard_tol = 0), "hard_tol")
  expect_error(
    metab_hgam(d, monotone = "hard", hard_lambda_init = 1e6,
               hard_lambda_max = 1e3),
    "hard_lambda_max"
  )
})

test_that("parent fractions outside [0, 1] or non-finite are rejected", {
  d <- gen_pf_data()
  for (bad in c(-0.1, 1.1)) {
    dd <- d
    dd$parentFraction[1] <- bad
    expect_error(metab_hgam(dd), "outside \\[0, 1\\]")
  }
  dd <- d
  dd$parentFraction[1] <- Inf
  expect_error(metab_hgam(dd), "finite")
})

test_that("rows with missing values are dropped with a warning", {
  d <- gen_pf_data()
  d$time[1] <- NA
  expect_warning(fit <- metab_hgam(d, monotone = "soft"), "Dropped 1 row")
  expect_equal(nrow(fit$model), nrow(d) - 1)
})

test_that("a formula not depending on time_var is rejected in monotone modes", {
  d <- gen_pf_data(violate = 1)
  d$t2 <- d$time
  # The derivative grid perturbs time_var, so a formula built on another time
  # column would otherwise report an increasing curve as monotone.
  expect_error(
    metab_hgam(d, formula = parentFraction ~ t2, monotone = "hard"),
    "does not depend on"
  )
})

test_that("monotone = 'soft' reduces positive slope versus unconstrained", {
  d <- gen_pf_data(violate = 1)
  fit_none <- metab_hgam(d, monotone = "none")
  fit_soft <- metab_hgam(d, monotone = "soft")
  expect_lt(max_slope(fit_soft, d), max_slope(fit_none, d))
})

test_that("monotone modes work with a custom three-level formula", {
  # violate = 1, not 0.2: at 0.2 the unconstrained fit is already monotone, so
  # the penalised branch never ran and this test asserted nothing.
  d <- gen_pf_data(violate = 1)
  expect_gt(max_slope(metab_hgam(d, monotone = "none"), d), 0)

  fit <- metab_hgam(
    d, monotone = "hard", hard_tol = 1e-5,
    formula = parentFraction ~ s(time, k = 8) +
      s(time, subject, bs = "fs", k = 5) +
      s(time, pet, bs = "fs", k = 5)
  )
  expect_s3_class(fit, "gam")
  expect_true(fit$metab_hgam$converged)
  expect_lt(fit$metab_hgam$max_slope, 1e-5)
})

test_that("already-monotone data returns quickly without modification", {
  d <- gen_pf_data(violate = 0)
  fit_soft <- metab_hgam(d, monotone = "soft")
  fit_none <- metab_hgam(d, monotone = "none")
  # When unconstrained fit is already monotone, soft mode should not change it
  expect_equal(stats::coef(fit_soft), stats::coef(fit_none), tolerance = 1e-6)
})

test_that("custom time_var, parentFraction_var and group_var names work", {
  d <- gen_pf_data()
  names(d)[names(d) == "time"]           <- "framemidpoint"
  names(d)[names(d) == "pet"]            <- "scan_id"
  names(d)[names(d) == "parentFraction"] <- "pf"
  fit <- metab_hgam(
    d, monotone = "soft",
    time_var = "framemidpoint", parentFraction_var = "pf",
    group_var = "scan_id"
  )
  expect_s3_class(fit, "gam")
})

test_that("predict and summary work on the returned object", {
  d <- gen_pf_data()
  fit <- metab_hgam(d, monotone = "soft")
  preds <- stats::predict(fit, newdata = d[1:5, ], type = "response")
  expect_length(preds, 5)
  expect_true(all(preds > 0 & preds < 1))
  expect_silent(summary(fit))
})

test_that("missing time column raises informative error", {
  d <- gen_pf_data()
  d$time <- NULL
  expect_error(
    metab_hgam(d, monotone = "soft"),
    "time_var"
  )
})

test_that("verbose = TRUE produces messages", {
  d <- gen_pf_data(violate = 0.2)
  expect_message(
    metab_hgam(d, monotone = "soft", verbose = TRUE),
    "positive slopes"
  )
})


# Helper: parent-fraction data whose beta precision falls across the scan.
gen_pf_hetero <- function(n_pet = 8, n_time = 14, th0 = 5, th1 = -3,
                          violate = 0, seed = 5) {
  set.seed(seed)
  d <- expand.grid(pet_i = seq_len(n_pet),
                   time  = seq(1, 60, length.out = n_time))
  pe <- stats::rnorm(n_pet, 0, 0.3)
  d$pet <- factor(sprintf("p%02d", d$pet_i))
  u  <- (d$time - 1) / 59
  mu <- stats::plogis(2 + pe[d$pet_i] - 0.07 * d$time +
                        violate * (d$time / 30)^2)
  theta_i <- exp(th0 + th1 * u)
  d$parentFraction <- pmin(pmax(
    stats::rbeta(nrow(d), mu * theta_i, (1 - mu) * theta_i), 1e-6), 1 - 1e-6)
  d[, c("pet", "time", "parentFraction")]
}


test_that("theta_time = TRUE recovers a falling beta precision", {
  d <- gen_pf_hetero(th0 = 5, th1 = -3)
  fit <- metab_hgam(d, theta_time = TRUE)

  expect_true(fit$metab_hgam$theta_time)
  th <- fit$metab_hgam$theta
  expect_named(th, c("intercept", "slope"))

  # The precision must be found to fall, and both parameters should land near
  # the truth. The tolerances are loose because this is one simulated sample.
  expect_lt(th[["slope"]], 0)
  expect_equal(th[["intercept"]], 5, tolerance = 0.5)
  expect_equal(th[["slope"]], -3, tolerance = 1.5)
})

test_that("theta_time = TRUE beats a constant precision on REML score", {
  d <- gen_pf_hetero(th0 = 5, th1 = -3)
  f_const <- metab_hgam(d, theta_time = FALSE)
  f_lin   <- metab_hgam(d, theta_time = TRUE)
  expect_lt(f_lin$gcv.ubre, f_const$gcv.ubre)
})

test_that("theta_time = TRUE composes with hard monotonicity", {
  d <- gen_pf_hetero(th0 = 5, th1 = -3, violate = 1)
  expect_gt(max_slope(metab_hgam(d, monotone = "none", theta_time = TRUE), d), 0)

  expect_no_warning(
    fit <- metab_hgam(d, monotone = "hard", theta_time = TRUE, hard_tol = 1e-5)
  )
  expect_true(fit$metab_hgam$converged)
  expect_lt(fit$metab_hgam$max_slope, 1e-5)
  expect_lt(th <- fit$metab_hgam$theta[["slope"]], 0)
})

test_that("theta_time records nothing when it is switched off", {
  d <- gen_pf_hetero()
  st <- metab_hgam(d, theta_time = FALSE)$metab_hgam
  expect_false(st$theta_time)
  expect_null(st$theta)
})

test_that("theta_time = TRUE is rejected for a non-beta family", {
  d <- gen_pf_hetero()
  expect_error(
    metab_hgam(d, family = stats::gaussian(), theta_time = TRUE),
    "needs a beta family"
  )
  expect_error(metab_hgam(d, theta_time = "yes"), "must be TRUE or FALSE")
})

test_that("the linear-precision family matches numeric theta derivatives", {
  # The wrapper maps betar's log-scale theta derivatives to two parameters by
  # multiplying by powers of t_std. Check that against finite differences of
  # the deviance so a future change to betar cannot break it silently.
  set.seed(7)
  n <- 200
  u <- sort(stats::runif(n))
  mu <- stats::plogis(1 - 1.5 * u)
  th <- c(3, -1.2)
  theta_i <- exp(th[1] + th[2] * u)
  y  <- stats::rbeta(n, mu * theta_i, (1 - mu) * theta_i)
  wt <- rep(1, n)

  fam <- betar_lintheta(t_std = u)
  dev <- function(tt) sum(fam$dev.resids(y, mu, wt, tt))
  r   <- fam$Dd(y, mu, th, wt, level = 2)

  h <- 1e-4
  for (j in 1:2) {
    e <- c(0, 0); e[j] <- h
    expect_equal(sum(r$Dth[, j]), (dev(th + e) - dev(th - e)) / (2 * h),
                 tolerance = 1e-5)
  }
  pairs <- list(c(1, 1), c(1, 2), c(2, 2))
  for (k in seq_along(pairs)) {
    a <- pairs[[k]][1]; b <- pairs[[k]][2]
    ea <- c(0, 0); ea[a] <- h
    eb <- c(0, 0); eb[b] <- h
    num <- (dev(th + ea + eb) - dev(th + ea - eb) -
              dev(th - ea + eb) + dev(th - ea - eb)) / (4 * h * h)
    expect_equal(sum(r$Dth2[, k]), num, tolerance = 1e-4)
  }
})

test_that("the linear-precision family rejects a misaligned covariate", {
  # t_std is captured in the family's environment, so a length mismatch would
  # otherwise assign the wrong precision to each observation silently.
  d <- gen_pf_hetero()
  fam <- betar_lintheta(t_std = rep(0.5, 3))
  expect_error(
    mgcv::gam(parentFraction ~ s(time, k = 5), family = fam, data = d),
    "incomplete rows|length"
  )
})


test_that("bd_addfit rejects a metab_hgam fit and points to bd_addfitted", {
  # metab_hgam models many measurements jointly, so it needs its grouping
  # variable to predict. bd_getdata() predicts from time alone, so storing the
  # fit with bd_addfit() would be accepted and then fail later.
  d <- gen_pf_data()
  fit <- metab_hgam(d)
  expect_error(bd_addfit(list(), fit, "parentFraction"), "bd_addfitted")
  expect_error(bd_addfit(list(), fit, "parentFraction"), "metab_hgam")
})


test_that("a failed penalty step returns the last good fit, not the unconstrained one", {
  # Inject a failure into mgcv::gam for large lambda. Before this was fixed, a
  # pass that achieved nothing returned the globally unconstrained fit while
  # reporting the previous pass's slope.
  d <- gen_pf_data(violate = 1)
  orig <- mgcv::gam
  patched <- function(...) {
    a <- list(...)
    if (!is.null(a$G)) {
      lsp0 <- a$G$lsp0
      if (length(lsp0) && exp(utils::tail(lsp0, 1)) >= 1e6) {
        stop("simulated mgcv step failure")
      }
    }
    do.call(orig, list(...))
  }
  suppressWarnings(utils::assignInNamespace("gam", patched, "mgcv"))
  on.exit(suppressWarnings(utils::assignInNamespace("gam", orig, "mgcv")),
          add = TRUE)

  expect_warning(
    fit <- metab_hgam(d, monotone = "hard", hard_tol = 1e-9,
                      hard_lambda_init = 1e4, hard_lambda_max = 1e8),
    "numerically unstable"
  )
  expect_false(fit$metab_hgam$converged)
  expect_false(is.null(fit$metab_hgam$error))

  # Not the unconstrained fit, and the reported slope must match the returned fit.
  fit_none <- metab_hgam(d, monotone = "none")
  expect_false(isTRUE(all.equal(stats::coef(fit), stats::coef(fit_none))))
  expect_equal(fit$metab_hgam$max_slope, max_slope(fit, d, n_grid = 60),
               tolerance = 1e-6)
})

test_that("soft mode warns when a penalty step fails", {
  d <- gen_pf_data(violate = 1)
  orig <- mgcv::gam
  patched <- function(...) {
    if (!is.null(list(...)$G)) stop("simulated mgcv step failure")
    do.call(orig, list(...))
  }
  suppressWarnings(utils::assignInNamespace("gam", patched, "mgcv"))
  on.exit(suppressWarnings(utils::assignInNamespace("gam", orig, "mgcv")),
          add = TRUE)

  expect_warning(fit <- metab_hgam(d, monotone = "soft"), "failed numerically")
  expect_false(fit$metab_hgam$converged)
})

test_that("the linear-precision family reports the fitted precision downstream", {
  d <- gen_pf_hetero(th0 = 5, th1 = -3)
  fit <- metab_hgam(d, theta_time = TRUE)

  # family$family must stay a single string; betar's postproc would paste both
  # theta values into it and make it length 2.
  expect_length(fit$family$family, 1)

  # qf and rd must use the fitted precision, not betar's untouched scalar.
  # A precision of 1 would put qf(0.1, mu = 0.7) near 0.15. A short call cannot
  # know which rows are meant, so it warns and recycles from the first row.
  expect_warning(q1 <- fit$family$qf(0.1, mu = 0.7), "Recycling")
  expect_gt(q1, 0.4)

  # Precision falls over the scan, so the interval must widen.
  qn <- fit$family$qf(rep(0.1, nrow(d)), mu = rep(0.7, nrow(d)))
  expect_lt(qn[length(qn)], qn[1])

  expect_silent(invisible(capture.output(summary(fit))))

  # gam.check plots, so send the output to a null device rather than leaving an
  # Rplots.pdf behind in the test directory.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(capture.output(mgcv::gam.check(fit)), NA)
})

test_that("rows with NA in a custom formula column are dropped", {
  d <- gen_pf_data(violate = 1)
  d$subject <- factor(substr(as.character(d$pet), 1, 6))
  d$subject[1] <- NA
  expect_warning(
    fit <- metab_hgam(d, monotone = "hard", hard_tol = 1e-5,
                      formula = parentFraction ~ s(time, k = 8) +
                        s(time, subject, bs = "fs", k = 5)),
    "Dropped 1 row"
  )
  expect_equal(nrow(fit$model), nrow(d) - 1)
  expect_true(fit$metab_hgam$converged)
})

test_that("theta_time is rejected for a family with a fixed precision", {
  d <- gen_pf_hetero()
  expect_error(
    metab_hgam(d, family = mgcv::betar(theta = 5), theta_time = TRUE),
    "fixed precision"
  )
})


test_that("a returned fit's precision matches its own deviance", {
  # Every refit shares one family object whose theta lives in a mutable
  # environment, so a stored earlier iterate used to report a later iterate's
  # theta and its deviance no longer matched its coefficients.
  d <- gen_pf_hetero(violate = 1)
  for (tt in c(FALSE, TRUE)) {
    fit <- metab_hgam(d, theta_time = tt, monotone = "hard", hard_tol = 1e-5)
    th <- fit$family$getTheta()
    sat <- fit$family$saturated.ll(
      fit$y, fit$prior.weights, if (tt) th else exp(th))
    recomputed <- 2 * sat$f +
      sum(fit$family$dev.resids(fit$y, fit$fitted.values,
                                fit$prior.weights, th))
    expect_equal(recomputed, fit$deviance, tolerance = 1e-6)
  }
})

test_that("theta derivatives are zero where the log precision is clamped", {
  u  <- c(0.05, 0.4, 0.9)
  y  <- c(0.21, 0.52, 0.81)
  mu <- c(0.3, 0.55, 0.72)
  wt <- c(0.7, 1.3, 2)
  fam <- betar_lintheta(t_std = u)
  dev <- function(z) sum(fam$dev.resids(y, mu, wt, z))
  h <- 1e-5

  # Fully clamped in both directions, and partially clamped.
  for (th in list(c(-21, 0), c(21, 0), c(15, 10))) {
    r <- fam$Dd(y, mu, th, wt, level = 2)
    num <- vapply(1:2, function(j) {
      e <- c(0, 0); e[j] <- h
      (dev(th + e) - dev(th - e)) / (2 * h)
    }, numeric(1))
    expect_equal(as.numeric(colSums(r$Dth)), num, tolerance = 1e-4)
  }
})

test_that("a precision pinned at its cap is flagged as unidentified", {
  set.seed(77)
  n <- 600
  t <- stats::runif(n)
  mu <- stats::plogis(0.2 - 0.4 * t)
  # True log precision 25, above the family's cap of 20.
  y <- stats::rbeta(n, mu * exp(25), (1 - mu) * exp(25))
  d <- data.frame(time = t, pet = factor("p1"),
                  parentFraction = pmin(pmax(y, 1e-14), 1 - 1e-14))

  # This regime is deliberately pathological, and mgcv's optimiser either
  # converges into saturation or cannot evaluate the likelihood, depending on
  # the platform's BLAS. Both outcomes are correct; a bare mgcv error is not.
  # Assert on whichever happens.
  w <- NULL
  res <- tryCatch(
    withCallingHandlers(
      metab_hgam(d, formula = parentFraction ~ time, theta_time = TRUE),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)

  if (inherits(res, "error")) {
    # It must name the cause and the remedy, not leak mgcv's internal message.
    expect_match(conditionMessage(res), "theta_time = FALSE")
  } else {
    expect_true(any(grepl("not identified", w)))
    expect_false(res$metab_hgam$theta_identified)
  }
})

test_that("hard mode reports success when the returned fit meets the tolerance", {
  # A pass can meet hard_tol and then fail on a later iteration. The fit that
  # comes back satisfies the constraint, so that is a success.
  set.seed(1)
  d <- data.frame(time = rep(seq(0, 1, length.out = 50), 2), pet = factor("p1"))
  d$parentFraction <- 0.3 + 0.4 * d$time + stats::rnorm(nrow(d), 0, 0.002)

  ctr <- 0
  orig <- mgcv::gam
  patched <- function(...) {
    ctr <<- ctr + 1
    if (ctr == 5) stop("forced step failure")
    do.call(orig, list(...))
  }
  suppressWarnings(utils::assignInNamespace("gam", patched, "mgcv"))
  on.exit(suppressWarnings(utils::assignInNamespace("gam", orig, "mgcv")),
          add = TRUE)

  fit <- metab_hgam(d, formula = parentFraction ~ time,
                    family = stats::gaussian(), monotone = "hard",
                    hard_tol = 1e-3, hard_lambda_init = 1e12,
                    hard_lambda_max = 1e12)
  expect_true(fit$metab_hgam$converged)
  expect_lt(fit$metab_hgam$max_slope, 1e-3)
})

test_that("formula variables outside data are handled", {
  set.seed(12)
  d <- expand.grid(pet = factor(paste0("p", 1:4)),
                   time = seq(0, 1, length.out = 20))
  z <- stats::rnorm(nrow(d))
  d$parentFraction <- stats::plogis(1 - d$time + 0.1 * z +
                                      stats::rnorm(nrow(d), 0, 0.1))
  z[7] <- NA

  # mgcv resolves formula variables through parent.frame(), so put z where it
  # can see it, as it would be at the top level.
  assign("z", z, envir = globalenv())
  on.exit(rm("z", envir = globalenv()), add = TRUE)

  # Rows must not be dropped from data here: z keeps its original length, so
  # mgcv has to drop the row itself or it fails on differing lengths.
  fit <- metab_hgam(d, formula = parentFraction ~ time + z)
  expect_equal(nrow(fit$model), nrow(d) - 1)

  # theta_time cannot align its covariate in that case, so it says so.
  expect_error(
    metab_hgam(d, formula = parentFraction ~ time + z, theta_time = TRUE),
    "column of 'data'"
  )
})

test_that("validation does not reject legitimate inputs", {
  set.seed(1)
  d <- data.frame(time = 0:20, pet = factor("p1"))
  d$parentFraction <- 0.8 - 0.01 * d$time + stats::rnorm(21, 0, 1e-4)

  # A small-scaled but genuinely time-dependent column.
  expect_error(
    suppressWarnings(metab_hgam(d, formula = parentFraction ~ I(1e-9 * time),
                                family = stats::gaussian(), monotone = "hard")),
    NA
  )
  # A positive tolerance below machine epsilon.
  expect_error(
    suppressWarnings(metab_hgam(d, formula = parentFraction ~ time,
                                family = stats::gaussian(), monotone = "hard",
                                hard_tol = 1e-20)),
    NA
  )
  # monotone = "none" builds no grid, so one distinct time is fine.
  single <- data.frame(time = rep(5, 20), pet = factor(1:20),
                       parentFraction = stats::rbeta(20, 10, 3))
  expect_error(
    metab_hgam(single, formula = parentFraction ~ 1, monotone = "none"),
    NA
  )
})
