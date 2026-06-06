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
  # Hard fit should be monotone within tolerance
  fit_hard <- metab_hgam(d, monotone = "hard", hard_tol = 1e-5)
  expect_lt(max_slope(fit_hard, d), 1e-3)
})

test_that("monotone = 'soft' reduces positive slope versus unconstrained", {
  d <- gen_pf_data(violate = 1)
  fit_none <- metab_hgam(d, monotone = "none")
  fit_soft <- metab_hgam(d, monotone = "soft")
  expect_lt(max_slope(fit_soft, d), max_slope(fit_none, d))
})

test_that("monotone modes work with a custom three-level formula", {
  d <- gen_pf_data(violate = 0.2)
  fit <- metab_hgam(
    d, monotone = "hard", hard_tol = 1e-5,
    formula = parentFraction ~ s(time, k = 8) +
      s(time, subject, bs = "fs", k = 5) +
      s(time, pet, bs = "fs", k = 5)
  )
  expect_s3_class(fit, "gam")
  expect_lt(max_slope(fit, d), 1e-3)
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
