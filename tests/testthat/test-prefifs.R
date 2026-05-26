context("test-prefifs.R")

data("pbr28")

meas <- 1
t_tac    <- pbr28$tacs[[meas]]$Times / 60
tac_pref <- pbr28$tacs[[meas]]$CBL
tac_tgt  <- pbr28$tacs[[meas]]$FC
weights  <- pbr28$tacs[[meas]]$Weights
bd       <- pbr28$procblood[[meas]]
t_blood  <- bd$Time / 60
blood    <- bd$Cbl_dispcorr

k2prime_val <- 0.05

set.seed(42)

# Fit the Feng+1TC smoother to the pseudo-reference TAC once and reuse the
# parameters across every test below (matches real workflows where many
# target ROIs share one pRef).
pref_fit <- feng_1tc_tac(t_tac, tac_pref, weights)
pref_par <- pref_fit$par

test_that("pRefIFS_1tcm runs and returns sensible parameters", {
  fit <- pRefIFS_1tcm(
    t_tac, pref_par, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val,
    weights = weights
  )

  expect_s3_class(fit, "pRefIFS_1tcm")
  expect_s3_class(fit, "kinfit")
  expect_true(all(c("K1", "k2", "Vt") %in% names(fit$par)))
  expect_true(all(is.finite(unlist(fit$par))))
  expect_gt(fit$par$K1, 0)
  expect_gt(fit$par$k2, 0)
  expect_gt(fit$par$Vt, 0)
  expect_identical(fit$k2prime, k2prime_val)
  expect_true(any(class(plot_kinfit(fit)) == "ggplot"))
  expect_true(any(class(plot(fit)) == "ggplot"))
})

test_that("pRefIFS_1tcm accepts a feng_1tc_tac fit object directly", {
  fit_a <- pRefIFS_1tcm(
    t_tac, pref_par, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val, weights = weights
  )
  fit_b <- pRefIFS_1tcm(
    t_tac, pref_fit, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val, weights = weights
  )
  expect_equal(fit_a$par$Vt, fit_b$par$Vt)
})

test_that("pref_par missing required columns errors clearly", {
  bad <- data.frame(A = 1, B = 1)  # missing C, alpha, beta, gamma, Ph1, Th1
  expect_error(
    pRefIFS_1tcm(t_tac, bad, tac_tgt,
                 t_blood = t_blood, blood = blood,
                 k2prime = k2prime_val, weights = weights),
    "missing required column"
  )
})

test_that("k2prime is required", {
  expect_error(
    pRefIFS_1tcm(t_tac, pref_par, tac_tgt,
                 t_blood = t_blood, blood = blood,
                 weights = weights),
    "k2prime must be supplied"
  )
  expect_error(
    pRefIFS_1tcm(t_tac, pref_par, tac_tgt,
                 t_blood = t_blood, blood = blood,
                 k2prime = c(0.05, 0.1),
                 weights = weights),
    "single positive numeric"
  )
  expect_error(
    pRefIFS_1tcm(t_tac, pref_par, tac_tgt,
                 t_blood = t_blood, blood = blood,
                 k2prime = -1,
                 weights = weights),
    "single positive numeric"
  )
})

test_that("analytical and symbolic derivatives are numerically equivalent", {
  coefs <- as.list(pref_par)
  if (is.null(coefs$t0)) coefs$t0 <- 0
  grid <- seq(0, max(t_tac), length.out = 6000)
  d_an <- feng_1tc_tac_deriv(grid, coefs, method = "analytical")
  d_sy <- feng_1tc_tac_deriv(grid, coefs, method = "symbolic")
  expect_lt(max(abs(d_an - d_sy)), 1e-8)
})

test_that("scaling is SF = auc_blood / auc_pref by construction", {
  # Note: SF is set BEFORE negative values are clamped to zero. Once clamping
  # removes negative-positive cancellation in the early window, the post-clamp
  # AIF AUC over [0, scale_time] generally exceeds auc_blood. The mathematical
  # invariant of the construction is SF = auc_blood / auc_pref, and that is
  # what the diagnostic fields should match.
  shape <- pRefIFS_shape(
    t_tac, pref_par,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val,
    scale_time = 5
  )
  sub <- shape$input$Time <= shape$scale_time
  b_interp <- stats::approx(t_blood, blood,
                            xout = shape$input$Time[sub], rule = 2)$y
  auc_blood_independent <- pracma::trapz(shape$input$Time[sub], b_interp)
  expect_equal(shape$auc_blood_early, auc_blood_independent, tolerance = 1e-10)
  expect_equal(shape$scale_factor,
               shape$auc_blood_early / shape$auc_pref_early,
               tolerance = 1e-10)
})

test_that("pRefIFS_shape produces the same input used by pRefIFS_1tcm", {
  # With the refactor the smoother is no longer re-fit, so two calls with
  # the same pref_par should be bitwise identical.
  fit <- pRefIFS_1tcm(
    t_tac, pref_par, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val, weights = weights
  )
  shape <- pRefIFS_shape(
    t_tac, pref_par,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val
  )
  expect_equal(fit$scale_factor, shape$scale_factor)
  expect_equal(fit$input$AIF, shape$input$AIF)
  expect_equal(fit$pRefIFS$pRefIFS_scaled, shape$pRefIFS$pRefIFS_scaled)
})

test_that("multstart_iter > 1 fit runs (regression: modelweights lookup)", {
  fit <- pRefIFS_1tcm(
    t_tac, pref_par, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val,
    weights = weights,
    multstart_iter = 3
  )
  expect_s3_class(fit, "pRefIFS_1tcm")
  expect_true(all(is.finite(unlist(fit$par))))
})

test_that("plot(fit) works when frameStartEnd is used (regression)", {
  fit <- pRefIFS_1tcm(
    t_tac, pref_par, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val, weights = weights,
    frameStartEnd = c(1, 25)
  )
  expect_true(any(class(plot(fit)) == "ggplot"))
})

test_that("timeStartEnd out of range errors clearly", {
  expect_error(
    pRefIFS_1tcm(t_tac, pref_par, tac_tgt,
                 t_blood = t_blood, blood = blood,
                 k2prime = k2prime_val, weights = weights,
                 timeStartEnd = c(200, 300)),
    "outside the t_tac range"
  )
})

test_that("NA in t_blood or blood is dropped with a warning", {
  expect_warning(
    pRefIFS_shape(t_tac, pref_par,
                  t_blood = c(0, 1, NA, 5, 10),
                  blood = c(0, 5, 10, 8, 4),
                  k2prime = k2prime_val),
    "NA/non-finite"
  )
})

test_that("Early flat-extrapolation is replaced with zero padding when t_blood[1] > 0", {
  # Drop the t=0 entry in the blood data, leaving t_blood[1] > 0 with
  # a non-zero peak value early on, then check SF matches the case where
  # we explicitly prepend a (0, 0) entry.
  bd_late_t <- bd$Time[bd$Time > 0] / 60
  bd_late_v <- bd$Cbl_dispcorr[bd$Time > 0]
  bd_late_v[1] <- 30   # simulate peak-like first sample
  shape_late <- pRefIFS_shape(t_tac, pref_par,
                              t_blood = bd_late_t, blood = bd_late_v,
                              k2prime = k2prime_val)
  shape_pad  <- pRefIFS_shape(t_tac, pref_par,
                              t_blood = c(0, bd_late_t),
                              blood = c(0, bd_late_v),
                              k2prime = k2prime_val)
  expect_equal(shape_late$scale_factor, shape_pad$scale_factor, tolerance = 1e-10)
})

test_that("auc_blood guard fires when early blood AUC is non-positive", {
  zero_blood <- rep(0, length(t_blood))
  expect_error(
    pRefIFS_shape(t_tac, pref_par,
                  t_blood = t_blood, blood = zero_blood,
                  k2prime = k2prime_val),
    "AUC of the supplied blood is non-positive"
  )
})

test_that("Output includes diagnostic AUC fields", {
  fit <- pRefIFS_1tcm(
    t_tac, pref_par, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val, weights = weights
  )
  expect_true(all(c("auc_blood_early", "auc_pref_early",
                    "auc_target_early", "frac_clamped") %in% names(fit)))
  expect_gt(fit$auc_blood_early, 0)
  expect_gt(fit$auc_pref_early, 0)
  expect_gt(fit$auc_target_early, 0)
  # By construction SF = auc_blood / auc_pref
  expect_equal(fit$scale_factor, fit$auc_blood_early / fit$auc_pref_early,
               tolerance = 1e-8)
})

test_that("Negative pRef-IFS values are clamped to zero in the AIF", {
  # Force the situation by passing k2prime much smaller than the smoother's Th1
  fit <- expect_message(
    pRefIFS_1tcm(t_tac, pref_par, tac_tgt,
                 t_blood = t_blood, blood = blood,
                 k2prime = 0.005, weights = weights),
    "clamped to 0"
  )
  expect_gte(min(fit$input$AIF), 0)
  expect_gt(fit$frac_clamped, 0)
})

test_that("coerce_pref_par informs when t0 is defaulted", {
  pref_no_t0 <- pref_par
  pref_no_t0$t0 <- NULL
  expect_message(
    coerce_pref_par(pref_no_t0),
    "defaulting to t0 = 0"
  )
})

test_that("Vt from pRefIFS_1tcm is in the same ballpark as onetcm with measured AIF", {
  fit <- pRefIFS_1tcm(
    t_tac, pref_par, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val, weights = weights
  )
  input_true <- blood_interp(
    t_blood, blood,
    t_blood, bd$Cpl_metabcorr,
    t_parentfrac = 1, parentfrac = 1
  )
  fit_true <- onetcm(t_tac, tac_tgt, input_true, weights,
                     inpshift = 0, vB = 0,
                     K1.start = 0.1, k2.start = 0.05)
  # pRef-IFS Vt should be within ~20% of the gold-standard 1TC Vt for this tracer/ROI
  expect_lt(abs(fit$par$Vt - fit_true$par$Vt) / fit_true$par$Vt, 0.2)
})
