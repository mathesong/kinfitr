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

test_that("pRefIFS_1tcm runs and returns sensible parameters", {
  fit <- pRefIFS_1tcm(
    t_tac, tac_pref, tac_tgt,
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

test_that("k2prime is required", {
  expect_error(
    pRefIFS_1tcm(t_tac, tac_pref, tac_tgt,
                 t_blood = t_blood, blood = blood,
                 weights = weights),
    "k2prime must be supplied"
  )
  expect_error(
    pRefIFS_1tcm(t_tac, tac_pref, tac_tgt,
                 t_blood = t_blood, blood = blood,
                 k2prime = c(0.05, 0.1),
                 weights = weights),
    "single positive numeric"
  )
  expect_error(
    pRefIFS_1tcm(t_tac, tac_pref, tac_tgt,
                 t_blood = t_blood, blood = blood,
                 k2prime = -1,
                 weights = weights),
    "single positive numeric"
  )
})

test_that("analytical and symbolic derivatives are numerically equivalent", {
  shape <- pRefIFS_shape(
    t_tac, tac_pref,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val,
    weights = weights
  )
  coefs <- as.list(coef(shape$pref_fit$fit))
  grid <- seq(0, max(t_tac), length.out = 6000)
  d_an <- feng_1tc_tac_deriv(grid, coefs, method = "analytical")
  d_sy <- feng_1tc_tac_deriv(grid, coefs, method = "symbolic")
  expect_lt(max(abs(d_an - d_sy)), 1e-8)
})

test_that("scaling matches early blood AUC by construction", {
  shape <- pRefIFS_shape(
    t_tac, tac_pref,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val,
    weights = weights,
    scale_time = 5
  )
  sub <- shape$input$Time <= shape$scale_time
  auc_aif <- pracma::trapz(shape$input$Time[sub], shape$input$AIF[sub])
  b_interp <- stats::approx(t_blood, blood,
                            xout = shape$input$Time[sub], rule = 2)$y
  auc_blood <- pracma::trapz(shape$input$Time[sub], b_interp)
  expect_lt(abs(auc_aif - auc_blood) / auc_blood, 1e-6)
})

test_that("pRefIFS_shape produces nearly identical input to pRefIFS_1tcm", {
  # The pRef smoother uses random multistart so two independent runs differ
  # slightly; check agreement to within tight relative tolerance rather than
  # bit-equality.
  fit <- pRefIFS_1tcm(
    t_tac, tac_pref, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val, weights = weights
  )
  shape <- pRefIFS_shape(
    t_tac, tac_pref,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val, weights = weights
  )
  expect_equal(fit$scale_factor, shape$scale_factor, tolerance = 1e-2)
  peak <- max(abs(shape$input$AIF))
  expect_lt(max(abs(fit$input$AIF - shape$input$AIF)) / peak, 5e-3)
  expect_lt(max(abs(fit$pRefIFS$pRefIFS_scaled - shape$pRefIFS$pRefIFS_scaled)) /
              peak, 5e-3)
})

test_that("multstart_iter > 1 fit runs (regression: modelweights lookup)", {
  fit <- pRefIFS_1tcm(
    t_tac, tac_pref, tac_tgt,
    t_blood = t_blood, blood = blood,
    k2prime = k2prime_val,
    weights = weights,
    multstart_iter = 3
  )
  expect_s3_class(fit, "pRefIFS_1tcm")
  expect_true(all(is.finite(unlist(fit$par))))
})

test_that("Vt from pRefIFS_1tcm is in the same ballpark as onetcm with measured AIF", {
  fit <- pRefIFS_1tcm(
    t_tac, tac_pref, tac_tgt,
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
