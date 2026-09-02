context("test-refreg")

data("simref")

set.seed(12345)
meas <- 2

reftac <- simref$tacs[[meas]]$Reference
roitac <- simref$tacs[[meas]]$ROI1
t_tac <- simref$tacs[[meas]]$Times
weights <- simref$tacs[[meas]]$Weights
dur <- simref$tacs[[meas]]$Duration

lowroi <- simref$tacs[[meas]]$ROI1
medroi <- simref$tacs[[meas]]$ROI2
highroi <- simref$tacs[[meas]]$ROI3


#### Reversible

# SRTM

test_that("srtm works", {
  srtmout <- srtm(t_tac, reftac, roitac, weights = weights)
  expect_gt(srtmout$par$bp, 1.5)
  expect_lt(srtmout$par$bp, 2.5)
  expect_true(any(class(plot(srtmout)) == "ggplot"))
})

test_that("srtm works with frameStartEnd", {
  srtmout <- srtm(t_tac, reftac, roitac,
    weights = weights,
    frameStartEnd = c(1, 33)
  )
  expect_gt(srtmout$par$bp, 1.5)
  expect_lt(srtmout$par$bp, 2.5)
  expect_true(any(class(plot(srtmout)) == "ggplot"))
})

test_that("srtm works with multstart", {
  srtmout <- srtm(t_tac, reftac, roitac,
    weights = weights,
    multstart_iter = 5
  )
  expect_gt(srtmout$par$bp, 1.5)
  expect_lt(srtmout$par$bp, 2.5)
  expect_true(any(class(plot(srtmout)) == "ggplot"))
})

test_that("srtm works with frameStartEnd and multstart", {
  srtmout <- srtm(t_tac, reftac, roitac,
    weights = weights,
    frameStartEnd = c(1, 33),
    multstart_iter = 5
  )
  expect_gt(srtmout$par$bp, 1.5)
  expect_lt(srtmout$par$bp, 2.5)
  expect_true(any(class(plot(srtmout)) == "ggplot"))
})

# SRTM2

test_that("srtm2 works", {
  srtm2out <- srtm2(t_tac, reftac, roitac, weights = weights)
  expect_gt(srtm2out$par$bp, 1.5)
  expect_lt(srtm2out$par$bp, 2.5)
  expect_true(any(class(plot(srtm2out)) == "ggplot"))
})

test_that("srtm2 works with frameStartEnd", {
  srtm2out <- srtm2(t_tac, reftac, roitac,
                  weights = weights,
                  frameStartEnd = c(1, 33)
  )
  expect_gt(srtm2out$par$bp, 1.5)
  expect_lt(srtm2out$par$bp, 2.5)
  expect_true(any(class(plot(srtm2out)) == "ggplot"))
})

test_that("srtm2 works with multstart", {
  srtm2out <- srtm2(t_tac, reftac, roitac,
                  weights = weights,
                  multstart_iter = 5
  )
  expect_gt(srtm2out$par$bp, 1.5)
  expect_lt(srtm2out$par$bp, 2.5)
  expect_true(any(class(plot(srtm2out)) == "ggplot"))
})

test_that("srtm2 works with frameStartEnd and multstart", {
  srtm2out <- srtm2(t_tac, reftac, roitac,
                  weights = weights,
                  frameStartEnd = c(1, 33),
                  multstart_iter = 5
  )
  expect_gt(srtm2out$par$bp, 1.5)
  expect_lt(srtm2out$par$bp, 2.5)
  expect_true(any(class(plot(srtm2out)) == "ggplot"))
})

test_that("srtm2 works with set k2prime", {
  srtm2out <- srtm2(t_tac, reftac, roitac, k2prime = 0.1,  weights = weights)
  expect_gt(srtm2out$par$bp, 1.5)
  expect_lt(srtm2out$par$bp, 2.5)
  expect_true(any(class(plot(srtm2out)) == "ggplot"))
})

test_that("srtm2 works with frameStartEnd and set k2prime", {
  srtm2out <- srtm2(t_tac, reftac, roitac, k2prime = 0.1,
                    weights = weights,
                    frameStartEnd = c(1, 33)
  )
  expect_gt(srtm2out$par$bp, 1.5)
  expect_lt(srtm2out$par$bp, 2.5)
  expect_true(any(class(plot(srtm2out)) == "ggplot"))
})

test_that("srtm2 works with multstart and set k2prime", {
  srtm2out <- srtm2(t_tac, reftac, roitac, k2prime = 0.1,
                    weights = weights,
                    multstart_iter = 5
  )
  expect_gt(srtm2out$par$bp, 1.5)
  expect_lt(srtm2out$par$bp, 2.5)
  expect_true(any(class(plot(srtm2out)) == "ggplot"))
})

test_that("srtm2 works with frameStartEnd and multstart and set k2prime", {
  srtm2out <- srtm2(t_tac, reftac, roitac, k2prime = 0.1,
                    weights = weights,
                    frameStartEnd = c(1, 33),
                    multstart_iter = 5
  )
  expect_gt(srtm2out$par$bp, 1.5)
  expect_lt(srtm2out$par$bp, 2.5)
  expect_true(any(class(plot(srtm2out)) == "ggplot"))
})


# FRTM

test_that("frtm works", {
  frtmout <- frtm(t_tac, reftac, roitac, weights = weights,
                  k4.upper=1e6) # Note: this model isn't really right here
  expect_gt(frtmout$par$bp, 1.5)
  expect_lt(frtmout$par$bp, 2.5)
  expect_true(any(class(plot(frtmout)) == "ggplot"))
})

test_that("frtm works with frameStartEnd", {
  frtmout <- frtm(t_tac, reftac, roitac,
                  weights = weights,
                  frameStartEnd = c(1, 33),
                  k4.upper=1e6 # Note: this model isn't really right here
  )
  expect_gt(frtmout$par$bp, 1.5)
  expect_lt(frtmout$par$bp, 2.5)
  expect_true(any(class(plot(frtmout)) == "ggplot"))
})

# test_that("frtm works with multstart", {
#   frtmout <- frtm(t_tac, reftac, roitac,
#                   weights = weights,
#                   multstart_iter = 2,
#                   k4.upper=1e6 # Note: this model isn't really right here
#   )
#   expect_gt(frtmout$par$bp, 1.5)
#   expect_lt(frtmout$par$bp, 2.5)
#   expect_true(any(class(plot(frtmout)) == "ggplot"))
# })
#
# test_that("frtm works with frameStartEnd and multstart", {
#   frtmout <- frtm(t_tac, reftac, roitac,
#                   weights = weights,
#                   frameStartEnd = c(1, 33),
#                   multstart_iter = 5
#   )
#   expect_gt(frtmout$par$bp, 1.5)
#   expect_lt(frtmout$par$bp, 2.5)
#   expect_true(any(class(plot(frtmout)) == "ggplot"))
# })


# refLogan

test_that("refLogan works", {
  refLoganout <- refLogan(t_tac, reftac, roitac, 0.1, 10, weights = weights)
  expect_gt(refLoganout$par$bp, 1.5)
  expect_lt(refLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refLoganout)) == "ggplot"))
})

test_that("refLogan works with frameStartEnd", {
  refLoganout <- refLogan(t_tac, reftac, roitac, 0.1, 10,
    weights = weights,
    frameStartEnd = c(1, 33)
  )
  expect_gt(refLoganout$par$bp, 1.5)
  expect_lt(refLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refLoganout)) == "ggplot"))
})

test_that("refLogan works with durations", {
  refLoganout <- refLogan(t_tac, reftac, roitac, 0.1, 10, weights = weights, dur=dur)
  expect_gt(refLoganout$par$bp, 1.5)
  expect_lt(refLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refLoganout)) == "ggplot"))
})

test_that("refLogan works with durations and frameStartEnd", {
  refLoganout <- refLogan(t_tac, reftac, roitac, 0.1, 10, weights = weights,
                          dur=dur, frameStartEnd = c(1, 33))
  expect_gt(refLoganout$par$bp, 1.5)
  expect_lt(refLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refLoganout)) == "ggplot"))
})

test_that("refLogan with time-based tstar works", {
  refLoganout <- refLogan(t_tac, reftac, roitac, 0.1, tstar = 30, tstar_type = "time", weights = weights)
  expect_gt(refLoganout$par$bp, 1.5)
  expect_lt(refLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refLoganout)) == "ggplot"))
})

test_that("refLogan tstarfinder works", {
  suppressWarnings(
    tstar <- refLogan_tstar(t_tac, reftac, lowroi, medroi, highroi,
      k2prime = 0.1
    )
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})


# refmlLogan

test_that("refmlLogan works", {
  refmlLoganout <- refmlLogan(t_tac, reftac, roitac, 0.1, 10, weights = weights)
  expect_gt(refmlLoganout$par$bp, 1.5)
  expect_lt(refmlLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refmlLoganout)) == "ggplot"))
})

test_that("refmlLogan works with frameStartEnd", {
  refmlLoganout <- refmlLogan(t_tac, reftac, roitac, 0.1, 10,
    weights = weights,
    frameStartEnd = c(1, 33)
  )
  expect_gt(refmlLoganout$par$bp, 1.5)
  expect_lt(refmlLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refmlLoganout)) == "ggplot"))
})

test_that("refmlLogan works with durations", {
  refmlLoganout <- refmlLogan(t_tac, reftac, roitac, 0.1, 10,
                              weights = weights, dur=dur)
  expect_gt(refmlLoganout$par$bp, 1.5)
  expect_lt(refmlLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refmlLoganout)) == "ggplot"))
})

test_that("refmlLogan works with durations and frameStartEnd", {
  refmlLoganout <- refmlLogan(t_tac, reftac, roitac, 0.1, 10,
                              weights = weights, dur=dur,
                              frameStartEnd = c(1, 33)
  )
  expect_gt(refmlLoganout$par$bp, 1.5)
  expect_lt(refmlLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refmlLoganout)) == "ggplot"))
})

test_that("refmlLogan with time-based tstar works", {
  refmlLoganout <- refmlLogan(t_tac, reftac, roitac, 0.1, tstar = 30, tstar_type = "time", weights = weights)
  expect_gt(refmlLoganout$par$bp, 1.5)
  expect_lt(refmlLoganout$par$bp, 2.5)
  expect_true(any(class(plot(refmlLoganout)) == "ggplot"))
})

test_that("refmlLogan tstarfinder works", {
  suppressWarnings(
    tstar <- refmlLogan_tstar(t_tac, reftac, lowroi, medroi, highroi,
      k2prime = 0.1
    )
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})


# MRTM1

test_that("mrtm1 works", {
  mrtm1out <- mrtm1(t_tac, reftac, roitac, weights = weights)
  expect_gt(mrtm1out$par$bp, 1.5)
  expect_lt(mrtm1out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm1out)) == "ggplot"))
})

test_that("mrtm1 works with frameStartEnd", {
  mrtm1out <- mrtm1(t_tac, reftac, roitac,
    weights = weights,
    frameStartEnd = c(1, 33)
  )
  expect_gt(mrtm1out$par$bp, 1.5)
  expect_lt(mrtm1out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm1out)) == "ggplot"))
})

test_that("mrtm1 works with tstar", {
  mrtm1out <- mrtm1(t_tac, reftac, roitac,
    weights = weights,
    tstar = 30
  )
  expect_gt(mrtm1out$par$bp, 1.5)
  expect_lt(mrtm1out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm1out)) == "ggplot"))
})

test_that("mrtm1 works with tstar and frameStartEnd", {
  mrtm1out <- mrtm1(t_tac, reftac, roitac,
    weights = weights,
    frameStartEnd = c(1, 33),
    tstar = 30
  )
  expect_gt(mrtm1out$par$bp, 1.5)
  expect_lt(mrtm1out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm1out)) == "ggplot"))
})

test_that("mrtm1 works with duration", {
  mrtm1out <- mrtm1(t_tac, reftac, roitac, weights = weights, dur=dur)
  expect_gt(mrtm1out$par$bp, 1.5)
  expect_lt(mrtm1out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm1out)) == "ggplot"))
})

test_that("mrtm1 works with duration and frameStartEnd", {
  mrtm1out <- mrtm1(t_tac, reftac, roitac, weights = weights, dur=dur,
                    frameStartEnd = c(1, 33))
  expect_gt(mrtm1out$par$bp, 1.5)
  expect_lt(mrtm1out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm1out)) == "ggplot"))
})

test_that("mrtm1 with time-based tstar works", {
  mrtm1out <- mrtm1(t_tac, reftac, roitac,
    weights = weights,
    tstar = 30, tstar_type = "time"
  )
  expect_gt(mrtm1out$par$bp, 1.5)
  expect_lt(mrtm1out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm1out)) == "ggplot"))
})

test_that("mrtm1 tstarfinder works", {
  suppressWarnings(
    tstar <- mrtm1_tstar(t_tac, reftac, lowroi, medroi, highroi)
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})


# MRTM2

# MRTM parameter definitions
#
# These check the identities that distinguish k2 from k2a, and R1 from DVR.
# Both were previously mislabelled: mrtm1() returned k2a as k2, and mrtm2()
# returned bp + 1 as R1.

test_that("mrtm1 recovers the parameters of a simulated SRTM TAC", {
  true_R1 <- 1.08
  true_k2 <- 0.146
  true_bp <- 1.85
  simtac <- srtm_model(t_tac, reftac, R1 = true_R1, k2 = true_k2, bp = true_bp)

  fit <- mrtm1(t_tac, reftac, simtac, weights = weights, dur = dur)

  expect_equal(fit$par$bp, true_bp, tolerance = 0.05)
  expect_equal(fit$par$R1, true_R1, tolerance = 0.05)
  expect_equal(fit$par$k2, true_k2, tolerance = 0.05)
  expect_equal(fit$par$k2prime, true_k2 / true_R1, tolerance = 0.05)
  expect_equal(fit$par$k2a, true_k2 / (1 + true_bp), tolerance = 0.05)
})

test_that("mrtm1 parameters are internally consistent", {
  fit <- mrtm1(t_tac, reftac, roitac, weights = weights, dur = dur)

  # k2 = R1 * k2prime by definition; returning k2a as k2 breaks this.
  expect_equal(fit$par$k2, fit$par$R1 * fit$par$k2prime)
  # k2a = k2 / (1 + BP).
  expect_equal(fit$par$k2a, fit$par$k2 / (1 + fit$par$bp))
  expect_lt(fit$par$k2a, fit$par$k2)
})

test_that("mrtm2 recovers the parameters of a simulated SRTM TAC", {
  true_R1 <- 1.08
  true_k2 <- 0.146
  true_bp <- 1.85
  simtac <- srtm_model(t_tac, reftac, R1 = true_R1, k2 = true_k2, bp = true_bp)

  fit <- mrtm2(t_tac, reftac, simtac,
    k2prime = true_k2 / true_R1, weights = weights, dur = dur
  )

  expect_equal(fit$par$bp, true_bp, tolerance = 0.05)
  expect_equal(fit$par$R1, true_R1, tolerance = 0.05)
  expect_equal(fit$par$k2, true_k2, tolerance = 0.05)
  expect_equal(fit$par$k2a, true_k2 / (1 + true_bp), tolerance = 0.05)
})

test_that("mrtm2 R1 is correct", {
  k2prime <- 0.1
  fit <- mrtm2(t_tac, reftac, roitac, k2prime, weights = weights, dur = dur)

  # R1 = Term1 / (-Term2) simplifies to bp + 1 for every input.
  expect_false(isTRUE(all.equal(fit$par$R1, fit$par$bp + 1)))
  expect_equal(fit$par$R1, fit$par$k2 / k2prime)
  expect_equal(fit$par$k2a, fit$par$k2 / (1 + fit$par$bp))
})

test_that("mrtm2 agrees with mrtm1 given mrtm1's own k2prime", {
  fit1 <- mrtm1(t_tac, reftac, roitac, weights = weights, dur = dur)
  fit2 <- mrtm2(t_tac, reftac, roitac,
    k2prime = fit1$par$k2prime, weights = weights, dur = dur
  )

  expect_equal(fit2$par$bp, fit1$par$bp, tolerance = 0.05)
  expect_equal(fit2$par$R1, fit1$par$R1, tolerance = 0.05)
  expect_equal(fit2$par$k2, fit1$par$k2, tolerance = 0.05)
})

test_that("mrtm2 works", {
  mrtm2out <- mrtm2(t_tac, reftac, roitac, 0.1, weights = weights)
  expect_gt(mrtm2out$par$bp, 1.5)
  expect_lt(mrtm2out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm2out)) == "ggplot"))
})

test_that("mrtm2 works with frameStartEnd", {
  mrtm2out <- mrtm2(t_tac, reftac, roitac, 0.1,
    weights = weights,
    frameStartEnd = c(1, 33)
  )
  expect_gt(mrtm2out$par$bp, 1.5)
  expect_lt(mrtm2out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm2out)) == "ggplot"))
})

test_that("mrtm2 works with tstar", {
  mrtm2out <- mrtm2(t_tac, reftac, roitac, 0.1,
    weights = weights,
    tstar = 10
  )
  expect_gt(mrtm2out$par$bp, 1.5)
  expect_lt(mrtm2out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm2out)) == "ggplot"))
})

test_that("mrtm2 works with tstar and frameStartEnd", {
  mrtm2out <- mrtm2(t_tac, reftac, roitac, 0.1,
    weights = weights,
    frameStartEnd = c(1, 33),
    tstar = 10
  )
  expect_gt(mrtm2out$par$bp, 1.5)
  expect_lt(mrtm2out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm2out)) == "ggplot"))
})

test_that("mrtm2 works with duration", {
  mrtm2out <- mrtm2(t_tac, reftac, roitac, 0.1, weights = weights, dur=dur)
  expect_gt(mrtm2out$par$bp, 1.5)
  expect_lt(mrtm2out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm2out)) == "ggplot"))
})

test_that("mrtm2 works with duration and frameStartEnd", {
  mrtm2out <- mrtm2(t_tac, reftac, roitac, 0.1, weights = weights,
                    dur=dur, frameStartEnd = c(1, 33))
  expect_gt(mrtm2out$par$bp, 1.5)
  expect_lt(mrtm2out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm2out)) == "ggplot"))
})

test_that("mrtm2 with time-based tstar works", {
  mrtm2out <- mrtm2(t_tac, reftac, roitac, 0.1,
    weights = weights,
    tstar = 30, tstar_type = "time"
  )
  expect_gt(mrtm2out$par$bp, 1.5)
  expect_lt(mrtm2out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm2out)) == "ggplot"))
})

test_that("mrtm2 tstarfinder works", {
  suppressWarnings(
    tstar <- mrtm2_tstar(t_tac, reftac, lowroi, medroi, highroi,
      k2prime = 0.1
    )
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})


# SRTM_V

input <- pbr28$input[[meas]]

newvals <- shift_timings(
  t_tac,
  roitac,
  input,
  inpshift = 0
)

bloodtac <- pracma::interp1(newvals$input$Time, newvals$input$Blood, t_tac)

test_that("srtm_v works", {
  srtm_vout <- srtm_v(t_tac, reftac, roitac, bloodtac, weights = weights,
                      vBt.lower = -5, vBr.lower = -5)
  # Note that using the crazy values of vBr and vBt avoids hitting limits
  expect_gt(srtm_vout$par$bp, 1.5)
  expect_lt(srtm_vout$par$bp, 2.5)
  expect_true(any(class(plot(srtm_vout)) == "ggplot"))
})

test_that("srtm_v works with vBr fitted", {
  srtm_vout <- srtm_v(t_tac, reftac, roitac, bloodtac,
    weights = weights, vBr = 0.05
  )
  expect_gt(srtm_vout$par$bp, 1.5)
  expect_lt(srtm_vout$par$bp, 2.5)
  expect_true(any(class(plot(srtm_vout)) == "ggplot"))
})

test_that("srtm_v works with frameStartEnd", {
  srtm_vout <- srtm_v(t_tac, reftac, roitac, bloodtac,
    weights = weights,
    frameStartEnd = c(1, 33),
    vBt.lower = -5, vBr.lower = -5,
    vBt.upper = 0.5, vBr.upper = 0.5)
  # Note that using the crazy values of vBr and vBt avoids hitting limits
  expect_gt(srtm_vout$par$bp, 1.5)
  expect_lt(srtm_vout$par$bp, 2.5)
  expect_true(any(class(plot(srtm_vout)) == "ggplot"))
})

test_that("srtm_v works with multstart", {
  srtm_vout <- srtm_v(t_tac, reftac, roitac, bloodtac,
    weights = weights,
    multstart_iter = 5,
    bp.lower = 0.5, bp.upper = 3,
    vBr = 0.05
  )
  expect_gt(srtm_vout$par$bp, 1.5)
  expect_lt(srtm_vout$par$bp, 2.5)
  expect_true(any(class(plot(srtm_vout)) == "ggplot"))
})

test_that("srtm_v works with frameStartEnd and multstart", {
  srtm_vout <- srtm_v(t_tac, reftac, roitac, bloodtac,
    weights = weights,
    frameStartEnd = c(1, 33),
    multstart_iter = 5,
    bp.lower = 0.5, bp.upper = 3,
    vBr = 0.05
  )
  expect_gt(srtm_vout$par$bp, 1.5)
  expect_lt(srtm_vout$par$bp, 2.5)
  expect_true(any(class(plot(srtm_vout)) == "ggplot"))
})



#### Irreversible

test_that("refPatlak works", {
  refPatlakout <- refPatlak(t_tac, reftac, roitac, 10, weights = weights)
  expect_gt(refPatlakout$par$K, -1)
  expect_lt(refPatlakout$par$K, 0.5)
  expect_true(any(class(plot(refPatlakout)) == "ggplot"))
})

test_that("refPatlak works with frameStartEnd", {
  refPatlakout <- refPatlak(t_tac, reftac, roitac, 10,
    weights = weights,
    frameStartEnd = c(1, 33)
  )
  expect_gt(refPatlakout$par$K, -1)
  expect_lt(refPatlakout$par$K, 0.5)
  expect_true(any(class(plot(refPatlakout)) == "ggplot"))
})

test_that("refPatlak works with duration", {
  refPatlakout <- refPatlak(t_tac, reftac, roitac, 10, weights = weights, dur=dur)
  expect_gt(refPatlakout$par$K, -1)
  expect_lt(refPatlakout$par$K, 0.5)
  expect_true(any(class(plot(refPatlakout)) == "ggplot"))
})

test_that("refPatlak works with duration and frameStartEnd", {
  refPatlakout <- refPatlak(t_tac, reftac, roitac, 10, weights = weights,
                            dur=dur, frameStartEnd = c(1, 33))
  expect_gt(refPatlakout$par$K, -1)
  expect_lt(refPatlakout$par$K, 0.5)
  expect_true(any(class(plot(refPatlakout)) == "ggplot"))
})

test_that("refPatlak with time-based tstar works", {
  refPatlakout <- refPatlak(t_tac, reftac, roitac, tstar = 30, tstar_type = "time", weights = weights)
  expect_gt(refPatlakout$par$K, -1)
  expect_lt(refPatlakout$par$K, 0.5)
  expect_true(any(class(plot(refPatlakout)) == "ggplot"))
})

test_that("refPatlak tstarfinder works", {
  suppressWarnings(
    tstar <- refPatlak_tstar(t_tac, reftac, lowroi, medroi, highroi)
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})


#### Reference TAC Fitting Functions

# spline_tac

test_that("spline_tac works", {
  fit <- spline_tac(t_tac, reftac, weights = weights)
  expect_true("spline_tac" %in% class(fit))
  expect_true("kinfit" %in% class(fit))
  expect_true(!is.null(fit$par$t0))
  expect_true(!is.null(fit$gam_fit))
  expect_gt(fit$par$t0, 0)
  expect_lt(fit$par$t0, 5)
  expect_equal(length(fit$tacs$TAC_fitted), length(reftac)) # tidyinput adds zero frame if needed
  expect_true(any(class(plot(fit)) == "ggplot"))
})

test_that("spline_tac works with frameStartEnd", {
  fit <- spline_tac(t_tac, reftac, weights = weights, frameStartEnd = c(1, 30))
  expect_true(!is.null(fit$par$t0))
  expect_equal(nrow(fit$tacs), 30)
  expect_true(any(class(plot(fit)) == "ggplot"))
})

test_that("spline_tac works with timeStartEnd", {
  fit <- spline_tac(t_tac, reftac, weights = weights, timeStartEnd = c(0, 50))
  expect_true(!is.null(fit$par$t0))
  expect_true(max(fit$tacs$Time) <= 50)
  expect_true(any(class(plot(fit)) == "ggplot"))
})

test_that("spline_tac predict method works at original times", {
  fit <- spline_tac(t_tac, reftac, weights = weights)
  pred <- predict(fit)
  expect_equal(length(pred), nrow(fit$tacs))
  expect_true(all(pred >= 0)) # All predictions should be non-negative
  expect_true(all(pred[fit$tacs$Time < fit$par$t0] == 0)) # Before t0 should be zero
})

test_that("spline_tac predict method works with newdata", {
  fit <- spline_tac(t_tac, reftac, weights = weights)
  new_times <- seq(0, max(t_tac), length.out = 100)
  pred <- predict(fit, newdata = list(t_tac = new_times))
  expect_equal(length(pred), 100)
  expect_true(all(pred >= 0)) # All predictions should be non-negative
  expect_true(all(pred[new_times < fit$par$t0] == 0)) # Before t0 should be zero
})

test_that("spline_tac fitted values are accessible", {
  fit <- spline_tac(t_tac, reftac, weights = weights)
  expect_true("TAC_fitted" %in% names(fit$tacs))
  expect_equal(length(fit$tacs$TAC_fitted), nrow(fit$tacs))
  expect_true(all(fit$tacs$TAC_fitted >= 0))
})

test_that("spline_tac weights are accessible", {
  fit <- spline_tac(t_tac, reftac, weights = weights)
  expect_true(!is.null(fit$weights))
  expect_equal(length(fit$weights), nrow(fit$tacs))
})

test_that("spline_tac k parameter controls wiggliness", {
  fit_smooth <- spline_tac(t_tac, reftac, weights = weights, k = 5)
  fit_flex <- spline_tac(t_tac, reftac, weights = weights, k = 15)

  # Check both fits work
  expect_true(!is.null(fit_smooth$gam_fit))
  expect_true(!is.null(fit_flex$gam_fit))

  # Flexible fit should have higher effective degrees of freedom
  edf_smooth <- sum(fit_smooth$gam_fit$edf)
  edf_flex <- sum(fit_flex$gam_fit$edf)
  expect_lt(edf_smooth, edf_flex)
})


# feng_1tc_tac

test_that("feng_1tc_tac works", {
  fit <- feng_1tc_tac(t_tac, reftac, weights = weights)
  expect_true("feng_1tc_tac" %in% class(fit))
  expect_true("kinfit" %in% class(fit))
  expect_true(!is.null(fit$par))
  expect_true("t0" %in% names(fit$par))
  expect_gt(fit$par$t0, 0)
  expect_lt(fit$par$t0, 5)
  expect_equal(length(fit$tacs$Reference_fitted), length(reftac)) # tidyinput adds zero frame if needed
  expect_true(any(class(plot(fit)) == "ggplot"))
})

test_that("feng_1tc_tac works with frameStartEnd", {
  fit <- feng_1tc_tac(t_tac, reftac, weights = weights, frameStartEnd = c(1, 30))
  expect_true(!is.null(fit$par$t0))
  expect_equal(nrow(fit$tacs), 30)
  expect_true(any(class(plot(fit)) == "ggplot"))
})

test_that("feng_1tc_tac works with timeStartEnd", {
  fit <- feng_1tc_tac(t_tac, reftac, weights = weights, timeStartEnd = c(0, 50))
  expect_true(!is.null(fit$par$t0))
  expect_true(max(fit$tacs$Time) <= 50)
  expect_true(any(class(plot(fit)) == "ggplot"))
})

test_that("feng_1tc_tac works without fitting t0", {
  fit <- feng_1tc_tac(t_tac, reftac, weights = weights, fit_t0 = FALSE)
  expect_true("feng_1tc_tac" %in% class(fit))
  expect_false("t0" %in% names(fit$par))
  expect_true(any(class(plot(fit)) == "ggplot"))
})

test_that("feng_1tc_tac fitted values are accessible", {
  fit <- feng_1tc_tac(t_tac, reftac, weights = weights)
  expect_true("Reference_fitted" %in% names(fit$tacs))
  expect_equal(length(fit$tacs$Reference_fitted), nrow(fit$tacs))
  expect_true(all(fit$tacs$Reference_fitted >= 0))
})


# SUVR

test_that("suvr works", {
  suvrout <- suvr(t_tac, reftac, roitac, dur = dur)

  expect_equal(suvrout$par$SUVR, suvrout$par$SUV_AUC / suvrout$par$SUV_ref_AUC)
  expect_gt(suvrout$par$SUVR, 1)
  expect_equal(suvrout$par$n_frames, length(t_tac))
  expect_true(all(suvrout$tacs$Included))
  expect_equal(suvrout$model, "suvr")
  expect_s3_class(suvrout, "kinfit")
})

test_that("suvr outcomes are consistent with SUV called directly", {
  suvrout <- suvr(t_tac, reftac, roitac, dur = dur, frameStartEnd = c(5, 15))

  suvout <- suv(roitac,
    t_tac = t_tac, dur = dur,
    frameStartEnd = c(5, 15)
  )

  expect_equal(suvrout$par$SUV_AUC, suvout$par$SUV_AUC)
  expect_equal(suvrout$par$SUV, suvout$par$SUV)
  expect_equal(suvrout$par$window_duration, sum(dur[5:15]))
  expect_equal(suvrout$par$n_frames, 11)
  # The mean is the integral over the duration which was integrated
  expect_equal(suvrout$par$SUV, suvrout$par$SUV_AUC / suvrout$par$window_duration)
  expect_equal(suvrout$par$SUV_ref, suvrout$par$SUV_ref_AUC / suvrout$par$window_duration)
})

test_that("the dose cancels in the SUVR but not in the SUV", {
  no_dose <- suvr(t_tac, reftac, roitac, dur = dur, timeStartEnd = c(20, 60))
  with_dose <- suvr(t_tac, reftac, roitac,
    dur = dur, timeStartEnd = c(20, 60),
    injRad = 150, bodymass = 85
  )

  expect_equal(no_dose$par$SUVR, with_dose$par$SUVR)

  expect_equal(no_dose$par$SUV_denominator, 1)
  expect_equal(with_dose$par$SUV_denominator, 150 / 85)
  expect_equal(with_dose$par$SUV, no_dose$par$SUV / (150 / 85))
  expect_equal(with_dose$par$SUV_AUC, no_dose$par$SUV_AUC / (150 / 85))
})

test_that("suvr time windows select whole frames by their midpoints", {
  suvrout <- suvr(t_tac, reftac, roitac, dur = dur, timeStartEnd = c(20, 60))

  expect_equal(suvrout$tacs$Included, t_tac >= 20 & t_tac <= 60)
  # The whole frames count, so the duration is theirs and not the 40 minutes
  # spanned by the request
  expect_equal(suvrout$par$window_duration, sum(dur[t_tac >= 20 & t_tac <= 60]))
})

test_that("suvr accepts one-sided and empty windows", {
  from_frame_5 <- suvr(t_tac, reftac, roitac, dur = dur, frameStartEnd = c(5, NA))
  expect_equal(from_frame_5$window$start, 5L)
  expect_equal(from_frame_5$window$end, length(t_tac))

  until_frame_5 <- suvr(t_tac, reftac, roitac, dur = dur, frameStartEnd = c(NA, 5))
  expect_equal(until_frame_5$window$start, 1L)
  expect_equal(until_frame_5$window$end, 5L)

  # c(0, 0) means no window at all
  all_frames <- suvr(t_tac, reftac, roitac, dur = dur, timeStartEnd = c(0, 0))
  expect_true(all(all_frames$tacs$Included))

  expect_error(
    suvr(t_tac, reftac, roitac, dur = dur, timeStartEnd = c(500, 600)),
    "No frames fall within"
  )
})

test_that("suvr works without durations", {
  suvrout <- suvr(t_tac, reftac, roitac)

  expect_gt(suvrout$par$SUVR, 1)
  expect_equal(suvrout$par$SUV_AUC, suv(roitac, t_tac = t_tac)$par$SUV_AUC)
})

test_that("suvr rejects mismatched inputs", {
  expect_error(suvr(t_tac, reftac[-1], roitac, dur = dur), "not equal")
  expect_error(suvr(t_tac, reftac, roitac, dur = dur[-1]), "not equal")
})

test_that("plot_suvrfit shades only the included frames", {
  suvrout <- suvr(t_tac, reftac, roitac, dur = dur, frameStartEnd = c(5, 15))

  suvrplot <- plot_suvrfit(suvrout, roiname = "ROI1")
  expect_s3_class(suvrplot, "ggplot")

  # The whole TAC is drawn, and the shaded layer carries the included frames of
  # both regions as one rectangle each
  expect_equal(nrow(suvrplot$data), 2 * length(t_tac))
  expect_s3_class(suvrplot$layers[[1]]$geom, "GeomRect")
  rects <- suvrplot$layers[[1]]$data
  expect_equal(nrow(rects), 2 * 11)
  expect_true(all(rects$Included))

  # The two sets of rectangle areas are the integrals whose ratio is the SUVR
  areas <- tapply((rects$xmax - rects$xmin) * rects$Radioactivity, rects$Region, sum)
  expect_equal(as.numeric(areas[["ROI1"]]), suvrout$par$SUV_AUC)
  expect_equal(as.numeric(areas[["Reference"]]), suvrout$par$SUV_ref_AUC)
  expect_equal(as.numeric(areas[["ROI1"]] / areas[["Reference"]]), suvrout$par$SUVR)

  # And plot() dispatches to it through the kinfit method
  expect_s3_class(plot(suvrout, roiname = "ROI1"), "ggplot")
})

test_that("plot_suvrfit shades under the curves without durations", {
  suvrout <- suvr(t_tac, reftac, roitac, frameStartEnd = c(5, 15))

  suvrplot <- plot_suvrfit(suvrout)
  expect_s3_class(suvrplot$layers[[1]]$geom, "GeomRibbon")
})
