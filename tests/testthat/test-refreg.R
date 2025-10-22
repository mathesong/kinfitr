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
