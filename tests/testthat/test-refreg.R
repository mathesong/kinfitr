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
    tstarIncludedFrames = 30
  )
  expect_gt(mrtm1out$par$bp, 1.5)
  expect_lt(mrtm1out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm1out)) == "ggplot"))
})

test_that("mrtm1 works with tstar and frameStartEnd", {
  mrtm1out <- mrtm1(t_tac, reftac, roitac,
    weights = weights,
    frameStartEnd = c(1, 33),
    tstarIncludedFrames = 30
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
    tstarIncludedFrames = 10
  )
  expect_gt(mrtm2out$par$bp, 1.5)
  expect_lt(mrtm2out$par$bp, 2.5)
  expect_true(any(class(plot(mrtm2out)) == "ggplot"))
})

test_that("mrtm2 works with tstar and frameStartEnd", {
  mrtm2out <- mrtm2(t_tac, reftac, roitac, 0.1,
    weights = weights,
    frameStartEnd = c(1, 33),
    tstarIncludedFrames = 10
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

test_that("refPatlak tstarfinder works", {
  suppressWarnings(
    tstar <- refPatlak_tstar(t_tac, reftac, lowroi, medroi, highroi)
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})
