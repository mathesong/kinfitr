context("test-artlinear.R")

data("pbr28")

set.seed(12345)
meas <- 2

tac <- pbr28$tacs[[meas]]$FC
t_tac <- pbr28$tacs[[meas]]$Times / 60
input <- pbr28$input[[meas]]
weights <- pbr28$tacs[[meas]]$Weights
inpshift <- 0.1438066
dur <- pbr28$tacs[[meas]]$Duration/60

lowroi <- pbr28$tacs[[meas]]$FC
medroi <- pbr28$tacs[[meas]]$CBL
highroi <- pbr28$tacs[[meas]]$THA


#### Reversible

# Logan

test_that("Loganplot works", {
  Loganout <- Loganplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift
  )
  expect_lt(Loganout$par$Vt, 3)
  expect_gt(Loganout$par$Vt, 2)
  expect_true(any(class(plot(Loganout)) == "ggplot"))
})

test_that("Loganplot with frameStartEnd works", {
  Loganout <- Loganplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift,
    frameStartEnd = c(1, 33)
  )
  expect_lt(Loganout$par$Vt, 3)
  expect_gt(Loganout$par$Vt, 2)
  expect_lt(max(Loganout$tacs$Time), max(t_tac))
  expect_true(any(class(plot(Loganout)) == "ggplot"))
})

test_that("Loganplot works with durations", {
  Loganout <- Loganplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift, dur=dur
  )
  expect_lt(Loganout$par$Vt, 3)
  expect_gt(Loganout$par$Vt, 2)
  expect_true(any(class(plot(Loganout)) == "ggplot"))
})

test_that("Loganplot works with durations and frameStartEnd", {
  Loganout <- Loganplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift, dur=dur,
    frameStartEnd = c(1, 33)
  )
  expect_lt(Loganout$par$Vt, 3)
  expect_gt(Loganout$par$Vt, 2)
  expect_true(any(class(plot(Loganout)) == "ggplot"))
})

test_that("Loganplot with time-based tstar works", {
  Loganout <- Loganplot(
    t_tac, tac, input, tstar = 30, tstar_type = "time", weights,
    inpshift = inpshift
  )
  expect_lt(Loganout$par$Vt, 3)
  expect_gt(Loganout$par$Vt, 2)
  expect_true(any(class(plot(Loganout)) == "ggplot"))
})

test_that("Loganplot tstarfinder works", {
  suppressWarnings(
    tstar <- Logan_tstar(t_tac, lowroi, medroi, highroi,
      input,
      inpshift = inpshift, vB = 0.05
    )
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})

# mlLogan

test_that("mlLoganplot works", {
  mlLoganout <- mlLoganplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift
  )
  expect_lt(mlLoganout$par$Vt, 3)
  expect_gt(mlLoganout$par$Vt, 2)
  expect_true(any(class(plot(mlLoganout)) == "ggplot"))
})

test_that("mlLoganplot with frameStartEnd works", {
  mlLoganout <- mlLoganplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift,
    frameStartEnd = c(1, 33)
  )
  expect_lt(mlLoganout$par$Vt, 3)
  expect_gt(mlLoganout$par$Vt, 2)
  expect_lt(max(mlLoganout$tacs$Time), max(t_tac))
  expect_true(any(class(plot(mlLoganout)) == "ggplot"))
})

test_that("mlLoganplot works with durations", {
  mlLoganout <- mlLoganplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift, dur=dur
  )
  expect_lt(mlLoganout$par$Vt, 3)
  expect_gt(mlLoganout$par$Vt, 2)
  expect_true(any(class(plot(mlLoganout)) == "ggplot"))
})

test_that("mlLoganplot works with durations and frameStartEnd", {
  mlLoganout <- mlLoganplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift, dur=dur,
    frameStartEnd = c(1, 33)
  )
  expect_lt(mlLoganout$par$Vt, 3)
  expect_gt(mlLoganout$par$Vt, 2)
  expect_true(any(class(plot(mlLoganout)) == "ggplot"))
})

test_that("mlLoganplot with time-based tstar works", {
  mlLoganout <- mlLoganplot(
    t_tac, tac, input, tstar = 30, tstar_type = "time", weights,
    inpshift = inpshift
  )
  expect_lt(mlLoganout$par$Vt, 3)
  expect_gt(mlLoganout$par$Vt, 2)
  expect_true(any(class(plot(mlLoganout)) == "ggplot"))
})

test_that("mlLoganplot tstarfinder works", {
  suppressWarnings(
    tstar <- mlLogan_tstar(t_tac, lowroi, medroi, highroi,
      input,
      inpshift = inpshift, vB = 0.05
    )
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})


# MA1

test_that("MA1 works", {
  ma1out <- ma1(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift
  )
  expect_lt(ma1out$par$Vt, 3)
  expect_gt(ma1out$par$Vt, 2)
  expect_true(any(class(plot(ma1out)) == "ggplot"))
})

test_that("MA1 with frameStartEnd works", {
  ma1out <- ma1(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift,
    frameStartEnd = c(1, 33)
  )
  expect_lt(ma1out$par$Vt, 3)
  expect_gt(ma1out$par$Vt, 2)
  expect_lt(max(ma1out$tacs$Time), max(t_tac))
  expect_true(any(class(plot(ma1out)) == "ggplot"))
})

test_that("MA1 works with durations", {
  ma1out <- ma1(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift, dur=dur
  )
  expect_lt(ma1out$par$Vt, 3)
  expect_gt(ma1out$par$Vt, 2)
  expect_true(any(class(plot(ma1out)) == "ggplot"))
})

test_that("MA1 works with durations and frameStartEnd", {
  ma1out <- ma1(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift, dur=dur,
    frameStartEnd = c(1, 33)
  )
  expect_lt(ma1out$par$Vt, 3)
  expect_gt(ma1out$par$Vt, 2)
  expect_true(any(class(plot(ma1out)) == "ggplot"))
})

test_that("MA1 with time-based tstar works", {
  ma1out <- ma1(
    t_tac, tac, input, tstar = 30, tstar_type = "time", weights,
    inpshift = inpshift
  )
  expect_lt(ma1out$par$Vt, 3)
  expect_gt(ma1out$par$Vt, 2)
  expect_true(any(class(plot(ma1out)) == "ggplot"))
})

test_that("MA1 tstarfinder works", {
  suppressWarnings(
    tstar <- ma1_tstar(t_tac, lowroi, medroi, highroi,
      input,
      inpshift = inpshift, vB = 0.05
    )
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})


# LEGA

test_that("LEGA works", {
  LEGAout <- LEGA(
    t_tac, dur, tac, input, 10, weights,
    inpshift = inpshift
  )
  expect_lt(LEGAout$par$Vt, 3)
  expect_gt(LEGAout$par$Vt, 2)
  expect_true(any(class(plot(LEGAout)) == "ggplot"))
})

test_that("LEGA with frameStartEnd works", {
  LEGAout <- LEGA(
    t_tac, dur, tac, input, 10, weights,
    inpshift = inpshift,
    frameStartEnd = c(1, 33)
  )
  expect_lt(LEGAout$par$Vt, 3)
  expect_gt(LEGAout$par$Vt, 2)
  expect_lt(max(LEGAout$tacs$Time), max(t_tac))
  expect_true(any(class(plot(LEGAout)) == "ggplot"))
})

test_that("LEGA with time-based tstar works", {
  LEGAout <- LEGA(
    t_tac, dur, tac, input, tstar = 30, tstar_type = "time", weights,
    inpshift = inpshift
  )
  expect_lt(LEGAout$par$Vt, 3)
  expect_gt(LEGAout$par$Vt, 2)
  expect_true(any(class(plot(LEGAout)) == "ggplot"))
})

test_that("LEGA works without weights", {
  LEGAout <- LEGA(
    t_tac, dur, tac, input, 10,
    inpshift = inpshift
  )
  expect_lt(LEGAout$par$Vt, 3)
  expect_gt(LEGAout$par$Vt, 2)
})

test_that("LEGA errors when dur is not supplied", {
  expect_error(
    LEGA(t_tac = t_tac, tac = tac, input = input, tstar = 10, weights = weights),
    "required"
  )
})

test_that("LEGA accepts Vt.start and gamma.start", {
  LEGAout <- LEGA(
    t_tac, dur, tac, input, 10, weights,
    inpshift = inpshift,
    Vt.start = 2.5, gamma.start = -40
  )
  expect_lt(LEGAout$par$Vt, 3)
  expect_gt(LEGAout$par$Vt, 2)
  # Supplying both starts and disabling the fallback bypasses the reference fit
  LEGAout2 <- LEGA(
    t_tac, dur, tac, input, 10, weights,
    inpshift = inpshift,
    Vt.start = 2.5, gamma.start = -40, fallback = "none"
  )
  expect_lt(LEGAout2$par$Vt, 3)
  expect_gt(LEGAout2$par$Vt, 2)
})

test_that("LEGA stability fallback uses the reference Vt", {
  # Force the fallback by providing an implausibly small reference Vt
  suppressMessages(
    LEGAout <- LEGA(
      t_tac, dur, tac, input, 10, weights,
      inpshift = inpshift, Vt.start = 0.5
    )
  )
  expect_true(LEGAout$fallback_used)
  expect_equal(LEGAout$par$Vt, 0.5)
})

test_that("LEGA fallback method can be selected", {
  fit_ma1 <- LEGA(t_tac, dur, tac, input, 10, weights, inpshift = inpshift, fallback = "MA1")
  fit_log <- LEGA(t_tac, dur, tac, input, 10, weights, inpshift = inpshift, fallback = "Logan")
  fit_none <- LEGA(t_tac, dur, tac, input, 10, weights, inpshift = inpshift, fallback = "none")
  # The healthy estimate is unchanged regardless of fallback method
  expect_false(fit_ma1$fallback_used)
  expect_equal(fit_ma1$par$Vt, fit_log$par$Vt, tolerance = 1e-6)
  expect_equal(fit_ma1$par$Vt, fit_none$par$Vt, tolerance = 1e-6)
})

test_that("LEGA reduces the negative bias of the Logan plot", {
  LEGAout <- LEGA(
    t_tac, dur, tac, input, 10, weights,
    inpshift = inpshift
  )
  Loganout <- Loganplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift, dur = dur
  )
  # LEGA is bias-free, the Logan plot is negatively biased, so LEGA >= Logan
  expect_gte(LEGAout$par$Vt, Loganout$par$Vt)
})

test_that("LEGA returns a finite standard error and intercept", {
  LEGAout <- LEGA(
    t_tac, dur, tac, input, 10, weights,
    inpshift = inpshift
  )
  expect_true(is.finite(LEGAout$par.se$Vt.se))
  expect_gt(LEGAout$par.se$Vt.se, 0)
  expect_true(is.finite(LEGAout$gamma))
  # fit object plugs into the shared kinfitr machinery
  expect_true(is.finite(maxpercres(LEGAout)))
  expect_true(is.finite(summary(LEGAout$fit)$r.squared))
})

test_that("LEGA fit object supports residuals and fitted", {
  LEGAout <- LEGA(
    t_tac, dur, tac, input, 10, weights,
    inpshift = inpshift
  )
  expect_equal(length(residuals(LEGAout$fit)), length(fitted(LEGAout$fit)))
  # fitted values are the noise-free tissue concentrations Z* (offset included)
  expect_equal(
    as.numeric(fitted(LEGAout$fit)),
    LEGAout$fitvals$Target_fitted
  )
})

test_that("LEGA tstarfinder works", {
  suppressWarnings(
    tstar <- LEGA_tstar(t_tac, dur, lowroi, medroi, highroi,
      input,
      inpshift = inpshift, vB = 0.05
    )
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})


# MA2

test_that("MA2 works", {
  ma2out <- ma2(
    t_tac, tac, input, weights,
    inpshift = inpshift
  )
  expect_lt(ma2out$par$Vt, 3)
  expect_gt(ma2out$par$Vt, 2)
  expect_true(any(class(plot(ma2out)) == "ggplot"))
})

test_that("MA2 with frameStartEnd works", {
  ma2out <- ma2(
    t_tac, tac, input, weights,
    inpshift = inpshift,
    frameStartEnd = c(1, 33)
  )
  expect_lt(ma2out$par$Vt, 3)
  expect_gt(ma2out$par$Vt, 2)
  expect_lt(max(ma2out$tacs$Time), max(t_tac))
  expect_true(any(class(plot(ma2out)) == "ggplot"))
})

test_that("MA2 works with durations", {
  ma2out <- ma2(
    t_tac, tac, input, weights,
    inpshift = inpshift, dur=dur
  )
  expect_lt(ma2out$par$Vt, 3)
  expect_gt(ma2out$par$Vt, 2)
  expect_true(any(class(plot(ma2out)) == "ggplot"))
})

test_that("MA2 with durations and frameStartEnd", {
  ma2out <- ma2(
    t_tac, tac, input, weights,
    inpshift = inpshift, dur=dur,
    frameStartEnd = c(1, 33)
  )
  expect_lt(ma2out$par$Vt, 3)
  expect_gt(ma2out$par$Vt, 2)
  expect_lt(max(ma2out$tacs$Time), max(t_tac))
  expect_true(any(class(plot(ma2out)) == "ggplot"))
})

# Linearised 2TCM including vB

test_that("lin2tcm works", {
  lin2tcmout <- lin2tcm(
    t_tac, tac, input, weights,
    inpshift = inpshift
  )
  expect_lt(lin2tcmout$par$Vt, 3)
  expect_gt(lin2tcmout$par$Vt, 2)
  expect_true(any(class(plot(lin2tcmout)) == "ggplot"))
})

test_that("lin2tcm with frameStartEnd works", {
  lin2tcmout <- lin2tcm(
    t_tac, tac, input, weights,
    inpshift = inpshift,
    frameStartEnd = c(1, 33)
  )
  expect_lt(lin2tcmout$par$Vt, 3)
  expect_gt(lin2tcmout$par$Vt, 2)
  expect_lt(max(lin2tcmout$tacs$Time), max(t_tac))
  expect_true(any(class(plot(lin2tcmout)) == "ggplot"))
})

test_that("lin2tcm works with durations", {
  lin2tcmout <- lin2tcm(
    t_tac, tac, input, weights,
    inpshift = inpshift, dur=dur
  )
  expect_lt(lin2tcmout$par$Vt, 3)
  expect_gt(lin2tcmout$par$Vt, 2)
  expect_true(any(class(plot(lin2tcmout)) == "ggplot"))
})

test_that("lin2tcm with durations and frameStartEnd", {
  lin2tcmout <- lin2tcm(
    t_tac, tac, input, weights,
    inpshift = inpshift, dur=dur,
    frameStartEnd = c(1, 33)
  )
  expect_lt(lin2tcmout$par$Vt, 3)
  expect_gt(lin2tcmout$par$Vt, 2)
  expect_lt(max(lin2tcmout$tacs$Time), max(t_tac))
  expect_true(any(class(plot(lin2tcmout)) == "ggplot"))
})

# Linearised 2TCM excluding vB

test_that("lin2tcm works", {
  lin2tcmout <- lin2tcm(
    t_tac, tac, input, weights,
    inpshift = inpshift, vB = 0.05,
  )
  expect_lt(lin2tcmout$par$Vt, 3)
  expect_gt(lin2tcmout$par$Vt, 2)
  expect_true(any(class(plot(lin2tcmout)) == "ggplot"))
})

test_that("lin2tcm with frameStartEnd works", {
  lin2tcmout <- lin2tcm(
    t_tac, tac, input, weights,
    inpshift = inpshift, vB = 0.05,
    frameStartEnd = c(1, 33)
  )
  expect_lt(lin2tcmout$par$Vt, 3)
  expect_gt(lin2tcmout$par$Vt, 2)
  expect_lt(max(lin2tcmout$tacs$Time), max(t_tac))
  expect_true(any(class(plot(lin2tcmout)) == "ggplot"))
})

test_that("lin2tcm works with durations", {
  lin2tcmout <- lin2tcm(
    t_tac, tac, input, weights, vB = 0.05,
    inpshift = inpshift, dur=dur
  )
  expect_lt(lin2tcmout$par$Vt, 3)
  expect_gt(lin2tcmout$par$Vt, 2)
  expect_true(any(class(plot(lin2tcmout)) == "ggplot"))
})

test_that("lin2tcm with durations and frameStartEnd", {
  lin2tcmout <- lin2tcm(
    t_tac, tac, input, weights, vB = 0.05,
    inpshift = inpshift, dur=dur,
    frameStartEnd = c(1, 33)
  )
  expect_lt(lin2tcmout$par$Vt, 3)
  expect_gt(lin2tcmout$par$Vt, 2)
  expect_lt(max(lin2tcmout$tacs$Time), max(t_tac))
  expect_true(any(class(plot(lin2tcmout)) == "ggplot"))
})


# Linearised 2TCM for inpshift profiling

test_that("lin2tcm inpshift profiling", {

  is1 <- lin2tcm_inpshiftProfile(t_tac, tac, input, weights)
  is2 <- lin2tcm_inpshiftProfile(t_tac, tac, input, weights, dur = dur)
  is3 <- lin2tcm_inpshiftProfile(t_tac, tac, input, weights, vB=0.05,
                          frameStartEnd = c(1,15))

  expect_true(any(class(is1) == "ggplot"))
  expect_true(any(class(is2) == "ggplot"))
  expect_true(any(class(is3) == "ggplot"))
})




#### Irreversible

# Patlak

test_that("Patlakplot works", {
  Patlakout <- Patlakplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift
  )
  expect_gt(Patlakout$par$K, 0)
  expect_lt(Patlakout$par$K, 0.015)
  expect_true(any(class(plot(Patlakout)) == "ggplot"))
})

test_that("Patlakplot with frameStartEnd works", {
  Patlakout <- Patlakplot(
    t_tac, tac, input, 10, weights,
    inpshift = inpshift,
    frameStartEnd = c(1, 33)
  )
  expect_gt(Patlakout$par$K, 0)
  expect_lt(Patlakout$par$K, 0.01)
  expect_lt(max(Patlakout$tacs$Time), max(t_tac))
  expect_true(any(class(plot(Patlakout)) == "ggplot"))
})

test_that("Patlakplot with time-based tstar works", {
  Patlakout <- Patlakplot(
    t_tac, tac, input, tstar = 30, tstar_type = "time", weights,
    inpshift = inpshift
  )
  expect_gt(Patlakout$par$K, 0)
  expect_lt(Patlakout$par$K, 0.015)
  expect_true(any(class(plot(Patlakout)) == "ggplot"))
})

test_that("Patlak tstarfinder works", {
  suppressWarnings(
    tstar <- Patlak_tstar(t_tac, lowroi, medroi, highroi,
      input,
      inpshift = inpshift, vB = 0.05
    )
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})
