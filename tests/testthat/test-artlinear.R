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

test_that("MA1 tstarfinder works", {
  suppressWarnings(
    tstar <- ma1_tstar(t_tac, lowroi, medroi, highroi,
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

test_that("Patlak tstarfinder works", {
  suppressWarnings(
    tstar <- Patlak_tstar(t_tac, lowroi, medroi, highroi,
      input,
      inpshift = inpshift, vB = 0.05
    )
  )
  expect_true(any(class(plot(tstar)) == "ggplot"))
})
