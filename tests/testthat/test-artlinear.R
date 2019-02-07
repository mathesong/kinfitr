context("test-artlinear.R")

data("pbr28")

meas <- 2

tac <- pbr28$tacs[[meas]]$FC
t_tac <- pbr28$tacs[[meas]]$Times / 60
input <- pbr28$input[[meas]]
weights <- pbr28$tacs[[meas]]$Weights
inpshift <- 0.1438066


#### Reversible

# Logan

test_that("Loganplot works", {
  Loganout <- Loganplot(
    t_tac, tac, input, 10, weights, inpshift = inpshift  )
  expect_lt(Loganout$par$Vt, 3)
  expect_gt(Loganout$par$Vt, 2)
  expect_true( any(class(plot(Loganout))=="ggplot") )
})

test_that("Loganplot with frameStartEnd works", {
  Loganout <- Loganplot(
    t_tac, tac, input, 10, weights, inpshift = inpshift,
    frameStartEnd = c(1,33))
  expect_lt(Loganout$par$Vt, 3)
  expect_gt(Loganout$par$Vt, 2)
  expect_lt(max(Loganout$tacs$Time), max(t_tac))
  expect_true( any(class(plot(Loganout))=="ggplot") )
})

# mlLogan

test_that("mlLoganplot works", {
  mlLoganout <- mlLoganplot(
    t_tac, tac, input, 10, weights, inpshift = inpshift  )
  expect_lt(mlLoganout$par$Vt, 3)
  expect_gt(mlLoganout$par$Vt, 2)
  expect_true( any(class(plot(mlLoganout))=="ggplot") )
})

test_that("mlLoganplot with frameStartEnd works", {
  mlLoganout <- mlLoganplot(
    t_tac, tac, input, 10, weights, inpshift = inpshift,
    frameStartEnd = c(1,33))
  expect_lt(mlLoganout$par$Vt, 3)
  expect_gt(mlLoganout$par$Vt, 2)
  expect_lt(max(mlLoganout$tacs$Time), max(t_tac))
  expect_true( any(class(plot(mlLoganout))=="ggplot") )
})


# MA1

test_that("MA1 works", {
  ma1out <- ma1(
    t_tac, tac, input, 10, weights, inpshift = inpshift  )
  expect_lt(ma1out$par$Vt, 3)
  expect_gt(ma1out$par$Vt, 2)
  expect_true( any(class(plot(ma1out))=="ggplot") )
})

test_that("MA1 with frameStartEnd works", {
  ma1out <- ma1(
    t_tac, tac, input, 10, weights, inpshift = inpshift,
    frameStartEnd = c(1,33))
  expect_lt(ma1out$par$Vt, 3)
  expect_gt(ma1out$par$Vt, 2)
  expect_lt(max(ma1out$tacs$Time), max(t_tac))
  expect_true( any(class(plot(ma1out))=="ggplot") )
})


# MA2

test_that("MA2 works", {
  ma2out <- ma2(
    t_tac, tac, input, weights, inpshift = inpshift  )
  expect_lt(ma2out$par$Vt, 3)
  expect_gt(ma2out$par$Vt, 2)
  expect_true( any(class(plot(ma2out))=="ggplot") )
})

test_that("MA2 with frameStartEnd works", {
  ma2out <- ma2(
    t_tac, tac, input, weights, inpshift = inpshift,
    frameStartEnd = c(1,33))
  expect_lt(ma2out$par$Vt, 3)
  expect_gt(ma2out$par$Vt, 2)
  expect_lt(max(ma2out$tacs$Time), max(t_tac))
  expect_true( any(class(plot(ma2out))=="ggplot") )
})




#### Irreversible

# Patlas

test_that("Patlakplot works", {
  Patlakout <- Patlakplot(
    t_tac, tac, input, 10, weights, inpshift = inpshift  )
  expect_gt(Patlakout$par$K, 0)
  expect_lt(Patlakout$par$K, 0.01)
  expect_true( any(class(plot(Patlakout))=="ggplot") )
})

test_that("Patlakplot with frameStartEnd works", {
  Patlakout <- Patlakplot(
    t_tac, tac, input, 10, weights, inpshift = inpshift,
    frameStartEnd = c(1,33))
  expect_gt(Patlakout$par$K, 0)
  expect_lt(Patlakout$par$K, 0.01)
  expect_lt(max(Patlakout$tacs$Time), max(t_tac))
  expect_true( any(class(plot(Patlakout))=="ggplot") )
})
