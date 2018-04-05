context("test-tcms.R")

data("pbr28")

meas <- 2

tac <- pbr28$tacs[[meas]]$FC
t_tac <- pbr28$tacs[[meas]]$Times / 60
input <- pbr28$input[[meas]]
weights <- pbr28$tacs[[meas]]$Weights
inpshift <- 0.1438066

# 1TCM

test_that("1TCM fitting delay works", {
  onetcmout <- onetcm(
    t_tac, tac, input, weights,
    K1.start = 0.1, k2.start = 0.05,
    K1.lower = 0.08, K1.upper = 0.12,
    k2.lower = 0.04, k2.upper = 0.06,
    inpshift.start = 0.15,
    inpshift.lower = 0.14,
    inpshift.upper = 0.16
  )
  expect_equal(round(onetcmout$par$Vt), 2)
})

test_that("1TCM with fitted delay works", {
  onetcmout <- onetcm(
    t_tac, tac, input, weights, inpshift = inpshift,
    K1.start = 0.1, k2.start = 0.05,
    K1.lower = 0.08, K1.upper = 0.12,
    k2.lower = 0.04, k2.upper = 0.06
  )
  expect_equal(round(onetcmout$par$Vt), 2)
})

test_that("1TCM with frameStartEnd works", {
  onetcmout <- onetcm(
    t_tac, tac, input, weights, inpshift = inpshift,
    K1.start = 0.1, k2.start = 0.05,
    K1.lower = 0.08, K1.upper = 0.12,
    k2.lower = 0.04, k2.upper = 0.06,
    frameStartEnd = c(1, round(0.75 * length(t_tac)))
  )
  expect_equal(round(onetcmout$par$Vt), 2)
  expect_lt(max(onetcmout$tacs$Time), max(t_tac))
})

test_that("1TCM fitting delay & multstart works", {
  onetcmout <- onetcm(
    t_tac, tac, input, weights,
    K1.start = 0.1, k2.start = 0.05,
    K1.lower = 0.08, K1.upper = 0.12,
    k2.lower = 0.04, k2.upper = 0.06,
    inpshift.start = 0.15,
    inpshift.lower = 0.14,
    inpshift.upper = 0.16,
    multstart_iter = 2
  )
  expect_equal(round(onetcmout$par$Vt), 2)
})

test_that("1TCM with fitted delay & multstart works", {
  onetcmout <- onetcm(
    t_tac, tac, input, weights, inpshift = inpshift,
    K1.start = 0.1, k2.start = 0.05,
    K1.lower = 0.08, K1.upper = 0.12,
    k2.lower = 0.04, k2.upper = 0.06,
    multstart_iter = 2
  )
  expect_equal(round(onetcmout$par$Vt), 2)
})


# 2TCM

test_that("2TCM fitting delay works", {
  twotcmout <- twotcm(
    t_tac, tac, input, weights,
    K1.start = 0.11, k2.start = 0.1,
    k3.start = 0.06, k4.start = 0.06,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.09, k2.upper = 0.10,
    k3.lower = 0.05, k3.upper = 0.07,
    k4.lower = 0.05, k4.upper = 0.06,
    inpshift.start = 0.15,
    inpshift.lower = 0.14,
    inpshift.upper = 0.16
  )
  expect_lt(twotcmout$par$Vt, 3)
  expect_gt(twotcmout$par$Vt, 2)
})

test_that("2TCM with fitted delay works", {
  twotcmout <- twotcm(
    t_tac, tac, input, weights, inpshift = inpshift,
    K1.start = 0.11, k2.start = 0.1,
    k3.start = 0.06, k4.start = 0.06,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.09, k2.upper = 0.10,
    k3.lower = 0.05, k3.upper = 0.07,
    k4.lower = 0.05, k4.upper = 0.06
  )
  expect_lt(twotcmout$par$Vt, 3)
  expect_gt(twotcmout$par$Vt, 2)
})

test_that("2TCM with frameStartEnd works", {
  twotcmout <- twotcm(
    t_tac, tac, input, weights, inpshift = inpshift,
    K1.start = 0.11, k2.start = 0.1,
    k3.start = 0.06, k4.start = 0.06,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.09, k2.upper = 0.10,
    k3.lower = 0.05, k3.upper = 0.07,
    k4.lower = 0.05, k4.upper = 0.06,
    frameStartEnd = c(1, round(0.75 * length(t_tac)))
  )
  expect_lt(twotcmout$par$Vt, 3)
  expect_gt(twotcmout$par$Vt, 2)
  expect_lt(max(twotcmout$tacs$Time), max(t_tac))
})

test_that("2TCM fitting delay & multstart works", {
  twotcmout <- twotcm(
    t_tac, tac, input, weights,
    K1.start = 0.11, k2.start = 0.1,
    k3.start = 0.06, k4.start = 0.06,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.09, k2.upper = 0.10,
    k3.lower = 0.05, k3.upper = 0.07,
    k4.lower = 0.05, k4.upper = 0.06,
    inpshift.start = 0.15,
    inpshift.lower = 0.14,
    inpshift.upper = 0.16,
    multstart_iter = 2
  )
  expect_lt(twotcmout$par$Vt, 3)
  expect_gt(twotcmout$par$Vt, 2)
})

test_that("2TCM with fitted delay & multstart works", {
  twotcmout <- twotcm(
    t_tac, tac, input, weights, inpshift = inpshift,
    K1.start = 0.11, k2.start = 0.1,
    k3.start = 0.06, k4.start = 0.06,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.09, k2.upper = 0.10,
    k3.lower = 0.05, k3.upper = 0.07,
    k4.lower = 0.05, k4.upper = 0.06,
    multstart_iter = 2
  )
  expect_lt(twotcmout$par$Vt, 3)
  expect_gt(twotcmout$par$Vt, 2)
})




# 2TCM1k

test_that("2TCM1k fitting delay works", {
  twotcm1kout <- twotcm1k(
    t_tac, tac, input, weights,
    K1.start = 0.11, k2.start = 0.14,
    k3.start = 0.16, k4.start = 0.13,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.12, k2.upper = 0.16,
    k3.lower = 0.14, k3.upper = 0.18,
    k4.lower = 0.12, k4.upper = 0.14,
    Kb.start = 0.12, Kb.upper = 0.14,
    Kb.lower = 0.1,
    inpshift.start = 0.15,
    inpshift.lower = 0.14,
    inpshift.upper = 0.16
  )
  expect_lt(twotcm1kout$par$Vt, 2)
  expect_gt(twotcm1kout$par$Vt, 1)
})

test_that("2TCM1k with fitted delay works", {
  twotcm1kout <- twotcm1k(
    t_tac, tac, input, weights, inpshift = inpshift,
    K1.start = 0.11, k2.start = 0.14,
    k3.start = 0.16, k4.start = 0.13,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.12, k2.upper = 0.16,
    k3.lower = 0.14, k3.upper = 0.18,
    k4.lower = 0.12, k4.upper = 0.14,
    Kb.start = 0.12, Kb.upper = 0.14,
    Kb.lower = 0.1
  )
  expect_lt(twotcm1kout$par$Vt, 2)
  expect_gt(twotcm1kout$par$Vt, 1)
})

test_that("2TCM1k with frameStartEnd works", {
  twotcm1kout <- twotcm1k(
    t_tac, tac, input, weights, inpshift = inpshift,
    K1.start = 0.11, k2.start = 0.14,
    k3.start = 0.16, k4.start = 0.13,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.12, k2.upper = 0.16,
    k3.lower = 0.14, k3.upper = 0.18,
    k4.lower = 0.12, k4.upper = 0.14,
    Kb.start = 0.12, Kb.upper = 0.14,
    Kb.lower = 0.1,
    frameStartEnd = c(1, round(0.75 * length(t_tac)))
  )
  expect_lt(twotcm1kout$par$Vt, 2)
  expect_gt(twotcm1kout$par$Vt, 1)
  expect_lt(max(twotcm1kout$tacs$Time), max(t_tac))
})

test_that("2TCM1k fitting delay & multstart works", {
  twotcm1kout <- twotcm1k(
    t_tac, tac, input, weights,
    K1.start = 0.11, k2.start = 0.14,
    k3.start = 0.16, k4.start = 0.13,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.12, k2.upper = 0.16,
    k3.lower = 0.14, k3.upper = 0.18,
    k4.lower = 0.12, k4.upper = 0.14,
    Kb.start = 0.12, Kb.upper = 0.14,
    Kb.lower = 0.1,
    inpshift.start = 0.15,
    inpshift.lower = 0.14,
    inpshift.upper = 0.16,
    multstart_iter = 2
  )
  expect_lt(twotcm1kout$par$Vt, 2)
  expect_gt(twotcm1kout$par$Vt, 1)
})

test_that("2TCM1k with fitted delay & multstart works", {
  twotcm1kout <- twotcm1k(
    t_tac, tac, input, weights, inpshift = inpshift,
    K1.start = 0.11, k2.start = 0.14,
    k3.start = 0.16, k4.start = 0.13,
    K1.lower = 0.09, K1.upper = 0.13,
    k2.lower = 0.12, k2.upper = 0.16,
    k3.lower = 0.14, k3.upper = 0.18,
    k4.lower = 0.12, k4.upper = 0.14,
    Kb.start = 0.12, Kb.upper = 0.14,
    Kb.lower = 0.1,
    multstart_iter = 2
  )
  expect_lt(twotcm1kout$par$Vt, 2)
  expect_gt(twotcm1kout$par$Vt, 1)
})
