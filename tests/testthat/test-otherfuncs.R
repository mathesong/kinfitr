context("test-otherfuncs")

data("pbr28")

meas <- 2

tac <- pbr28$tacs[[meas]]$STR
t_tac <- pbr28$tacs[[meas]]$Times / 60
dur_tac <- pbr28$tacs[[1]]$Duration / 60
weights <- pbr28$tacs[[meas]]$Weights

# SUV

test_that("SUV works with trapz", {
  SUVout <- SUV(tac, t_tac = t_tac, injRad = 150, bodymass = 85)
  expect_gt(SUVout$par$intSUV, 5000)
  expect_lt(SUVout$par$intSUV, 8000)
})

test_that("SUV works with durations", {
  SUVout <- SUV(tac, dur_tac = dur_tac, injRad = 150, bodymass = 85)
  expect_gt(SUVout$par$intSUV, 5000)
  expect_lt(SUVout$par$intSUV, 8000)
})

test_that("SUV works with trapz and frameStartEnd", {
  SUVout <- SUV(tac, t_tac = t_tac, injRad = 150,
                bodymass = 85, frameStartEnd = c(1,33))
  expect_gt(SUVout$par$intSUV, 4000)
  expect_lt(SUVout$par$intSUV, 6000)
})

test_that("SUV works with durations and frameStartEnd", {
  SUVout <- SUV(tac, dur_tac = dur_tac, injRad = 150,
                bodymass = 85, frameStartEnd = c(1,33))
  expect_gt(SUVout$par$intSUV, 4000)
  expect_lt(SUVout$par$intSUV, 6000)
})
