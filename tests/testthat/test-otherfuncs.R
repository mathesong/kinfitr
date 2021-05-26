context("test-otherfuncs")

data("pbr28")

meas <- 2

tac <- pbr28$tacs[[meas]]$STR
t_tac <- pbr28$tacs[[meas]]$Times / 60
dur_tac <- pbr28$tacs[[meas]]$Duration / 60
weights <- pbr28$tacs[[meas]]$Weights

# SUV

test_that("SUV works with trapz", {
  SUVout <- SUV(tac, t_tac = t_tac, injRad = 150, bodymass = 85)
  expect_gt(SUVout$par$intSUV, 200)
  expect_lt(SUVout$par$intSUV, 500)
})

test_that("SUV works with durations", {
  SUVout <- SUV(tac, dur_tac = dur_tac, injRad = 150, bodymass = 85)
  expect_gt(SUVout$par$intSUV, 150)
  expect_lt(SUVout$par$intSUV, 500)
})

test_that("SUV works with trapz and frameStartEnd", {
  SUVout <- SUV(tac,
    t_tac = t_tac, injRad = 150,
    bodymass = 85, frameStartEnd = c(1, 33)
  )
  expect_gt(SUVout$par$intSUV, 150)
  expect_lt(SUVout$par$intSUV, 500)
})

test_that("SUV works with durations and frameStartEnd", {
  SUVout <- SUV(tac,
    dur_tac = dur_tac, injRad = 150,
    bodymass = 85, frameStartEnd = c(1, 33)
  )
  expect_gt(SUVout$par$intSUV, 150)
  expect_lt(SUVout$par$intSUV, 500)
})

test_that("kBq to nCi works", {
  kBq_nCi <- round(unit_convert(1, "kBq", "nCi"))
  expect_equal(kBq_nCi, 27)
})

test_that("nCi to Bq works", {
  nCi_Bq <- unit_convert(1, "nCi", "Bq")
  expect_equal(nCi_Bq, 37)
})

s1 <- pbr28$tacs[[1]]

test_that("weights_create default works", {
  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=FALSE)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)

  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=TRUE)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)
})

test_that("weights_create other options work", {
  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=FALSE,
    method = 1)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)

  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=FALSE,
    method = 2)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)

  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=FALSE,
    method = 3)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)

  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=FALSE,
    method = 4)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)

  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=FALSE,
    method = 5)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)

  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=FALSE,
    method = 6)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)

  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=FALSE,
    method = 7)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)

  w <- weights_create(
    s1$StartTime/60,
    (s1$StartTime + s1$Duration)/60,
    radioisotope = "C11",
    tac = s1$WB, minweight_risetopeak=FALSE,
    method = 8)

  expect_true(all(is.numeric(w)))
  expect_true(sum(is.na(w)) == 0)
})
