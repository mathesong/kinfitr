context("test-otherfuncs")

data("pbr28")

meas <- 2

tac <- pbr28$tacs[[meas]]$STR
t_tac <- pbr28$tacs[[meas]]$Times / 60
dur <- pbr28$tacs[[meas]]$Duration / 60
weights <- pbr28$tacs[[meas]]$Weights

# suv

test_that("suv works with trapz", {
  suvout <- suv(tac, t_tac = t_tac, injRad = 150, bodymass = 85)
  expect_gt(suvout$par$SUV_AUC, 200)
  expect_lt(suvout$par$SUV_AUC, 500)
  expect_equal(suvout$par$SUV, suvout$par$SUV_AUC / suvout$par$window_duration)
  expect_equal(suvout$par$SUV_denominator, 150 / 85)
  expect_equal(suvout$par$n_frames, length(tac))
})

test_that("suv works with durations", {
  suvout <- suv(tac, dur = dur, injRad = 150, bodymass = 85)
  expect_gt(suvout$par$SUV_AUC, 150)
  expect_lt(suvout$par$SUV_AUC, 500)
  expect_equal(suvout$par$window_duration, sum(dur))
  expect_equal(suvout$par$SUV, suvout$par$SUV_AUC / sum(dur))
})

test_that("suv works with trapz and frameStartEnd", {
  suvout <- suv(tac,
    t_tac = t_tac, injRad = 150,
    bodymass = 85, frameStartEnd = c(1, 33)
  )
  expect_gt(suvout$par$SUV_AUC, 150)
  expect_lt(suvout$par$SUV_AUC, 500)
  expect_equal(suvout$par$n_frames, 33)
})

test_that("suv works with durations and frameStartEnd", {
  suvout <- suv(tac,
    dur = dur, injRad = 150,
    bodymass = 85, frameStartEnd = c(1, 33)
  )
  expect_gt(suvout$par$SUV_AUC, 150)
  expect_lt(suvout$par$SUV_AUC, 500)
  expect_equal(suvout$par$window_duration, sum(dur[1:33]))
  expect_equal(suvout$par$n_frames, 33)
})

test_that("suv without a dose is in units of radioactivity concentration", {
  suvout <- suv(tac, dur = dur)
  expect_equal(suvout$par$SUV_denominator, 1)
  expect_equal(suvout$par$SUV_AUC, sum(tac * dur))
})

test_that("plot_suvfit works", {
  suvout <- suv(tac, t_tac = t_tac, dur = dur, injRad = 150, bodymass = 85)

  suvplot <- plot_suvfit(suvout, roiname = "FC")
  expect_s3_class(suvplot, "ggplot")

  # The dashed line is the mean SUV which was reported
  hline <- suvplot$layers[[4]]
  expect_s3_class(hline$geom, "GeomHline")
  expect_equal(hline$data$yintercept, suvout$par$SUV)

  # Durations were given, so each included frame is drawn as a rectangle whose
  # area is its contribution, and the areas sum to the reported integral
  expect_s3_class(suvplot$layers[[1]]$geom, "GeomRect")
  rects <- suvplot$layers[[1]]$data
  expect_equal(sum((rects$xmax - rects$xmin) * rects$SUV), suvout$par$SUV_AUC)

  # With a dose, the values shown are SUV
  expect_equal(suvplot$labels$y, "SUV")

  # plot() dispatches to it through the kinfit method
  expect_s3_class(plot(suvout), "ggplot")
  # The model name is what plot_kinfit() dispatches on, to plot_suvfit()
  expect_equal(suvout$model, "suv")
  expect_s3_class(suvout, "suv")
  expect_s3_class(suvout, "kinfit")
})

test_that("a trapezoidal window is integrated over its own frames alone", {
  # tidyinput_art() subsets to the window and then prepends a frame at time zero.
  # Integrating the result would add the whole area between injection and the
  # start of the window, and measure the duration from zero.
  windowed <- suv(tac, t_tac = t_tac, injRad = 150, bodymass = 85,
                  frameStartEnd = c(33, 38))

  expect_equal(
    windowed$par$SUV_AUC,
    pracma::trapz(t_tac[33:38], tac[33:38]) / (150 / 85)
  )
  expect_equal(windowed$par$window_duration, t_tac[38] - t_tac[33])
  expect_equal(windowed$par$SUV, windowed$par$SUV_AUC / windowed$par$window_duration)
  expect_equal(windowed$par$n_frames, 6)

  # The window is late in the acquisition, so its duration cannot reach back to
  # the origin
  expect_lt(windowed$par$window_duration, max(t_tac))
})

test_that("an unwindowed trapezoidal TAC is still integrated from injection", {
  # The frame at time zero belongs here, and this behaviour is unchanged
  whole <- suv(tac, t_tac = t_tac, injRad = 150, bodymass = 85)

  expect_equal(
    whole$par$SUV_AUC,
    pracma::trapz(c(0, t_tac), c(0, tac)) / (150 / 85)
  )
  expect_equal(whole$par$window_duration, max(t_tac))
  expect_equal(whole$par$n_frames, length(tac))
})

test_that("a trapezoidal window starting at the first frame integrates from zero", {
  windowed <- suv(tac, t_tac = t_tac, injRad = 150, bodymass = 85,
                  frameStartEnd = c(1, 10))

  expect_equal(
    windowed$par$SUV_AUC,
    pracma::trapz(c(0, t_tac[1:10]), c(0, tac[1:10])) / (150 / 85)
  )
  expect_equal(windowed$par$window_duration, t_tac[10])
})

test_that("suv returns the whole TAC marked with the frames it integrated", {
  all_frames <- suv(tac, t_tac = t_tac, dur = dur, injRad = 150, bodymass = 85)
  windowed <- suv(tac,
    t_tac = t_tac, dur = dur, injRad = 150,
    bodymass = 85, frameStartEnd = c(10, 20)
  )

  # Both carry every measured frame, not just the window
  expect_equal(nrow(all_frames$tacs), length(tac))
  expect_equal(nrow(windowed$tacs), length(tac))
  expect_equal(windowed$tacs$Target, tac)

  expect_true(all(all_frames$tacs$Included))
  expect_equal(windowed$tacs$Included, seq_along(tac) %in% 10:20)
  expect_equal(windowed$par$n_frames, 11)

  # And the integral is over the marked frames alone
  inc <- windowed$tacs$Included
  expect_equal(
    sum(windowed$tacs$SUV[inc] * windowed$tacs$Duration[inc]),
    windowed$par$SUV_AUC
  )
  expect_equal(windowed$par$window_duration, sum(dur[10:20]))
})

test_that("plot_suvfit shades the window within the whole curve", {
  suvout <- suv(tac,
    t_tac = t_tac, dur = dur, injRad = 150,
    bodymass = 85, frameStartEnd = c(10, 20)
  )

  suvplot <- plot_suvfit(suvout)

  # The curve is the whole TAC; only the window is shaded
  expect_equal(nrow(suvplot$data), length(tac))
  rects <- suvplot$layers[[1]]$data
  expect_equal(nrow(rects), 11)
  expect_true(all(rects$Included))

  # Each rectangle spans its own frame, so they run from the start of the first
  # included frame to the end of the last, and abut one another
  expect_equal(rects$xmin, t_tac[10:20] - dur[10:20] / 2)
  expect_equal(rects$xmax, t_tac[10:20] + dur[10:20] / 2)
  expect_equal(max(rects$xmax) - min(rects$xmin), suvout$par$window_duration)
})

test_that("plot_suvfit shades under the curve when integration is trapezoidal", {
  # No durations, so suv() integrates with trapz and the area under the curve is
  # what was taken, not a set of frame rectangles
  suvout <- suv(tac, t_tac = t_tac, injRad = 150, bodymass = 85,
                frameStartEnd = c(10, 20))

  suvplot <- plot_suvfit(suvout)
  expect_s3_class(suvplot$layers[[1]]$geom, "GeomRibbon")
  expect_equal(nrow(suvplot$layers[[1]]$data), 11)
})

test_that("plot_suvfit copes with durations alone, and names its units", {
  suvout <- suv(tac, dur = dur)

  suvplot <- plot_suvfit(suvout)
  expect_s3_class(suvplot, "ggplot")

  # No times were given, so the frames are placed end to end
  expect_equal(suvplot$labels$x, "Cumulative frame duration (min)")
  # And without a dose these are concentrations, not SUV
  expect_equal(suvplot$labels$y, "Radioactivity")
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
