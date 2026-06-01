context("test-nested.R")

data("pbr28")

meas <- 2

# Prepare long-format data from pbr28 (selecting a few regions)
tac_wide <- pbr28$tacs[[meas]]
t_tac_single <- tac_wide$Times / 60
weights_single <- tac_wide$Weights
input <- pbr28$input[[meas]]

# Select 3 regions for testing
selected_regions <- c("FC", "TC", "THA")

long_data <- do.call(rbind, lapply(selected_regions, function(r) {
  data.frame(
    t_tac = t_tac_single,
    tac = tac_wide[[r]],
    region = r,
    weights = weights_single,
    stringsAsFactors = FALSE
  )
}))

set.seed(42)

# --- nested_1tcm_delay tests ---

test_that("nested_1tcm_delay works", {
  out <- nested_1tcm_delay(
    long_data$t_tac, long_data$tac, long_data$region, input,
    vB = 0.05,
    timeStartEnd = c(0, 5)
  )

  # Check structure
  expect_true("nested_1tcm_delay" %in% class(out))
  expect_true("kinfit" %in% class(out))
  expect_true(all(c("par", "par.se", "outer_pars", "fit", "tacs") %in% names(out)))

  # Check outer_pars is character vector
  expect_true(is.character(out$outer_pars))
  expect_equal(out$outer_pars, "inpshift")

  # Check par column order: region first
  expect_equal(names(out$par)[1], "region")

  # Check per-region results
  expect_equal(nrow(out$par), length(selected_regions))
  expect_true(all(out$par$Vt > 0))
  expect_true(all(out$par$Vt < 20))

  # Check vB stored
  expect_equal(out$vB, 0.05)

  # Check plotting
  expect_true(any(class(plot_kinfit(out)) == "ggplot"))

  # Check predict
  preds <- predict(out)
  expect_true(is.numeric(preds))
  expect_equal(length(preds), nrow(out$tacs))
})

# --- nested_2tcm_delay tests ---

test_that("nested_2tcm_delay works", {
  out <- nested_2tcm_delay(
    long_data$t_tac, long_data$tac, long_data$region, input,
    vB = 0.05,
    timeStartEnd = c(0, 5)
  )

  expect_true("nested_2tcm_delay" %in% class(out))
  expect_true(is.character(out$outer_pars))
  expect_equal(out$outer_pars, "inpshift")
  expect_equal(names(out$par)[1], "region")
  expect_equal(nrow(out$par), length(selected_regions))
  expect_true(all(out$par$Vt > 0))
  expect_true(any(class(plot_kinfit(out)) == "ggplot"))

  # Check predict
  preds <- predict(out)
  expect_true(is.numeric(preds))
  expect_equal(length(preds), nrow(out$tacs))
})

# --- nested_2tcm tests ---

test_that("nested_2tcm with shared Vnd works", {
  out <- nested_2tcm(
    long_data$t_tac, long_data$tac, long_data$region, input,
    shared = "Vnd",
    inpshift = 0.08, vB = 0.05,
    weights = long_data$weights
  )

  expect_true("nested_2tcm" %in% class(out))
  expect_equal(nrow(out$par), length(selected_regions))
  expect_true(is.character(out$outer_pars))
  expect_equal(out$outer_pars, "Vnd")
  expect_equal(out$shared, "Vnd")
  expect_equal(names(out$par)[1], "region")

  # k4 should be per-region (different values possible)
  expect_true("k4" %in% names(out$par))

  # Check plotting
  expect_true(any(class(plot_kinfit(out)) == "ggplot"))

  # Check predict
  preds <- predict(out)
  expect_true(is.numeric(preds))
  expect_equal(length(preds), nrow(out$tacs))
})

test_that("nested_2tcm with shared k4 works", {
  out <- nested_2tcm(
    long_data$t_tac, long_data$tac, long_data$region, input,
    shared = "k4",
    inpshift = 0.14, vB = 0.05,
    weights = long_data$weights
  )

  expect_true(is.character(out$outer_pars))
  expect_equal(out$outer_pars, "k4")
  expect_equal(out$shared, "k4")

  # Vnd should be per-region
  expect_true("Vnd" %in% names(out$par))
})

test_that("nested_2tcm with shared Vnd and k4 works", {
  out <- nested_2tcm(
    long_data$t_tac, long_data$tac, long_data$region, input,
    shared = "Vnd_k4",
    inpshift = 0.14, vB = 0.05,
    weights = long_data$weights
  )

  expect_true(is.character(out$outer_pars))
  expect_true(all(c("Vnd", "k4") %in% out$outer_pars))
  expect_equal(out$shared, "Vnd_k4")
})

# --- nested_srtm tests ---

test_that("nested_srtm works", {
  data("simref")

  # Build long-format data from simref
  sim_t <- simref$tacs[[2]]$Times
  sim_ref <- simref$tacs[[2]]$Reference
  sim_weights <- simref$tacs[[2]]$Weights

  sim_regions <- c("ROI1", "ROI2", "ROI3")
  sim_long <- do.call(rbind, lapply(sim_regions, function(r) {
    data.frame(
      t_tac = sim_t,
      roitac = simref$tacs[[2]][[r]],
      reftac = sim_ref,
      region = r,
      weights = sim_weights,
      stringsAsFactors = FALSE
    )
  }))

  out <- nested_srtm(
    sim_long$t_tac, sim_long$roitac, sim_long$reftac, sim_long$region,
    weights = sim_long$weights
  )

  expect_true("nested_srtm" %in% class(out))
  expect_equal(nrow(out$par), length(sim_regions))
  expect_true(is.character(out$outer_pars))
  expect_equal(out$outer_pars, "k2prime")
  expect_equal(names(out$par)[1], "region")

  # k2prime from fit object
  expect_true(out$fit$par[["k2prime"]] > 0)
  expect_true(out$fit$par[["k2prime"]] < 1)

  expect_true(any(class(plot_kinfit(out)) == "ggplot"))

  # Check predict
  preds <- predict(out)
  expect_true(is.numeric(preds))
  expect_equal(length(preds), nrow(out$tacs))
})


# --- roiweights with prepended zero point (regression) ---

test_that("per-observation roiweights work when a zero frame is prepended", {
  # Drop the t=0 frame so tidyinput_long must prepend one per region. A
  # per-observation roiweights vector sized to the input TACs must still be
  # accepted (previously it ended up one-per-region too short after the zero
  # point was added).
  nz <- long_data$t_tac > 0
  ld <- long_data[nz, ]

  roiw <- as.numeric(factor(ld$region, levels = selected_regions))  # one per obs

  expect_false(any(ld$t_tac == 0))

  out <- nested_2tcm(
    ld$t_tac, ld$tac, ld$region, input,
    vB = 0.05, roiweights = roiw
  )

  expect_true("nested_2tcm" %in% class(out))
  expect_equal(nrow(out$par), length(selected_regions))
  # roiweights collapse to one (normalised) value per region, in region order
  expect_equal(length(out$roiweights), length(selected_regions))
  expect_equal(names(out$roiweights), selected_regions)
  expect_equal(max(out$roiweights), 1)
})
