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


# --- Review fixes: argument validation and per-region handling ---

test_that("vector multstart_iter is accepted by the nested models", {
  out <- nested_1tcm_delay(
    long_data$t_tac, long_data$tac, long_data$region, input,
    timeStartEnd = c(0, 5), multstart_iter = c(3, 3)
  )

  expect_equal(nrow(out$par), length(selected_regions))
  expect_true(all(out$par$Vt > 0))
})

test_that("multstart_iter of the wrong length is rejected", {
  expect_error(
    nested_1tcm_delay(
      long_data$t_tac, long_data$tac, long_data$region, input,
      timeStartEnd = c(0, 5), multstart_iter = c(3, 3, 3)
    ),
    "one value per fitted parameter"
  )
})

test_that("named roiweights missing a region is rejected", {
  expect_error(
    nested_1tcm_delay(
      long_data$t_tac, long_data$tac, long_data$region, input,
      timeStartEnd = c(0, 5), roiweights = c(FC = 1, TC = 2)
    ),
    "missing a value for the region"
  )
})

test_that("nesting a single region warns", {
  one_region <- long_data[long_data$region == "FC", ]

  expect_warning(
    nested_1tcm_delay(
      one_region$t_tac, one_region$tac, one_region$region, input,
      timeStartEnd = c(0, 5)
    ),
    "only one region"
  )
})

test_that("frame weights are applied per region rather than shared", {
  # Downweight the early frames of THA only: THA should move, and the
  # returned weights should be the full stacked vector.
  w_flat <- rep(1, length(t_tac_single))
  w_tha <- w_flat
  w_tha[1:5] <- 0

  out_even <- nested_1tcm_delay(
    long_data$t_tac, long_data$tac, long_data$region, input,
    timeStartEnd = c(0, 5), weights = rep(w_flat, 3)
  )
  out_uneven <- nested_1tcm_delay(
    long_data$t_tac, long_data$tac, long_data$region, input,
    timeStartEnd = c(0, 5), weights = c(w_flat, w_flat, w_tha)
  )

  tha <- which(out_even$par$region == "THA")
  expect_false(isTRUE(all.equal(out_even$par$K1[tha],
                                out_uneven$par$K1[tha])))

  expect_equal(length(out_uneven$weights), nrow(out_uneven$tacs))
  expect_equal(out_uneven$tacs$Weights, out_uneven$weights)
  expect_true(any(out_uneven$weights == 0 &
                    out_uneven$tacs$Region == "THA"))
})

# --- Standard errors ---

test_that("nested_1tcm_delay recovers the delay SE of the unnested model", {
  # With a single region the profile objective is the 1TCM delay likelihood,
  # so the Hessian-based SE should agree closely with the one nls gives.
  one_region <- long_data[long_data$region == "FC", ]

  nested <- suppressWarnings(nested_1tcm_delay(
    one_region$t_tac, one_region$tac, one_region$region, input,
    timeStartEnd = c(0, 5)
  ))
  unnested <- onetcm(t_tac_single, tac_wide$FC, input, vB = 0.05,
                     timeStartEnd = c(0, 5))

  expect_equal(nested$par$inpshift[1], unnested$par$inpshift,
               tolerance = 1e-3)
  expect_equal(nested$par.se$inpshift.se[1], unnested$par.se$inpshift,
               tolerance = 0.05)
})

test_that("nested models report an SE for the shared parameter", {
  out <- nested_1tcm_delay(
    long_data$t_tac, long_data$tac, long_data$region, input,
    timeStartEnd = c(0, 5)
  )

  expect_true("inpshift.se" %in% names(out$par.se))
  expect_true(all(is.finite(out$par.se$inpshift.se)))
  # Shared parameter, so one value for all regions
  expect_equal(length(unique(out$par.se$inpshift.se)), 1)
})

test_that("nested_2tcm gives derived SEs in every shared mode", {
  derived <- c("Vt.se", "BPnd.se", "k2.se", "k3.se")

  for (sh in c("Vnd", "k4", "Vnd_k4")) {
    out <- nested_2tcm(
      long_data$t_tac, long_data$tac, long_data$region, input,
      shared = sh, inpshift = 0.08, vB = 0.05,
      weights = long_data$weights
    )

    expect_true(all(derived %in% names(out$par.se)), info = sh)
    expect_true(all(is.finite(as.matrix(out$par.se[, derived]))), info = sh)

    # SEs for the shared parameters themselves
    for (pname in out$outer_pars) {
      expect_true(paste0(pname, ".se") %in% names(out$par.se), info = sh)
    }
  }
})

# --- Plot pagination ---

test_that("faceted plots past three regions print as a paginated set", {
  more_regions <- c("FC", "TC", "THA", "STR", "CBL")
  wide_data <- do.call(rbind, lapply(more_regions, function(r) {
    data.frame(
      t_tac = t_tac_single,
      tac = tac_wide[[r]],
      region = r,
      stringsAsFactors = FALSE
    )
  }))

  out <- nested_2tcm(
    wide_data$t_tac, wide_data$tac, wide_data$region, input,
    shared = "Vnd", inpshift = 0.08, vB = 0.05
  )

  p <- plot(out)
  expect_s3_class(p, "nested_kinfit_plots")
  expect_equal(length(p), 2)

  # Each page is a ggplot, and the set prints without error
  expect_true(all(vapply(p, function(x) inherits(x, "ggplot"), logical(1))))
  pdf(NULL)
  expect_silent(print(p))
  dev.off()
})
