# Parsing a study must not be all-or-nothing: one malformed measurement should
# not take down every other measurement in the study.

make_bids_study <- function(complete = "01", no_sidecar = character(0),
                            empty_sidecar = character(0)) {

  root <- tempfile("kinfitr-bids-")

  write_pet <- function(sub, json) {
    dir.create(file.path(root, paste0("sub-", sub), "pet"), recursive = TRUE)
    stem <- file.path(root, paste0("sub-", sub), "pet", paste0("sub-", sub, "_pet"))
    file.create(paste0(stem, ".nii.gz"))
    if (!is.null(json)) {
      writeLines(json, paste0(stem, ".json"))
    }
  }

  for (s in complete) {
    write_pet(s, '{"FrameTimesStart":[0,60,120],"FrameDuration":[60,60,60]}')
  }
  for (s in no_sidecar) {
    write_pet(s, NULL)
  }
  for (s in empty_sidecar) {
    write_pet(s, '{"TracerRadionuclide":"[18F]"}')
  }

  root
}

test_that("a measurement missing its _pet.json does not abort the whole parse", {

  root <- make_bids_study(complete = c("01", "03"), no_sidecar = "02")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  # Previously this failed outright with "object 'dur' not found": the tibble()
  # call referenced `dur`, which was never created because the sidecar filter
  # returned no rows. Being mapped over every measurement, that single failure
  # took the entire study down with it.
  expect_warning(measurements <- bids_parse_study(root),
                 "No _pet.json sidecar found")

  # The two good measurements survive; only the malformed one is dropped
  expect_equal(nrow(measurements), 2)
  expect_equal(sort(measurements$sub), c("01", "03"))
})

test_that("the missing-sidecar warning names the measurement", {

  root <- make_bids_study(complete = "01", no_sidecar = "02")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  expect_warning(bids_parse_study(root), "sub-02")
})

test_that("a _pet.json without frame timing is reported, not crashed on", {

  root <- make_bids_study(complete = "01", empty_sidecar = "02")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  # A sidecar that exists but carries no FrameTimesStart/FrameDuration reaches
  # the same tibble() and fails the same way.
  expect_warning(measurements <- bids_parse_study(root),
                 "missing FrameTimesStart and FrameDuration")

  expect_equal(nrow(measurements), 1)
  expect_equal(measurements$sub, "01")
})

test_that("a fully specified study parses its frame times", {

  root <- make_bids_study(complete = "01")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  measurements <- bids_parse_study(root)

  expect_equal(nrow(measurements), 1)

  # Times are converted from seconds to minutes, and time is the frame midpoint
  tactimes <- measurements$tactimes[[1]]
  expect_equal(tactimes$start, c(0, 1, 2))
  expect_equal(tactimes$dur, c(1, 1, 1))
  expect_equal(tactimes$time, c(0.5, 1.5, 2.5))
})

test_that("entity labels differing only in case are warned about", {

  root <- tempfile("kinfitr-case-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  # Same tracer, two spellings. BIDS labels are case-sensitive, so these become
  # two distinct measurements -- valid parsing of an invalid dataset.
  for (spec in list(c("01", "PF974"), c("02", "pf974"))) {
    dir.create(file.path(root, paste0("sub-", spec[1]), "pet"), recursive = TRUE)
    stem <- file.path(root, paste0("sub-", spec[1]), "pet",
                      paste0("sub-", spec[1], "_trc-", spec[2], "_pet"))
    file.create(paste0(stem, ".nii.gz"))
    writeLines('{"FrameTimesStart":[0,60],"FrameDuration":[60,60]}',
               paste0(stem, ".json"))
  }

  expect_warning(bids_parse_files(root), "differing only in case")
  expect_warning(bids_parse_files(root), "trc-PF974")
  expect_warning(bids_parse_files(root), "bids-validator")

  # The replacement parser checks too, comparing entity by entity
  expect_warning(bids_parse_filenames(root), "trc-PF974")

  # It warns; it must never refuse the dataset
  files <- suppressWarnings(bids_parse_files(root))
  expect_true(nrow(files) > 0)
})

test_that("a consistently cased study raises no case warning", {

  root <- make_bids_study(complete = c("01", "02"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  # Not expect_silent: bids_parse_files() is deprecated and says so. What must
  # be absent is the case warning specifically.
  warnings <- character(0)
  withCallingHandlers(
    bids_parse_filenames(root),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  expect_false(any(grepl("differing only in case", warnings)))
})
