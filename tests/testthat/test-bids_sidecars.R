# bids_resolve_sidecars() merges the sidecars applying to one data file.

sidecar_study <- function(files) {
  root <- tempfile("kinfitr-sidecars-")
  for (f in names(files)) {
    dir.create(file.path(root, dirname(f)), recursive = TRUE, showWarnings = FALSE)
    writeLines(files[[f]], file.path(root, f))
  }
  root
}

test_that("a nearer sidecar overrides a further one, field by field", {

  root <- sidecar_study(list(
    "task-rest_pet.json" =
      '{"TracerRadionuclide":"[18F]","InjectedRadioactivity":100,"Manufacturer":"GE"}',
    "sub-01/ses-test/pet/sub-01_ses-test_pet.json" =
      '{"InjectedRadioactivity":183416.77}'))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  resolved <- suppressWarnings(bids_resolve_sidecars(
    "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz",
    c("task-rest_pet.json", "sub-01/ses-test/pet/sub-01_ses-test_pet.json"),
    sidecar_root = root))

  # The nearer file wins for the field it sets...
  expect_equal(resolved$values$InjectedRadioactivity, 183416.77)
  # ...and the root's other fields are inherited rather than lost
  expect_equal(resolved$values$TracerRadionuclide, "[18F]")
  expect_equal(resolved$values$Manufacturer, "GE")
})

test_that("resolution reports where each field came from", {

  root <- sidecar_study(list(
    "task-rest_pet.json" = '{"TracerRadionuclide":"[18F]","InjectedRadioactivity":100}',
    "sub-01/ses-test/pet/sub-01_ses-test_pet.json" = '{"InjectedRadioactivity":183416.77}'))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  resolved <- suppressWarnings(bids_resolve_sidecars(
    "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz",
    c("task-rest_pet.json", "sub-01/ses-test/pet/sub-01_ses-test_pet.json"),
    sidecar_root = root))

  expect_equal(resolved$provenance[["TracerRadionuclide"]], "task-rest_pet.json")
  expect_equal(resolved$provenance[["InjectedRadioactivity"]],
               "sub-01/ses-test/pet/sub-01_ses-test_pet.json")
})

test_that("two sidecars in one directory are refused", {

  root <- sidecar_study(list(
    "task-rest_pet.json" = '{"InjectedRadioactivity":100}',
    "task-memory_pet.json" = '{"InjectedRadioactivity":200}'))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  # Both are lenient-applicable to a task-less image, and nothing distinguishes
  # them. Choosing by enumeration order would make the injected dose depend on
  # how the filesystem lists files.
  expect_error(
    suppressWarnings(bids_resolve_sidecars(
      "sub-01/pet/sub-01_pet.nii.gz",
      c("task-rest_pet.json", "task-memory_pet.json"), sidecar_root = root)),
    "same directory")
})

test_that("task and rec are lenient, trc and run are not", {

  data_file <- "sub-01/pet/sub-01_pet.nii.gz"

  # A sidecar naming task or rec still applies to an image naming neither
  expect_true(is.na(bids_sidecar_applies("task-rest_pet.json", data_file,
                                         c("task", "rec"))))
  expect_true(is.na(bids_sidecar_applies("rec-acdyn_pet.json", data_file,
                                         c("task", "rec"))))

  # trc does not: attaching the wrong tracer's sidecar would supply the wrong
  # half-life and injected dose without any sign of it
  expect_match(bids_sidecar_applies("trc-pf974_pet.json", data_file,
                                    c("task", "rec")),
               "names trc-pf974")
  expect_match(bids_sidecar_applies("run-01_pet.json", data_file,
                                    c("task", "rec")),
               "names run-01")
})

test_that("a sidecar naming a different value never applies", {

  data_file <- "sub-01/pet/sub-01_task-memory_pet.nii.gz"

  # Leniency covers an entity the data file omits, not one it contradicts
  expect_match(bids_sidecar_applies("task-rest_pet.json", data_file,
                                    c("task", "rec")),
               "but the data file names task-memory")
})

test_that("suffix and directory ancestry are required", {

  expect_match(bids_sidecar_applies("sub-01/pet/sub-01_blood.json",
                                    "sub-01/pet/sub-01_pet.nii.gz", "task"),
               "suffix differs")

  # A sibling directory is not an ancestor
  expect_match(bids_sidecar_applies("sub-02/pet/sub-02_pet.json",
                                    "sub-01/pet/sub-01_pet.nii.gz", "task"),
               "not in a parent directory")
})

test_that("skipped sidecars are named with their reason", {

  root <- sidecar_study(list("trc-pf974_pet.json" = '{"InjectedRadioactivity":100}'))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  # A sidecar that quietly fails to apply is how metadata goes missing without
  # anything looking wrong
  expect_warning(
    bids_resolve_sidecars("sub-01/pet/sub-01_pet.nii.gz",
                          "trc-pf974_pet.json", sidecar_root = root),
    "names trc-pf974")
})

test_that("resolution is independent of the order sidecars are supplied", {

  root <- sidecar_study(list(
    "task-rest_pet.json" = '{"InjectedRadioactivity":100}',
    "sub-01/ses-test/pet/sub-01_ses-test_pet.json" = '{"InjectedRadioactivity":183416.77}'))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  data_file <- "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz"
  paths <- c("task-rest_pet.json", "sub-01/ses-test/pet/sub-01_ses-test_pet.json")

  forwards <- suppressWarnings(
    bids_resolve_sidecars(data_file, paths, sidecar_root = root))
  backwards <- suppressWarnings(
    bids_resolve_sidecars(data_file, rev(paths), sidecar_root = root))

  expect_equal(forwards$values, backwards$values)
  expect_equal(forwards$provenance, backwards$provenance)
})

test_that("no sidecars yields an empty result rather than an error", {

  resolved <- bids_resolve_sidecars("sub-01/pet/sub-01_pet.nii.gz", character(0))

  expect_equal(length(resolved$values), 0)
  expect_equal(length(resolved$provenance), 0)
})
