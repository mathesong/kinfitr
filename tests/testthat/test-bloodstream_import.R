# The bloodstream importers read a derivatives folder, which may hold several
# analyses and files belonging to the analysis rather than to any acquisition.

bloodstream_folder <- function(files) {
  root <- tempfile("kinfitr-bloodstream-")
  for (f in names(files)) {
    dir.create(file.path(root, dirname(f)), recursive = TRUE, showWarnings = FALSE)
    writeLines(files[[f]], file.path(root, f))
  }
  root
}

input_json <- paste0(
  '{"time":{"Units":"s"},"whole_blood_radioactivity":{"Units":"kBq"},',
  '"plasma_radioactivity":{"Units":"kBq"},',
  '"metabolite_parent_fraction":{"Units":"unitless"}}')

input_tsv <- paste(
  "time\twhole_blood_radioactivity\tplasma_radioactivity\tmetabolite_parent_fraction",
  "0\t0\t0\t1", "60\t10\t12\t0.9", sep = "\n")

test_that("a config belonging to the analysis is not read as a measurement", {

  # A bloodstream folder holds its own bloodstream_config.json beside the
  # per-measurement AIF fits. Only the latter carry fitted parameters, and
  # reading the former as though it were one fails.
  root <- bloodstream_folder(list(
    "Primary_Analysis/bloodstream_config.json" = '{"Subsets":[],"Model":"x"}',
    # Shaped as bloodstream writes it: Method and Units are length-one arrays,
    # Parameters an array of one object.
    "Primary_Analysis/sub-01/ses-01/pet/sub-01_ses-01_config.json" =
      paste0('{"AIF":{"Method":["Fit Individually: FengConv"],',
             '"Fit":{"Parameters":[{"A":54526.4,"alpha":0.2029,"B":6712.6,',
             '"beta":0.0103,"C":3367.1,"gamma":0.0008,"t0":26.0,"ti":7.2}],',
             '"Units":{"time":["s"],"AIF":["Bq"]}}}}')))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  pars <- suppressWarnings(bloodstream_import_aifpars(root))

  # One row for the one measurement; the analysis-level config contributes none
  expect_equal(nrow(pars), 1)
  expect_equal(pars$sub, "01")
  expect_equal(pars$ses, "01")
  expect_true("AIFpars" %in% colnames(pars))

  # Times are converted from seconds to minutes on the way out
  expect_equal(pars$AIFpars[[1]]$t0, 26.0 / 60, tolerance = 1e-6)
})

test_that("an input function without its sidecar is named, not crashed on", {

  root <- bloodstream_folder(list(
    "Primary_Analysis/sub-01/ses-01/pet/sub-01_ses-01_inputfunction.tsv" = input_tsv,
    "Primary_Analysis/sub-01/ses-01/pet/sub-01_ses-01_inputfunction.json" = input_json,
    # sub-02 has no json: the units and availability flags live there
    "Primary_Analysis/sub-02/ses-01/pet/sub-02_ses-01_inputfunction.tsv" = input_tsv))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  expect_warning(imported <- bloodstream_import_inputfunctions(root),
                 "complete tsv/json pair")
  expect_warning(bloodstream_import_inputfunctions(root), "sub-02")

  # The usable one still comes through
  expect_equal(nrow(imported), 1)
  expect_equal(imported$sub, "01")
})

test_that("the same acquisition in two analyses is imported once each", {

  root <- bloodstream_folder(list(
    "Primary_Analysis/sub-01/ses-01/pet/sub-01_ses-01_inputfunction.tsv" = input_tsv,
    "Primary_Analysis/sub-01/ses-01/pet/sub-01_ses-01_inputfunction.json" = input_json,
    "tutorial/sub-01/ses-01/pet/sub-01_ses-01_inputfunction.tsv" = input_tsv,
    "tutorial/sub-01/ses-01/pet/sub-01_ses-01_inputfunction.json" = input_json))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  imported <- suppressWarnings(bloodstream_import_inputfunctions(root))

  # Two analyses of one acquisition: two rows, not four from a collided pairing
  expect_equal(nrow(imported), 2)
  expect_true(all(imported$sub == "01"))
})

test_that("importing yields the entity columns callers join on", {

  root <- bloodstream_folder(list(
    "Primary_Analysis/sub-01/ses-01/pet/sub-01_ses-01_inputfunction.tsv" = input_tsv,
    "Primary_Analysis/sub-01/ses-01/pet/sub-01_ses-01_inputfunction.json" = input_json))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  imported <- suppressWarnings(bloodstream_import_inputfunctions(root))

  # petfit natural-joins its TAC data against this, so the entity columns are
  # the contract, and the internal keys must not leak into it
  expect_true(all(c("sub", "ses", "measurement", "input") %in% colnames(imported)))
  expect_false(any(c("source_key", "artifact_key", "analysis_scope_key") %in%
                     colnames(imported)))

  # Entities are what the filenames carry, not defaults
  expect_equal(imported$sub, "01")
  expect_equal(imported$ses, "01")
  expect_false("acq" %in% colnames(imported))
})

test_that("an importer with only incomplete input functions returns empty", {

  # Every input function lacking its sidecar: the incomplete rows are warned
  # about and dropped, and what remains must come back as a zero-row table
  # rather than a subscript error.
  root <- bloodstream_folder(list(
    "Primary_Analysis/sub-01/ses-01/pet/sub-01_ses-01_inputfunction.tsv" =
      input_tsv))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  expect_warning(imported <- bloodstream_import_inputfunctions(root),
                 "complete tsv/json pair")

  expect_equal(nrow(imported), 0)
  expect_true("input" %in% colnames(imported))
})
