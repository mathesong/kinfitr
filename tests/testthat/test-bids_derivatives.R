# bids_parse_derivatives() does filename grouping over a derivatives folder.
# A derivatives folder holds no PET images, so there is nothing to enumerate and
# the raw parser's entity grid does not apply to it.

make_derivatives <- function(files) {
  root <- tempfile("kinfitr-derivs-")
  for (f in files) {
    dir.create(file.path(root, dirname(f)), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(root, f))
  }
  root
}

# The shape petfit actually writes: session in the directory, not the filename.
petfit_shaped <- c(
  "desc-petfitoptions_config.json",
  "model-2TCM_desc-model1_kinpar.tsv",
  "reports/model1_report.html",
  "sub-01/ses-test/pet/sub-01_desc-combinedregions_tacs.tsv",
  "sub-01/ses-test/pet/sub-01_inputfunction.tsv",
  "sub-01/ses-test/pet/sub-01_inputfunction.json",
  "sub-02/pet/sub-02_desc-combinedregions_tacs.tsv"
)

test_that("bids_parse_derivatives lets the directory complete the filename", {

  root <- make_derivatives(petfit_shaped)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)

  # sub-01 stores its session in the path only. A filename-only key would fail
  # to join for exactly the subjects that have sessions.
  expect_equal(unique(d$source_key[d$sub == "01" & !is.na(d$sub)]),
               "sub-01_ses-test")

  # sub-02 has no session anywhere, and must not acquire one. This is the shape
  # that made a subject vanish when run alongside others: the raw parser gave it
  # its neighbours' session, and it then matched nothing.
  expect_equal(unique(d$source_key[d$sub == "02" & !is.na(d$sub)]),
               "sub-02")
})

test_that("bids_parse_derivatives ignores reports and other non-data files", {

  root <- make_derivatives(petfit_shaped)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)

  # Reports carry no entities, so every one of them would key identically
  expect_false(any(grepl("\\.html$", d$path)))
  expect_equal(nrow(d), 6)
})

test_that("bids_parse_derivatives scopes group-level files instead of joining them", {

  root <- make_derivatives(petfit_shaped)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)
  group <- d[is.na(d$source_key), ]

  # A file with no sub is a group result. It gets no source_key, so it cannot be
  # joined to an acquisition by accident and masquerade as a per-subject result.
  expect_equal(sort(group$path),
               c("desc-petfitoptions_config.json",
                 "model-2TCM_desc-model1_kinpar.tsv"))
  expect_true(all(!is.na(group$analysis_scope_key)))
  expect_equal(unique(group$analysis_scope_key), basename(root))

  # Everything else is attributable
  expect_true(all(is.na(d$analysis_scope_key[!is.na(d$source_key)])))
})

test_that("a tsv and its json sidecar are one artifact", {

  root <- make_derivatives(petfit_shaped)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)
  inputfunction <- d[d$suffix == "inputfunction" & !is.na(d$suffix), ]

  # Two files, one product: a tsv should nearly always have a json beside it
  expect_equal(nrow(inputfunction), 2)
  expect_equal(length(unique(inputfunction$artifact_key)), 1)
  expect_setequal(inputfunction$extension, c("tsv", "json"))
})

test_that("distinct products get distinct artifact keys", {

  root <- make_derivatives(c(
    "sub-01/pet/sub-01_model-1TCM_desc-model1_kinpar.tsv",
    "sub-01/pet/sub-01_model-2TCM_desc-model1_kinpar.tsv",
    "sub-01/pet/sub-01_model-2TCM_desc-model1fitted_tacs.tsv",
    # An entity the parser has no special knowledge of must still separate files
    "sub-01/pet/sub-01_model-2TCM_space-MNI_desc-model1_kinpar.tsv"
  ))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)

  expect_equal(nrow(d), 4)
  expect_equal(length(unique(d$artifact_key)), 4)
  expect_equal(length(unique(d$source_key)), 1)

  # The unfamiliar entity is carried into the key rather than dropped
  expect_true(any(grepl("space-MNI", d$artifact_key)))
})

test_that("an entity a file does not name imposes no constraint on it", {

  root <- make_derivatives(c(
    "sub-01/ses-test/pet/sub-01_trc-A_desc-model1_kinpar.tsv",
    "sub-01/ses-test/pet/sub-01_trc-B_desc-model1_kinpar.tsv",
    # Weights carry no trc, and apply to both of the results above
    "sub-01/ses-test/pet/sub-01_desc-weights_weights.tsv"
  ))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)
  weights <- d[d$suffix == "weights", ]

  # The tracer-specific results key with their tracer...
  expect_setequal(
    d$source_key[d$suffix == "kinpar"],
    c("sub-01_ses-test_trc-A", "sub-01_ses-test_trc-B"))

  # ...while the weights key with only what they actually say. Matching that
  # against acquisitions is a subset test, so it reaches both tracers rather
  # than being attributed to one of them or refused as ambiguous.
  expect_equal(weights$source_key, "sub-01_ses-test")
  expect_true(is.na(weights$trc))
})

test_that("bids_parse_derivatives errors when filename and directory disagree", {

  root <- make_derivatives(
    "sub-01/ses-test/pet/sub-01_ses-retest_desc-combinedregions_tacs.tsv")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  expect_error(bids_parse_derivatives(root), "different values for the same entity")
})

test_that("bids_parse_derivatives handles empty and missing folders", {

  root <- make_derivatives("reports/model1_report.html")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  # Only non-data files present
  expect_equal(nrow(bids_parse_derivatives(root)), 0)

  expect_error(bids_parse_derivatives(file.path(root, "nope")),
               "Derivatives folder not found")
})

test_that("path_absolute points at the file that was found", {

  root <- make_derivatives("sub-01/ses-test/pet/sub-01_inputfunction.json")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  writeLines('{"time": {"Units": "s"}}',
             file.path(root, "sub-01/ses-test/pet/sub-01_inputfunction.json"))

  d <- bids_parse_derivatives(root)

  # path is relative to the folder, path_absolute resolves on its own. These
  # came out equal once, because tibble() evaluates its arguments in its own
  # scope and `path` there meant the column, not the argument.
  expect_true(file.exists(d$path_absolute))
  expect_false(identical(d$path, d$path_absolute))
  expect_equal(d$path, "sub-01/ses-test/pet/sub-01_inputfunction.json")

  # It is a real file, so it can actually be read
  expect_equal(jsonlite::fromJSON(d$path_absolute)$time$Units, "s")
})

test_that("the same acquisition in two analyses gives two artifacts", {

  # A derivatives tree can hold several analyses side by side. Without the path
  # scope in the key, one acquisition's input function in two analyses is a
  # single artifact, and pairing its tsv with its json then multiplies rows.
  root <- make_derivatives(c(
    "Primary_Analysis/sub-01/ses-test/pet/sub-01_ses-test_inputfunction.tsv",
    "Primary_Analysis/sub-01/ses-test/pet/sub-01_ses-test_inputfunction.json",
    "tutorial/sub-01/ses-test/pet/sub-01_ses-test_inputfunction.tsv",
    "tutorial/sub-01/ses-test/pet/sub-01_ses-test_inputfunction.json"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)

  expect_equal(nrow(d), 4)
  expect_equal(length(unique(d$artifact_key)), 2)

  # Same acquisition either way: the scope separates the products, not the source
  expect_equal(length(unique(d$source_key)), 1)
  expect_true(any(grepl("^Primary_Analysis_", d$artifact_key)))
  expect_true(any(grepl("^tutorial_", d$artifact_key)))
})

test_that("group-level files in different analyses stay distinct", {

  root <- make_derivatives(c(
    "Primary_Analysis/bloodstream_config.json",
    "tutorial/bloodstream_config.json"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)

  expect_equal(nrow(d), 2)
  expect_equal(length(unique(d$artifact_key)), 2)
  expect_true(all(is.na(d$source_key)))
})

test_that("files directly under the root keep an unscoped key", {

  # petfit calls this on a single analysis folder, where there is no scope
  # directory. Those keys must not gain a prefix.
  root <- make_derivatives(c(
    "desc-petfitoptions_config.json",
    "sub-01/ses-test/pet/sub-01_desc-combinedregions_tacs.tsv"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)

  expect_equal(d$artifact_key[d$suffix == "tacs"],
               "sub-01_ses-test_desc-combinedregions_tacs")
  expect_equal(d$artifact_key[d$suffix == "config"],
               paste0(basename(root), "_desc-petfitoptions_config"))
})

test_that("analysis_scope_key names each analysis, not the parse root", {

  # Parsing a parent of several analyses: each group-level file's scope is the
  # analysis holding it. One shared scope would make "which analysis?"
  # unanswerable exactly where it matters.
  root <- make_derivatives(c(
    "analysisA/bloodstream_config.json",
    "analysisB/bloodstream_config.json"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- bids_parse_derivatives(root)

  expect_setequal(d$analysis_scope_key, c("analysisA", "analysisB"))

  # Directly under the root there is no analysis directory, and the root's own
  # name remains the scope -- petfit parses single analysis folders this way.
  single <- make_derivatives("bloodstream_config.json")
  on.exit(unlink(single, recursive = TRUE), add = TRUE)
  expect_equal(bids_parse_derivatives(single)$analysis_scope_key,
               basename(single))
})
