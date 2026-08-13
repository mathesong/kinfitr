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
