# bids_parse_filenames() enumerates the acquisitions that exist on disk.
#
# The parser it replaces built a grid of every entity value against every other
# and kept the combinations matching a file, which invents acquisitions: on the
# fixture below it returns 192 measurements where there are 12.

ragged_fixture <- function() {

  root <- tempfile("kinfitr-ragged-")
  mk <- function(p) {
    dir.create(file.path(root, dirname(p)), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(root, p))
  }

  # Two reconstructions sharing one blood file, which names no rec
  mk("sub-01/ses-test/pet/sub-01_ses-test_rec-A_pet.nii.gz")
  mk("sub-01/ses-test/pet/sub-01_ses-test_rec-B_pet.nii.gz")
  mk("sub-01/ses-test/pet/sub-01_ses-test_recording-manual_blood.tsv")
  # Two runs sharing one blood file, which names no run
  mk("sub-02/ses-test/pet/sub-02_ses-test_run-01_pet.nii.gz")
  mk("sub-02/ses-test/pet/sub-02_ses-test_run-02_pet.nii.gz")
  mk("sub-02/ses-test/pet/sub-02_ses-test_recording-manual_blood.tsv")
  # No session anywhere
  mk("sub-03/anat/sub-03_T1w.nii.gz")
  mk("sub-03/pet/sub-03_pet.nii.gz")
  # Two sessions, blood in only one of them
  mk("sub-04/ses-test/pet/sub-04_ses-test_pet.nii.gz")
  mk("sub-04/ses-test/pet/sub-04_ses-test_recording-manual_blood.tsv")
  mk("sub-04/ses-retest/pet/sub-04_ses-retest_pet.nii.gz")
  # A study-level inheritance sidecar, which is not an acquisition
  mk("task-rest_pet.json")

  root
}

test_that("bids_parse_filenames enumerates acquisitions, not combinations", {

  root <- ragged_fixture()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  m <- suppressWarnings(bids_parse_filenames(root))

  # Seven images, and a study-level sidecar that is not one of them
  expect_equal(nrow(m), 7)

  # Nothing is invented. Absent entities stay absent rather than being filled
  # in from whatever the neighbouring subjects happen to use.
  expect_false("task" %in% colnames(m))
  expect_false("acq" %in% colnames(m))
})

test_that("bids_parse_filenames leaves a sessionless subject sessionless", {

  root <- ragged_fixture()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  m <- suppressWarnings(bids_parse_filenames(root))
  sub03 <- m[m$sub == "03", ]

  # Exactly one acquisition, with no session. Giving it its neighbours' sessions
  # is what made a subject vanish when run alongside others.
  expect_equal(nrow(sub03), 1)
  expect_true(is.na(sub03$ses))
})

test_that("a file reaches every acquisition it does not exclude", {

  root <- ragged_fixture()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  m <- suppressWarnings(bids_parse_filenames(root))
  blood_of <- function(rows) {
    vapply(m$filedata[rows], function(f) sum(f$measurement == "blood"), integer(1))
  }

  # One blood file, no rec entity, two reconstructions: it belongs to both
  expect_equal(unname(blood_of(m$sub == "01")), c(1L, 1L))

  # Likewise across runs
  expect_equal(unname(blood_of(m$sub == "02")), c(1L, 1L))
})

test_that("blood does not leak between sessions", {

  root <- ragged_fixture()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  m <- suppressWarnings(bids_parse_filenames(root))
  sub04 <- m[m$sub == "04", ]
  blood <- vapply(sub04$filedata,
                  function(f) sum(f$measurement == "blood"), integer(1))

  # Blood was collected in ses-test only, and names that session, so it must not
  # reach ses-retest
  expect_equal(sum(blood), 1)
  expect_equal(blood[sub04$ses == "test"], 1L, ignore_attr = TRUE)
  expect_equal(blood[sub04$ses == "retest"], 0L, ignore_attr = TRUE)
})

test_that("only _pet files inside pet/ count as acquisitions", {

  root <- ragged_fixture()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  m <- suppressWarnings(bids_parse_filenames(root))

  # task-rest_pet.json sits at the study root: it is an inheritance sidecar, not
  # an acquisition, and would otherwise appear as a measurement of its own
  expect_false(any(is.na(m$sub)))
})

test_that("an image and its sidecar are one acquisition", {

  root <- tempfile("kinfitr-pair-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/pet/sub-01_pet.nii.gz"))
  file.create(file.path(root, "sub-01/pet/sub-01_pet.json"))

  expect_equal(nrow(suppressWarnings(bids_parse_filenames(root))), 1)
})

test_that("a _pet.json with no image is still an acquisition", {

  # An acquisition whose image is not found is still a measurement, so blood
  # arriving before the image still has something to attach to.
  root <- tempfile("kinfitr-jsononly-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/pet/sub-01_pet.json"))
  file.create(file.path(root, "sub-01/pet/sub-01_recording-manual_blood.tsv"))

  m <- suppressWarnings(bids_parse_filenames(root))
  expect_equal(nrow(m), 1)
  expect_equal(sum(m$filedata[[1]]$measurement == "blood"), 1)
})

test_that("two images describing one acquisition is an error", {

  root <- tempfile("kinfitr-dupe-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/pet/sub-01_pet.nii.gz"))
  file.create(file.path(root, "sub-01/pet/sub-01_pet.nii"))

  # A duplicate is the validator's business, but silently keeping whichever
  # comes first would be ours.
  expect_error(suppressWarnings(bids_parse_filenames(root)),
               "More than one PET image")
})

test_that("bids_parse_files keeps its old behaviour and says it is deprecated", {

  root <- ragged_fixture()
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  expect_warning(bids_parse_files(root), "deprecated")

  # The old name is retained unchanged, grid and all, so existing callers keep
  # working until they move over
  old <- suppressWarnings(bids_parse_files(root))
  new <- suppressWarnings(bids_parse_filenames(root))
  expect_gt(nrow(old), nrow(new))
})

test_that("a shared inheritance sidecar in pet/ is not an acquisition", {

  # One sub-01_pet.json serving run-01 and run-02 is metadata for both runs
  # (BIDS inheritance), not a third, run-less acquisition.
  root <- tempfile("kinfitr-shared-sidecar-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/pet/sub-01_run-01_pet.nii.gz"))
  file.create(file.path(root, "sub-01/pet/sub-01_run-02_pet.nii.gz"))
  file.create(file.path(root, "sub-01/pet/sub-01_pet.json"))

  m <- suppressWarnings(bids_parse_filenames(root))

  expect_equal(nrow(m), 2)
  expect_setequal(m$run, c("01", "02"))

  # The sidecar reaches both runs rather than becoming a measurement of its own
  for (fd in m$filedata) {
    expect_true(any(grepl("sub-01_pet\\.json$", fd$path)))
  }
})

test_that("a json-only acquisition survives beside unrelated images", {

  # The shared-sidecar rule must not swallow D18's case: a _pet.json that
  # applies to no other acquisition is still a measurement of its own.
  root <- tempfile("kinfitr-jsononly-beside-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/pet"), recursive = TRUE)
  dir.create(file.path(root, "sub-02/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/pet/sub-01_pet.nii.gz"))
  file.create(file.path(root, "sub-02/pet/sub-02_pet.json"))

  m <- suppressWarnings(bids_parse_filenames(root))

  expect_equal(nrow(m), 2)
  expect_setequal(m$sub, c("01", "02"))
})

test_that("a sidecar naming an entity the image omits does not apply", {

  # A sidecar specific enough to name an entity applies only to data that names
  # it too. task-rest_pet.json does not describe an image with no task, any
  # more than trc-abc_pet.json describes one with no tracer -- attaching the
  # wrong tracer's sidecar would silently supply the wrong half-life and dose.
  root <- tempfile("kinfitr-strict-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/pet/sub-01_pet.nii.gz"))
  writeLines('{"FrameTimesStart":[0,60],"FrameDuration":[60,60]}',
             file.path(root, "task-rest_pet.json"))
  writeLines('{"FrameTimesStart":[0,60],"FrameDuration":[60,60]}',
             file.path(root, "trc-abc_pet.json"))

  m <- bids_parse_filenames(root)
  expect_equal(nrow(m), 1)
  expect_false(any(grepl("task-rest_pet\\.json$", m$filedata[[1]]$path)))
  expect_false(any(grepl("trc-abc_pet\\.json$", m$filedata[[1]]$path)))

  # With its frame times only in inapplicable sidecars, the acquisition is
  # excluded -- with a warning, not silently
  expect_warning(s <- bids_parse_study(root), "No _pet.json sidecar")
  expect_equal(nrow(s), 0)
})

test_that("a study-level sidecar naming no extra entities is inherited", {

  # The other direction is ordinary inheritance and always applies: an image
  # accepts a sidecar that names fewer entities than it does.
  root <- tempfile("kinfitr-inherit-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/pet/sub-01_task-rest_pet.nii.gz"))
  writeLines('{"FrameTimesStart":[0,60],"FrameDuration":[60,60]}',
             file.path(root, "task-rest_pet.json"))

  m <- bids_parse_filenames(root)
  expect_equal(nrow(m), 1)
  expect_true(any(grepl("task-rest_pet\\.json$", m$filedata[[1]]$path)))

  s <- bids_parse_study(root)
  expect_equal(nrow(s), 1)
  expect_equal(nrow(s$tactimes[[1]]), 2)
})
