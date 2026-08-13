# Characterization tests for BIDS filename parsing.
#
# Nothing here asserts new behaviour — these pin the contract that petfit's
# subsetting validation now depends on. The most consequential part is that an
# absent entity produces no column at all, which is why callers must cope with a
# missing column rather than expecting NA.

test_that("bids_filename_attributes reads every entity in a full filename", {

  attrs <- bids_filename_attributes(
    "sub-01_ses-test_trc-pf974_rec-acdyn_task-rest_run-1_pet.nii.gz")

  expect_equal(attrs$sub, "01")
  expect_equal(attrs$ses, "test")
  expect_equal(attrs$trc, "pf974")
  expect_equal(attrs$rec, "acdyn")
  expect_equal(attrs$task, "rest")
  expect_equal(attrs$run, "1")
  expect_equal(attrs$measurement, "pet")
})

test_that("bids_filename_attributes omits absent entities entirely", {

  # The function ends in tidyr::spread(), so a filename that carries no `ses`
  # yields a table with no `ses` column — not ses = NA. petfit fills absent
  # entities with NA in its own combined TACs table, so the two layers differ,
  # and anything filtering this output must check the column exists first.

  attrs <- bids_filename_attributes("sub-01_pet.nii.gz")

  expect_equal(attrs$sub, "01")
  expect_false("ses" %in% colnames(attrs))
  expect_false("trc" %in% colnames(attrs))
  expect_false("task" %in% colnames(attrs))
  expect_false("run" %in% colnames(attrs))
})

test_that("bids_filename_attributes completes filename entities from the path", {

  # petfit's own derivatives put ses in the directory but not the filename, so
  # path entities have to be picked up for those files to be identifiable.
  attrs <- bids_filename_attributes("sub-01/ses-test/pet/sub-01_pet.nii.gz")

  expect_equal(attrs$sub, "01")
  expect_equal(attrs$ses, "test")
})

test_that("bids_filename_attributes prefers the filename over the directory", {

  # Deduplication keeps the last occurrence, and the filename comes last. Where
  # a path and a filename disagree, the filename is the more specific claim.
  attrs <- bids_filename_attributes("sub-01/ses-test/pet/sub-01_ses-retest_pet.nii.gz")

  expect_equal(attrs$ses, "retest")
})

test_that("bids_filename_attributes reads the measurement suffix", {

  expect_equal(
    bids_filename_attributes(
      "sub-01/ses-01/pet/sub-01_ses-01_recording-continuous_blood.json")$measurement,
    "blood")

  expect_equal(
    bids_filename_attributes(
      "sub-01/ses-01/pet/sub-01_ses-01_recording-continuous_blood.json")$recording,
    "continuous")
})
