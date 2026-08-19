# bids_associate_blood() decides which blood files belong to one acquisition.

blood_filedata <- function(paths, key = "sub-01_ses-test") {

  fd <- tibble::tibble(
    path = paths,
    path_absolute = file.path("/study", paths),
    extension = ifelse(grepl("\\.json$", paths), "json", "tsv"),
    measurement = ifelse(grepl("_blood\\.", paths), "blood", "pet"),
    recording = stringr::str_match(paths, "recording-([a-zA-Z0-9]+)")[, 2]
  )
  attr(fd, "pet_key") <- key
  fd
}

test_that("bids_associate_blood pairs a tsv with its json per recording", {

  fd <- blood_filedata(c(
    "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-manual_blood.tsv",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-manual_blood.json"))

  blood <- bids_associate_blood(fd)

  expect_equal(nrow(blood), 1)
  expect_equal(blood$recording, "manual")
  expect_true(grepl("blood\\.tsv$", blood$tsv))
  expect_true(grepl("blood\\.json$", blood$json))
})

test_that("several recordings for one acquisition are valid", {

  # BIDS explicitly permits manual samples alongside an autosampler
  fd <- blood_filedata(c(
    "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-manual_blood.tsv",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-manual_blood.json",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-autosampler_blood.tsv",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-autosampler_blood.json"))

  blood <- bids_associate_blood(fd)

  expect_equal(nrow(blood), 2)
  expect_setequal(blood$recording, c("autosampler", "manual"))
})

test_that("two blood files claiming one recording is an error", {

  fd <- blood_filedata(c(
    "sub-06/ses-test/pet/sub-06_ses-test_recording-manual_blood.tsv",
    "sub-06/ses-test/pet/sub-06_ses-test_recording-manual_blood.json",
    "sub-06/ses-test/pet/sub-06_ses-test_rec-A_recording-manual_blood.tsv",
    "sub-06/ses-test/pet/sub-06_ses-test_rec-A_recording-manual_blood.json"),
    key = "sub-06_ses-test_rec-A")

  # Both claim to be the manual recording for this acquisition, and there is
  # nothing to choose between them.
  expect_error(bids_associate_blood(fd), "claims recording-manual")
  expect_error(bids_associate_blood(fd), "sub-06_ses-test_rec-A")
})

test_that("blood outside pet/ is refused", {

  fd <- blood_filedata(c(
    "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz",
    "sub-01/ses-test/blood/sub-01_ses-test_recording-manual_blood.tsv"))

  expect_error(bids_associate_blood(fd), "outside a pet/ directory")
})

test_that("blood without a recording entity is refused", {

  fd <- blood_filedata(c(
    "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz",
    "sub-01/ses-test/pet/sub-01_ses-test_blood.tsv"))

  # recording is what tells one recording from another; without it there is no
  # way to know whether two files are the same measurement or different ones.
  expect_error(bids_associate_blood(fd), "no recording entity")
})

test_that("an acquisition with no blood is recorded, not silently empty", {

  fd <- blood_filedata("sub-05/ses-test/pet/sub-05_ses-test_pet.nii.gz")

  blood <- bids_associate_blood(fd)

  # Not an error and not an exclusion: reference tissue models need no blood at
  # all. There is simply no blood to associate.
  expect_equal(nrow(blood), 0)
  # But it says why, so it cannot be confused with blood that failed to attach
  expect_equal(attr(blood, "excluded"), "no blood files")
})

test_that("an incomplete pair yields no blood, and says which half is missing", {

  fd <- blood_filedata(c(
    "sub-02/ses-test/pet/sub-02_ses-test_pet.nii.gz",
    "sub-02/ses-test/pet/sub-02_ses-test_recording-manual_blood.tsv"))

  expect_warning(blood <- bids_associate_blood(fd), "missing json")

  expect_equal(nrow(blood), 0)
  expect_match(attr(blood, "excluded"), "missing json")
  expect_match(attr(blood, "excluded"), "recording-manual")
})

test_that("bids_parse_filenames attaches the acquisition's identity to its files", {

  root <- tempfile("kinfitr-attrs-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/ses-test/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz"))

  m <- suppressWarnings(bids_parse_filenames(root))

  # Carried as attributes rather than columns or an extra argument, so
  # bids_create_blooddata() keeps taking a bare table
  expect_equal(attr(m$filedata[[1]], "pet_key"), "sub-01_ses-test")
  expect_equal(normalizePath(attr(m$filedata[[1]], "study_root")),
               normalizePath(root))
})

test_that("blood files differing by a selector entity do not collide", {

  # Two recording-manual files for one subject, one per run. They belong to
  # different acquisitions, so they are never compared with each other:
  # filedata holds a single acquisition's files.
  root <- tempfile("kinfitr-runs-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  pet <- file.path(root, "sub-01/ses-test/pet")
  dir.create(pet, recursive = TRUE)
  for (run in c("01", "02")) {
    file.create(file.path(pet, sprintf("sub-01_ses-test_run-%s_pet.nii.gz", run)))
    file.create(file.path(pet, sprintf("sub-01_ses-test_run-%s_recording-manual_blood.tsv", run)))
    file.create(file.path(pet, sprintf("sub-01_ses-test_run-%s_recording-manual_blood.json", run)))
  }

  m <- suppressWarnings(bids_parse_filenames(root))
  expect_equal(nrow(m), 2)

  blood <- lapply(m$filedata, bids_associate_blood)

  expect_equal(vapply(blood, nrow, integer(1)), c(1L, 1L))
  expect_true(grepl("run-01", blood[[1]]$tsv))
  expect_true(grepl("run-02", blood[[2]]$tsv))
})

test_that("an acquisition without blood still parses and is returned", {

  # The excluded attribute describes the blood association, not the
  # acquisition: nothing drops a scan for having no blood.
  root <- tempfile("kinfitr-noblood-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  pet <- file.path(root, "sub-01/pet")
  dir.create(pet, recursive = TRUE)
  file.create(file.path(pet, "sub-01_pet.nii.gz"))
  writeLines('{"FrameTimesStart":[0,60],"FrameDuration":[60,60]}',
             file.path(pet, "sub-01_pet.json"))

  m <- suppressWarnings(bids_parse_filenames(root))
  expect_equal(nrow(m), 1)

  study <- suppressWarnings(bids_parse_study(root))
  expect_equal(nrow(study), 1)
})

test_that("blood present but unusable is warned about, not silently dropped", {

  # A manual tsv with no json cannot be read, but it must not be
  # indistinguishable from an acquisition that never had blood.
  fd <- blood_filedata(c(
    "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-manual_blood.tsv"))

  expect_warning(res <- bids_parse_blood(fd), "complete tsv/json pair")
  expect_warning(bids_parse_blood(fd), "sub-01_ses-test")
  expect_true(identical(res, NA))
})

test_that("an acquisition with no blood at all passes in silence", {

  fd <- blood_filedata("sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz")

  expect_no_warning(res <- bids_parse_blood(fd))
  expect_true(identical(res, NA))
})

test_that("a json naming entities its tsv does not is not its sidecar", {

  # recording alone is not enough to pair the files: a rec-A json beside a tsv
  # that names no rec may describe different data, and its units and
  # availability flags would then be wrong for this tsv.
  fd <- blood_filedata(c(
    "sub-01/pet/sub-01_rec-A_pet.nii.gz",
    "sub-01/pet/sub-01_recording-manual_blood.tsv",
    "sub-01/pet/sub-01_rec-A_recording-manual_blood.json"),
    key = "sub-01_rec-A")

  warnings <- character(0)
  blood <- withCallingHandlers(
    bids_associate_blood(fd),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })

  expect_match(warnings[1], "not treated as that tsv's sidecar")
  # The mismatched pair is then incomplete, not silently accepted
  expect_match(warnings[2], "missing json")
  expect_equal(nrow(blood), 0)
})

test_that("an incomplete recording is reported even when another is complete", {

  # A manual pair plus an autosampler tsv missing its json: the manual
  # recording is returned, and the autosampler one must not vanish in silence.
  fd <- blood_filedata(c(
    "sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-manual_blood.tsv",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-manual_blood.json",
    "sub-01/ses-test/pet/sub-01_ses-test_recording-autosampler_blood.tsv"))

  expect_warning(blood <- bids_associate_blood(fd),
                 "recording-autosampler missing json")

  expect_equal(blood$recording, "manual")
})

test_that("blood in another directory's pet/ is refused", {

  # Subject-level blood under sub-01/pet/ names no session, so the subset rule
  # alone would let it attach to the session-level acquisition. BIDS places
  # blood alongside the PET data it belongs to: same directory, not merely
  # some pet/ directory.
  root <- tempfile("kinfitr-bloodscope-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/ses-01/pet"), recursive = TRUE)
  dir.create(file.path(root, "sub-01/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/ses-01/pet/sub-01_ses-01_pet.nii.gz"))
  file.create(file.path(root, "sub-01/pet/sub-01_recording-manual_blood.tsv"))
  file.create(file.path(root, "sub-01/pet/sub-01_recording-manual_blood.json"))

  m <- suppressWarnings(bids_parse_filenames(root))
  expect_error(bids_associate_blood(m$filedata[[1]]),
               "not beside its PET data")
})

test_that("an inherited sidecar does not widen the blood anchor", {

  # A subject-level sub-01_pet.json attaches to the session acquisition as
  # inherited metadata, but its directory is not where the acquisition lives.
  # Blood there must still be refused: the anchor is where the acquisition's
  # own files sit, not every directory a sidecar was inherited from.
  root <- tempfile("kinfitr-anchor-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "sub-01/ses-01/pet"), recursive = TRUE)
  dir.create(file.path(root, "sub-01/pet"), recursive = TRUE)
  file.create(file.path(root, "sub-01/ses-01/pet/sub-01_ses-01_pet.nii.gz"))
  file.create(file.path(root, "sub-01/pet/sub-01_pet.json"))
  file.create(file.path(root, "sub-01/pet/sub-01_recording-manual_blood.tsv"))
  file.create(file.path(root, "sub-01/pet/sub-01_recording-manual_blood.json"))

  m <- suppressWarnings(bids_parse_filenames(root))
  expect_equal(nrow(m), 1)
  expect_equal(attr(m$filedata[[1]], "pet_dir"), "sub-01/ses-01/pet")
  expect_error(bids_associate_blood(m$filedata[[1]]),
               "not beside its PET data")
})
