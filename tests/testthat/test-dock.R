context("test-dock.R")

test_that("dock_find_tracer finds raclopride", {
  skip_if_offline()

  result <- dock_find_tracer("RAC")

  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) > 0)
  expect_true("11C_raclopride" %in% result$tracer)
  expect_true("tracer" %in% names(result))
  expect_true("target" %in% names(result))
})

test_that("dock_query_tracer returns data for 11C_raclopride", {
  skip_if_offline()

  result <- dock_query_tracer("11C_raclopride")

  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) > 0)
  expect_true(all(c("tracer", "target", "method", "reference_region",
                     "parameter", "region", "mean", "sd",
                     "author", "year", "journal", "doi",
                     "n_subjects") %in% names(result)))
  expect_true(all(result$tracer == "11C_raclopride"))
  expect_true(any(!is.na(result$mean)))
})
