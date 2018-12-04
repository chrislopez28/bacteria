context("Cleaning Functions")

test_that("column check is case sensitive", {
  expect_false(col_check(data.frame(SampleDate = 1), "sampledate"))
})
