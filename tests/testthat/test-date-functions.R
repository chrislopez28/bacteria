context("Date Functions")

testdf <- data.frame(SampleDate = as.Date(c("2017-01-01", "2016-02-02", "2018-03-03")),
                     Result = c(100, 200, 300))

test_that("first_date returns earliest SampleDate", {
  expect_equal(first_date(testdf), as.Date("2016-02-02"))
})


test_that("last_date returns last SampleDate", {
  expect_equal(last_date(testdf), as.Date("2018-03-03"))
})
