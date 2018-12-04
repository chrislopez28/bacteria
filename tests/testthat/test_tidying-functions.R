context("Tidying Functions")

# TODO: Add test for tidying functions

bcols <- c("StationCode", "SampleDate", "WeatherCondition",
           "AnalyteName", "Result", "Samples", "ResQualCode", "MDL", "RL")
bnames <- c("E. coli", "Coliform, Fecal", "Coliform, Total", "Enterococcus")

untidy <- data.frame(
  StationCode = "Station",
  SampleDate = as.Date(c(rep("2011-01-04", 4), rep("2011-02-05", 4))),
  WeatherCondition = "Wet",
  AnalyteName = rep(c("E. coli", "Coliform, Fecal", "Enterococcus", "Coliform, Total"), 2),
  Result = 101:108, Samples = 1, ResQualCode = "=", MDL = 10, RL = NA)

tidy <- data.frame(
  StationCode = "Station", SampleDate = as.Date(c("2011-01-04", "2011-02-05")),
  WeatherCondition = "Wet",
  ecoli = c(101, 105), ecoli_n = 1, ecoli_qual = "=", ecoli_mdl = 10,
  ecoli_rl = NA,
  fecal_coliform = c(102, 106), fc_n = 1, fc_qual = "=", fc_mdl = 10, fc_rl = NA,
  total_coliform = c(104, 108), tc_n = 1, tc_qual = "=", tc_mdl = 10, tc_rl = NA,
  enterococcus = c(103, 107), ent_n = 1, ent_qual = "=", ent_mdl = 10, ent_rl = NA
)

test_that("tidy_bacteria correclty tidies input data frame", {
  expect_equal(tidy_bacteria(untidy), tidy)
})
