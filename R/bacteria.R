#' Column Headers Check
#'
#' Checks if a data frame includes specified column names.
#'
#' @param df a data frame of water quality data
#' @param x a character vector of column names
#' @return a logical value that is TRUE if all elements of names are columns in the data frame
#' @examples
#' col_check(df, c("StationCode", "SampleDate", "WeatherCondition"))
#' col_check(df, required_columns)
#' @export

col_check <- function(df, x){
  stopifnot(is.character(x))

  nametest <- x %in% names(df)

  if (FALSE %in% nametest) return(FALSE) else return(TRUE)
}


average_results_daily <- function(df){
  df_new <- df %>%
    group_by(StationCode, SampleDate, AnalyteName) %>%
    dplyr::summarize(Result = mean(Result), Samples = n(),
              WeatherCondition = calc_mode(WeatherCondition),
              ResQualCode = ifelse(length(ResQualCode) > 1, "avg", ResQualCode),
              MDL = max(MDL),
              RL = max(RL)) %>%
    as.data.frame

  # df_old <- df %>%
  #   select(StationCode, SampleDate, AnalyteName, WeatherCondition, ResQualCode, MDL, RL)
  #
  # df <- semi_join(x = df_new, y = df_old, by = c("StationCode", "SampleDate", "AnalyteName")) %>%
  #   as.data.frame()

  df_new
}

#' Calculate Mode
#'
#' Compute the sample mode
#'
#' @param x a vector
#' @return This function returns a length-one object of the same type as x. If there are multiple
#' modes, this function returns the mode with the value that appears first in vector x.
#' @examples
#' calc_mode(c(1, 2, 3, 3, 4, 5))
#' calc_mode(c("ab", "ab", "ac", "ad"))
#' @export

calc_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' Calculate Geometric Means
#'
#' Expands the input data frame by adding columns for 42-day or 30-day rolling
#' geometric means and counts of samples used to calculate the geometric mean.
#'
#' @param df a data frame of indicator bacteria monitoring data that has been
#' tidied by tidy_bacteria(df) and has daily rows between the first and last
#' dates in the "SampleDate" column using expand_dates(df)
#'
#' @param six_week a logical value indicating whether the calculation should use
#'  a six-week (i.e. 42 day) period for the geometric mean calculation. If TRUE,
#'  the function will calculate 42-day geometric means on every Sunday. If
#'  FALSE, the function will calculate 30-day geometric means on each date that
#'  there is a result for the associated constituent (either E. coli, Fecal
#'  Coliform, Total Coliform, or Enterococcus)
#'
#' @return This function returns a data frame with additional columns for geometric means and sample counts.
#' The geometric mean columns are ecoli_geomean, fc_geomean, tc_geomean, and ent_geomean. The sample count
#' columns are ecoli_geo_count, fc_geo_count, tc_geo_count, and ent_geo_count.
#' @export

bact_geomeans <- function(df, six_week = TRUE, ...){
  # TODO: Check if df has appropriate columns / Use colcheck()
  # TODO: Check if df has consecutive SampleDates from first row to last row / Create consecutivecheck()

  dt <- tibble::frame_data(
    ~result, ~geomean, ~geocount,
    "ecoli", "ecoli_geomean", "ecoli_geo_count",
    "fecal_coliform", "fc_geomean", "fc_geo_count",
    "total_coliform", "tc_geomean", "tc_geo_count",
    "enterococcus", "ent_geomean", "ent_geo_count"
  )

  #
  # dt <- tibble::frame_data(
  #   ~result, ~qual, ~mdl, ~rl, ~WQO_ss, ~WQO_ss_val, ~exceed_WQO_ss, ~geomean, ~geocount,  ~WQO_geo, ~WQO_geo_val, ~exceed_WQO_geo, ~other,
  #   "ecoli", "ecoli_qual", "ecoli_mdl", "ecoli_rl", "ecoli_WQO_ss", NA, "exceed_ecoli_WQO_ss", "ecoli_geomean", "ecoli_geo_count", "ecoli_WQO_geo", NA, "exceed_ecoli_WQO_geo", NA,
  #   "fecal_coliform", "fc_qual", "fc_mdl", "fc_rl", "fc_WQO_ss", 400, "exceed_fc_WQO_ss", "fc_geomean", "fc_geo_count", "fc_WQO_geo", 200, "exceed_fc_WQO_geo", NA,
  #   "total_coliform", "tc_qual", "tc_mdl", "tc_rl", "tc_WQO_ss", 10000, "exceed_tc_WQO_ss", "tc_geomean", "tc_geo_count", "tc_WQO_geo", 1000, "exceed_tc_WQO_geo", NA,
  #   "enterococcus", "ent_qual", "ent_mdl", "ent_rl", "ent_WQO_ss", 104, "exceed_ent_WQO_ss", "ent_geomean", "ent_geo_count", "ent_WQO_geo", 35, "exceed_ent_WQO_geo", NA
  # )


  if (six_week){

    for (i in seq_along(dt$result)){
      df <- df %>% dplyr::mutate(!!as.name(dt$geomean[[i]]) := dplyr::if_else(lubridate::wday(SampleDate, label = TRUE) == "Sun",
                                                                              zoo::rollapply(!!as.name(dt$result[[i]]), 42, function(x) EnvStats::geoMean(x, na.rm = TRUE), fill = NA, align = "right"),
                                                                              as.double(NA)),
                                 !!as.name(dt$geocount[[i]]) := dplyr::if_else(lubridate::wday(SampleDate, label = TRUE) == "Sun",
                                                                               zoo::rollapply(!!as.name(dt$result[[i]]), 42, function(x) length(x[!is.na(x)]), fill = NA, align = "right"),
                                                                               as.integer(NA)))
    }

  } else {

    for (i in seq_along(dt$result)){
      df <- df %>% dplyr::mutate(!!as.name(dt$geomean[[i]]) := dplyr::if_else(Data_Row == TRUE,
                                                                              zoo::rollapply(!!as.name(dt$result[[i]]), 30, function(x) EnvStats::geoMean(x, na.rm = TRUE), fill = NA, align = "right"),
                                                                              as.double(NA)),
                                 !!as.name(dt$geocount[[i]]) := dplyr::if_else(Data_Row == TRUE,
                                                                               zoo::rollapply(!!as.name(dt$result[[i]]), 30, function(x) length(x[!is.na(x)]), fill = NA, align = "right"),
                                                                               as.integer(NA)))
    }

  }

  return(df)
}


#' Check Geometric Means Against Limitations
#'
#' Expands the input data frame to include columns that list relevant geometric mean water quality
#' objectives and compliance with the objective.
#'
#' @param df a data frame of indicator bacteria monitoring data that has been tidied by
#' tidy_bacteria(df), has daily rows between the first and last SampleDate using expand_dates(df)
#' and includes geometric mean calculations using bact_geomeans(df)
#' @param six_week a logical value indicating whether 42-day geomeans were calculated every Sunday
#' periods or 30-day geomeans were calculated every day with a result.
#' @return an expanded data frame with columns for relevant water quality objectives and compliance
#' @export

check_geolimits <- function(df, BU = "REC-1", water_type = "marine", six_week = TRUE, ...){
  # TODO: Check if df has appropriate columns / Use colcheck()
  # TODO: Check if df has consecutive SampleDates from first row to last row / Create consecutivecheck()

  limits <- grab_limits_geo(BU, water_type)

  dt <- tibble::frame_data(
    ~geomean, ~geocount,  ~WQO_geo, ~WQO_geo_val, ~exceed_WQO_geo,
    "ecoli_geomean", "ecoli_geo_count", "ecoli_WQO_geo", limits["ecoli_WQO"], "exceed_ecoli_WQO_geo",
    "fc_geomean", "fc_geo_count", "fc_WQO_geo", limits["fc_WQO"], "exceed_fc_WQO_geo",
    "tc_geomean", "tc_geo_count", "tc_WQO_geo", limits["tc_WQO"], "exceed_tc_WQO_geo",
    "ent_geomean", "ent_geo_count", "ent_WQO_geo", limits["ent_WQO"], "exceed_ent_WQO_geo"
  )

  if (six_week){
    for (i in seq_along(dt$geomean)){
      df <- df %>% dplyr::mutate(!!as.name(dt$WQO_geo[[i]]) := dplyr::if_else(lubridate::wday(SampleDate, label = TRUE) == "Sun", dt$WQO_geo_val[[i]], as.double(NA))) %>%
        dplyr::mutate(!!as.name(dt$exceed_WQO_geo[[i]]) :=
                        dplyr::if_else(!!as.name(dt$geomean[[i]]) > !!as.name(dt$WQO_geo[[i]]) & !!as.name(dt$geocount[[i]]) >= 5, TRUE, FALSE))
    }
  } else {

    for (i in seq_along(dt$geomean)){
      df <- df %>% dplyr::mutate(!!as.name(dt$WQO_geo[[i]]) := dplyr::if_else(Data_Row == TRUE, dt$WQO_geo_val[[i]], as.double(NA))) %>%
        dplyr::mutate(!!as.name(dt$exceed_WQO_geo[[i]]) :=
                        dplyr::if_else(!!as.name(dt$geomean[[i]]) > !!as.name(dt$WQO_geo[[i]]) & !!as.name(dt$geocount[[i]]) >= 5, TRUE, FALSE))
    }

  }
  return(df)
}


#' Check Single Sample Results Against Limitations
#'
#' Expands the input data frame to include columns that list relevant single sample water quality
#' objectives and compliance with the objective.
#'
#' @param df a data frame of indicator bacteria monitoring data that has been tidied by
#' tidy_bacteria(df)
#' @return an expanded data frame with columns for relevant water quality objectives and compliance
#' @export

check_sslimits <- function(df, BU = "REC-1", water_type = "marine"){
  # TODO: Check if df has appropriate columns / Use colcheck()
  # TODO: Check if df has consecutive SampleDates from first row to last row / Create consecutivecheck()

  limits <- grab_limits_ss(BU, water_type)

  dt <- tibble::frame_data(
    ~result, ~qual, ~mdl, ~rl, ~WQO_ss, ~WQO_ss_val, ~exceed_WQO_ss,
    "ecoli", "ecoli_qual", "ecoli_mdl", "ecoli_rl", "ecoli_WQO_ss", limits["ecoli_WQO"], "exceed_ecoli_WQO_ss",
    "fecal_coliform", "fc_qual", "fc_mdl", "fc_rl", "fc_WQO_ss", limits["fc_WQO"], "exceed_fc_WQO_ss",
    "total_coliform", "tc_qual", "tc_mdl", "tc_rl", "tc_WQO_ss", limits["tc_WQO"], "exceed_tc_WQO_ss",
    "enterococcus", "ent_qual", "ent_mdl", "ent_rl", "ent_WQO_ss", limits["ent_WQO"], "exceed_ent_WQO_ss"
  )

  for (i in seq_along(dt$result)){
    df <- df %>% dplyr::mutate(!!as.name(dt$WQO_ss[[i]]) := dplyr::if_else(!is.na(!!as.name(dt$result[[i]])),
                                                             dt$WQO_ss_val[[i]],
                                                             as.double(NA)),
                        !!as.name(dt$exceed_WQO_ss[[i]]) := dplyr::if_else(!!as.name(dt$result[[i]]) > !!as.name(dt$WQO_ss[[i]]),
                                                                    TRUE,
                                                                    FALSE))
  }
  return(df)
}


#' Convert Weather Condition
#'
#' Convert WeatherCondition entries of "Dry" to "Winter Dry" and "Summer Dry" depending on the date.
#' "Winter Dry" are dry weather days occuring on November 1 through March 31. "Summer Dry" are dry
#' weather days occurring on April 1 through October 31.
#'
#' @param df a data frame with WeatherCondition and SampleDate columns
#' @return a data frame
#' @export

convertWeather <- function(df){
  df <- df %>%
    dplyr::mutate(WeatherCondition = dplyr::if_else(WeatherCondition == "Dry",
                                      dplyr::if_else(month(SampleDate) %in% c(11, 12, 1, 2, 3), "Winter Dry", "Summer Dry"),
                                      WeatherCondition))
  return(df)
}

#' Check Fecal-To-Total Ratio
#'
#' Inserts columns that calculate the fecal coliform to total coliform ratio, show the single sample
#' objective for total coliform if applicable, and show compliance with applicable objectives.
#'
#' @param df a tidy data frame
#' @param water_type a string of either "fresh" or "marine' representing the water type
#'     of the water body
#' @return a data frame

check_fecal_to_total <- function(df, ...){
  df %>% dplyr::mutate(fc_to_tc = fecal_coliform/total_coliform,
                  tc_WQO_ss_2 = dplyr::if_else(fc_to_tc > 0.1, 1000, as.double(NA)),
                  exceed_tc_WQO_ss_2 = dplyr::if_else(total_coliform > tc_WQO_ss_2, TRUE, FALSE))
}

exceed_ss <- function(df){
  df <- df %>%
    dplyr::mutate(exceed_day = dplyr::if_else(exceed_ecoli_WQO_ss == TRUE | exceed_fc_WQO_ss == TRUE | exceed_tc_WQO_ss == TRUE | exceed_ent_WQO_ss == TRUE | exceed_tc_WQO_ss_2 == TRUE,
                                TRUE, FALSE)) %>%
    dplyr::mutate(exceed_day = exceed_day & !is.na(exceed_day))
  return(df)
}

set_start <- function(dt, april_start){
  stopifnot(lubridate::is.Date(dt))
  start <- dt
  lubridate::day(start) <- 1

  if (april_start == FALSE){
    lubridate::year(start) <- lubridate::year(dt)
    lubridate::month(start) <- 11
    lubridate::day(start) <- 1
  } else if (april_start == TRUE){
    lubridate::year(start) <- lubridate::year(dt)
    lubridate::month(start) <- 4
    lubridate::day(start) <- 1
  }

  if (start > dt){
    lubridate::year(start) <- (lubridate::year(start) - 1)
  }

  return(start)
}

set_end <- function(dt, april_start){
  stopifnot(lubridate::is.Date(dt))
  end <- dt
  lubridate::day(end) <- 1

  if (april_start == FALSE){
    lubridate::year(end) <- lubridate::year(dt)
    lubridate::month(end) <- 10
    lubridate::day(end) <- 31
  } else if (april_start == TRUE){
    lubridate::year(end) <- lubridate::year(dt)
    lubridate::month(end) <- 3
    lubridate::day(end) <- 31
  }

  if (end < dt){
    lubridate::year(end) <- (lubridate::year(end) + 1)
  }

  return(end)
}


#' Calculate Annual Exceedances
#'
#' Calculates the number of annual exceedances
#'
#' @param tidy_df a df
#' @return a data frame of results
annual_exceedances <- function(tidy_df, station, start_date, end_date, april_start = FALSE){
  start <- set_start(start_date, april_start)
  end <- set_end(end_date, april_start)

  winter_start <- start
  winter_end <- start + months(5) - lubridate::days(1)
  summer_start <- start + months(5)
  summer_end <- start + lubridate::years(1) - lubridate::days(1)

  period <- length(lubridate::year(start):(lubridate::year(end) - lubridate::years(1)))

  result <- data.frame(StationCode = character(period), year_start = integer(period), year_end = integer(period),
                       winter_dry = integer(period), winter_dry_n = integer(period),
                       summer_dry = integer(period), summer_dry_n = integer(period),
                       wet = integer(period), wet_n = integer(period), stringsAsFactors = FALSE,
                       ec_geomean = integer(period), fc_geomean = integer(period),
                       tc_geomean = integer(period), ent_geomean = integer(period),
                       ec_geocalc_n = integer(period), fc_geocalc_n = integer(period),
                       tc_geocalc_n = integer(period), ent_geocalc_n = integer(period))

  for (i in seq_along(result$year_start)){
    result$StationCode[[i]] <- station
    result$year_start[[i]] <- winter_start
    result$year_end[[i]] <- summer_end

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= winter_end) %>%
      dplyr::filter(WeatherCondition == "Dry") %>%
      dplyr::summarize(exceeds = sum(exceed_day, na.rm = TRUE), samples = n())

    result$winter_dry[[i]] <- temp[["exceeds"]][1]
    result$winter_dry_n[[i]] <- temp[["samples"]][1]

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= summer_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::filter(WeatherCondition == "Dry") %>%
      dplyr::summarize(exceeds = sum(exceed_day, na.rm = TRUE), samples = n())

    result$summer_dry[[i]] <- temp[["exceeds"]][1]
    result$summer_dry_n[[i]] <- temp[["samples"]][1]

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::filter(WeatherCondition == "Wet") %>%
      dplyr::summarize(exceeds = sum(exceed_day, na.rm = TRUE), samples = n())

    result$wet[[i]] <- temp[["exceeds"]][1]
    result$wet_n[[i]] <- temp[["samples"]][1]

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::summarize(ec.geo = sum(exceed_ecoli_WQO_geo, na.rm = TRUE),
                       fc.geo = sum(exceed_fc_WQO_geo, na.rm = TRUE),
                       tc.geo = sum(exceed_tc_WQO_geo, na.rm = TRUE),
                       ent.geo = sum(exceed_ent_WQO_geo, na.rm = TRUE))

    result$ec_geomean[[i]] <- temp[["ec.geo"]][1]
    result$fc_geomean[[i]] <- temp[["fc.geo"]][1]
    result$tc_geomean[[i]] <- temp[["tc.geo"]][1]
    result$ent_geomean[[i]] <- temp[["ent.geo"]][1]



    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::filter(ecoli_geo_count >= 5) %>%
      dplyr::summarize(ec.geo.n = n())

    result$ec_geocalc_n[[i]] <- temp[["ec.geo.n"]][[1]]

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::filter(fc_geo_count >= 5) %>%
      dplyr::summarize(fc.geo.n = n())

    result$fc_geocalc_n[[i]] <- temp[["fc.geo.n"]][[1]]

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::filter(tc_geo_count >= 5) %>%
      dplyr::summarize(tc.geo.n = n())

    result$tc_geocalc_n[[i]] <- temp[["tc.geo.n"]][[1]]

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::filter(ent_geo_count >= 5) %>%
      dplyr::summarize(ent.geo.n = n())

    result$ent_geocalc_n[[i]] <- temp[["ent.geo.n"]][[1]]

    winter_start <- winter_start + lubridate::years(1)
    winter_end <- winter_end + lubridate::years(1)
    summer_start <- summer_start + lubridate::years(1)
    summer_end <- summer_end + lubridate::years(1)
  }
  result$year_start <- as.Date(result$year_start, origin = "1970-01-01")
  result$year_end <- as.Date(result$year_end, origin = "1970-01-01")

  return(result)
}


#' Calculate Geomean Annual Exceedances
#'
#' Calculates the number of annual exceedances
#'
#' @param tidy_df a df
#' @return a data frame of results

annual_geo_exceedances <- function(tidy_df, station, start_date, end_date, april_start = FALSE){
  start <- set_start(start_date, april_start)
  end <- set_end(end_date, april_start)

  winter_start <- start
  winter_end <- start + months(5) - lubridate::days(1)
  summer_start <- start + months(5)
  summer_end <- start + years(1) - lubridate::days(1)

  period <- length(year(start):(year(end) - years(1)))

  result <- data.frame(StationCode = character(period), year_name = integer(period),
                       ecoli_geo_exceeds = integer(period), ecoli_n = integer(period),
                       tc_geo_exceed = integer(period), tc_n = integer(period),
                       fc_geo_exceed = integer(period), fc_n = integer(period),
                       ent_geo_exceed = integer(period), ent_n = integer(period),
                       stringsAsFactors = FALSE)

  for (i in seq_along(result$year_name)){
    result$StationCode[[i]] <- station
    result$year_name[[i]] <- year(winter_start)

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::mutate(ecoli_sample = ifelse(is.na(ecoli), FALSE, TRUE)) %>%
      dplyr::summarize(exceeds = sum(exceed_ecoli_WQO_geo, na.rm = TRUE), samples = sum(ecoli_sample, na.rm = TRUE))

    result$ecoli_geo_exceeds[[i]] <- temp[["exceeds"]][1]
    result$ecoli_n[[i]] <- temp[["samples"]][1]

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::mutate(tc_sample = ifelse(is.na(ecoli), FALSE, TRUE)) %>%
      dplyr::summarize(exceeds = sum(exceed_tc_WQO_geo, na.rm = TRUE), samples = sum(tc_sample, na.rm = TRUE))

    result$tc_geo_exceed[[i]] <- temp[["exceeds"]][1]
    result$tc_n[[i]] <- temp[["samples"]][1]

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::mutate(fc_sample = ifelse(is.na(ecoli), FALSE, TRUE)) %>%
      dplyr::summarize(exceeds = sum(exceed_fc_WQO_geo, na.rm = TRUE), samples = sum(fc_sample, na.rm = TRUE))

    result$fc_geo_exceed[[i]] <- temp[["exceeds"]][1]
    result$fc_n[[i]] <- temp[["samples"]][1]

    temp <- tidy_df %>%
      dplyr::filter(SampleDate >= winter_start) %>%
      dplyr::filter(SampleDate <= summer_end) %>%
      dplyr::mutate(ent_sample = ifelse(is.na(ecoli), FALSE, TRUE)) %>%
      dplyr::summarize(exceeds = sum(exceed_ent_WQO_geo, na.rm = TRUE), samples = sum(ent_sample, na.rm = TRUE))

    result$ent_geo_exceed[[i]] <- temp[["exceeds"]][1]
    result$ent_n[[i]] <- temp[["samples"]][1]

    winter_start <- winter_start + years(1)
    winter_end <- winter_end + years(1)
    summer_start <- summer_start + years(1)
    summer_end <- summer_end + years(1)
  }
  result
}


order_bacteria_columns <- function(df){

  # df %>%
  #   select(StationCode, SampleDate, WeatherCondition,
  #          ecoli, ecoli_n, ecoli_qual, ecoli_mdl, ecoli_rl, ecoli_WQO_ss, exceed_ecoli_WQO_ss,
  #          ecoli_geomean, ecoli_geo_count, ecoli_WQO_geo, exceed_ecoli_WQO_geo,
  #          fecal_coliform, fc_n, fc_qual, fc_mdl, fc_rl, fc_WQO_ss, exceed_fc_WQO_ss,
  #          fc_geomean, fc_geo_count, fc_WQO_geo, exceed_fc_WQO_geo,
  #          total_coliform, tc_n, tc_qual, tc_mdl, tc_rl, tc_WQO_ss, exceed_tc_WQO_ss,
  #          tc_geomean, tc_geo_count, tc_WQO_geo, exceed_tc_WQO_geo,
  #          enterococcus, ent_n, ent_qual, ent_mdl, ent_rl, ent_WQO_ss, exceed_ent_WQO_ss,
  #          ent_geomean, ent_geo_count, ent_WQO_geo, exceed_ent_WQO_geo,
  #          fc_to_tc, tc_WQO_ss_2, exceed_tc_WQO_ss_2,
  #          exceed_day, Data_Row,
  #          everything()
  #          )

  variables <- c(
    "StationCode", "SampleDate", "WeatherCondition",
    # E. coli
    "ecoli", "ecoli_n",  "ecoli_qual", "ecoli_mdl", "ecoli_rl", "ecoli_WQO_ss", "exceed_ecoli_WQO_ss",
    "ecoli_geomean", "ecoli_geo_count", "ecoli_WQO_geo", "exceed_ecoli_WQO_geo",
    # Fecal Coliform
    "fecal_coliform", "fc_n", "fc_qual", "fc_mdl", "fc_rl", "fc_WQO_ss", "exceed_fc_WQO_ss",
    "fc_geomean", "fc_geo_count", "fc_WQO_geo", "exceed_fc_WQO_geo",
    # Total Coliform
    "total_coliform", "tc_n", "tc_qual", "tc_mdl", "tc_rl", "tc_WQO_ss", "exceed_tc_WQO_ss",
    "tc_geomean", "tc_geo_count", "tc_WQO_geo", "exceed_tc_WQO_geo",
    # Enterococcus
    "enterococcus", "ent_n", "ent_qual", "ent_mdl", "ent_rl", "ent_WQO_ss", "exceed_ent_WQO_ss",
    "ent_geomean", "ent_geo_count", "ent_WQO_geo", "exceed_ent_WQO_geo",
    # Fecal to Total Ratio / Total Coliform Limit
    "fc_to_tc", "tc_WQO_ss_2", "exceed_tc_WQO_ss_2",
    # Additional Info
    "exceed_day", "Data_Row")

  df %>%
    select(one_of(variables), everything())

  # df <- df[, names(df) %in% c(
  #   c("StationCode", "SampleDate", "WeatherCondition"),
  #   c("ecoli", "ecoli_n",  "ecoli_qual", "ecoli_mdl", "ecoli_rl", "ecoli_WQO_ss", "exceed_ecoli_WQO_ss", "ecoli_geomean", "ecoli_geo_count", "ecoli_WQO_geo", "exceed_ecoli_WQO_geo"),
  #   c("fecal_coliform", "fc_n", "fc_qual", "fc_mdl", "fc_rl", "fc_WQO_ss", "exceed_fc_WQO_ss",  "fc_geomean", "fc_geo_count", "fc_WQO_geo", "exceed_fc_WQO_geo"),
  #   c("total_coliform", "tc_n", "tc_qual", "tc_mdl", "tc_rl", "tc_WQO_ss", "exceed_tc_WQO_ss",  "tc_geomean", "tc_geo_count", "tc_WQO_geo", "exceed_tc_WQO_geo"),
  #   c("enterococcus", "ent_n", "ent_qual", "ent_mdl", "ent_rl", "ent_WQO_ss", "exceed_ent_WQO_ss", "ent_geomean", "ent_geo_count", "ent_WQO_geo", "exceed_ent_WQO_geo"),
  #   c("fc_to_tc", "tc_WQO_ss_2", "exceed_tc_WQO_ss_2"),
  #   c("exceed_day", "Data_Row")
  # )]
  #
  # df
}


update_fecal <- function(df, sub_ecoli_for_fecal = FALSE, ...){
  if (sub_ecoli_for_fecal){

    temp <- df %>%
      dplyr::mutate(counter = ifelse(is.na(fecal_coliform) & !is.na(ecoli), TRUE, FALSE))

    if (sum(temp$counter, na.rm = TRUE) > 0) {
      print(paste0("Substituting ", sum(temp$counter, na.rm = TRUE), " E. coli results for Fecal Coliform."))
    }

    df <- df %>%
      dplyr::mutate(fc_qual = ifelse(is.na(fecal_coliform) & !is.na(ecoli),
                                     "ECsub", fc_qual)) %>%
      dplyr::mutate(fecal_coliform = ifelse(is.na(fecal_coliform) & !is.na(ecoli),
                                            ecoli, fecal_coliform))
  }
  df
}




collapse_bact_data <- function(tidy_df){
  if (is.null(tidy_df)) return(tidy_df)

  tidy_df <- tidy_df %>%
    dplyr::filter(Data_Row == TRUE | (!is.na(ecoli_geomean)) | (!is.na(fc_geomean)) | (!is.na(tc_geomean)) | (!is.na(ent_geomean)))
}

bact_heatmap <- function(tidy_df, title, subtitle){
  # Code Adapted from:
  # http://margintale.blogspot.in/2012/04/ggplot2-time-series-heatmaps.html
  df <- tidy_df %>% filter(Data_Row == TRUE)

  df <- df %>%
    dplyr::mutate(year = lubridate::year(SampleDate),
           month = lubridate::month(SampleDate),
           monthf = lubridate::month(SampleDate, label = TRUE, abbr = TRUE),
           weekday = lubridate::wday(SampleDate),
           weekdayf = lubridate::wday(SampleDate, label = TRUE, abbr = TRUE),
           week = lubridate::week(SampleDate),
           day = lubridate::day(SampleDate))

  df <- df[df$SampleDate >= first_date(tidy_df), ]  # filter reqd years

  df$monthf <- factor(df$month,levels = as.character(1:12),
                      labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul",
                                 "Aug","Sep","Oct","Nov","Dec"), ordered = TRUE)

  df$weekdayf <- factor(df$weekday,
                        levels=rev(1:7),
                        labels = rev(c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")),
                        ordered = TRUE)

  # Create Month Week
  df$yearmonth <- zoo::as.yearmon(df$SampleDate)
  df$yearmonthf <- factor(df$yearmonth)
  df <- plyr::ddply(df, .(yearmonthf), transform, monthweek=1+week-min(week))  # compute week number of month
  df <- df[, c("year", "yearmonthf", "monthf", "week", "monthweek", "weekdayf", "enterococcus", "exceed_day", "WeatherCondition", "Data_Row")]

  df <- df %>% dplyr::mutate(exceed_day = ifelse(exceed_day, "Exceedance", "No Exceedance"))

  fillcolor <- c("red", "green")
  plot_title = title
  sub_title = subtitle

  # Plot
  ggplot2::ggplot(df) +
    ggplot2::geom_tile(aes(monthweek, weekdayf, fill = exceed_day),
                       color = "black") +
    ggplot2::facet_grid(year~monthf) +
    ggplot2::scale_fill_manual(name = "", values = fillcolor) +
    ggplot2::labs(x="Week of Month",
         y="",
         title = plot_title,
         subtitle = sub_title,
         fill="") +
    ggplot2::geom_point(data = filter(df, Data_Row == TRUE),
                        aes(monthweek, weekdayf, shape = WeatherCondition)) +
    ggplot2::scale_shape_manual(values = c(0, 16))
}





