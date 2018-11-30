#' Get First Date
#'
#' Returns the first date found in the SampleDate column of a data frame.
#'
#' @param df a data frame with a column of dates named "SampleDate"
#' @return This function returns a length-one Date object

first_date <- function(df){
  first <- df %>%
    dplyr::arrange(SampleDate) %>%
    as.data.frame() %>%
    dplyr::slice(1L)

  return (first[["SampleDate"]])
}

#' Get Last Date
#'
#' Returns the last date found in the SampleDate column of a data frame.
#'
#' @param df a data frame with a column of dates named "SampleDate"
#' @return This function returns a length-one Date object

last_date <- function(df){
  last <- df %>%
    dplyr::arrange(SampleDate) %>%
    as.data.frame() %>%
    dplyr::slice(n())

  return(last[["SampleDate"]])
}


#' Insert Daily Rows
#'
#' Inserts a row for each day between the first and last dates found in the SampleDate
#' column of a data frame.
#'
#' @param df a data frame of monitoring data that has been tidied by tidy_bacteria
#' @return This functions returns a data frame with the same columns as df.
#' @export

expand_dates <- function(df){
  if (length(unique(df$StationCode)) > 1) stop(">1 StationCode in data frame. This function only works for data frames with data from a single monitoring station.")

  station <- unique(df$StationCode)

  start_dt <- first_date(df)
  end_dt <- last_date(df)

  obs <- data.frame(SampleDate = seq(from = start_dt, to = end_dt, by = "1 day"), StationCode = station, stringsAsFactors = FALSE)

  df <- df %>%
    mutate(Data_Row = TRUE)

  joined <- dplyr::left_join(x = obs, y = df, by = c("SampleDate", "StationCode"))
  joined
}
