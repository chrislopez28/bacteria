#' Analyze Multiple Stations
#'
#' Runs bacteria checker on multiple stations
#'
#' @param df a data frame
#' @param sites a character vector with names of the sites
#' @param BU a character vector with the highest beneficial use for the associated site
#' @param water_type a string of either "fresh" or "marine' representing the water type for the site
#' @return a list of data frames that correspond to each site
bact_check <- function(df, sites, BU, water_type, ...){

  df <- df %>%
    average_results_daily()

  df <- df %>%
    tidy_bacteria()

  df <- df %>%
    replace_nd() %>%
    update_fecal(sub_ecoli_for_fecal, ...)

  analysis_sites <- data.frame(sites, BU, water_type, stringsAsFactors = FALSE)
  names(analysis_sites) <- c("StationCode", "BU", "water_type")

  analysis_sites <- analysis_sites[(analysis_sites$StationCode %in% unique(df$StationCode)), ]

  out <- vector("list", length(analysis_sites$StationCode))
  results <- vector("list", length(analysis_sites$StationCode))

  for (i in seq_along(analysis_sites$StationCode)){
    out[[i]] <- df %>%
      dplyr::filter(StationCode == analysis_sites$StationCode[[i]])

    out[[i]] <- expand_dates(out[[i]])
    out[[i]] <- bact_geomeans(out[[i]], ...) %>%
      check_geolimits(BU = analysis_sites$BU[[i]], water_type = analysis_sites$water_type[[i]], ...) %>%
      check_sslimits(BU = analysis_sites$BU[[i]], water_type = analysis_sites$water_type[[i]])

    if (analysis_sites$water_type[[i]] == "marine"){
      out[[i]] <- out[[i]] %>% check_fecal_to_total(...)
    } else if (analysis_sites$water_type[[i]] == "fresh"){
      out[[i]] <- out[[i]] %>% mutate("fc_to_tc" = NA, "tc_WQO_ss_2" = NA, "exceed_tc_WQO_ss_2" = NA)
    }

    out[[i]] <- out[[i]] %>%
      exceed_ss() %>%
      order_bacteria_columns()
  }
  out
}

#' Calculate Annual Exceedances
#'
#' Runs bacteria checker on multiple stations and calculates the number of annual exceedances
#'
#' @param df a data frame
#' @param sites a character vector with names of the sites
#' @param BU a character vector with the highest beneficial use for the associated site
#' @param water_type a string of either "fresh" or "marine' representing the water type for the site
#' @return a list of data frames that correspond to each site
bact_ann_exceeds <- function(df, sites, BU, water_type, ...){

  analysis_sites <- data.frame(sites, BU, water_type, stringsAsFactors = FALSE)
  names(analysis_sites) <- c("StationCode", "BU", "water_type")

  analysis_sites <- analysis_sites[(analysis_sites$StationCode %in% unique(df$StationCode)), ]

  out <- bact_check(df, sites, BU, water_type, ...)
  results <- vector("list", length(analysis_sites$StationCode))

  for (i in seq_along(analysis_sites$StationCode)){
    station <- as.character(analysis_sites$StationCode[[i]])
    start <- first_date(out[[i]])
    end <- last_date(out[[i]])

    print(i)
    print(start)
    print(end)

    results[[i]] <- annual_exceedances(out[[i]], station = station, start_date = start, end_date = end)
  }
  results
}
