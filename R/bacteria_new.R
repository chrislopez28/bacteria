


grab_limits_geo_new <- function(BU, water_type){
  # TODO: warnings
  # TODO: SHELL limits

  if (water_type == "marine"){
    if(BU == "REC-1"){
      limits <- c(NA, NA, NA, 30)
    }
  }

  if (water_type == "fresh"){
    if (BU == "REC-1"){
      limits <- c(100, NA, NA, NA)
    } else if (BU == "LREC-1"){
      limits <- c(126, NA, NA, NA)
    } else if (BU == "REC-2"){
      limits <- c(NA, 2000, NA, NA)
    }
  }

  names(limits) <- c("ecoli_WQO", "fc_WQO", "tc_WQO", "ent_WQO")
  return(limits)
}



check_geolimits_new <- function(df, BU = "REC-1", water_type = "marine"){
  # TODO: Check if df has appropriate columns / Use colcheck()
  # TODO: Check if df has consecutive SampleDates from first row to last row / Create consecutivecheck()

  limits <- grab_limits_geo_new(BU, water_type)

  dt <- tibble::frame_data(
    ~geomean, ~geocount,  ~WQO_geo, ~WQO_geo_val, ~exceed_WQO_geo,
    "ecoli_geomean", "ecoli_geo_count", "ecoli_WQO_geo", limits["ecoli_WQO"], "exceed_ecoli_WQO_geo",
    "fc_geomean", "fc_geo_count", "fc_WQO_geo", limits["fc_WQO"], "exceed_fc_WQO_geo",
    "tc_geomean", "tc_geo_count", "tc_WQO_geo", limits["tc_WQO"], "exceed_tc_WQO_geo",
    "ent_geomean", "ent_geo_count", "ent_WQO_geo", limits["ent_WQO"], "exceed_ent_WQO_geo"
  )

  for (i in seq_along(dt$geomean)){
    df <- df %>% dplyr::mutate(!!as.name(dt$WQO_geo[[i]]) := dplyr::if_else(lubridate::wday(SampleDate, label = TRUE) == "Sun", dt$WQO_geo_val[[i]], as.double(NA))) %>%
      dplyr::mutate(!!as.name(dt$exceed_WQO_geo[[i]]) :=
                      dplyr::if_else(!!as.name(dt$geomean[[i]]) > !!as.name(dt$WQO_geo[[i]]) & !!as.name(dt$geocount[[i]]) >= 5, TRUE, FALSE))
  }
  return(df)
}

check_sslimits_new <- function(df, BU = "REC-1", water_type = "marine"){
  # TODO: Check if df has appropriate columns / Use colcheck()
  # TODO: Check if df has consecutive SampleDates from first row to last row / Create consecutivecheck()

  limits <- grab_limits_ss_new(BU, water_type)

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

grab_limits_ss_new <- function(BU, water_type){
  # TODO: warnings
  # TODO: SHELL limits

  if (water_type == "marine"){
    if(BU == "REC-1"){
      limits <- c(NA, NA, NA, 110)
    }
  }

  if (water_type == "fresh"){
    if (BU == "REC-1"){
      limits <- c(320, NA, NA, NA)
    } else if (BU == "LREC-1"){
      limits <- c(576, NA, NA, NA)
    } else if (BU == "REC-2"){
      limits <- c(NA, 4000, NA, NA)
    }
  }

  names(limits) <- c("ecoli_WQO", "fc_WQO", "tc_WQO", "ent_WQO")
  return(limits)
}

bact_check_new <- function(df, sites, BU, water_type, ...){

  df <- df %>%
    average_results_daily()

  df <- tidy_bacteria(df)
  df <- replace_nd(df)

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
      check_geolimits_new(BU = analysis_sites$BU[[i]], water_type = analysis_sites$water_type[[i]], ...) %>%
      check_sslimits_new(BU = analysis_sites$BU[[i]], water_type = analysis_sites$water_type[[i]])

    out[[i]] <- out[[i]] %>% mutate("fc_to_tc" = NA,
                                      "tc_WQO_ss_2" = NA,
                                      "exceed_tc_WQO_ss_2" = NA)

    out[[i]] <- exceed_ss(out[[i]])
    out[[i]] <- order_bacteria_columns(out[[i]])

    #results[[i]] <- convertWeather(out[[i]])

    #results[[i]] <- results[[i]] %>% filter(Data_Row == TRUE)

    #results[[i]] <- results[[i]] %>%
    #  group_by(StationCode, WeatherCondition)

    #results[[i]] %>% summarize(exceedances = sum(exceed_day, na.rm = TRUE), n = n())

  }
  out
}

bact_ann_exceeds_new <- function(df, sites, BU, water_type, ...){

  analysis_sites <- data.frame(sites, BU, water_type, stringsAsFactors = FALSE)
  names(analysis_sites) <- c("StationCode", "BU", "water_type")

  analysis_sites <- analysis_sites[(analysis_sites$StationCode %in% unique(df$StationCode)), ]

  out <- bact_check_new(df, sites, BU, water_type)
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

stv_check <- function(df, sites, BU, water_type){

  df <- df %>%
    filter(Data_Row %in% c(TRUE)) %>%
    mutate(year = year(SampleDate),
           month = month(SampleDate))

  by_month <- group_by(df, WeatherCondition, StationCode, year, month)

  dplyr::summarize(by_month,
            above_STV = sum(exceed_day, na.rm = TRUE),
            samples = n(),
            percent_above = above_STV/samples,
            exceed_STV_WQO = ifelse(percent_above > 0.1, TRUE, FALSE))
}
