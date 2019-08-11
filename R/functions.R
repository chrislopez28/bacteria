
test_bact_check <- function(df, BU, water_type) {
  df <- df %>%
    expand_dates() %>%
    bact_geomeans() %>%
    check_geolimits(BU, water_type, ...) %>%
    check_sslimits(BU, water_type, ...)

  # If water_type = "marine" --> Check fecal-total ratio; Add limits if necessary.

  if (water_type == "marine"){
    df <- df %>% check_fecal_to_total(...)
  } else if (analysis_sites$water_type == "fresh"){
    df$fc_to_tc <- NA
    df$tc_WQO_ss_2 <- NA
    df$exceed_tc_WQO_ss_2 <- NA
  }

  df<- df %>%
    exceed_ss() %>%
    order_bacteria_columns()

  df
}


test_check_geolimits <- function(df, BU = "REC-1", water_type = "marine", six_week = TRUE, ...){
  # TODO: Check if df has appropriate columns / Use colcheck()
  # TODO: Check if df has consecutive SampleDates from first row to last row / Create consecutivecheck()

  limits <- grab_limits_geo(BU, water_type)

  dt <- data.frame(
    geomean = c("ecoli_geomean", "fc_geomean", "tc_geomean", "ent_geomean"),
    geocount = c("ecoli_geo_count", "fc_geo_count", "tc_geo_count", "ent_geo_count"),
    WQO_geo = c("ecoli_WQO_geo", "fc_WQO_geo", "tc_WQO_geo", "ent_WQO_geo"),
    WQO_geo_val = c(limits["ecoli_WQ"], limits["fc_WQO"], limits["tc_WQO"], limits["ent_WQO"]),
    exceed_WQO_geo = c("exceed_ecoli_WQO_geo", "exceed_fc_WQO_geo", "exceed_tc_WQO_geo",
                        "exceed_ent_WQO_geo")
  )

  # dt <- tibble::frame_data(
  #   ~geomean, ~geocount,  ~WQO_geo, ~WQO_geo_val, ~exceed_WQO_geo,
  #   "ecoli_geomean", "ecoli_geo_count", "ecoli_WQO_geo", limits["ecoli_WQO"], "exceed_ecoli_WQO_geo",
  #   "fc_geomean", "fc_geo_count", "fc_WQO_geo", limits["fc_WQO"], "exceed_fc_WQO_geo",
  #   "tc_geomean", "tc_geo_count", "tc_WQO_geo", limits["tc_WQO"], "exceed_tc_WQO_geo",
  #   "ent_geomean", "ent_geo_count", "ent_WQO_geo", limits["ent_WQO"], "exceed_ent_WQO_geo"
  # )

  if (six_week){
    for (i in seq_along(dt$geomean)){
      df <- df %>%
        dplyr::mutate(!!as.name(dt$WQO_geo[[i]]) :=
                        dplyr::if_else(lubridate::wday(SampleDate, label = TRUE) == "Sun",
                                       dt$WQO_geo_val[[i]], as.double(NA))) %>%
        dplyr::mutate(!!as.name(dt$exceed_WQO_geo[[i]]) :=
                        dplyr::if_else(!!as.name(dt$geomean[[i]]) > !!as.name(dt$WQO_geo[[i]]) &
                                         !!as.name(dt$geocount[[i]]) >= 5, TRUE, FALSE))
    }
  } else {

    for (i in seq_along(dt$geomean)){
      df <- df %>%
        dplyr::mutate(!!as.name(dt$WQO_geo[[i]]) :=
                        dplyr::if_else(Data_Row == TRUE, dt$WQO_geo_val[[i]], as.double(NA))) %>%
        dplyr::mutate(!!as.name(dt$exceed_WQO_geo[[i]]) :=
                        dplyr::if_else(!!as.name(dt$geomean[[i]]) > !!as.name(dt$WQO_geo[[i]]) &
                                         !!as.name(dt$geocount[[i]]) >= 5, TRUE, FALSE))
    }

  }
  return(df)
}
