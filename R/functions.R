
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
