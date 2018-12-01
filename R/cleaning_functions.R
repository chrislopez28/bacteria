#' Tidy Bacteria Data
#'
#' Tidies water quality data so that bacteria data is grouped together.
#'
#' @param df a data frame of monitoring data.
#' @return a data frame of monitoring data that selects indicator bacteria results and reorganizes
#' results by sampling date.
#' @examples
#' tidy_bacteria(smb_beachdata_15-16)
#' tidy_bacteria(CEDEN_export)
#' @section Warning:
#' For this code to work, the input data frame should in an analytical record format
#' (e.g. CEDEN) that includes columns for the following:
#'
#' StationCode, SampleDate, WeatherCode, AnalyteName, Result, ResQualCode, MDL, and RL
#'
#' Indiciator bacteria names in the "AnalyteName" field should be labeled as followed:
#'
#' "E. coli", "Enterococcus", "Coliform, Fecal", and "Coliform, Total"
#' @export

tidy_bacteria <- function(df){
  tidycols <- c("StationCode", "SampleDate", "WeatherCondition", "AnalyteName", "Result", "Samples", "ResQualCode", "MDL", "RL")
  indicators <- c("E. coli", "Coliform, Fecal", "Coliform, Total", "Enterococcus")

  if (col_check(df, tidycols) == FALSE) {
    stop("Missing columns. Data frame should contain the following variables: StationCode, SampleDate, WeatherCondition, AnalyteName, Result, Samples, ResQualCode, MDL, and RL")
  }

  df <- df %>%
    select(StationCode, SampleDate, WeatherCondition, AnalyteName, Result, Samples, ResQualCode, MDL, RL) %>%
    filter(AnalyteName %in% indicators)

  ecoli <- df %>% filter(AnalyteName == "E. coli")
  fecal_coli <- df %>% filter(AnalyteName == "Coliform, Fecal")
  total_coli <- df %>% filter(AnalyteName == "Coliform, Total")
  entero <- df %>% filter(AnalyteName == "Enterococcus")

  bact <- tibble::frame_data(
    ~variable, ~name, ~result, ~samples, ~qual, ~mdl, ~rl,
    ecoli, "E. coli", "ecoli",  "ecoli_n", "ecoli_qual", "ecoli_mdl", "ecoli_rl",
    fecal_coli, "Coliform, Fecal",  "fecal_coliform", "fc_n", "fc_qual", "fc_mdl", "fc_rl",
    total_coli, "Coliform, Total",  "total_coliform", "tc_n", "tc_qual", "tc_mdl", "tc_rl",
    entero, "Enterococcus", "enterococcus", "ent_n", "ent_qual", "ent_mdl", "ent_rl"
  )

  #Split Data Frames by Constituent
  for (i in seq_along(bact$variable)){
    bact$variable[[i]] <- bact$variable[[i]] %>%
      dplyr::filter(AnalyteName == bact$name[[i]]) %>%
      dplyr::mutate(!!bact$result[[i]] := Result,
                    !!bact$samples[[i]] := Samples,
                    !!bact$qual[[i]] := ResQualCode,
                    !!bact$mdl[[i]] := MDL,
                    !!bact$rl[[i]] := RL)

    bact$variable[[i]] <- bact$variable[[i]] %>%
      dplyr::select(-AnalyteName, -Result, -Samples, -ResQualCode, -MDL, -RL)
  }

  #Join Tables Together
  join_list = c("StationCode", "SampleDate", "WeatherCondition")

  df_join <- dplyr::full_join(bact$variable[[1]], bact$variable[[2]], by = join_list)
  df_join <- dplyr::full_join(df_join, bact$variable[[3]], by = join_list)
  df_join <- dplyr::full_join(df_join, bact$variable[[4]], by = join_list)

  return(df_join)
}


#' Replace Non Detects
#'
#' Replace non detects in the data to allow for the calculation of geometric means.
#'
#' @param df a data frame of monitoring data that has been tidied by tidy_bacteria
#' @param assume_mdl a logical
#' @param ecoli_mdl a numeric
#' @param fc_mdl a numeric
#' @param tc_mdl a numeric
#' @param ent_mdl a numeric
#' @return a data frame
#' @export

replace_nd <- function(df, assume_mdl = TRUE, ecoli_mdl = 10, fc_mdl = 10, tc_mdl = 10, ent_mdl = 10){
  # TODO: Check if df is in FORMAT B

  dt <- tibble::frame_data(
    ~result, ~qual, ~mdl, ~rl, ~assume,
    "ecoli", "ecoli_qual", "ecoli_mdl", "ecoli_rl", ecoli_mdl,
    "fecal_coliform", "fc_qual", "fc_mdl", "fc_rl", fc_mdl,
    "total_coliform", "tc_qual", "tc_mdl", "tc_rl", tc_mdl,
    "enterococcus", "ent_qual", "ent_mdl", "ent_rl", ent_mdl
  )

  for (i in seq_along(dt$result)){

    ### Replace blank results that have ND/DNQ/< qualifier with MDL or RL
    temp <- df %>%
      dplyr::mutate(counter := if_else(is.na(!!as.name(dt$result[[i]])) & !!as.name(dt$qual[[i]]) %in% c("ND", "DNQ", "<"),
                                       TRUE, FALSE))

    if (sum(temp$counter, na.rm = TRUE) > 0) {
      print(paste0("Replacing ", sum(temp$counter, na.rm = TRUE), " ", dt$result[[i]], " results with MDL or RL."))
    }
    df <- df %>%
      dplyr::mutate(!!as.name(dt$result[[i]]) := if_else(is.na(!!as.name(dt$result[[i]])) & !!as.name(dt$qual[[i]]) %in% c("ND", "DNQ", "<"),
                                                         !!as.name(dt$mdl[[i]]), !!as.name(dt$result[[i]]))) %>%
      dplyr::mutate(!!as.name(dt$result[[i]]) := if_else(is.na(!!as.name(dt$result[[i]])) & !!as.name(dt$qual[[i]]) %in% c("ND", "DNQ", "<"),
                                                         !!as.name(dt$rl[[i]]), !!as.name(dt$result[[i]])))

    if (assume_mdl == TRUE){

      ### Replace blank results that have ND/DNQ/< qualifier with assumed MDL
      temp2 <- df %>%
        dplyr::mutate(counter := if_else(is.na(!!as.name(dt$result[[i]])) & !!as.name(dt$qual[[i]]) %in% c("ND", "DNQ", "<"),
                                         TRUE, FALSE))
      if (sum(temp2$counter, na.rm = TRUE) > 0) {
        print(paste0("Replacing ", sum(temp2$counter, na.rm = TRUE), " ", dt$result[[i]], " results with assumed MDL, since no MDL or RL specified."))
      }
      df <- df %>%
        dplyr::mutate(!!as.name(dt$result[[i]]) := if_else(is.na(!!as.name(dt$result[[i]])) & !!as.name(dt$qual[[i]]) %in% c("ND", "DNQ", "<"),
                                                           dt$assume[[i]], !!as.name(dt$result[[i]])))

      ### Replace results < 1 with assumed MDL
      temp3 <- df %>%
        dplyr::mutate(counter := if_else((!!as.name(dt$result[[i]]) < 1),
                                         TRUE, FALSE))
      if (sum(temp3$counter, na.rm = TRUE) > 0) {
        print(paste0("Replacing ", sum(temp3$counter, na.rm = TRUE), " ", dt$result[[i]], " results with assumed MDL, since result < 1."))
      }
      df <- df %>%
        dplyr::mutate(!!as.name(dt$result[[i]]) := if_else((!!as.name(dt$result[[i]]) < 1),
                                                           dt$assume[[i]], !!as.name(dt$result[[i]])))
    }
  }
  return(df)
}
