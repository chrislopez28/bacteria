#' Grab Single Sample Limitations
#'
#' Retrieves single sample water quality objectives for indicator bacteria.
#'
#' @param BU a string of either "REC-1", "LREC-1", or "REC-2" representing the highest
#'     beneficial use of the water body
#' @param water_type a string of either "fresh" or "marine' representing the water type
#'     of the water body
#' @return a numeric vector with the single sample objectives for e. coli, fecal coliform,
#' total coliform, and enterococcus

grab_limits_ss <- function(BU, water_type){
  # TODO: warnings
  # TODO: SHELL limits

  if (water_type == "marine"){
    if(BU == "REC-1"){
      limits <- c(NA, 400, 10000, 104)
    }
  }

  if (water_type == "fresh"){
    if (BU == "REC-1"){
      limits <- c(235, NA, NA, NA)
    } else if (BU == "LREC-1"){
      limits <- c(576, NA, NA, NA)
    } else if (BU == "REC-2"){
      limits <- c(NA, 4000, NA, NA)
    }
  }

  names(limits) <- c("ecoli_WQO", "fc_WQO", "tc_WQO", "ent_WQO")
  return(limits)
}

#' Grab Geomteric Mean Limitations
#'
#' Retrieves geometric mean water quality objectives for indicator bacteria.
#'
#' @param BU a string of either "REC-1", "LREC-1", or "REC-2" representing the highest
#'     beneficial use of the water body
#' @param water_type a string of either "fresh" or "marine' representing the water type
#'     of the water body
#' @return a numeric vector with the geometric mean objectives for e. coli, fecal coliform,
#' total coliform, and enterococcus

grab_limits_geo <- function(BU, water_type){
  # TODO: warnings
  # TODO: SHELL limits for total coliform: 30-day median = 70;
  # >10% samples over 30-day period cannot exceed 230 (5-tube decimal dilution) or 330 (3-tube decimal dilution)

  if (water_type == "marine"){
    if(BU == "REC-1"){
      limits <- c(NA, 200, 1000, 35)
    }
  }

  if (water_type == "fresh"){
    if (BU == "REC-1"){
      limits <- c(126, NA, NA, NA)
    } else if (BU == "LREC-1"){
      limits <- c(126, NA, NA, NA)
    } else if (BU == "REC-2"){
      limits <- c(NA, 2000, NA, NA)
    }
  }

  names(limits) <- c("ecoli_WQO", "fc_WQO", "tc_WQO", "ent_WQO")
  return(limits)
}
