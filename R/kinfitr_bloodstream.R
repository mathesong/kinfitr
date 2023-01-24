#' Import bloodstream input functions
#'
#' Imports input function data from bloodstream for a BIDS dataset, and converts all the units to minutes and kBq.
#'
#' @param blstream_folder The specific bloodstream derivative folder (for the whole study) from which the data should be extracted
#'
#' @return The BIDS attributes and kinfitr input objects
#' @export
#'
#' @examples
#' \dontrun{
#' bloodstream_import_inputfunction("../OpenNeuro/ds004230/derivatives/bloodstream2022-11-24_id-9AAO/")
#' }
bloodstream_import_inputfunctions <- function(blstream_folder) {

  blstream_files <- bids_parse_files(blstream_folder)
  blstream_files <- tidyr::unnest(blstream_files,
                                  cols = filedata)
  blstream_files <- dplyr::filter(blstream_files,
                                  measurement == "inputfunction")

  # JSON Data
  blstream_json <- dplyr::filter(blstream_files, extension=="json")
  blstream_json$jsondata <- purrr::map(blstream_json$path_absolute,
                                       jsonlite::fromJSON)
  blstream_json <- dplyr::select(blstream_json, -path_absolute, -path, -extension)

  # TSV Data
  blstream_tsv <- dplyr::filter(blstream_files, extension=="tsv")
  blstream_tsv$tsvdata <- purrr::map(blstream_tsv$path_absolute,
                                       ~read.delim(.x, sep="\t"))
  blstream_tsv <- dplyr::select(blstream_tsv, -path_absolute, -path, -extension)

  # Combine
  suppressMessages(
    blstream_data <- dplyr::full_join(blstream_json, blstream_tsv)
  )

  # Unit Conversions

  for(i in 1:nrow(blstream_data)) {

    # time
    if("time" %in% colnames(blstream_data$tsvdata[[i]])) {
      # correct if units are seconds
      if(blstream_data$jsondata[[i]]$time$Units == "s") {
        blstream_data$tsvdata[[i]]$time <- blstream_data$tsvdata[[i]]$time / 60
      }
    }

    # whole_blood_radioactivity
    if("whole_blood_radioactivity" %in% colnames(blstream_data$tsvdata[[i]])) {

      blstream_data$tsvdata[[i]]$whole_blood_radioactivity <- unit_convert(
        blstream_data$tsvdata[[i]]$whole_blood_radioactivity,
        from_units = blstream_data$jsondata[[i]]$whole_blood_radioactivity$Units,
        to_units = "kBq"
      )

    }

    # plasma_radioactivity
    if("plasma_radioactivity" %in% colnames(blstream_data$tsvdata[[i]])) {

      blstream_data$tsvdata[[i]]$plasma_radioactivity <- unit_convert(
        blstream_data$tsvdata[[i]]$plasma_radioactivity,
        from_units = blstream_data$jsondata[[i]]$plasma_radioactivity$Units,
        to_units = "kBq"
      )

    }

    # AIF
    if("AIF" %in% colnames(blstream_data$tsvdata[[i]])) {

      blstream_data$tsvdata[[i]]$AIF <- unit_convert(
        blstream_data$tsvdata[[i]]$AIF,
        from_units = blstream_data$jsondata[[i]]$AIF$Units,
        to_units = "kBq"
      )

    }

  }

  blstream_data <- dplyr::mutate(blstream_data,
                         input = purrr::map(tsvdata,
                                     ~blood_interp(.x$time,
                                                   .x$whole_blood_radioactivity,
                                                   .x$time,
                                                   .x$plasma_radioactivity,
                                                   .x$time,
                                                   .x$metabolite_parent_fraction,
                                                   .x$time,
                                                   .x$AIF)))

  blstream_data <- dplyr::select(blstream_data, -jsondata, -tsvdata)

  return(blstream_data)

}



#' Import bloodstream AIF fitted parameters
#'
#' Imports bloodstream AIF fitted parameters for a BIDS dataset, and converts all the units to minutes and kBq.
#'
#' @param blstream_folder The specific bloodstream derivative folder (for the whole study) from which the data should be extracted
#'
#' @return The BIDS attributes and fitted AIF parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' bloodstream_import_aifpars("../OpenNeuro/ds004230/derivatives/bloodstream2022-11-24_id-9AAO/")
#' }
bloodstream_import_aifpars <- function(blstream_folder) {

  blstream_files <- bids_parse_files(blstream_folder)
  blstream_files <- tidyr::unnest(blstream_files,
                                  cols = filedata)
  blstream_files <- dplyr::filter(blstream_files,
                                  measurement == "config")

  # JSON Data
  blstream_json <- dplyr::filter(blstream_files, extension=="json")
  blstream_json$jsondata <- purrr::map(blstream_json$path_absolute,
                                       jsonlite::fromJSON)
  blstream_json <- dplyr::select(blstream_json, -path_absolute, -path, -extension)

  # AIF Parameters
  blstream_aif <- dplyr::mutate(blstream_json,
                                AIFmodel = purrr::map_chr(jsondata,
                                                          c("AIF", "Method")))

  blstream_aif <- dplyr::mutate(blstream_aif,
                                Unit_time = purrr::map_chr(jsondata,
                                                           c("AIF", "Fit",
                                                             "Units", "time")))

  blstream_aif <- dplyr::mutate(blstream_aif,
                                Unit_AIF = purrr::map_chr(jsondata,
                                                          c("AIF", "Fit",
                                                            "Units", "AIF")))

  extract_parameters <- function(data) { data$AIF$Fit$Parameters }
  possibly_extract_pars <- purrr::possibly(extract_parameters,
                                           otherwise = NA)

  blstream_aif <- dplyr::mutate(blstream_aif,
                                AIFpars = purrr::map(jsondata,
                                                     possibly_extract_pars))


  blstream_aif <- dplyr::select(blstream_aif, -measurement, -jsondata)

  blstream_aif$AIFpars <- purrr::pmap(
    list( blstream_aif$AIFmodel, blstream_aif$AIFpars,
          blstream_aif$Unit_time, blstream_aif$Unit_AIF ),
    bloodstream_parameter_fix)

  blstream_aif <- dplyr::select(blstream_aif, -AIFmodel, -Unit_time,
                          -Unit_AIF)

  # blstream_aif <- tidyr::unnest(blstream_aif, cols=c(AIFpars))

  return(blstream_aif)

}


bloodstream_parameter_fix <- function(AIFmodel, AIFpars, Unit_time, Unit_AIF) {

  triexp <- AIFmodel == "Fit Individually: Linear Rise, Triexponential Decay"
  feng <- AIFmodel == "Fit Individually: Feng"
  fengconv <- AIFmodel == "Fit Individually: FengConv"


  ###### Tri-exponential ######
  if(triexp) {

    if(Unit_time == "s") {
      AIFpars$alpha <- AIFpars$alpha*60
      AIFpars$beta  <- AIFpars$beta *60
      AIFpars$gamma <- AIFpars$gamma*60

      AIFpars$t0       <- AIFpars$t0/60
      AIFpars$peaktime <- AIFpars$peaktime/60
    }

    if(Unit_AIF != "kBq") {
      conversion <- unit_convert(1,
                                 from_units=Unit_AIF,
                                 to_units = "kBq")
      AIFpars$A <- AIFpars$A * conversion
      AIFpars$B <- AIFpars$B * conversion
      AIFpars$C <- AIFpars$C * conversion

      AIFpars$peakval <- AIFpars$peakval * conversion
    }
  }

  ###### Feng ######
  if(feng) {

    if(Unit_time == "s") {
      AIFpars$alpha    <- AIFpars$alpha*60
      AIFpars$beta     <- AIFpars$beta *60
      AIFpars$gamma    <- AIFpars$gamma*60
      AIFpars$A        <- AIFpars$A*60

      AIFpars$t0       <- AIFpars$t0/60
    }

    if(Unit_AIF != "kBq") {
      conversion <- unit_convert(1,
                                 from_units=Unit_AIF,
                                 to_units = "kBq")
      AIFpars$A <- AIFpars$A * conversion
      AIFpars$B <- AIFpars$B * conversion
      AIFpars$C <- AIFpars$C * conversion
    }
  }

  ###### FengConv ######
  if(fengconv) {

    if(Unit_time == "s") {
      AIFpars$alpha    <- AIFpars$alpha*60
      AIFpars$beta     <- AIFpars$beta *60
      AIFpars$gamma    <- AIFpars$gamma*60
      AIFpars$A        <- AIFpars$A*60

      AIFpars$t0       <- AIFpars$t0/60
      AIFpars$ti       <- AIFpars$ti/60
    }

    if(Unit_AIF != "kBq") {
      conversion <- unit_convert(1,
                                 from_units=Unit_AIF,
                                 to_units = "kBq")
      AIFpars$A <- AIFpars$A * conversion
      AIFpars$B <- AIFpars$B * conversion
      AIFpars$C <- AIFpars$C * conversion
    }
  }

  return(AIFpars)

}
