#' Import bloodstream input functions
#'
#' Imports input function data from bloodstream for a BIDS dataset.
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
