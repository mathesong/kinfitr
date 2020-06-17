#' Extract the filenames from a PET BIDS study
#'
#' This function returns a data frame of the files for each measurement in a
#' PET BIDS study, its path, and what it is.
#'
#' @param studypath The BIDS study path for the main study.
#'
#' @return A data frame of the files, what they are, and where they are
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' filedata <- bids_parse_files(studypath)
#' }
bids_parse_files <- function(studypath) {

  extensions= paste(c('*.nii.gz', "*.tsv", "*.json"), collapse="|")

  files <- fs::dir_info(studypath, recurse = T, type = "file",
                        glob = extensions)

  files$path_relative <- fs::path_rel(files$path, studypath)

  attributes <- dplyr::select(files, path_absolute=path, path=path_relative) %>%
    dplyr::mutate(extension  = fs::path_ext(path)) %>%
    dplyr::mutate(attr  = purrr::map(path, bids_filename_attributes)) %>%
    tidyr::unnest(cols=all_of("attr")) %>%
    dplyr::mutate(extension = ifelse( extension=="gz" &
                                      grepl(".nii.gz", path),
                                      "nii.gz", extension))

  if(!("sub" %in% colnames(attributes))) {
    attributes$sub <- "01"
  }

  if(!("ses" %in% colnames(attributes))) {
    attributes$ses <- "01"
  }

  if(!("task" %in% colnames(attributes))) {
    attributes$task <- "rest"
  }

  if(!("acq" %in% colnames(attributes))) {
    attributes$acq <- "acq"
  }

  measurements <- attributes %>%
    dplyr::group_by(sub, ses, task, acq) %>%
    tidyr::nest() %>%
    dplyr::rename(filedata = data) %>%
    dplyr::ungroup()

  return(measurements)




}

#' Extract BIDS attributes from filenames and file paths
#'
#' Get all the stuff into a table
#'
#' @param filename The filename and filepath of the relevant file. Recommended
#' to provide the relative path within the studypath.
#'
#' @return A tibble containing all the data
#' @export
#'
#' @examples
#' bids_filename_attributes(
#'    "sub-01/ses-01/pet/sub-01_ses-01_recording-continuous_blood.json")
bids_filename_attributes <- function(filename) {

  attr <- stringr::str_match_all(filename, "([a-z0-9]*-[a-z0-9]*)[/_]")[[1]]
  attr_val <- stringr::str_match(attr[,2], "([a-z0-9]*)-([a-z0-9]*)")
  attr_vals <- tibble::tibble(
    attribute = attr_val[,2],
    value = attr_val[,3]
  )

  attr_vals <- attr_vals[!duplicated(
    attr_vals$attribute, fromLast = T),]

  attr_vals$measurement <- stringr::str_match(filename, "\\_([a-z]*)\\.")[,2]

  tidyr::spread(attr_vals, attribute, value)
}


#' Extract blood data from BIDS study folder
#'
#' Extracts the Data section of the blooddata object
#'
#' @param filedata A table of the file data, created using bids_parse_files()
#'
#' @return Data from these files: the Data section of the blooddata object
#' @export
#'
#' @examples
#' \dontrun{
#' bd_dat <- bids_parse_blood(filedata)
#' }
bids_parse_blood <- function(filedata) {

  if(!("blood" %in% filedata$measurement)) {
    return(NA)
  }

  ### Get the filenames ###

  json_pet <- filedata %>%
    dplyr::filter(measurement=="pet" & extension=="json") %>%
    dplyr::pull(path_absolute)

  json_blood_discrete <- filedata %>%
    dplyr::filter(measurement=="blood" &
                    extension=="json" &
                    recording=="discrete") %>%
    dplyr::pull(path_absolute)

  tsv_blood_discrete <- filedata %>%
    dplyr::filter(measurement=="blood" &
                    extension=="tsv" &
                    recording=="discrete") %>%
    dplyr::pull(path_absolute)

  json_blood_cont <- filedata %>%
    dplyr::filter(measurement=="blood" &
                    extension=="json" &
                    recording=="continuous") %>%
    dplyr::pull(path_absolute)

  tsv_blood_cont <- filedata %>%
    dplyr::filter(measurement=="blood" &
                    extension=="tsv" &
                    recording=="continuous") %>%
    dplyr::pull(path_absolute)

  json_blood_discrete <- filedata %>%
    dplyr::filter(measurement=="blood" &
                    extension=="json" &
                    recording=="discrete") %>%
    dplyr::pull(path_absolute)





  ### Get the data ###

  jsondat_pet <- jsonlite::fromJSON(json_pet)

  # Read data

  ## Discrete
  if( jsondat_pet$PlasmaAvail | jsondat_pet$DiscreteBloodAvail ) {

    # Checks
    if( length(json_blood_discrete) == 0 ) {
      stop("No discrete blood JSON found")
    }
    if( length(tsv_blood_discrete) == 0 ) {
      stop("No discrete blood TSV file found")
    }

    jsondat_blood_discrete <- jsonlite::fromJSON(json_blood_discrete)
    tsvdat_blood_discrete  <- read.delim(tsv_blood_discrete, sep = "\t") %>%
      dplyr::filter(!is.na(time))
  }

  ## Continuous
  if( jsondat_pet$ContinuousBloodAvail ) {

    # Checks
    if( length(json_blood_cont) == 0 ) {
      stop("No continuous blood JSON found")
    }
    if( length(tsv_blood_cont) == 0 ) {
      stop("No continuous blood TSV file found")
    }

    jsondat_blood_cont <- jsonlite::fromJSON(json_blood_cont)
    tsvdat_blood_cont  <- read.delim(tsv_blood_cont, sep = "\t") %>%
      dplyr::filter(!is.na(time))
  }

  ## Get Metabolite
  if( jsondat_pet$MetaboliteAvail ) {

    pf <- dplyr::select(tsvdat_blood_discrete, time,
                           starts_with("metabolite_")) %>%
      dplyr::filter(!is.na(metabolite_parent_fraction))

    pf_desc <- jsondat_blood_discrete[names(pf)]

    names(pf) <- gsub("metabolite_", "", names(pf))
    names(pf) <- gsub('\\_(\\w?)', '\\U\\1', names(pf), perl = TRUE)

    names(pf_desc) <- gsub("metabolite_", "", names(pf_desc))
    names(pf_desc) <- gsub('\\_(\\w?)', '\\U\\1', names(pf_desc), perl = TRUE)
  }

  ## Get Plasma
  if( jsondat_pet$PlasmaAvail ) {

    plasma <- dplyr::select(tsvdat_blood_discrete, time,
                            activity = plasma_radioactivity) %>%
      dplyr::filter(!is.na(activity))

    plasma_desc <- jsondat_blood_discrete[c("time", "plasma_radioactivity")]
    names(plasma_desc)[2] <- "activity"
  }

  ## Get Discrete Whole Blood
  if( jsondat_pet$DiscreteBloodAvail ) {

    blood_discrete <- dplyr::select(tsvdat_blood_discrete, time,
                            activity = whole_blood_radioactivity) %>%
      dplyr::filter(!is.na(activity))

    blood_discrete_desc <- jsondat_blood_discrete[c("time",
                                                  "whole_blood_radioactivity")]
    names(blood_discrete_desc)[2] <- "activity"
  }

  ## Get Continuous Whole Blood
  if( jsondat_pet$ContinuousBloodAvail ) {

    blood_cont <- dplyr::select(tsvdat_blood_cont, time,
                                    activity = whole_blood_radioactivity) %>%
      dplyr::filter(!is.na(activity))

    blood_cont_desc <- jsondat_blood_cont[c("time","whole_blood_radioactivity")]
    names(blood_cont_desc)[2] <- "activity"
  }

  ## Edge cases

  ### No plasma, but only whole blood: use blood instead of plasma
  if( !jsondat_pet$PlasmaAvail & jsondat_pet$DiscreteBloodAvail ) {

    plasma <- dplyr::select(tsvdat_blood_discrete, time,
                            activity = whole_blood_radioactivity) %>%
      dplyr::filter(!is.na(activity))

    plasma_desc <- jsondat_blood_discrete[c("time","whole_blood_radioactivity")]
    names(plasma_desc)[2] <- "activity"
    plasma_desc$activity$Description <- paste(
      "Whole blood used as no plasma available.",
      plasma_desc$activity$Description)
  }

  ### No metabolite, but blood/plasma: metab=1
  if( !jsondat_pet$MetaboliteAvail &
      (jsondat_pet$DiscreteBloodAvail | jsondat_pet$PlasmaAvail) ) {

    pf <- dplyr::select(tsvdat_blood_discrete, time,
                            parentFraction = 1)

    pf_desc <- list(jsondat_blood_discrete[c("time")])
    pf_desc$parentFraction <- list(Description =
                            "All set to 1 because no metabolite data available",
                            Units = "unitless")
  }

  ### No whole blood, only plasma: use plasma as blood
  if( jsondat_pet$PlasmaAvail & !jsondat_pet$DiscreteBloodAvail ) {

    blood_discrete <- dplyr::select(tsvdat_blood_discrete, time,
                            activity = plasma_radioactivity) %>%
      dplyr::filter(!is.na(activity))

    blood_discrete_desc <- jsondat_blood_discrete[c("time",
                                                    "plasma_radioactivity")]
    names(blood_discrete_desc)[2] <- "activity"
    blood_discrete_desc$activity$Description <- paste(
      "Plasma used as no whole blood available.",
      blood_discrete_desc$activity$Description)


  }




  ### Arrange the Data ###

  MetaboliteData <- jsondat_pet[grep("Metabolite", names(jsondat_pet))]
  names(MetaboliteData) <- gsub("Metabolite", "", names(MetaboliteData))
  MetaboliteData$Values <- pf
  MetaboliteData <- c(MetaboliteData, pf_desc)

  PlasmaData <- jsondat_pet[grep("Plasma", names(jsondat_pet))]
  names(PlasmaData) <- gsub("Plasma", "", names(PlasmaData))
  PlasmaData$Values <- plasma
  PlasmaData <- c(PlasmaData, plasma_desc)

  DBloodData <- jsondat_pet[grep("DiscreteBlood", names(jsondat_pet))]
  names(DBloodData) <- gsub("DiscreteBlood", "", names(DBloodData))
  DBloodData$Values <- blood_discrete
  DBloodData <- c(DBloodData, blood_discrete_desc)

  CBloodData <- jsondat_pet[grep("ContinuousBlood", names(jsondat_pet))]
  names(CBloodData) <- gsub("ContinuousBlood", "", names(CBloodData))
  if( CBloodData$Avail ) {
    CBloodData$Values <- blood_cont
    CBloodData <- c(CBloodData, blood_cont_desc)
  }


  ### Unit conversions ###

  # Metabolite

  if( MetaboliteData$time$Units != "s" ) {
    if( MetaboliteData$time$Units == "min" ) {
      MetaboliteData$Values$time <- MetaboliteData$Values$time * 60
      MetaboliteData$time$Units <- "s"
    } else {
      stop(paste("Unrecognised time units for metabolite:",
                 MetaboliteData$time$Units))
    }
  }

  if( MetaboliteData$parentFraction$Units != "unitless" ) {
      stop(paste("Unrecognised parentFraction units for Metabolite:",
                 MetaboliteData$parentFraction$Units))
  }


  # Plasma

  if( PlasmaData$time$Units != "s" ) {
    if( PlasmaData$time$Units == "min" ) {
      PlasmaData$Values$time <- PlasmaData$Values$time * 60
      PlasmaData$time$Units <- "s"
    } else {
      stop(paste("Unrecognised time units for Plasma:",
                 PlasmaData$time$Units))
    }
  }

  plasmarad <- get_units_radioactivity(PlasmaData$activity$Units)

  if( plasmarad$rad != "kBq" ) {
    PlasmaData$Values$activity <- unit_convert(PlasmaData$Values$activity,
                                               plasmarad$rad, "kBq")
  }

  if( !(plasmarad$vol %in% c("ml", "cc", "mL") )) {
    stop(paste("Unrecognised activity units for Plasma:",
               PlasmaData$activity$Units,
               ": units should be in kBq/ml"))
  }

  PlasmaData$activity$Units <- "kBq/ml"



  # Blood Discrete

  if( DBloodData$time$Units != "s" ) {
    if( DBloodData$time$Units == "min" ) {
      DBloodData$Values$time <- DBloodData$Values$time * 60
      DBloodData$time$Units <- "s"
    } else {
      stop(paste("Unrecognised time units for BloodDiscrete:",
                 DBloodData$time$Units))
    }
  }

  dbloodrad <- get_units_radioactivity(DBloodData$activity$Units)

  if( dbloodrad$rad != "kBq" ) {
    DBloodData$Values$activity <- unit_convert(DBloodData$Values$activity,
                                               dbloodrad$rad, "kBq")
  }

  if( !(dbloodrad$vol %in% c("ml", "cc", "mL") )) {
    stop(paste("Unrecognised activity units for BloodDiscrete:",
               DBloodData$activity$Units,
               ": units should be in kBq/ml"))
  }

  DBloodData$activity$Units <- "kBq/ml"


  # Blood Continuous

  if( CBloodData$Avail ) {


    if( CBloodData$time$Units != "s" ) {
      if( CBloodData$time$Units == "min" ) {
        CBloodData$Values$time <- CBloodData$Values$time * 60
        CBloodData$time$Units <- "s"
      } else {
        stop(paste("Unrecognised time units for BloodContinuous:",
                   CBloodData$time$Units))
      }
    }

    cbloodrad <- get_units_radioactivity(CBloodData$activity$Units)

    if( cbloodrad$rad != "kBq" ) {
      CBloodData$Values$activity <- unit_convert(CBloodData$Values$activity,
                                                 cbloodrad$rad, "kBq")
    }

    if( !(cbloodrad$vol %in% c("ml", "cc", "mL") )) {
      stop(paste("Unrecognised activity units for BloodContinuous:",
                 CBloodData$activity$Units,
                 ": units should be in kBq/ml"))
    }

    CBloodData$activity$Units <- "kBq/ml"

  }




  ### Output ###

  bids_data <- list(
    Blood = list(
      Discrete = DBloodData,
      Continuous = CBloodData),
    Plasma = PlasmaData,
    Metabolite = MetaboliteData
  )

  return(bids_data)

}



#' Create a blooddata object from BIDS data
#'
#' This function takes a set of files from one measurement of a BIDS-structured
#' study, reads in the data, and creates a blooddata file.
#'
#' @param filedata A table of the file data, created using bids_parse_files()
#'
#' @return A blooddata object
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' bd <- bids_create_blooddata(filedata)
#' }
bids_create_blooddata <- function(filedata) {

  bids_data <- bids_parse_blood(filedata)

  blooddata <- list(
    Data = bids_data,
    Models = list(
      Blood = list(Method = "interp", Data = NULL),
      BPR = list(Method = "interp", Data = NULL),
      parentFraction = list(Method = "interp", Data = NULL),
      AIF = list(Method = "interp", Data = NULL)
    ),
    TimeShift = 0
  )

  blooddata <- rapply(blooddata, function(x) ifelse(x == "true", TRUE, x),
                      how = "replace")
  blooddata <- rapply(blooddata, function(x) ifelse(x == "false", FALSE, x),
                      how = "replace")

  class(blooddata) <- "blooddata"

  return(blooddata)

}

bids_parse_pet <- function(filedata) {

  if(!("pet" %in% filedata$measurement)) {
    return(NA)
  }

  ### Get the filenames ###

  json_pet <- filedata %>%
    dplyr::filter(measurement=="pet" & extension=="json") %>%
    dplyr::pull(path_absolute)

  ### Get the data ###

  jsondat_pet <- jsonlite::fromJSON(json_pet)

  tacdata <- tibble::tibble(
    start = jsondat_pet$FrameTimesStart,
    dur = jsondat_pet$FrameDuration,
    time = start + 0.5*dur
  ) %>%
    dplyr::mutate_all(~./60)

  return(tacdata)

}

#' Parse the contents of a PET BIDS study
#'
#' This function parses a PET BIDS study, and returns a nested tibble with the
#' relevant information nested appropriately.
#'
#' @param studypath The BIDS study path for the main study.
#'
#' @return A data frame of the files, what they are, and where they are
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' filedata <- bids_parse_files(studypath)
#' }
bids_parse_study <- function(studypath) {

  measurements <- bids_parse_files(studypath) %>%
    dplyr::group_by(ses, sub, task, acq)


  measurements$blooddata <- purrr::map(measurements$filedata,
                                       bids_create_blooddata)

  measurements$tactimes <- purrr::map(measurements$filedata,
                                     bids_parse_pet)

  measurements <- dplyr::filter(measurements, !is.na(tactimes))

  return(measurements)

}
