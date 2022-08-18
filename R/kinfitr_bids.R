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
    dplyr::filter(!stringr::str_detect(path, "^derivatives/")) %>%
    dplyr::mutate(attr  = purrr::map(path, bids_filename_attributes)) %>%
    tidyr::unnest(cols=all_of("attr")) %>%
    dplyr::mutate(extension = ifelse( extension=="gz" &
                                      grepl(".nii.gz", path),
                                      "nii.gz", extension)) %>%
    dplyr::filter(!is.na(measurement))

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

  if(!("run" %in% colnames(attributes))) {
    attributes$run <- 1
  }

  if(!("rec" %in% colnames(attributes))) {
    attributes$rec <- "rec"
  }

  # Test
  # attributes <- bind_rows(attributes, attributes[1,])
  # attributes$sub[7] <- NA

  # Inheritance

  # Here I split the files off where information isn't provided to apply to
  # multiple files by inheritance

  attributes_complete <- attributes %>%
    dplyr::filter(extension=="json" || extension=="tsv") %>%
    dplyr::filter(!is.na(sub) &
                  !is.na(ses) &
                  !is.na(task) &
                  !is.na(acq) &
                  !is.na(run) &
                  !is.na(rec))

  # Remove completely empty columns
  attributes_complete <- attributes_complete[,
                            !colMeans(is.na(attributes_complete))==1]

  attributes_inherit <- attributes %>%
    dplyr::filter(extension=="json" || extension=="tsv") %>%
    dplyr::filter(is.na(sub) |
                    is.na(ses) |
                    is.na(task) |
                    is.na(acq) |
                    is.na(run) |
                    is.na(rec))


  # Note to future self: The idea with the inheritance here is to assign things
  # which have multiple matches to the corresponding fields. So, blood data,
  # for instance can be matched to multiple reconstruction types of the PET
  # data. To accomplish this, I remove measurement and extension from the
  # matching. However, I feel like the inheritance might have originally been
  # for something else. If problems occur, check here.

  if(nrow(attributes_inherit) > 0) {
    ## https://stackoverflow.com/questions/50483890/dplyr-join-na-match-to-any

    attributes_inherit <- suppressMessages(
      attributes_inherit %>%
        split(seq(nrow(.))) %>%
        purrr::map_dfr(
          ~purrr::modify_if(.,is.na,~NULL) %>%
              dplyr::inner_join(.,
                dplyr::select(attributes_complete,
                                -path_absolute,
                                -path,
                                -extension,
                                -measurement)))) %>%
      dplyr::distinct()
  }

  attributes <- dplyr::bind_rows(attributes_complete,
                                   attributes_inherit)





  measurements <- attributes %>%
    dplyr::group_by(sub, ses, task, acq, run, rec) %>%
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

  attr <- stringr::str_match_all(filename, "([a-z0-9]*-[a-zA-Z0-9]*)[/_]")[[1]]
  attr_val <- stringr::str_match(attr[,2], "([a-z0-9]*)-([a-zA-Z0-9]*)")
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
                    recording=="manual") %>%
    dplyr::pull(path_absolute)

  tsv_blood_discrete <- filedata %>%
    dplyr::filter(measurement=="blood" &
                    extension=="tsv" &
                    recording=="manual") %>%
    dplyr::pull(path_absolute)

  json_blood_cont <- filedata %>%
    dplyr::filter(measurement=="blood" &
                    extension=="json" &
                    recording=="autosampler") %>%
    dplyr::pull(path_absolute)

  tsv_blood_cont <- filedata %>%
    dplyr::filter(measurement=="blood" &
                    extension=="tsv" &
                    recording=="autosampler") %>%
    dplyr::pull(path_absolute)





  ### Get the data ###

  jsondat_blood_discrete <- jsonlite::fromJSON(json_blood_discrete)

  # Read data

  ## Discrete
  if( jsondat_blood_discrete$PlasmaAvail | jsondat_blood_discrete$WholeBloodAvail ) {

    # Checks
    if( length(json_blood_discrete) == 0 ) {
      stop("No manual blood JSON found")
    }
    if( length(tsv_blood_discrete) == 0 ) {
      stop("No manual blood TSV file found")
    }

    tsvdat_blood_discrete  <- read.delim(tsv_blood_discrete, sep = "\t") %>%
      dplyr::filter(!is.na(time)) %>%
      dplyr::mutate(dplyr::across(everything(), ~ifelse(.x=="n/a",
                                                      yes = NA,
                                                      no = .x))) %>%
      dplyr::mutate(dplyr::across(everything(), as.numeric)) %>%
      dplyr::as_tibble()
  }

  ## Continuous
  if( length(json_blood_cont) > 0 ) {

    # Checks
    if( length(json_blood_cont) == 0 ) {
      stop("No autosampler blood JSON found")
    }
    if( length(tsv_blood_cont) == 0 ) {
      stop("No autosampler blood TSV file found")
    }

    jsondat_blood_cont <- jsonlite::fromJSON(json_blood_cont)
    tsvdat_blood_cont  <- read.delim(tsv_blood_cont, sep = "\t") %>%
      dplyr::filter(!is.na(time))
  }

  ## Get Metabolite
  if( jsondat_blood_discrete$MetaboliteAvail ) {

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
  if( jsondat_blood_discrete$PlasmaAvail ) {

    plasma <- dplyr::select(tsvdat_blood_discrete, time,
                            activity = plasma_radioactivity) %>%
      dplyr::filter(!is.na(activity))

    plasma_desc <- jsondat_blood_discrete[c("time", "plasma_radioactivity")]
    names(plasma_desc)[2] <- "activity"
  }

  ## Get Discrete Whole Blood
  if( jsondat_blood_discrete$WholeBloodAvail ) {

    blood_discrete <- dplyr::select(tsvdat_blood_discrete, time,
                            activity = whole_blood_radioactivity) %>%
      dplyr::filter(!is.na(activity))

    blood_discrete_desc <- jsondat_blood_discrete[c("time",
                                                  "whole_blood_radioactivity")]
    names(blood_discrete_desc)[2] <- "activity"
  }

  ## Get Continuous Whole Blood
  if( length(json_blood_cont) > 0 ) {

    blood_cont <- dplyr::select(tsvdat_blood_cont, time,
                                    activity = whole_blood_radioactivity) %>%
      dplyr::filter(!is.na(activity))

    blood_cont_desc <- jsondat_blood_cont[c("time","whole_blood_radioactivity")]
    names(blood_cont_desc)[2] <- "activity"
  }

  ## Edge cases

  ### No plasma, but only whole blood: use blood instead of plasma
  if( !jsondat_blood_discrete$PlasmaAvail &
      jsondat_blood_discrete$WholeBloodAvail ) {

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
  if( !jsondat_blood_discrete$MetaboliteAvail &
      (jsondat_blood_discrete$WholeBloodAvail | jsondat_blood_discrete$PlasmaAvail) ) {

    pf <- dplyr::select(tsvdat_blood_discrete, time,
                            parentFraction = 1)

    pf_desc <- list(jsondat_blood_discrete[c("time")])
    pf_desc$parentFraction <- list(Description =
                            "All set to 1 because no metabolite data available",
                            Units = "arbitrary")
  }

  ### No whole blood, only plasma: use plasma as blood
  if( jsondat_blood_discrete$PlasmaAvail & !jsondat_blood_discrete$WholeBloodAvail ) {

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

  MetaboliteData <- jsondat_blood_discrete[grep("Metabolite", names(jsondat_blood_discrete))]
  names(MetaboliteData) <- gsub("Metabolite", "", names(MetaboliteData))
  MetaboliteData$Values <- pf
  MetaboliteData <- c(MetaboliteData, pf_desc)

  PlasmaData <- jsondat_blood_discrete[grep("Plasma", names(jsondat_blood_discrete))]
  names(PlasmaData) <- gsub("Plasma", "", names(PlasmaData))
  PlasmaData$Values <- plasma
  PlasmaData <- c(PlasmaData, plasma_desc)

  DBloodData <- jsondat_blood_discrete[grep("WholeBlood", names(jsondat_blood_discrete))]
  names(DBloodData) <- gsub("WholeBlood", "", names(DBloodData))
  DBloodData$Values <- blood_discrete
  DBloodData <- c(DBloodData, blood_discrete_desc)

  CBloodData <- json_blood_cont[grep("WholeBlood", names(json_blood_cont))]
  names(CBloodData) <- gsub("WholeBlood", "", names(CBloodData))
  if( length(json_blood_cont) > 0 ) {
    CBloodData$Values <- blood_cont
    CBloodData <- c(CBloodData, blood_cont_desc)
    CBloodData$Avail = TRUE
  } else {
    CBloodData <- as.list(CBloodData)
    CBloodData$Avail = FALSE
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

  # This check seems irrelevant: I've removed it
  # if( MetaboliteData$parentFraction$Units != "arbitrary" ) {
  #     stop(paste("Unrecognised parentFraction units for Metabolite:",
  #                MetaboliteData$parentFraction$Units))
  # }

  # Checking the scaling between 0 and 1, or percentages
  if( !all(MetaboliteData$Values$parentFraction >= 0 &
          MetaboliteData$Values$parentFraction <= 1) ) {

    # Not all between 0 and 1

    ### If they're (probably) scaled between 0 and 100, no error
    if( all(MetaboliteData$Values$parentFraction >= 0 &
            MetaboliteData$Values$parentFraction <= 100) &
        max(MetaboliteData$Values$parentFraction > 50) ) {

      MetaboliteData$Values$parentFraction <- MetaboliteData$Values$parentFraction / 100

      warning(paste0("It seems your parent fraction values are scaled between 0",
                    " and 100. They have been rescaled between 0 and 1 by diving",
                     " by 100. Please check whether this is correct."))
    } else {

      ### Error if they're not 0-1 or 0-100
      stop("Parent fraction values should be scaled between 0 and 1.")

    }
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

bids_parse_pettimes <- function(filedata) {

  if(!("pet" %in% filedata$measurement)) {
    return(NA)
  }

  ### Get the filenames ###

  json_pet <- filedata %>%
    dplyr::filter(measurement=="pet" & extension=="json") %>%
    dplyr::mutate(jsondat_pet = purrr::map(
      path_absolute, jsonlite::fromJSON
    ))

  ### Extract the data ###

  jsondat_pet <- purrr::flatten(json_pet$jsondat_pet)

  tacdata <- tibble::tibble(
    start = jsondat_pet$FrameTimesStart,
    dur = jsondat_pet$FrameDuration,
    time = start + 0.5*dur
  ) %>%
    dplyr::mutate_all(~./60)

  return(tacdata)

}

bids_parse_petinfo <- function(filedata) {

  if(!("pet" %in% filedata$measurement)) {
    return(NA)
  }

  ### Get the filenames ###

  json_pet <- filedata %>%
    dplyr::filter(measurement=="pet" & extension=="json") %>%
    dplyr::mutate(jsondat_pet = purrr::map(
      path_absolute, jsonlite::fromJSON
    ))

  ### Extract the data ###

  jsondat_pet <- purrr::flatten(json_pet$jsondat_pet)

  return(jsondat_pet)

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
#' studydata <- bids_parse_study(studypath)
#' }
bids_parse_study <- function(studypath) {

  measurements <- bids_parse_files(studypath) %>%
    dplyr::group_by(ses, sub, task, acq)


  measurements$petinfo <- purrr::map(measurements$filedata,
                                      bids_parse_petinfo)

  measurements$blooddata <- purrr::map(measurements$filedata,
                                       bids_create_blooddata)

  measurements$tactimes <- purrr::map(measurements$filedata,
                                     bids_parse_pettimes)

  measurements <- dplyr::filter(measurements, !is.na(tactimes))

  return(measurements)

}


# join_by_nonmissing <- function(tbl_complete, tbl_incomplete) {
#
#   tbl_incomp_sep <- tbl_incomplete %>%
#     dplyr::mutate(idx = 1:dplyr::n()) %>%
#     dplyr::group_by(idx) %>%
#     tidyr::nest(matchdata=-idx) %>%
#     dplyr::mutate(matchdata = purrr::map(matchdata,
#                                          ~Filter(function(x) !all(is.na(x)),
#                                                  .x)))
#
#   tbl_incomp_sep_joined <- suppressMessages( tbl_incomp_sep %>%
#     dplyr::mutate(matchdata = purrr::map(matchdata,
#                                          ~dplyr::left_join(.x, tbl_complete))) )
#
#   tbl_incomp_sep_joined %>%
#     dplyr::ungroup() %>%
#     dplyr::select(-idx) %>%
#     tidyr::unnest(matchdata)
# }

