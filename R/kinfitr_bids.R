# Warn about entity labels that differ only in case.
#
# BIDS labels are case-sensitive, so trc-PF974 and trc-pf974 are two different
# tracers as far as any tool is concerned -- which is almost never what was
# meant, and produces confusing duplicate measurements.
#
# Warn, never error: a collision simply yields two distinct measurements, so
# refusing the dataset would only turn away work that can be done correctly.
bids_warn_case_collisions <- function(attributes) {

  entity_cols <- setdiff(colnames(attributes),
                         c("path", "path_absolute", "extension", "measurement"))

  for (entity in entity_cols) {

    values <- unique(attributes[[entity]])
    values <- values[!is.na(values)]
    if (length(values) < 2) next

    collided <- split(values, tolower(values))
    collided <- collided[vapply(collided, length, integer(1)) > 1]
    if (length(collided) == 0) next

    groups <- vapply(collided, function(v) {
      paste(paste0(entity, "-", v), collapse = " / ")
    }, character(1))

    warning("Entity labels differing only in case: ",
            paste(groups, collapse = "; "),
            ". BIDS labels are case-sensitive, so these are treated as ",
            "distinct measurements. If that was not intended, run the BIDS ",
            "validator (https://bids-standard.github.io/bids-validator/) and ",
            "correct the filenames.",
            call. = FALSE)
  }

  invisible(attributes)
}

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

  extensions = paste(c('*.nii.gz', "*.tsv", "*.json"), collapse="|")

  files <- fs::dir_info(studypath, recurse = T, type = "file",
                        glob = extensions)

  files$path_relative <- fs::path_rel(files$path, studypath)

  attributes <- dplyr::select(files, path_absolute=path, path=path_relative) %>%
    dplyr::mutate(extension  = fs::path_ext(path)) %>%
    dplyr::filter(!stringr::str_detect(path, "^derivatives/")) %>%
    dplyr::filter(!stringr::str_detect(path, "^code/")) %>%
    dplyr::filter(!stringr::str_detect(path, "^phenotype/")) %>%
    dplyr::filter(!stringr::str_detect(path, "^sourcedata/")) %>%
    dplyr::mutate(attr  = purrr::map(path, bids_filename_attributes)) %>%
    tidyr::unnest(cols=all_of("attr")) %>%
    dplyr::mutate(extension = ifelse( extension=="gz" &
                                      grepl(".nii.gz", path),
                                      "nii.gz", extension)) %>%
    dplyr::filter(!is.na(measurement))

  # Check before the defaults below are filled in, so that synthetic values
  # such as trc = "trc" cannot themselves collide.
  bids_warn_case_collisions(attributes)

  if(!("sub" %in% colnames(attributes))) {
    attributes$sub <- "01"
  }

  if(!("ses" %in% colnames(attributes))) {
    attributes$ses <- "01"
  }

  if(!("task" %in% colnames(attributes))) {
    attributes$task <- "rest"
  }

  if(!("trc" %in% colnames(attributes))) {
    attributes$trc <- "trc"
  }

  # This seems to have actually been removed: this should be deprecated!
  if(!("acq" %in% colnames(attributes))) {
    attributes$acq <- "acq"
  }

  if(!("run" %in% colnames(attributes))) {
    attributes$run <- "01"
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


  attributes_ultra <- tidyr::expand_grid(
    sub = unique(na.omit(attributes$sub)),
    ses = unique(na.omit(attributes$ses)),
    task = unique(na.omit(attributes$task)),
    trc = unique(na.omit(attributes$trc)),
    acq = unique(na.omit(attributes$acq)),
    run = unique(na.omit(attributes$run)),
    rec = unique(na.omit(attributes$rec))
  )

  attributes_completed <- suppressMessages(
    attributes %>%
      split(seq(nrow(.))) %>%
      purrr::map(., ~.x %>%
                   dplyr::select(dplyr::where(~!all(is.na(.x)))) %>%
                   dplyr::inner_join(attributes_ultra)) %>%
      dplyr::bind_rows() %>%
      dplyr::distinct())

  # attributes_complete <- attributes %>%
  #   dplyr::filter(extension=="json" | extension=="tsv" | extension=="nii.gz" ) %>%
  #   dplyr::filter(!is.na(sub) &
  #                 !is.na(ses) &
  #                 !is.na(task) &
  #                 !is.na(trc) &
  #                 !is.na(acq) &
  #                 !is.na(run) &
  #                 !is.na(rec)) %>%
  #   dplyr::select(dplyr::where(~!all(is.na(.x)))) # remove totally empty
  #
  # attributes_inherit <- attributes %>%
  #   dplyr::filter(extension=="json" | extension=="tsv" | extension=="nii.gz" ) %>%
  #   dplyr::filter(is.na(sub) |
  #                 is.na(ses) |
  #                 is.na(task) |
  #                 is.na(trc) |
  #                 is.na(acq) |
  #                 is.na(run) |
  #                 is.na(rec))
  #
  #
  # # Sometimes, the attributes are incompletely provided for the PET images,
  # # resulting in the PET files ending up in attributes_inherit. So this completes
  # # their attributes, and brings them back to attributes_complete
  #
  # if(nrow(attributes_inherit) > 0) {
  #
  #   attributes_inherit <- attributes_inherit %>%
  #     dplyr::mutate(attr_ic_id = 1:dplyr::n())  # create identifier column
  #
  #   attributes_incompletepet <- attributes_inherit %>%
  #     dplyr::filter(measurement=="pet") %>%
  #     dplyr::filter(extension=="json" | extension=="tsv" | extension=="nii.gz") %>%
  #     dplyr::mutate(
  #       ses = bids_complete_attribute(ses, "ses", "01"),
  #       task = bids_complete_attribute(task, "task", "rest"),
  #       trc = bids_complete_attribute(trc, "trc", "trc"),
  #       acq = bids_complete_attribute(acq, "acq", "acq"),
  #       run = bids_complete_attribute(run, "run", "01"),
  #       rec = bids_complete_attribute(rec, "rec", "rec")
  #     ) %>%
  #     dplyr::select(dplyr::where(~!all(is.na(.x)))) # remove totally empty
  #
  #   attributes_inherit <- suppressMessages(
  #     attributes_inherit %>%
  #       dplyr::anti_join(attributes_incompletepet %>%
  #                          dplyr::select(attr_ic_id))) %>%
  #     dplyr::select(dplyr::where(~!all(is.na(.x)))) %>%  # remove totally empty
  #     dplyr::select(-attr_ic_id)
  #
  #   attributes_incompletepet <- attributes_incompletepet %>%
  #     dplyr::select(-attr_ic_id)
  #
  #   attributes_complete <- attributes_complete %>%
  #     dplyr::bind_rows(attributes_incompletepet)
  #
  # }
  #
  # # Note to future self: The idea with the inheritance here is to assign things
  # # which have multiple matches to the corresponding fields. So, blood data,
  # # for instance can be matched to multiple reconstruction types of the PET
  # # data. To accomplish this, I remove measurement and extension from the
  # # matching. However, I feel like the inheritance might have originally been
  # # for something else. If problems occur, check here.
  #
  # if(nrow(attributes_inherit) > 0) {
  #   ## https://stackoverflow.com/questions/50483890/dplyr-join-na-match-to-any
  #
  #   attributes_inherit_completed <- suppressMessages(
  #     attributes_inherit %>%
  #       split(seq(nrow(.))) %>%
  #       purrr::map_dfr(
  #         ~purrr::modify_if(.,is.na,~NULL) %>%
  #             dplyr::inner_join(.,
  #               dplyr::select(attributes_complete,
  #                               -path_absolute,
  #                               -path,
  #                               -extension,
  #                               -measurement)))) %>%
  #     dplyr::distinct()
  #
  #   attributes <- dplyr::bind_rows(attributes_complete,
  #                                  attributes_inherit_completed)
  #
  # } else {
  #
  #   attributes <- attributes_complete
  #
  # }
  #






  measurements <- attributes_completed %>%
    dplyr::group_by(sub, ses, task, acq, trc, run, rec) %>%
    tidyr::nest() %>%
    dplyr::rename(filedata = data) %>%
    dplyr::ungroup()

  return(measurements)




}

# bids_complete_attribute <- function(attr_vec, attr_name, attr_replacement) {
#
#   if( sum( is.na(attr_vec)) > 0 ) {
#
#     warning(paste0(
#       "The ", attr_name, " attribute is incomplete, and has been automatically ",
#       "completed for some PET measurements missing this attribute. ",
#       "It has been replaced with \"", attr_replacement, "\". ",
#       "Please check that this does not conflict with anything."
#     ))
#
#     attr_vec[is.na(attr_vec)] <- attr_replacement
#   }
#
#   return(attr_vec)
#
# }


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

  # "+" is a legal character in a BIDS entity label. Without it here the entity
  # does not merely lose the label -- the whole key-value pair fails to match and
  # is dropped, so "task-A+B" parses as no task at all rather than as task "A".
  attr <- stringr::str_match_all(filename, "([a-z0-9]*-[a-zA-Z0-9+]*)[/_]")[[1]]
  attr_val <- stringr::str_match(attr[,2], "([a-z0-9]*)-([a-zA-Z0-9+]*)")
  attr_vals <- tibble::tibble(
    attribute = attr_val[,2],
    value = attr_val[,3]
  )

  attr_vals <- attr_vals[!duplicated(
    attr_vals$attribute, fromLast = T),]

  attr_vals$measurement <- stringr::str_match(filename, "\\_([a-zA-Z1-9]*)\\.")[,2]

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

  ### No metabolite, but blood/plasma: metab=1
  if( !jsondat_blood_discrete$MetaboliteAvail &
      (jsondat_blood_discrete$WholeBloodAvail | jsondat_blood_discrete$PlasmaAvail) ) {

    pf <- dplyr::select(tsvdat_blood_discrete, time) %>%
      dplyr::mutate(parentFraction = 1)

    pf_desc <- list(time = jsondat_blood_discrete$time)
    pf_desc$parentFraction <- list(Description =
                                     "All set to 1 because no metabolite data available",
                                   Units = "arbitrary")
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


  if( length(json_blood_cont) > 0 ) {
    CBloodData <- jsondat_blood_cont[grep("WholeBlood", names(jsondat_blood_cont))]
    names(CBloodData) <- gsub("WholeBlood", "", names(CBloodData))

    CBloodData$Values <- blood_cont
    CBloodData <- c(CBloodData, blood_cont_desc)
    CBloodData$Avail = TRUE
  } else {
    CBloodData <- list()
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
      warning(paste0("Parent fraction values should be scaled between 0 and 1.\n",
                     "Anomalies detected in ", basename(json_blood_discrete), "\n",
                     "Parent Fraction ranges from ",
                     min(MetaboliteData$Values$parentFraction)," to ",
                     max(MetaboliteData$Values$parentFraction)," .\n\n"))

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

  if(length(bids_data) == 1) {
    return(list(NA))
  }

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

# A short label naming a measurement, for use in messages. Prefers the PET
# file's relative path, which is what a user would go looking for.
bids_describe_measurement <- function(filedata) {

  paths <- filedata$path[filedata$measurement == "pet"]
  if (length(paths) == 0) {
    paths <- filedata$path
  }
  if (length(paths) == 0) {
    return("this measurement")
  }

  paths[1]
}

bids_parse_pettimes <- function(filedata) {

  if(!("pet" %in% filedata$measurement)) {
    return(NA)
  }

  ### Get the filenames ###

  json_pet <- filedata %>%
    dplyr::filter(measurement=="pet" & extension=="json")

  # Frame times come from the _pet.json sidecar. Report a measurement that has
  # none and return NA so it is dropped, rather than failing: this function is
  # mapped over every measurement in the study, so an error here used to abort
  # the whole parse -- one missing sidecar took down every other measurement
  # too, with the opaque message "object 'dur' not found".
  if (nrow(json_pet) == 0) {
    warning("No _pet.json sidecar found for ", bids_describe_measurement(filedata),
            ". Frame times cannot be determined, so this measurement is ",
            "excluded. Every PET image needs an accompanying _pet.json.",
            call. = FALSE)
    return(NA)
  }

  json_pet <- json_pet %>%
    dplyr::mutate(jsondat_pet = purrr::map(
      path_absolute, jsonlite::fromJSON
    ))

  ### Extract the data ###

  jsondat_pet <- purrr::flatten(json_pet$jsondat_pet)

  missing_fields <- setdiff(c("FrameTimesStart", "FrameDuration"),
                            names(jsondat_pet))
  if (length(missing_fields) > 0) {
    warning("The _pet.json for ", bids_describe_measurement(filedata),
            " is missing ", paste(missing_fields, collapse = " and "),
            ". Frame times cannot be determined, so this measurement is ",
            "excluded.", call. = FALSE)
    return(NA)
  }

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



# --- Derivative parsing -------------------------------------------------------
#
# A derivatives folder is not a normal BIDS folder: with no PET images, there are
# no acquisitions to enumerate, and the raw parser's inheritance and entity
# grid do not apply to it. Using bids_parse_files() on one produces fabricated
# entities -- a subject with no session directory comes back carrying ses-test.
# These functions do filename grouping and nothing else.

# The BIDS PET filename template, and thus the entities that identify which
# acquisition a derivative belongs to. Deliberately hard-coded rather than
# derived: an entity added here changes what counts as the same measurement.
bids_selector_entities <- function() {
  c("sub", "ses", "task", "trc", "rec", "run")
}

# Every "key-label" pair in a path or filename, in order of appearance.
bids_path_entities <- function(x) {
  m <- stringr::str_match_all(x, "(?:^|[_/])([a-z0-9]+)-([a-zA-Z0-9+]+)")[[1]]
  if (nrow(m) == 0) {
    return(stats::setNames(character(0), character(0)))
  }
  stats::setNames(m[, 3], m[, 2])
}

# The trailing token of a filename stem: "kinpar" in
# sub-01_model-2TCM_desc-model1_kinpar.tsv. A stem whose last token is itself an
# entity pair has no suffix.
bids_filename_suffix <- function(basename) {
  stem <- stringr::str_remove(basename, "\\.(tsv|json|nii\\.gz|nii)$")
  last <- stringr::str_extract(stem, "[^_]+$")
  if (is.na(last) || stringr::str_detect(last, "-")) NA_character_ else last
}

bids_compose_key <- function(entities, keys) {
  keys <- keys[keys %in% names(entities)]
  if (length(keys) == 0) {
    return(NA_character_)
  }
  paste(paste0(keys, "-", entities[keys]), collapse = "_")
}

#' Parse the derivative files of a PET BIDS derivatives folder
#'
#' @description Returns one row per derivative file, with the keys needed to
#'   join it to the acquisition it came from and to group the files that make up
#'   a single product.
#'
#'   This is deliberately a narrower contract than BIDS Derivatives. It groups
#'   by filename and does no inheritance, no entity grid and no acquisition
#'   enumeration, because a derivatives folder contains no PET images to
#'   enumerate. Running the raw parser over one instead invents entities: a
#'   subject stored without a session directory comes back carrying whichever
#'   session its neighbours happen to use, and then fails to join.
#'
#'   Three keys are returned per file:
#'
#'   \describe{
#'     \item{`source_key`}{*Which acquisition?* The selector entities
#'       (`sub`, `ses`, `task`, `trc`, `rec`, `run`) the file carries, with the
#'       directory supplying any the filename omits. `NA` for group-level files.}
#'     \item{`artifact_key`}{*Which product?* The source key plus every
#'       remaining entity plus the suffix. A `.tsv` and its `.json` sidecar
#'       describe one product and so share an artifact key.}
#'     \item{`analysis_scope_key`}{*Which analysis?* Set for group-level files,
#'       `NA` otherwise.}
#'   }
#'
#'   Files with no `sub` are group-level results. They get no `source_key`, so
#'   they cannot be joined to an acquisition by accident.
#'
#'   `.html` reports and any other non-data file are ignored: they are output
#'   for people to read, not artifacts to join on, and they carry no entities to
#'   tell them apart.
#'
#' @param path Path to the derivatives folder, e.g. a petfit analysis folder.
#'
#' @return A tibble with one row per derivative file: `path`, `path_absolute`,
#'   `extension`, `suffix`, the three keys above, and a column per entity found.
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' derivatives <- bids_parse_derivatives("derivatives/petfit/analysis3")
#' }
bids_parse_derivatives <- function(path) {

  if (!dir.exists(path)) {
    stop("Derivatives folder not found: ", path, call. = FALSE)
  }

  relative <- list.files(path, recursive = TRUE)

  # Data files only. Reports and other non-data files carry no entities, so
  # every one of them would key identically to the rest.
  relative <- relative[stringr::str_detect(relative, "\\.(tsv|json|nii\\.gz|nii)$")]

  scope <- basename(normalizePath(path, mustWork = FALSE))
  selectors <- bids_selector_entities()

  if (length(relative) == 0) {
    return(tibble::tibble(
      path = character(0), path_absolute = character(0),
      extension = character(0), suffix = character(0),
      source_key = character(0), artifact_key = character(0),
      analysis_scope_key = character(0)
    ))
  }

  parsed <- purrr::map(relative, function(rel) {

    base <- basename(rel)
    file_entities <- bids_path_entities(base)
    dir_entities <- bids_path_entities(paste0(dirname(rel), "/"))

    # A filename that names an entity outranks the directory it sits in.
    conflicting <- intersect(names(file_entities), names(dir_entities))
    conflicting <- conflicting[file_entities[conflicting] != dir_entities[conflicting]]

    # The directory completes what the filename omits. petfit's own derivatives
    # put ses in the path but not the filename, so a filename-only key would
    # fail to join for exactly the subjects that have sessions.
    #
    # What neither the filename nor the directory says stays unsaid. An entity a
    # file does not name imposes no constraint on it: a weights file called
    # sub-01_desc-weights_weights.tsv sitting beside trc-A and trc-B results
    # applies to both. So source_key is whatever the file does specify, and
    # matching it against acquisitions is a subset test rather than string
    # equality. This function could not decide that in any case -- a derivatives
    # folder holds no PET images, so it never sees the acquisitions to match.
    completed <- file_entities
    for (entity in names(dir_entities)) {
      if (!entity %in% names(completed)) {
        completed[entity] <- dir_entities[entity]
      }
    }

    list(
      rel = rel,
      entities = completed,
      conflicting = conflicting,
      from_directory = setdiff(names(dir_entities), names(file_entities)),
      extension = stringr::str_extract(base, "(nii\\.gz|tsv|json|nii)$"),
      suffix = bids_filename_suffix(base)
    )
  })

  conflicts <- purrr::keep(parsed, ~ length(.x$conflicting) > 0)
  if (length(conflicts) > 0) {
    detail <- purrr::map_chr(conflicts, function(p) {
      paste0("  ", p$rel, " (disagrees on: ",
             paste(p$conflicting, collapse = ", "), ")")
    })
    stop("A filename and the directory holding it name different values for ",
         "the same entity, so the acquisition it belongs to is ambiguous:\n",
         paste(detail, collapse = "\n"), call. = FALSE)
  }

  out <- purrr::map_dfr(parsed, function(p) {

    entities <- p$entities
    is_group <- !("sub" %in% names(entities))

    source_key <- if (is_group) {
      NA_character_
    } else {
      bids_compose_key(entities, selectors)
    }

    # Every remaining entity, never a hand-maintained list of the ones we happen
    # to recognise: an unfamiliar but valid entity would otherwise collapse two
    # distinct files into one identity.
    qualifiers <- sort(setdiff(names(entities), selectors))

    artifact_key <- paste(stats::na.omit(c(
      if (is_group) scope else source_key,
      bids_compose_key(entities, qualifiers),
      p$suffix
    )), collapse = "_")

    row <- tibble::tibble(
      path = p$rel,
      path_absolute = file.path(normalizePath(path, mustWork = FALSE), p$rel),
      extension = p$extension,
      suffix = p$suffix,
      source_key = source_key,
      artifact_key = artifact_key,
      analysis_scope_key = if (is_group) scope else NA_character_
    )

    for (entity in names(entities)) {
      row[[entity]] <- unname(entities[entity])
    }

    row
  })


  out
}

