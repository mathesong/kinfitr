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

  warning(
    "bids_parse_files() is deprecated in favour of bids_parse_filenames(), ",
    "and will be removed in 2027.",
    call. = FALSE)


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

  # Which blood belongs to this acquisition, and is it resolvable? Selecting the
  # files here rather than filtering inline gets the checks that filtering never
  # did: that blood sits in pet/, that it carries a recording entity, and that
  # two files do not both claim to be the same recording. Filtering also matched
  # only the two labels named below, silently ignoring any other.
  # bids_associate_blood() refuses blood outside pet/, blood with no recording
  # entity, and two files both claiming one recording. Report and skip: this
  # runs once per acquisition, so raising would end the whole study parse.
  blood <- tryCatch(
    bids_associate_blood(filedata),
    error = function(e) {
      warning(conditionMessage(e),
              "\nThe blood for this acquisition was skipped. The acquisition ",
              "itself is unaffected, as is the rest of the study.",
              call. = FALSE)
      NULL
    })

  if (is.null(blood) || nrow(blood) == 0) {
    # An acquisition with no blood files at all is ordinary and passes in
    # silence. Blood that is present but unusable has already been warned
    # about, per recording, by bids_associate_blood() itself.
    return(NA)
  }

  pick <- function(recording, extension) {
    row <- blood[blood$recording == recording, , drop = FALSE]
    if (nrow(row) == 0) return(character(0))
    row[[extension]]
  }

  json_blood_discrete <- pick("manual", "json")
  tsv_blood_discrete  <- pick("manual", "tsv")
  json_blood_cont     <- pick("autosampler", "json")
  tsv_blood_cont      <- pick("autosampler", "tsv")

  # The reader below understands manual and autosampler recordings. Any other
  # label is carried by bids_associate_blood() but cannot be interpreted here,
  # so say so rather than dropping it without comment.
  unhandled <- setdiff(blood$recording, c("manual", "autosampler"))
  if (length(unhandled) > 0) {
    warning("Blood recording", if (length(unhandled) > 1) "s" else "", " ",
            paste0("recording-", unhandled, collapse = ", "),
            " for ", attr(filedata, "pet_key") %||% "this acquisition",
            " could not be read: only manual and autosampler are understood.",
            call. = FALSE)
  }

  if (length(json_blood_discrete) == 0) {
    warning("No manual blood recording for ",
            attr(filedata, "pet_key") %||% "this acquisition",
            "; its blood could not be read.", call. = FALSE)
    return(NA)
  }





  ### Get the data ###

  jsondat_blood_discrete <- bids_read_sidecar(json_blood_discrete, filedata,
                                             "manual blood sidecar")

  if (is.null(jsondat_blood_discrete)) {
    return(NA)
  }

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
# file's relative path.
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

# Read a JSON sidecar, returning NULL and reporting against the acquisition if it
# cannot be parsed. Each caller runs once per acquisition, so raising would end
# the whole study parse.
bids_read_sidecar <- function(path, filedata, what) {
  tryCatch(jsonlite::fromJSON(path), error = function(e) {
    warning("The ", what, " for ", bids_describe_measurement(filedata),
            " could not be read: ", conditionMessage(e),
            call. = FALSE)
    NULL
  })
}

bids_parse_pettimes <- function(filedata) {

  if(!("pet" %in% filedata$measurement)) {
    return(NA)
  }

  ### Get the filenames ###

  json_pet <- filedata %>%
    dplyr::filter(measurement=="pet" & extension=="json")

  # Frame times come from the _pet.json sidecar. Report an acquisition that has
  # none and return NA so it is dropped: this is mapped over every measurement,
  # so raising would end the whole parse.
  if (nrow(json_pet) == 0) {
    warning("No _pet.json sidecar found for ", bids_describe_measurement(filedata),
            ". Frame times cannot be determined, so this measurement is ",
            "excluded. Every PET image needs an accompanying _pet.json.",
            call. = FALSE)
    return(NA)
  }

  ### Extract the data ###

  # Merge the applicable sidecars nearest-last rather than flattening them, so a
  # field set close to the data file wins and the result does not depend on the
  # order the files were listed in.
  data_file <- filedata$path[filedata$measurement == "pet" &
                               filedata$extension %in% c("nii", "nii.gz")]
  if (length(data_file) == 0) data_file <- json_pet$path
  resolved <- bids_resolve_sidecars(data_file[1], json_pet$path,
                                    sidecar_root = attr(filedata, "study_root"))
  jsondat_pet <- resolved$values

  if (length(jsondat_pet) == 0) {
    return(NA)
  }

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

  json_pet <- dplyr::filter(filedata,
                            measurement == "pet" & extension == "json")

  if (nrow(json_pet) == 0) {
    warning("No _pet.json sidecar found for ",
            bids_describe_measurement(filedata),
            ". This acquisition has no metadata; the rest of the study is ",
            "unaffected.", call. = FALSE)
    return(NA)
  }

  ### Extract the data ###

  # Merge the applicable sidecars nearest-last rather than flattening them, so a
  # field set close to the data file wins and the result does not depend on the
  # order the files were listed in.
  data_file <- filedata$path[filedata$measurement == "pet" &
                               filedata$extension %in% c("nii", "nii.gz")]
  if (length(data_file) == 0) data_file <- json_pet$path
  resolved <- bids_resolve_sidecars(data_file[1], json_pet$path,
                                    sidecar_root = attr(filedata, "study_root"))
  jsondat_pet <- resolved$values

  if (length(jsondat_pet) == 0) {
    return(NA)
  }

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

  # One row per acquisition actually on disk. Grouping is by whichever selector
  # entities this study uses -- naming them literally would fail on a study that
  # has no task or no acq, which the previous parser hid by fabricating both.
  measurements <- bids_parse_filenames(studypath)

  grouping <- intersect(bids_selector_entities(), colnames(measurements))
  measurements <- dplyr::group_by(measurements,
                                  dplyr::across(dplyr::all_of(grouping)))


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

  # Resolved once, outside the row loop: tibble() evaluates its arguments in its
  # own scope, so referring to `path` there would pick up the path column rather
  # than this argument.
  root <- normalizePath(path, mustWork = FALSE)
  scope <- basename(root)
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

    # The directory completes what the filename omits: petfit's derivatives put
    # ses in the path but not the filename, and a filename-only key would not
    # join for the subjects that have sessions.
    #
    # An entity a file does not name imposes no constraint on it: a weights file
    # called sub-01_desc-weights_weights.tsv sitting beside trc-A and trc-B
    # results applies to both. source_key is therefore whatever the file itself
    # specifies, and matching it against acquisitions is a subset test rather
    # than string equality. A derivatives folder holds no PET images, so that
    # match cannot be resolved here.
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

    # Every remaining entity, so that an unfamiliar but valid one still
    # separates two files rather than collapsing them into one identity.
    qualifiers <- sort(setdiff(names(entities), selectors))

    # The directories between the parse root and the subject folder. A
    # derivatives tree can hold several analyses of the same acquisition, and
    # without this they would be one artifact rather than several.
    parts <- strsplit(p$rel, "/", fixed = TRUE)[[1]]
    subject_at <- which(startsWith(parts, "sub-"))
    path_scope <- if (length(subject_at) > 0) {
      if (subject_at[1] > 1) paste(parts[seq_len(subject_at[1] - 1)], collapse = "/") else NA_character_
    } else {
      if (length(parts) > 1) paste(utils::head(parts, -1), collapse = "/") else NA_character_
    }

    artifact_key <- paste(stats::na.omit(c(
      if (!is.na(path_scope)) path_scope else if (is_group) scope else NA_character_,
      if (is_group) NA_character_ else source_key,
      bids_compose_key(entities, qualifiers),
      p$suffix
    )), collapse = "_")

    row <- tibble::tibble(
      path = p$rel,
      path_absolute = file.path(root, p$rel),
      extension = p$extension,
      suffix = p$suffix,
      source_key = source_key,
      artifact_key = artifact_key,
      # The analysis a group-level file belongs to is the directory holding it
      # within the parse root, not the root itself: parsing a parent of
      # several analyses must not give their configs one shared scope.
      analysis_scope_key = if (is_group) {
        if (!is.na(path_scope)) path_scope else scope
      } else {
        NA_character_
      }
    )

    for (entity in names(entities)) {
      row[[entity]] <- unname(entities[entity])
    }

    row
  })


  out
}



# --- Measurement enumeration --------------------------------------------------

# Does a file belong to a given acquisition?
#
# The subset rule: every selector entity the *file* names must match the
# acquisition's value for it. A selector the file leaves out imposes no
# constraint, so blood recorded once for a subject reaches both of that
# subject's reconstructions. A file naming a selector the acquisition does not
# have at all does not belong to it -- blood claiming rec-A cannot attach to an
# acquisition with no rec.
#
# `lenient` relaxes that last condition for the entities it names, and is set
# only for metadata sidecars (Part D: sidecar -> data is the one lenient
# direction). A root-level task-rest_pet.json is meant to apply to images that
# name no task; matched strictly it applies to nothing, and the frame times it
# carries become unreachable. Data files always pass character(0): blood
# claiming a rec no acquisition has stays unattached rather than guessing.
bids_file_belongs_to <- function(file_entities, pet_entities, selectors,
                                 lenient = character(0)) {

  claimed <- intersect(names(file_entities), selectors)

  for (entity in claimed) {
    if (!entity %in% names(pet_entities)) {
      if (entity %in% lenient) next
      return(FALSE)
    }
    if (!identical(unname(file_entities[entity]),
                   unname(pet_entities[entity]))) return(FALSE)
  }

  TRUE
}

# The entities a sidecar may name that its data file does not. Kept in one
# place so enumeration-time attachment and bids_resolve_sidecars() cannot
# drift apart. trc and run stay strict: attaching the wrong tracer's sidecar
# silently supplies the wrong half-life and injected dose.
bids_lenient_entities <- function() {
  c("task", "rec")
}

#' Extract the measurements of a PET BIDS study from its filenames
#'
#' @description Returns one row per PET acquisition found on disk, with the
#'   files belonging to it nested in `filedata`.
#'
#'   An acquisition's entities are exactly those its filename and directory
#'   carry. Nothing is filled in: a study whose files name no `trc` has no `trc`
#'   column, rather than a fabricated one. This matters because the previous
#'   behaviour built a grid of every entity value against every other and kept
#'   the combinations that matched a file, which invents acquisitions that do
#'   not exist -- on a deliberately ragged 12-acquisition fixture it returns 192
#'   measurements, giving subjects sessions and tracers belonging to their
#'   neighbours.
#'
#'   An acquisition is any `_pet.nii`, `_pet.nii.gz` or `_pet.json` inside a
#'   `pet/` directory; the image and its sidecar count once. Requiring `pet/`
#'   excludes study-level inheritance sidecars such as `task-rest_pet.json`. A
#'   `_pet.json` that applies as an inheritance sidecar to another acquisition
#'   -- `sub-01_pet.json` beside `run-01` and `run-02` images -- is metadata for
#'   those acquisitions, not one of its own. One that applies to no other
#'   acquisition does count on its own, so a `_pet.json` whose image is not
#'   found still yields a measurement -- most commonly one whose blood is
#'   collected before the image is available.
#'
#'   Files are attached by the subset rule: every selector entity a file names
#'   must match, and one it omits imposes no constraint. Blood recorded once per
#'   subject therefore reaches every reconstruction of it, while blood naming a
#'   `rec` no acquisition has attaches to nothing. Metadata sidecars are matched
#'   leniently for `task` and `rec`, so a study-level `task-rest_pet.json`
#'   reaches images that name no task; data files are always matched strictly.
#'
#' @param studypath The BIDS study path for the main study.
#'
#' @return A tibble with one row per acquisition: a column per entity the study
#'   actually uses, and `filedata` holding that acquisition's files.
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' measurements <- bids_parse_filenames(studypath)
#' }
bids_parse_filenames <- function(studypath) {

  selectors <- bids_selector_entities()
  study_root <- normalizePath(studypath, mustWork = FALSE)

  extensions <- paste(c("*.nii.gz", "*.nii", "*.tsv", "*.json"), collapse = "|")
  files <- fs::dir_info(studypath, recurse = TRUE, type = "file",
                        glob = extensions)

  filedata <- tibble::tibble(
    path_absolute = as.character(files$path),
    path = as.character(fs::path_rel(files$path, studypath))
  )

  filedata <- filedata[!stringr::str_detect(
    filedata$path, "^(derivatives|code|phenotype|sourcedata)/"), , drop = FALSE]

  if (nrow(filedata) == 0) {
    return(tibble::tibble())
  }

  filedata$extension <- ifelse(
    stringr::str_detect(filedata$path, "\\.nii\\.gz$"), "nii.gz",
    fs::path_ext(filedata$path))
  filedata$measurement <- purrr::map_chr(basename(filedata$path),
                                         ~ bids_filename_suffix(.x) %||% NA_character_)
  filedata <- filedata[!is.na(filedata$measurement), , drop = FALSE]

  file_entities <- purrr::map(filedata$path, bids_path_entities)

  # filedata carries the entities describing the file rather than the
  # acquisition: recording, desc, seg and the like. Consumers filter on them --
  # bids_parse_blood() needs `recording`. The selectors are on the measurement
  # row itself.
  qualifiers <- setdiff(unique(unlist(lapply(file_entities, names))), selectors)
  for (entity in qualifiers) {
    filedata[[entity]] <- purrr::map_chr(file_entities, function(e) {
      if (entity %in% names(e)) unname(e[entity]) else NA_character_
    })
  }

  # One column per entity, so labels are compared entity by entity.
  all_entities <- unique(unlist(lapply(file_entities, names)))
  entity_table <- tibble::tibble(.rows = length(file_entities))
  for (entity in all_entities) {
    entity_table[[entity]] <- purrr::map_chr(file_entities, function(e) {
      if (entity %in% names(e)) unname(e[entity]) else NA_character_
    })
  }
  bids_warn_case_collisions(entity_table)

  # An acquisition is a _pet.* file inside pet/. The image and its sidecar are
  # the same acquisition, so they collapse to one entry.
  is_pet <- filedata$measurement == "pet" &
    stringr::str_detect(filedata$path, "(^|/)pet/") &
    filedata$extension %in% c("nii", "nii.gz", "json")

  if (!any(is_pet)) {
    return(tibble::tibble())
  }

  pet_entities <- file_entities[is_pet]
  pet_keys <- purrr::map_chr(pet_entities, function(e) {
    bids_compose_key(e, selectors[selectors %in% names(e)])
  })

  # Two images for one acquisition cannot be told apart, so refuse rather than
  # keeping whichever comes first.
  images <- filedata$extension[is_pet] %in% c("nii", "nii.gz")
  duplicated_images <- unique(pet_keys[images][duplicated(pet_keys[images])])
  if (length(duplicated_images) > 0) {
    offenders <- filedata$path[is_pet][images][pet_keys[images] %in% duplicated_images]
    stop("More than one PET image describes the same acquisition:\n",
         paste0("  ", offenders, collapse = "\n"),
         "\nEach acquisition needs entities that tell it apart from the others.",
         call. = FALSE)
  }

  keep <- !duplicated(pet_keys)
  measurements_entities <- pet_entities[keep]
  measurement_keys <- pet_keys[keep]
  has_image <- measurement_keys %in% pet_keys[images]

  # A _pet.json inside pet/ is not necessarily an acquisition: BIDS inheritance
  # lets one sidecar serve several images, so sub-01_pet.json beside run-01 and
  # run-02 describes both runs rather than a third, run-less scan. A json-only
  # entry is therefore an acquisition only if it does not apply as a sidecar to
  # another acquisition -- one with an image, or a json-only one it is less
  # specific than. A _pet.json matching no image at all still counts: an
  # acquisition whose image is not found (most commonly, blood collected
  # before the image is available) remains a measurement.
  selector_count <- vapply(measurements_entities, function(e) {
    length(intersect(names(e), selectors))
  }, integer(1))

  is_inherited_sidecar <- vapply(seq_along(measurements_entities), function(i) {
    if (has_image[i]) return(FALSE)
    for (j in seq_along(measurements_entities)) {
      if (j == i) next
      if (!has_image[j] && selector_count[j] <= selector_count[i]) next
      if (bids_file_belongs_to(measurements_entities[[i]],
                               measurements_entities[[j]], selectors,
                               lenient = bids_lenient_entities())) {
        return(TRUE)
      }
    }
    FALSE
  }, logical(1))

  measurements_entities <- measurements_entities[!is_inherited_sidecar]

  used <- unique(unlist(lapply(measurements_entities, names)))
  used <- selectors[selectors %in% used]

  out <- purrr::map_dfr(seq_along(measurements_entities), function(i) {

    entities <- measurements_entities[[i]]

    row <- tibble::tibble(.rows = 1)
    for (entity in used) {
      row[[entity]] <- if (entity %in% names(entities)) {
        unname(entities[entity])
      } else {
        NA_character_
      }
    }

    # Metadata sidecars attach leniently (a task-rest_pet.json reaches an image
    # that names no task); data files attach strictly. A _blood.json follows
    # its tsv rather than getting leniency of its own, so a pair never splits.
    belongs <- purrr::map_lgl(seq_along(file_entities), function(k) {
      lenient <- if (filedata$measurement[k] == "pet" &&
                     filedata$extension[k] == "json") {
        bids_lenient_entities()
      } else {
        character(0)
      }
      bids_file_belongs_to(file_entities[[k]], entities, selectors, lenient)
    })

    # The acquisition's identity travels with its files as attributes, so
    # bids_create_blooddata() keeps taking a bare table. A hand-built one works
    # too; it simply has no key to name itself with in messages.
    acquisition_files <- filedata[belongs, , drop = FALSE]
    attr(acquisition_files, "pet_key") <-
      bids_compose_key(entities, selectors[selectors %in% names(entities)])
    attr(acquisition_files, "study_root") <- study_root

    row$filedata <- list(acquisition_files)

    row
  })

  out
}


# --- Blood association --------------------------------------------------------

#' Associate blood files with a PET acquisition
#'
#' @description Works out which blood files belong to one acquisition, and
#'   refuses the combinations that cannot be resolved.
#'
#'   BIDS requires blood to live beside the PET data it belongs to -- the same
#'   `pet/` directory as the acquisition's own `_pet.*` files, not merely some
#'   `pet/` directory -- and to carry a `recording` entity, so both are
#'   enforced here. Several recordings for one acquisition are normal and
#'   expected -- manual samples alongside an autosampler is explicitly
#'   permitted -- but two complete file pairs claiming the *same* recording
#'   cannot both be it, and a json pairs with a tsv only if every entity the
#'   json names is named by the tsv with the same value.
#'
#'   For an acquisition with no blood, the
#'   returned table is simply empty, and its `excluded` attribute says what was
#'   left out of the *blood association*, so that "this scan has no blood" can
#'   be told apart from "this scan's blood failed to attach".
#'
#' @param filedata The files belonging to one acquisition, as nested by
#'   [bids_parse_filenames()]. Its `pet_key` attribute, if present, is used to
#'   name the acquisition in messages.
#'
#' @return A tibble with one row per recording: `recording`, `tsv`, `json`. Zero
#'   rows when the acquisition has no usable blood, carrying an `excluded`
#'   attribute explaining which case it is. That refers to the blood, not the
#'   acquisition, which remains perfectly usable without it.
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' blood <- bids_associate_blood(measurements$filedata[[1]])
#' }
bids_associate_blood <- function(filedata) {

  key <- attr(filedata, "pet_key")
  if (is.null(key) || is.na(key)) key <- "this acquisition"

  empty <- function(reason) {
    out <- tibble::tibble(recording = character(0),
                          tsv = character(0),
                          json = character(0))
    attr(out, "excluded") <- reason
    out
  }

  if (is.null(filedata) || nrow(filedata) == 0 ||
      !("blood" %in% filedata$measurement)) {
    return(empty("no blood files"))
  }

  blood <- filedata[filedata$measurement == "blood", , drop = FALSE]

  # Spec-mandated: blood belongs in the acquisition's pet/ directory.
  misplaced <- blood$path[!stringr::str_detect(blood$path, "(^|/)pet/")]
  if (length(misplaced) > 0) {
    stop("Blood files for ", key, " are outside a pet/ directory:\n",
         paste0("  ", misplaced, collapse = "\n"),
         "\nBIDS places blood recordings alongside the PET data they belong to.",
         call. = FALSE)
  }

  # And not merely any pet/ directory: the one holding this acquisition's own
  # _pet.* files. Without this, subject-level blood under sub-01/pet/ attaches
  # to an acquisition under sub-01/ses-01/pet/, silently crossing scope.
  pet_dirs <- unique(dirname(filedata$path[filedata$measurement == "pet"]))
  pet_dirs <- pet_dirs[stringr::str_detect(pet_dirs, "(^|/)pet$")]
  if (length(pet_dirs) > 0) {
    astray <- blood$path[!dirname(blood$path) %in% pet_dirs]
    if (length(astray) > 0) {
      stop("Blood files for ", key, " are not beside its PET data (",
           paste(pet_dirs, collapse = ", "), "):\n",
           paste0("  ", astray, collapse = "\n"),
           "\nBIDS places blood recordings alongside the PET data they belong ",
           "to, in the same directory.", call. = FALSE)
    }
  }

  # Spec-mandated: recording distinguishes one blood recording from another.
  if (!"recording" %in% colnames(blood)) {
    blood$recording <- NA_character_
  }
  unlabelled <- blood$path[is.na(blood$recording)]
  if (length(unlabelled) > 0) {
    stop("Blood files for ", key, " carry no recording entity:\n",
         paste0("  ", unlabelled, collapse = "\n"),
         "\nWithout it there is no way to tell one recording from another.",
         call. = FALSE)
  }

  recordings <- sort(unique(blood$recording))

  pairs <- purrr::map_dfr(recordings, function(recording) {

    this <- blood[blood$recording == recording, , drop = FALSE]
    tsv <- this$path_absolute[this$extension == "tsv"]
    json <- this$path_absolute[this$extension == "json"]

    # Several recordings per acquisition are fine, and so are files that differ
    # by a selector entity: run-01 and run-02 blood belong to different
    # acquisitions and are never compared here, because filedata holds one
    # acquisition's files. What cannot be resolved is two files claiming the
    # same recording of the *same* acquisition.
    if (length(tsv) > 1 || length(json) > 1) {
      stop("More than one blood file pair for ", key,
           " claims recording-", recording, ":\n",
           paste0("  ", this$path, collapse = "\n"),
           "\nGive them distinct recording labels.", call. = FALSE)
    }

    # A sidecar describes the tsv it names: every entity the json carries must
    # be carried by the tsv with the same value (naming fewer is fine). A
    # rec-A json beside a tsv that names no rec may describe different data,
    # so recording alone is not enough to pair them.
    if (length(tsv) == 1 && length(json) == 1) {
      tsv_entities <- bids_path_entities(basename(this$path[this$extension == "tsv"]))
      json_entities <- bids_path_entities(basename(this$path[this$extension == "json"]))
      extra <- setdiff(names(json_entities), names(tsv_entities))
      shared <- intersect(names(json_entities), names(tsv_entities))
      mismatched <- shared[json_entities[shared] != tsv_entities[shared]]
      disagreeing <- c(extra, mismatched)
      if (length(disagreeing) > 0) {
        warning("The blood json for ", key, ", recording-", recording,
                " names ", paste(disagreeing, collapse = ", "),
                " differently from its tsv, so it is not treated as that ",
                "tsv's sidecar:\n",
                paste0("  ", this$path, collapse = "\n"),
                call. = FALSE)
        json <- character(0)
      }
    }

    tibble::tibble(
      recording = recording,
      tsv = if (length(tsv)) tsv else NA_character_,
      json = if (length(json)) json else NA_character_
    )
  })

  # Every incomplete recording is reported, whether or not another recording is
  # complete: an autosampler tsv missing its json must not vanish in silence
  # just because the manual samples are fine.
  incomplete <- pairs[is.na(pairs$tsv) | is.na(pairs$json), , drop = FALSE]
  if (nrow(incomplete) > 0) {
    warning("Blood for ", key,
            " without a complete tsv/json pair, skipped:\n",
            paste0("  recording-", incomplete$recording,
                   ifelse(is.na(incomplete$tsv), " missing tsv", " missing json"),
                   collapse = "\n"), call. = FALSE)
  }

  complete <- pairs[!is.na(pairs$tsv) & !is.na(pairs$json), , drop = FALSE]

  if (nrow(complete) == 0) {
    return(empty(paste0(
      "blood files present but no complete tsv/json pair (",
      paste0("recording-", pairs$recording,
             ifelse(is.na(pairs$tsv), " missing tsv", " missing json"),
             collapse = "; "), ")")))
  }

  complete
}


# --- Sidecar resolution -------------------------------------------------------

# Is a sidecar applicable to a data file?
#
# Three conditions: the same suffix, a directory that is an ancestor of the data
# file's, and every entity the sidecar names also named by the data file with the
# same value.
#
# The last condition is relaxed for the entities in `lenient`. A root-level
# task-rest_pet.json is meant to apply to PET files that name no task, and
# without leniency it applies to nothing at all. trc and run stay strict: a
# sidecar specific enough to name the tracer should match an image that names it
# too, and attaching the wrong one would silently supply the wrong half-life and
# injected dose.
bids_sidecar_applies <- function(sidecar_path, data_path, lenient) {

  if (!identical(bids_filename_suffix(basename(sidecar_path)),
                 bids_filename_suffix(basename(data_path)))) {
    return("suffix differs")
  }

  sidecar_dir <- dirname(sidecar_path)
  data_dir <- dirname(data_path)
  ancestor <- sidecar_dir == "." || sidecar_dir == data_dir ||
    startsWith(paste0(data_dir, "/"), paste0(sidecar_dir, "/"))
  if (!ancestor) {
    return("not in a parent directory of the data file")
  }

  sidecar_entities <- bids_path_entities(sidecar_path)
  data_entities <- bids_path_entities(data_path)

  for (entity in names(sidecar_entities)) {
    if (!entity %in% names(data_entities)) {
      if (entity %in% lenient) next
      return(paste0("names ", entity, "-", sidecar_entities[entity],
                    ", which the data file does not"))
    }
    if (!identical(unname(sidecar_entities[entity]),
                   unname(data_entities[entity]))) {
      return(paste0("names ", entity, "-", sidecar_entities[entity],
                    " but the data file names ", entity, "-",
                    data_entities[entity]))
    }
  }

  NA_character_
}

#' Resolve the metadata sidecars applying to a data file
#'
#' @description Merge the JSON sidecars that apply to one data file, nearest
#'   directory last, and report where each field came from.
#'
#'   Merging is key by key from the study root inwards, so a field set close to
#'   the data file overrides the same field set further out, and fields set only
#'   at the root are inherited. Two applicable sidecars in the same directory are
#'   refused: nothing distinguishes them, and picking by enumeration order makes
#'   the frame times depend on how the filesystem happens to list files.
#'
#'   Every sidecar that does not apply is reported, with the reason.
#'
#' @param data_path Path of the data file, relative to the study root.
#' @param sidecar_paths Candidate sidecar paths, relative to the same root.
#' @param sidecar_root Optional directory to prepend when reading the files.
#' @param lenient Entities a sidecar may name that the data file does not.
#'   Defaults to `task` and `rec`.
#'
#' @return A list with `values`, the merged fields, and `provenance`, naming the
#'   sidecar each field came from.
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' bids_resolve_sidecars("sub-01/ses-test/pet/sub-01_ses-test_pet.nii.gz",
#'                       c("task-rest_pet.json",
#'                         "sub-01/ses-test/pet/sub-01_ses-test_pet.json"))
#' }
bids_resolve_sidecars <- function(data_path, sidecar_paths, sidecar_root = NULL,
                                  lenient = bids_lenient_entities()) {

  empty <- list(values = list(), provenance = character(0))

  if (length(sidecar_paths) == 0) {
    return(empty)
  }

  reasons <- vapply(sidecar_paths, bids_sidecar_applies, character(1),
                    data_path, lenient, USE.NAMES = FALSE)

  skipped <- !is.na(reasons)
  if (any(skipped)) {
    warning("Sidecars not applied to ", data_path, ":\n",
            paste0("  ", sidecar_paths[skipped], " -- ", reasons[skipped],
                   collapse = "\n"), call. = FALSE)
  }

  applicable <- sidecar_paths[!skipped]
  if (length(applicable) == 0) {
    return(empty)
  }

  depth <- lengths(strsplit(applicable, "/", fixed = TRUE))

  duplicated_depth <- unique(depth[duplicated(depth)])
  if (length(duplicated_depth) > 0) {
    clashing <- applicable[depth %in% duplicated_depth]
    stop("More than one sidecar applies to ", data_path,
         " from the same directory:\n",
         paste0("  ", clashing, collapse = "\n"),
         "\nNothing distinguishes them, so which one wins would depend on the ",
         "order the files are listed in.", call. = FALSE)
  }

  applicable <- applicable[order(depth)]

  values <- list()
  provenance <- character(0)

  for (sidecar in applicable) {
    full <- if (is.null(sidecar_root)) sidecar else file.path(sidecar_root, sidecar)
    fields <- tryCatch(jsonlite::fromJSON(full), error = function(e) {
      warning("The sidecar ", sidecar, " could not be read: ",
              conditionMessage(e), call. = FALSE)
      NULL
    })
    if (is.null(fields) || length(fields) == 0) next

    for (key in names(fields)) {
      values[[key]] <- fields[[key]]
      provenance[key] <- sidecar
    }
  }

  list(values = values, provenance = provenance)
}
