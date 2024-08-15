#' Create a blooddata object from data vectors
#'
#' This function creates a blooddata object from data vectors. This function
#' creates data objects which contain the minimum amount of data for modelling.
#' Leave blank fields which do not exist.
#'
#' Ideally, I recommend storing your data according to the BIDS specification
#' and using the bids_create_blooddata() command instead. Please refer to the PET
#' BIDS standard for further details about the inputs.
#'
#' @param Blood.Discrete.Values.time Sample times in seconds.
#' @param Blood.Discrete.Values.activity Measured radioactivity in kBq/ml.
#' @param Plasma.Values.time In seconds. Sample times in seconds.
#' @param Plasma.Values.activity Measured radioactivity in kBq/ml.
#' @param Metabolite.Values.time Sample times in seconds.
#' @param Metabolite.Values.parentFraction Measured fraction (unitless).
#' @param Blood.Continuous.Values.time In seconds.
#' @param Blood.Continuous.Values.activity in kBq/ml.
#' @param Blood.Continuous.DispersionConstant External dispersion time constant
#'   resulting from tubing in seconds.
#' @param Blood.Continuous.DispersionCorrected Boolean flag specifying whether
#'   the continuous blood data have been dispersion-corrected.
#' @param TimeShift The extent to which all the times in the data should be
#'   shifted (in seconds). Defaults to 0.
#'
#' @return a blooddata object
#'
#' @examples
#'
#' blooddata <- pbr28$blooddata[[1]]
#'
#' blooddata2 <- create_blooddata_components(
#'    Blood.Discrete.Values.time =
#'      blooddata$Data$Blood$Discrete$Values$time,
#'    Blood.Discrete.Values.activity =
#'      blooddata$Data$Blood$Discrete$Values$activity,
#'    Plasma.Values.time =
#'      blooddata$Data$Plasma$Values$time,
#'    Plasma.Values.activity =
#'      blooddata$Data$Plasma$Values$activity,
#'    Metabolite.Values.time =
#'      blooddata$Data$Metabolite$Values$time,
#'    Metabolite.Values.parentFraction =
#'      blooddata$Data$Metabolite$Values$parentFraction,
#'    Blood.Continuous.Values.time =
#'      blooddata$Data$Blood$Continuous$Values$time,
#'    Blood.Continuous.Values.activity =
#'      blooddata$Data$Blood$Continuous$Values$activity,
#'    Blood.Continuous.DispersionConstant =
#'      blooddata$Data$Blood$Continuous$DispersionConstant,
#'    Blood.Continuous.DispersionCorrected = FALSE,
#'    TimeShift = 0)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
create_blooddata_components <- function(
                                        Blood.Discrete.Values.time = NULL,
                                        Blood.Discrete.Values.activity = NULL,
                                        Plasma.Values.time = NULL,
                                        Plasma.Values.activity = NULL,
                                        Metabolite.Values.time = NULL,
                                        Metabolite.Values.parentFraction = NULL,
                                        Blood.Continuous.Values.time = NULL,
                                        Blood.Continuous.Values.activity = NULL,
                                        Blood.Continuous.DispersionConstant = NULL,
                                        Blood.Continuous.DispersionCorrected = TRUE,
                                        TimeShift = 0) {




  # Blood

  Blood <- list()

  # Continuous blood

  if( !is.null(Blood.Continuous.Values.time) ) {

    Blood$Continuous <- list(
      Values = tibble::tibble(
        time = Blood.Continuous.Values.time,
        activity = Blood.Continuous.Values.activity
      ),
      DispersionConstant = Blood.Continuous.DispersionConstant,
      DispersionCorrected = Blood.Continuous.DispersionCorrected,
      Avail = TRUE,
      time = list(
        Description = "No description",
        Units = "s"
      ),
      activity = list(
        Description = "No description",
        Units = "kBq/ml"
      )
    )
  } else {
    Blood$Continuous <- list(
      Avail = FALSE)
  }


  # Discrete blood

  if( !is.null(Blood.Discrete.Values.time) ) {

    Blood$Discrete <- list(
      Values = tibble::tibble(
        time = Blood.Discrete.Values.time,
        activity = Blood.Discrete.Values.activity
      ),
      Avail=TRUE,
      time = list(
        Description = "No description",
        Units = "s"
      ),
      activity = list(
        Description = "No description",
        Units = "kBq/ml"
      )
    )
  } else {
    Blood$Discrete <- list(
      Avail = FALSE)
  }

  # Plasma

  if( !is.null(Plasma.Values.time) ) {

  Plasma <- list(
    Values = tibble::tibble(
      time = Plasma.Values.time,
      activity = Plasma.Values.activity
    ),
    Avail = TRUE,
    time = list(
      Description = "No description",
      Units = "s"
    ),
    activity = list(
      Description = "No description",
      Units = "kBq/ml"
    )
  )
  } else {
    Plasma <- list(
      Avail = FALSE
    )
  }

  # Metabolite
  if( !is.null(Metabolite.Values.time) ) {

    Metabolite <- list(
      Values = tibble::tibble(
        time = Metabolite.Values.time,
        parentFraction = Metabolite.Values.parentFraction
      ),
      Avail = TRUE,
      time = list(
        Description = "No description",
        Units = "s"
      ),
      parentFraction = list(
        Description = "No description",
        Units = "unitless"
      )
    )
  } else {
    Metabolite <- list(
      Avail = FALSE
    )
  }


  ## Edge cases

  ### No plasma, but only whole blood: use blood instead of plasma
  if( !Plasma$Avail & Blood$Discrete$Avail ) {

    Plasma <- Blood$Discrete

    Plasma$activity$Description <- paste(
      "Whole blood used as no plasma available.")
  }

  ### No metabolite, but blood/plasma: metab=1
  if( !Metabolite$Avail ) {

    Metabolite$Values <- tibble::tibble(
      time = Plasma$Values$time,
      parentFraction = rep(1, nrow(Plasma$Values))
    )

    Metabolite$Avail = TRUE
    Metabolite$time <- Plasma$time
    Metabolite$parentFraction <- list(
      Description = "All set to 1 because no metabolite data available",
      Units = "unitless"
    )
  }

  ### No whole blood, only plasma: use plasma as blood
  if( Plasma$Avail & !Blood$Discrete$Avail ) {

    Blood$Discrete <- Plasma

    Blood$Discrete$activity$Description <- paste(
      "Plasma used as no whole blood available."
    )
  }


  blooddata <- list(
    Data = list(
      Blood = Blood,
      Plasma = Plasma,
      Metabolite = Metabolite
    ),
    Models = list(
      Blood = list(Method = "interp", Data = NULL),
      BPR = list(Method = "interp", Data = NULL),
      parentFraction = list(Method = "interp", Data = NULL),
      AIF = list(Method = "interp", Data = NULL)
    ),
    TimeShift = TimeShift
  )

  blooddata <- rapply(blooddata, function(x) ifelse(x == "true", TRUE, x), how = "replace")
  blooddata <- rapply(blooddata, function(x) ifelse(x == "false", FALSE, x), how = "replace")

  class(blooddata) <- "blooddata"

  return(blooddata)
}


#' Create a blooddata object from BIDS data
#'
#' Deprecated. This function creates a blooddata object from JSON data structured according to
#' the old PET BIDS standard.
#'
#' @param bids_data The filename of a PET BIDS json sidecar, or a list
#'   containing the information contained within a PET BIDS json sidecar.
#' @param TimeShift The extent to which all the times in the data should be
#'   shifted (in seconds). Defaults to 0.
#'
#' @return A blooddata object
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' a <- create_blooddata_bids("bids_sidecar.json")
#' }
create_blooddata_bids <- function(bids_data, TimeShift = 0) {

  message( paste("Deprecation warning: this function is ",
                 "based on the old PET BIDS standard, and will be phased out",
                 "of future releases.") )

  if (!is.list(bids_data)) {
    if (file.exists(bids_data) &
      tools::file_ext(bids_data) == "json") {
      bids_data <- jsonlite::fromJSON(bids_data)
    }
  }

  tibblify_bidsjson <- function(list) {
    list$Data$Values <- tibble::as_tibble(list$Data$Values,
                                          .name_repair = "unique")
    colnames(list$Data$Values) <- list$Data$Labels
    return(list)
  }


  # Blood
  Blood <- list()

  # Continuous blood
  Blood$Continuous <- tibblify_bidsjson(bids_data$Blood$Continuous)

  # Discrete blood
  Blood$Discrete <- tibblify_bidsjson(bids_data$Blood$Discrete)

  # Plasma
  Plasma <- tibblify_bidsjson(bids_data$Plasma)

  # Metabolite
  Metabolite <- tibblify_bidsjson(bids_data$Metabolite)

  blooddata <- list(
    Data = list(
      Blood = Blood,
      Plasma = Plasma,
      Metabolite = Metabolite
    ),
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


#' Add a model fit to a blooddata object
#'
#' This function adds a fit object which can predict values to the blooddata
#' object
#'
#' @param blooddata A blooddata object created using one of the
#'   create_blooddata_* functions.
#' @param fit A model fit object, for which the predict function can be used for
#'   new data.
#' @param modeltype The function which the model predicts. One of the following:
#'   Blood for models of how the blood data should be described, BPR for models
#'   of the blood-to-plasma ratio, parentFraction for models of metabolism, and
#'   AIF for models of the arterial input function.
#'
#' @return A blooddata object with the model inserted.
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' blooddata <- bd_addfit(blooddata, hillfit, "parentFraction")
#' }
bd_addfit <- function(blooddata, fit, modeltype = c(
                           "Blood",
                           "BPR",
                           "parentFraction",
                           "AIF"
                         )) {

  # Verify fit object
  if (!(length(as.numeric(predict(fit))) > 1)) {
    stop("The fit object should be capable of producing
         predicted values with the predict command.")
  }

  if (length(modeltype) > 1) {
    stop("modeltype should be defined for the model as one of
         the following: Blood, BPR, parentFraction or AIF")
  }
  modeltype <- match.arg(modeltype, c(
    "Blood",
    "BPR",
    "parentFraction",
    "AIF"
  ))

  blooddata$Models[[modeltype]] <- list(Method = "fit", Data = fit)

  return(blooddata)
}

#' Add model fit parameters to a blooddata object
#'
#' This function adds the fitted parameters of a known model to the blooddata object.
#'
#' @param blooddata A blooddata object created using one of the create_blooddata_* functions.
#' @param modelname The name of the model function which produces fitted outcomes for given parameters.
#' @param fitpars The fitted parameters as a named list, data.frame or tibble.
#' @param modeltype The function which the model predicts. One of the following:
#'   Blood for models of how the blood data should be described, BPR for models
#'   of the blood-to-plasma ratio, parentFraction for models of metabolism, and
#'   AIF for models of the arterial input function.
#'
#' @return A blooddata object with the model inserted.
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' blooddata <- bd_addfitpars(blooddata, hillfit_model, hillfit$pars, "parentFraction")
#' }
bd_addfitpars <- function(blooddata, modelname, fitpars,
                             modeltype = c("Blood", "BPR", "parentFraction", "AIF")) {


  # Verify fit type
  if (!exists(modelname,
    mode = "function"
  )) {
    stop("The modelname should match the model name of
         the corresponding model function.")
  }

  if (length(modeltype) > 1) {
    stop("modeltype should be defined for the model as one of
         the following: Blood, BPR, parentFraction or AIF")
  }
  modeltype <- match.arg(modeltype, c(
    "Blood",
    "BPR",
    "parentFraction",
    "AIF"
  ))



  blooddata$Models[[modeltype]] <- list(
    Method = "fitpars",
    Data = list(
      Model = modelname,
      Pars = fitpars
    )
  )

  return(blooddata)
}


#' Add fitted values from another fit to a blooddata object
#'
#' This function adds the predicted values of an external function to a blooddata object. This allows one to use functions outside of kinfitr to estimate blood values.
#'
#' @param blooddata A blooddata object created using one of the create_blooddata_* functions.
#' @param time The times of the predictions at a sufficiently detailed level of interpolation.
#' @param predicted The predicted values for the given set of times.
#' @param modeltype The function which the model predicts. One of the following:
#'   Blood for models of how the blood data should be described, BPR for models
#'   of the blood-to-plasma ratio, parentFraction for models of metabolism, and
#'   AIF for models of the arterial input function.
#'
#' @return A blooddata object with the model fits inserted.
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' blooddata <- bd_addfitpars(
#'   blooddata, matlabout$time,
#'   matlabout$predicted, "AIF"
#' )
#' }
bd_addfitted <- function(blooddata, time, predicted,
                            modeltype = c("Blood", "BPR", "parentFraction", "AIF")) {
  if (length(modeltype) > 1) {
    stop("modeltype should be defined for the model as one of
         the following: Blood, BPR, parentFraction or AIF")
  }
  modeltype <- match.arg(modeltype, c(
    "Blood",
    "BPR",
    "parentFraction",
    "AIF"
  ))

  fitted <- tibble::tibble(time = time,
                           predicted = as.numeric(predicted) )

  blooddata$Models[[modeltype]] <- list(Method = "fitted", Data = fitted)

  return(blooddata)
}



# bd_getblood <- function(blooddata) {
#
#   blood_discrete <- blooddata$Data$Blood$Discrete$Data$Values
#
#   blood_discrete$time <- blood_discrete$sampleStartTime +
#     0.5*blood_discrete$sampleDuration
#
#   blood_discrete <- dplyr::arrange(blood_discrete, time)
#
#   blood_continuous <- blooddata$Data$Blood$Continuous$Data$Values
#
#   blood <- dplyr::bind_rows(blood_discrete, blood_continuous)
#   blood <- dplyr::arrange(blood, time)
#   blood$Method <- ifelse(is.na(blood$sampleDuration),
#                          yes = "Continuous",
#                          no = "Discrete")
#
#   return(blood)
#
# }
#
# bd_getplasma <- function(blooddata) {
#
#   plasma <- blooddata$Data$Plasma$Data$Values
#
#   plasma$time <- plasma$sampleStartTime +
#     0.5*plasma$sampleDuration
#
#   plasma <- dplyr::arrange(plasma, time)
#
#   return(plasma)
#
# }
#
# bd_getbpr <- function(blooddata) {
#
#   blood <- bd_getblood(blooddata)
#   plasma <- bd_getplasma(blooddata)
#
#   commonvalues <- intersect(plasma$time, blood_discrete$time)
#
#   bprvec <- blood_discrete$activity[blood_discrete$time %in% commonvalues] /
#     plasma$activity[plasma$time %in% commonvalues]
#
#
#   bpr <- tibble::tibble(time = commonvalues, bpr = bprvec)
#
#   return(bpr)
#
# }
#
# bd_getpf <- function(blooddata) {
#
#   pf <- blooddata$Data$Metabolite$Data$Values
#   pf$time <- pf$sampleStartTime + 0.5*pf$sampleDuration
#
#   return(pf)
#
# }
#
#
# bd_getaif <- function(blooddata) {
#
#   aif <- bd_getblood(blooddata)
#   bpr <- bd_getbpr(blooddata)
#
#   aif <- dplyr::rename(aif, blood = activity)
#   aif <- dplyr::mutate(aif,
#                        plasma = blood * (1/bpr),
#                        aif = plasma * parentFraction)
#
#   if(output == "AIF") {
#     return(aif)
#   }
#
# }









#' Tidy time information for blooddata operations
#'
#' Just a little tidier function.
#'
#' @param blooddata A blooddata object.
#' @param startTime The starting time for the interpolation. Defaults to zero. If, after application of the TimeShift value in the blooddata object, the startTime is still after zero, it will be set to zero.
#' @param stopTime The end time for the interpolation. Defaults to the maximum measured time.
#' @param interpPoints The number of points to interpolate over between the start and stop times. Defaults to 6000.
#'
#' @return The startTime, stopTime and interptimes after the relevant checks
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' bd_tidy_times(blooddata)
#' }
bd_tidy_times <- function(blooddata,
                              startTime = 0,
                              stopTime = NULL,
                              interpPoints = 6000) {

  # Check that it actually is a blooddata object, even if old
  if( class(blooddata) != "blooddata" ) {
    stop("The input data does not appear to be a blooddata object. Please make
           sure that the input data is correct.")
  }

  # Compatibility with old PET BIDS

  if(is.null(blooddata$Data$Plasma$Avail)) {
    stop("It appears that you are using a blooddata object created using an
         older version of kinfitr, based on an older version of the PET BIDS
         standard. Ideally, please create your blooddata objects afresh.
         Otherwise, please use the update_blooddata_bids command to keep only
         the bare essentials.")
  }


  if (is.null(stopTime)) {
    stopTime <- max(c(
        blooddata$Data$Blood$Discrete$Values$time,
        blooddata$Data$Plasma$Values$time,
        blooddata$Data$Metabolite$Values$time,
        blooddata$Data$Blood$Continuous$Values$time
    ),
    na.rm = T
    )
  }

  if (startTime > blooddata$TimeShift * (-1)) {
    startTime <- 0
  }

  interptime <- seq(startTime, stopTime, length.out = interpPoints)

  out <- list(startTime = startTime,
              stopTime = stopTime,
              interptime = interptime)

  return(out)

}

bd_extract_blood <- function(blooddata,
                             what = c("raw", "pred", "interp"),
                             startTime = 0,
                             stopTime = NULL,
                             interpPoints = 6000) {

  # Housekeeping
  what <- match.arg(what, c("raw", "pred", "interp"))

  tidied <- bd_tidy_times(blooddata,
                          startTime,
                          stopTime,
                          interpPoints)

  startTime <- tidied$startTime
  stopTime <- tidied$stopTime
  interptime <- tidied$interptime

  # Preparing

  blood_discrete <- blooddata$Data$Blood$Discrete$Values

  blood_discrete <- dplyr::filter(blood_discrete, !is.na(activity)) %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(Method = "Discrete")

  if(blooddata$Data$Blood$Continuous$Avail) {
    blood_continuous <- blooddata$Data$Blood$Continuous$Values %>%
      dplyr::filter(!is.na(activity)) %>%
      dplyr::mutate(Method = "Continuous")

    blood <- dplyr::bind_rows(blood_discrete, blood_continuous)
  } else {
    blood <- blood_discrete
  }

  blood <- dplyr::filter(blood, !is.na(activity)) %>%
    dplyr::arrange(blood, time)

  if (what == "raw") {
    return(blood)
  }

  ## Interp
  if (blooddata$Models$Blood$Method == "interp") {

    blood_for_interp <- blood

    ### Comment - for interpolation of blood data, I remove the discrete samples
    ### from the first half of the overlap between discrete and continuous
    ### samples. This is where the peak is, and where the timing of discrete
    ### samples can mess stuff up.

    if ( blooddata$Data$Blood$Continuous$Avail ) {
      if (max(blood_continuous$time) > min(blood_discrete$time)) {
        overlap_start <- min(blood_discrete$time)
        overlap_stop <- max(blood_continuous$time)
        overlap_time <- overlap_stop - overlap_start

        blood_for_interp$keep <- ifelse(
          blood_for_interp$Method == "Discrete" &
            blood_for_interp$time <
            (overlap_start + 0.5 * overlap_time),
          yes = FALSE, no = TRUE
        )

        blood_for_interp <- dplyr::filter(blood_for_interp, keep == TRUE)
      }
    }

    suppressWarnings(
      i_blood <- tibble::tibble(
        time = interptime,
        activity = interpends(blood_for_interp$time,
                              blood_for_interp$activity,
                              interptime,
                              method = "linear",
                              yzero = 0
        )
      )
    )
  }

  ## Fit
  if (blooddata$Models$Blood$Method == "fit") {
    i_blood <- tibble::tibble(
      time = interptime,
      activity = as.numeric(
        predict(blooddata$Models$Blood$Data,
                         newdata = list(time = interptime))
      )
    )

    blood$activity <- as.numeric(
      predict(blooddata$Models$Blood$Data,
                              newdata = list(time = blood$time))
    )
  }

  ## Fit pars
  if (blooddata$Models$Blood$Method == "fitpars") {
    modelname <- blooddata$Models$Blood$Data$Model
    Pars <- append(
      list(time = interptime),
      as.list(blooddata$Models$Blood$Data$Pars)
    )

    i_blood <- tibble::tibble(
      time = interptime,
      activity = do.call(
        what = modelname,
        args = Pars
      )
    )

    Pars <- append(
      list(time = blood$time),
      as.list(blooddata$Models$Blood$Data$Pars)
    )

    blood$activity <- do.call(
      what = modelname,
      args = Pars
    )
  }

  ## Fitted
  if (blooddata$Models$Blood$Method == "fitted") {
    i_blood <- tibble::tibble(
      time = interptime,
      activity = interpends(blooddata$Models$Blood$Data$time,
                            blooddata$Models$Blood$Data$predicted,
                            interptime,
                            method = "linear",
                            yzero = 0
      )
    )

    blood$activity <- interpends(blooddata$Models$Blood$Data$time,
                                 blooddata$Models$Blood$Data$predicted,
                                 blood$time,
                                 method = "linear"
    )
  }

  # Return

  if( what=="pred" ) {
    return(blood)
  }

  if( what=="interp" ) {
    return(i_blood)
  }
}





bd_extract_bpr <- function(blooddata,
                           what = c("raw", "pred", "interp"),
                           startTime = 0,
                           stopTime = NULL,
                           interpPoints = 6000,
                           bpr_peakfrac_cutoff = 0.001) {

  # Housekeeping
  what <- match.arg(what, c("raw", "pred", "interp"))

  tidied <- bd_tidy_times(blooddata,
                          startTime,
                          stopTime,
                          interpPoints)

  startTime <- tidied$startTime
  stopTime <- tidied$stopTime
  interptime <- tidied$interptime

  # Preparing

  plasma <- blooddata$Data$Plasma$Values %>%
    dplyr::filter(!is.na(activity)) %>%
    dplyr::arrange(time) %>%
    dplyr::filter(!duplicated(time)) %>%
    dplyr::filter(activity > (bpr_peakfrac_cutoff*max(activity)))

  blood <- bd_extract_blood(blooddata, startTime, stopTime, what = "pred")


  blood_discrete <- blood %>%
    dplyr::filter(!is.na(activity)) %>%
    dplyr::filter(Method=="Discrete") %>%
    dplyr::filter(!duplicated(time)) %>%
    dplyr::filter(activity > (bpr_peakfrac_cutoff*max(activity)))

  commonvalues <- intersect(plasma$time, blood_discrete$time)

  bprvec <- blood_discrete$activity[blood_discrete$time %in% commonvalues] /
    plasma$activity[plasma$time %in% commonvalues]

  if( is.nan(bprvec[1]) ) {
    bprvec[1] <- bprvec[2]
  }


  bpr <- tibble::tibble(time = commonvalues, bpr = bprvec) %>%
    dplyr::filter(!is.infinite(bpr))

  # # Remove massive outliers
  # bpr_median <- median(bpr$bpr)
  # bpr_mad <- mad(bpr$bpr)
  # bpr_min <- bpr_median - 100*bpr_mad
  # bpr_max <- bpr_median + 100*bpr_mad
  #
  # if( any(bpr$bpr > bpr_max) || any(bpr$bpr < bpr_min) ) {
  #
  #   warning("BPR outlier(s) excluded: >100*MAD")
  #
  #   bpr <- dplyr::filter(bpr, bpr < bpr_max)
  #   bpr <- dplyr::filter(bpr, bpr > bpr_min)
  # }

  if (what == "raw") {
    return(bpr)
  }



  ## Interp
  if (blooddata$Models$BPR$Method == "interp") {
    i_bpr <- tibble::tibble(
      time = interptime,
      bpr = interpends(bpr$time,
        bpr$bpr,
        interptime,
        method = "linear"
      )
    )
    blood$bpr <- interpends(bpr$time, bpr$bpr, blood$time,
      method = "linear"
    )
  }

  ## Fit
  if (blooddata$Models$BPR$Method == "fit") {
    i_bpr <- tibble::tibble(
      time = interptime,
      bpr = as.numeric(
        predict(blooddata$Models$BPR$Data,
        newdata = list(time = interptime))
      )
    )

    blood$bpr <- as.numeric(
      predict(blooddata$Models$BPR$Data,
      newdata = list(time = blood$time))
    )
  }

  ## Fit pars
  if (blooddata$Models$BPR$Method == "fitpars") {
    modelname <- blooddata$Models$BPR$Data$Model
    Pars <- append(
      list(time = interptime),
      as.list(blooddata$Models$BPR$Data$Pars)
    )

    i_bpr <- tibble::tibble(
      time = interptime,
      activity = do.call(
        what = modelname,
        args = Pars
      )
    )
    Pars <- append(
      list(time = blood$time),
      as.list(blooddata$Models$BPR$Data$Pars)
    )

    blood$bpr <- do.call(
      what = modelname,
      args = Pars
    )
  }

  ## Fitted
  if (blooddata$Models$BPR$Method == "fitted") {
    i_bpr <- tibble::tibble(
      time = interptime,
      bpr = interpends(blooddata$Models$BPR$Data$time,
        blooddata$Models$BPR$Data$predicted,
        interptime,
        method = "linear"
      )
    )

    blood$bpr <- interpends(blooddata$Models$BPR$Data$time,
      blooddata$Models$BPR$Data$predicted,
      blood$time,
      method = "linear"
    )
  }

  # Return

  if( what=="pred" ) {
    return(blood)
  }

  if( what=="interp" ) {
    return(i_bpr)
  }

}

bd_extract_pf <- function(blooddata,
                          what = c("raw", "pred", "interp"),
                           startTime = 0,
                           stopTime = NULL,
                           interpPoints = 6000) {

  # Housekeeping
  what <- match.arg(what, c("raw", "pred", "interp"))

  tidied <- bd_tidy_times(blooddata,
                          startTime,
                          stopTime,
                          interpPoints)

  startTime <- tidied$startTime
  stopTime <- tidied$stopTime
  interptime <- tidied$interptime

  # Preparing

  pf <- blooddata$Data$Metabolite$Values %>%
    dplyr::filter(!is.na(parentFraction))

  if (what == "raw") {
    return(pf)
  }

  # For pred and interp, we want interpolation to start at 1 at time 0
  if(!(0 %in% pf$time)) {
    pf <- rbind(c(0,1),
                pf)
  }

  blood <- bd_extract_bpr(blooddata, startTime, stopTime, what = "pred")

  ## Interp
  if (blooddata$Models$parentFraction$Method == "interp") {
    i_pf <- tibble::tibble(
      time = interptime,
      parentFraction = interpends(pf$time,
        pf$parentFraction,
        interptime,
        method = "linear",
        yzero = 1
      )
    )
    blood$parentFraction <- interpends(pf$time, pf$parentFraction,
      blood$time,
      method = "linear"
    )
  }

  ## Fit
  if (blooddata$Models$parentFraction$Method == "fit") {
    i_pf <- tibble::tibble(
      time = interptime,
      parentFraction = as.numeric(
        predict(blooddata$Models$parentFraction$Data,
        newdata = list(time = interptime))
      )
    )

    blood$parentFraction <- as.numeric(
      predict(blooddata$Models$parentFraction$Data,
      newdata = list(time = blood$time))
    )
  }

  ## Fit pars
  if (blooddata$Models$parentFraction$Method == "fitpars") {
    modelname <- blooddata$Models$parentFraction$Data$Model
    Pars <- append(
      list(time = interptime),
      as.list(blooddata$Models$parentFraction$Data$Pars)
    )

    i_pf <- tibble::tibble(
      time = interptime,
      parentFraction = do.call(
        what = modelname,
        args = Pars
      )
    )

    Pars <- append(
      list(time = blood$time),
      as.list(blooddata$Models$parentFraction$Data$Pars)
    )

    blood$parentFraction <- do.call(
      what = modelname,
      args = Pars
    )
  }

  ## Fitted
  if (blooddata$Models$parentFraction$Method == "fitted") {
    i_pf <- tibble::tibble(
      time = interptime,
      parentFraction = interpends(
        blooddata$Models$parentFraction$Data$time,
        blooddata$Models$parentFraction$Data$predicted,
        interptime,
        method = "linear",
        yzero = 1
      )
    )

    blood$parentFraction <- interpends(
      blooddata$Models$parentFraction$Data$time,
      blooddata$Models$parentFraction$Data$predicted,
      blood$time,
      method = "linear",
      yzero = 1
    )
  }

  # Return

  if( what=="pred" ) {
    return(blood)
  }

  if( what=="interp" ) {
    return(i_pf)
  }

}

bd_extract_aif <- function(blooddata,
                           what = c("raw", "pred", "interp"),
                          startTime = 0,
                          stopTime = NULL,
                          interpPoints = 6000) {

  # Housekeeping
  what <- match.arg(what, c("raw", "pred", "interp"))

  tidied <- bd_tidy_times(blooddata,
                          startTime,
                          stopTime,
                          interpPoints)

  startTime <- tidied$startTime
  stopTime <- tidied$stopTime
  interptime <- tidied$interptime

  # Preparing

  ### Note: here, the plasm_uncor is generated from the blood and bpr
  ### curves. But, if we've fitted the bpr, we're actually not using
  ### the true raw plasma data when possible.
  aif <- bd_extract_pf(blooddata, startTime, stopTime, what = "pred") %>%
    dplyr::rename(blood = activity) %>%
    dplyr::mutate(plasma_uncor = blood * (1 / bpr),
                  aif = plasma_uncor * parentFraction
    )

  ### Here the raw plasma data is extracted to replace the calculated
  ### values where raw plasma measurements were actually taken.
  bpr <- bd_extract_bpr(blooddata, startTime, stopTime, what = "raw")
  blood <- bd_extract_blood(blooddata, startTime, stopTime, what = "raw") %>%
    dplyr::filter(Method=="Discrete")

  rawplasma <- dplyr::left_join(bpr, blood, by="time") %>%
    dplyr::mutate(measplasma = activity * (1 / bpr)) %>%
    dplyr::select(time, measplasma, Method)

  ### Now the calculated raw plasma values are replaced with the measured ones.
  aif <- dplyr::left_join(aif, rawplasma, by=c("time", "Method")) %>%
    dplyr::mutate(plasma_uncor = ifelse( is.na(measplasma),
                                         yes = plasma_uncor,
                                         no = measplasma)) %>%
    dplyr::select(-measplasma) %>%
    dplyr::mutate(aif = plasma_uncor * parentFraction)

  if (what == "raw") {
    aif <- dplyr::mutate(aif, aif = plasma_uncor * parentFraction) # raw values
    return(aif)
  }


  ## Interp
  if (blooddata$Models$AIF$Method == "interp") {
    aif_for_interp <- aif

    aif_discrete <- dplyr::filter(aif, Method=="Discrete")

    aif_continuous <- dplyr::filter(aif, Method=="Continuous")

    ### Comment - for interpolation of blood data, I remove the discrete samples
    ### from the first half of the overlap between discrete and continuous
    ### samples. This is where the peak is, and where the timing of discrete
    ### samples can mess stuff up.

    if (nrow(aif_continuous) > 0) {
      if (max(aif_continuous$time) > min(aif_discrete$time)) {
        overlap_start <- min(aif_discrete$time)
        overlap_stop <- max(aif_continuous$time)
        overlap_time <- overlap_stop - overlap_start

        aif_for_interp$keep <- ifelse(
          aif_for_interp$Method == "Discrete" &
            aif_for_interp$time <
              (overlap_start + 0.5 * overlap_time),
          yes = FALSE, no = TRUE
        )

        aif_for_interp <- dplyr::filter(aif_for_interp, keep == TRUE)
      }
    }

    suppressWarnings(
      i_aif <- tibble::tibble(
        time = interptime,
        aif = interpends(aif_for_interp$time,
          aif_for_interp$aif,
          interptime,
          method = "linear",
          yzero = 0
        )
      )
    )
  }

  ## Fit
  if (blooddata$Models$AIF$Method == "fit") {
    i_aif <- tibble::tibble(
      time = interptime,
      aif = as.numeric(
        predict(blooddata$Models$AIF$Data,
        newdata = list(time = interptime))
      )
    )

    aif$aif = as.numeric(
      predict(blooddata$Models$AIF$Data,
              newdata = list(time = aif$time)))
  }

  ## Fit pars
  if (blooddata$Models$AIF$Method == "fitpars") {
    modelname <- blooddata$Models$AIF$Data$Model


    # Interp
    Pars <- append(
      list(time = interptime),
      as.list(blooddata$Models$AIF$Data$Pars)
    )

    i_aif <- tibble::tibble(
      time = interptime,
      activity = do.call(
        what = modelname,
        args = Pars
      )
    )

    # Original times
    Pars <- append(
      list(time = aif$time),
      as.list(blooddata$Models$AIF$Data$Pars)
    )

    aif$aif <- do.call(
        what = modelname,
        args = Pars)
  }

  ## Fitted
  if (blooddata$Models$AIF$Method == "fitted") {
    i_aif <- tibble::tibble(
      time = interptime,
      aif = interpends(
        blooddata$Models$AIF$Data$time,
        blooddata$Models$AIF$Data$predicted,
        interptime,
        method = "linear",
        yzero = 0
      )
    )

    aif$aif <- interpends(
      blooddata$Models$AIF$Data$time,
      blooddata$Models$AIF$Data$predicted,
      aif$time,
      method = "linear",
      yzero = 0
    )
  }

  # Return

  if( what=="pred" ) {
    return(aif)
  }

  if( what=="interp" ) {
    return(i_aif)
  }

}

#' Create an input object from a blooddata object
#'
#' Get an input object for kinetic modelling
#'
#'
#' @param blooddata The blooddata object.
#' @param startTime The starting time for the interpolation. Defaults to zero.
#'   If, after application of the TimeShift value in the blooddata object, the
#'   startTime is still after zero, it will be set to zero.
#' @param stopTime The end time for the interpolation. Defaults to the maximum
#'   measured time.
#' @param interpPoints The number of points to interpolate over between the
#'   start and stop times. Defaults to 6000.
#'
#' @return An input object
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' bd_create_input(blooddata)
#' }
bd_create_input <- function(blooddata,
                             startTime = 0,
                             stopTime = NULL,
                             interpPoints = 6000) {

  # Housekeeping
  tidied <- bd_tidy_times(blooddata,
                          startTime,
                          stopTime,
                          interpPoints)

  startTime <- tidied$startTime
  stopTime <- tidied$stopTime
  interptime <- tidied$interptime

  # Preparing

  i_blood <- bd_extract_blood(blooddata, what = "interp",
                              startTime,stopTime,interpPoints)

  i_bpr <- bd_extract_bpr(blooddata, what = "interp",
                          startTime,stopTime,interpPoints)

  i_pf <- bd_extract_pf(blooddata, what = "interp",
                        startTime,stopTime,interpPoints)

  i_aif <- suppressMessages(
    bd_extract_aif(blooddata, what = "interp",
                   startTime,stopTime,interpPoints)
  )

  input <- tibble::tibble(
    Time = (interptime + blooddata$TimeShift) / 60,
    Blood = i_blood$activity,
    Plasma = i_blood$activity / i_bpr$bpr,
    ParentFraction = i_pf$parentFraction,
    AIF = i_aif$aif
  )

  class(input) <- c("interpblood", class(input))

  return(input)
}


#' Extract data from blooddata objects
#'
#' Pull out data for modelling purposes.
#'
#' @param blooddata The blooddata object.
#' @param output What curve should be extracted. It can be "Blood", "BPR"
#'   (whole-blood-to-plasma-ratio), "parentFraction" or "AIF".
#' @param startTime The starting time for the interpolation. Defaults to zero.
#'   If, after application of the TimeShift value in the blooddata object, the
#'   startTime is still after zero, it will be set to zero.
#' @param stopTime The end time for the interpolation. Defaults to the maximum
#'   measured time.
#' @param interpPoints The number of points to interpolate over between the
#'   start and stop times. Defaults to 6000.
#' @param what What kind of data should be extracted? The "raw" data is the
#'   unaltered data for modelling, the "pred" data is the predicted data at the
#'   times of the original samples for assessing correspondence, and the
#'   "interp" values are the values interpolated through the interptimes.
#' @param bpr_peakfrac_cutoff What is the lowest permissible value of the blood
#'   or plasma curves as a fraction of the peak that should be included in the
#'   raw BPR extraction for modelling. This is to avoid ratios which tend
#'   towards 0 or Inf at early time points when blood and plasma are very low.
#'
#' @return A tibble with the relevant data
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' bd_extract(blooddata, output="parentFraction", what="raw")
#' }
bd_extract <- function(blooddata,
                       output = c("Blood", "BPR", "parentFraction", "AIF"),
                       startTime = 0,
                       stopTime = NULL,
                       interpPoints = 6000,
                       what = c("raw", "pred", "interp"),
                       bpr_peakfrac_cutoff = 0.001) {

  # Housekeeping
  tidied <- bd_tidy_times(blooddata,
                          startTime,
                          stopTime,
                          interpPoints)

  startTime <- tidied$startTime
  stopTime <- tidied$stopTime
  interptime <- tidied$interptime

  output <- match.arg(output, c("Blood", "BPR", "parentFraction", "AIF"))
  what <- match.arg(what, c("raw", "pred", "interp"))

  if(output=="Blood") {
    out <- bd_extract_blood(blooddata, startTime,stopTime,interpPoints,
                            what = what)
  } else if(output=="BPR") {
    out <- bd_extract_bpr(blooddata, startTime,stopTime,interpPoints,
                            what = what, bpr_peakfrac_cutoff)
  } else if(output=="parentFraction") {
    out <- bd_extract_pf(blooddata, startTime,stopTime,interpPoints,
                            what = what)
  } else if(output=="AIF") {
    out <- bd_extract_aif(blooddata, startTime,stopTime,interpPoints,
                            what = what)
  }

  return(out)

}

#' Create an input object from a blooddata object.
#'
#' DEPRECATION WARNING: this function will be slowly phased out of future
#' releases as it is based on the old PET BIDS standard. in favour of
#' bd_extract() and bd_create_input(). blooddata objects can be updated to the
#' new format using update_blooddata(), or (better) by creating it afresh from
#' the source.
#'
#' This function extracts data from blooddata objects to create either an
#' input object for kinetic modelling, or sets of values for modelling of
#' blood-related curves.
#'
#' @param blooddata A blooddata object, to which all the desired fits have been applied and added.
#' @param startTime The starting time for the interpolation. Defaults to zero. If, after application of the TimeShift value in the blooddata object, the startTime is still after zero, it will be set to zero.
#' @param stopTime The end time for the interpolation. Defaults to the maximum measured time.
#' @param interpPoints The number of points to interpolate over between the start and stop times.
#' @param output The output. This defaults to an "input" object, which can be used in a kinetic model fit. But if set to "Blood", "BPR", "parentFraction" or "AIF", it yields the appropriate input for the function which will be used to model these.
#'
#' @return A tibble containing the output specified.
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' bd_getdata(blooddata)
#' bd_getdata(blooddata, output = "parentFraction")
#' }
bd_getdata <- function(blooddata,
                         startTime = 0,
                         stopTime = NULL,
                         interpPoints = 6000,
                         output = c(
                           "input",
                           "Blood",
                           "BPR",
                           "parentFraction",
                           "AIF"
                         )) {

  message("DEPRECATION WARNING: based on recent changes to the PET BIDS standard,
          this function will slowly be deprecated, in favour of bd_extract() and
          bd_create_input(). blooddata objects can be updated to the new format
          using update_blooddata(), or (better) by creating it afresh from
          the source.")


  if (is.null(stopTime)) {
    stopTime <- max(c(
      with(
        blooddata$Data$Blood$Discrete$Data$Values,
        sampleStartTime + sampleDuration
      ),
      with(
        blooddata$Data$Plasma$Data$Values,
        sampleStartTime + sampleDuration
      ),
      with(
        blooddata$Data$Metabolite$Data$Values,
        sampleStartTime + sampleDuration
      ),
      blooddata$Data$Blood$Continuous$Data$Values$time
    ),
    na.rm = T
    )
  }

  if (startTime > blooddata$TimeShift * (-1)) {
    startTime <- 0
  }

  interptime <- seq(startTime, stopTime, length.out = interpPoints)

  output <- match.arg(output, c(
    "input",
    "Blood",
    "BPR",
    "parentFraction",
    "AIF"
  ))

  # Blood

  blood_discrete <- blooddata$Data$Blood$Discrete$Data$Values
  blood_discrete$time <- blood_discrete$sampleStartTime +
    0.5 * blood_discrete$sampleDuration
  blood_discrete <- dplyr::filter(blood_discrete, !is.na(activity))
  blood_discrete <- dplyr::arrange(blood_discrete, time)


  blood_continuous <- blooddata$Data$Blood$Continuous$Data$Values
  blood_continuous <- dplyr::filter(blood_continuous, !is.na(activity))

  blood <- dplyr::bind_rows(blood_discrete, blood_continuous)
  blood <- dplyr::filter(blood, !is.na(activity))
  blood <- dplyr::arrange(blood, time)
  blood$Method <- ifelse(is.na(blood$sampleDuration),
                         yes = "Continuous",
                         no = "Discrete"
  )

  if (output == "Blood") {
    return(blood)
  }

  ## Interp
  if (blooddata$Models$Blood$Method == "interp") {
    blood_for_interp <- blood

    ### Comment - for interpolation of blood data, I remove the discrete samples
    ### from the first half of the overlap between discrete and continuous
    ### samples. This is where the peak is, and where the timing of discrete
    ### samples can mess stuff up.

    if (nrow(blood_continuous) > 0) {
      if (max(blood_continuous$time) > min(blood_discrete$time)) {
        overlap_start <- min(blood_discrete$time)
        overlap_stop <- max(blood_continuous$time)
        overlap_time <- overlap_stop - overlap_start

        blood_for_interp$keep <- ifelse(
          blood_for_interp$Method == "Discrete" &
            blood_for_interp$time <
            (overlap_start + 0.5 * overlap_time),
          yes = FALSE, no = TRUE
        )

        blood_for_interp <- dplyr::filter(blood_for_interp, keep == TRUE)
      }
    }

    suppressWarnings(
      i_blood <- tibble::tibble(
        time = interptime,
        activity = interpends(blood_for_interp$time,
                              blood_for_interp$activity,
                              interptime,
                              method = "linear",
                              yzero = 0
        )
      )
    )
  }

  ## Fit
  if (blooddata$Models$Blood$Method == "fit") {
    i_blood <- tibble::tibble(
      time = interptime,
      activity = as.numeric(
        predict(blooddata$Models$Blood$Data,
                         newdata = list(time = interptime))
      )
    )

    blood$activity <- as.numeric(
      predict(blooddata$Models$Blood$Data,
                              newdata = list(time = blood$time))
    )
  }

  ## Fit pars
  if (blooddata$Models$Blood$Method == "fitpars") {
    modelname <- blooddata$Models$Blood$Data$Model
    Pars <- append(
      list(time = interptime),
      as.list(blooddata$Models$Blood$Data$Pars)
    )

    i_blood <- tibble::tibble(
      time = interptime,
      activity = do.call(
        what = modelname,
        args = Pars
      )
    )

    Pars <- append(
      list(time = blood$time),
      as.list(blooddata$Models$Blood$Data$Pars)
    )

    blood$activity <- do.call(
      what = modelname,
      args = Pars
    )
  }

  ## Fitted
  if (blooddata$Models$Blood$Method == "fitted") {
    i_blood <- tibble::tibble(
      time = interptime,
      activity = interpends(blooddata$Models$Blood$Data$time,
                            blooddata$Models$Blood$Data$predicted,
                            interptime,
                            method = "linear",
                            yzero = 0
      )
    )

    blood$activity <- interpends(blooddata$Models$Blood$Data$time,
                                 blooddata$Models$Blood$Data$predicted,
                                 blood$time,
                                 method = "linear"
    )
  }

  out <- list()

  # Blood-to-Plasma Ratio

  plasma <- blooddata$Data$Plasma$Data$Values
  plasma <- dplyr::filter(plasma, !is.na(activity))

  plasma$time <- plasma$sampleStartTime +
    0.5 * plasma$sampleDuration
  plasma <- dplyr::arrange(plasma, time)

  commonvalues <- intersect(plasma$time, blood_discrete$time)

  bprvec <- blood_discrete$activity[blood_discrete$time %in% commonvalues] /
    plasma$activity[plasma$time %in% commonvalues]

  if( is.nan(bprvec[1]) ) {
    bprvec[1] <- bprvec[2]
  }

  bpr <- tibble::tibble(time = commonvalues, bpr = bprvec)

  if (output == "BPR") {
    return(bpr)
  }

  ## Interp
  if (blooddata$Models$BPR$Method == "interp") {
    i_bpr <- tibble::tibble(
      time = interptime,
      bpr = interpends(bpr$time,
                       bpr$bpr,
                       interptime,
                       method = "linear"
      )
    )
    blood$bpr <- interpends(bpr$time, bpr$bpr, blood$time,
                            method = "linear"
    )
  }

  ## Fit
  if (blooddata$Models$BPR$Method == "fit") {
    i_bpr <- tibble::tibble(
      time = interptime,
      bpr = as.numeric(predict(blooddata$Models$BPR$Data,
                    newdata = list(time = interptime))
      )
    )

    blood$bpr <- as.numeric(predict(blooddata$Models$BPR$Data,
                         newdata = list(time = blood$time))
    )
  }

  ## Fit pars
  if (blooddata$Models$BPR$Method == "fitpars") {
    modelname <- blooddata$Models$BPR$Data$Model
    Pars <- append(
      list(time = interptime),
      as.list(blooddata$Models$BPR$Data$Pars)
    )

    i_bpr <- tibble::tibble(
      time = interptime,
      activity = do.call(
        what = modelname,
        args = Pars
      )
    )
    Pars <- append(
      list(time = blood$time),
      as.list(blooddata$Models$BPR$Data$Pars)
    )

    blood$bpr <- do.call(
      what = modelname,
      args = Pars
    )
  }

  ## Fitted
  if (blooddata$Models$BPR$Method == "fitted") {
    i_bpr <- tibble::tibble(
      time = interptime,
      bpr = interpends(blooddata$Models$BPR$Data$time,
                       blooddata$Models$BPR$Data$predicted,
                       interptime,
                       method = "linear"
      )
    )

    blood$bpr <- interpends(blooddata$Models$BPR$Data$time,
                            blooddata$Models$BPR$Data$predicted,
                            blood$time,
                            method = "linear"
    )
  }


  # Parent Fraction

  pf <- blooddata$Data$Metabolite$Data$Values
  pf <- dplyr::filter(pf, !is.na(parentFraction))
  pf$time <- pf$sampleStartTime + 0.5 * pf$sampleDuration

  if (output == "parentFraction") {
    return(pf)
  }

  ## Interp
  if (blooddata$Models$parentFraction$Method == "interp") {
    i_pf <- tibble::tibble(
      time = interptime,
      parentFraction = interpends(pf$time,
                                  pf$parentFraction,
                                  interptime,
                                  method = "linear",
                                  yzero = 1
      )
    )
    blood$parentFraction <- interpends(pf$time, pf$parentFraction,
                                       blood$time,
                                       method = "linear"
    )
  }

  ## Fit
  if (blooddata$Models$parentFraction$Method == "fit") {
    i_pf <- tibble::tibble(
      time = interptime,
      parentFraction = as.numeric(
        predict(blooddata$Models$parentFraction$Data,
                               newdata = list(time = interptime))
      )
    )

    blood$parentFraction <- as.numeric(
      predict(blooddata$Models$parentFraction$Data,
                                    newdata = list(time = blood$time))
    )
  }

  ## Fit pars
  if (blooddata$Models$parentFraction$Method == "fitpars") {
    modelname <- blooddata$Models$parentFraction$Data$Model
    Pars <- append(
      list(time = interptime),
      as.list(blooddata$Models$parentFraction$Data$Pars)
    )

    i_pf <- tibble::tibble(
      time = interptime,
      parentFraction = do.call(
        what = modelname,
        args = Pars
      )
    )

    Pars <- append(
      list(time = blood$time),
      as.list(blooddata$Models$parentFraction$Data$Pars)
    )

    blood$parentFraction <- do.call(
      what = modelname,
      args = Pars
    )
  }

  ## Fitted
  if (blooddata$Models$parentFraction$Method == "fitted") {
    i_pf <- tibble::tibble(
      time = interptime,
      parentFraction = interpends(
        blooddata$Models$parentFraction$Data$time,
        blooddata$Models$parentFraction$Data$predicted,
        interptime,
        method = "linear",
        yzero = 1
      )
    )

    blood$parentFraction <- interpends(
      blooddata$Models$parentFraction$Data$time,
      blooddata$Models$parentFraction$Data$predicted,
      blood$time,
      method = "linear",
      yzero = 1
    )
  }

  # AIF

  aif <- blood
  aif <- dplyr::rename(aif, blood = activity)
  aif <- dplyr::mutate(aif,
                       plasma_uncor = blood * (1 / bpr),
                       aif = plasma_uncor * parentFraction
  )

  if (output == "AIF") {
    return(aif)
  }


  ## Interp
  if (blooddata$Models$AIF$Method == "interp") {
    aif_for_interp <- aif

    ### Comment - for interpolation of blood data, I remove the discrete samples
    ### from the first half of the overlap between discrete and continuous
    ### samples. This is where the peak is, and where the timing of discrete
    ### samples can mess stuff up.

    if (nrow(blood_continuous) > 0) {
      if (max(blood_continuous$time) > min(blood_discrete$time)) {
        overlap_start <- min(blood_discrete$time)
        overlap_stop <- max(blood_continuous$time)
        overlap_time <- overlap_stop - overlap_start

        aif_for_interp$keep <- ifelse(
          aif_for_interp$Method == "Discrete" &
            aif_for_interp$time <
            (overlap_start + 0.5 * overlap_time),
          yes = FALSE, no = TRUE
        )

        aif_for_interp <- dplyr::filter(aif_for_interp, keep == TRUE)
      }
    }

    suppressWarnings(
      i_aif <- tibble::tibble(
        time = interptime,
        aif = interpends(aif_for_interp$time,
                         aif_for_interp$aif,
                         interptime,
                         method = "linear",
                         yzero = 0
        )
      )
    )
  }

  ## Fit
  if (blooddata$Models$AIF$Method == "fit") {
    i_aif <- tibble::tibble(
      time = interptime,
      aif = as.numeric(
        predict(blooddata$Models$AIF$Data,
                    newdata = list(time = interptime))
      )
    )
  }

  ## Fit pars
  if (blooddata$Models$AIF$Method == "fitpars") {
    modelname <- blooddata$Models$AIF$Data$Model
    Pars <- append(
      list(time = interptime),
      as.list(blooddata$Models$AIF$Data$Pars)
    )

    i_aif <- tibble::tibble(
      time = interptime,
      activity = do.call(
        what = modelname,
        args = Pars
      )
    )
  }

  ## Fitted
  if (blooddata$Models$AIF$Method == "fitted") {
    i_aif <- tibble::tibble(
      time = interptime,
      aif = interpends(
        blooddata$Models$AIF$Data$time,
        blooddata$Models$AIF$Data$predicted,
        interptime,
        method = "linear",
        yzero = 0
      )
    )
  }

  if (output == "input") {
    input <- tibble::tibble(
      Time = (interptime + blooddata$TimeShift) / 60,
      Blood = i_blood$activity,
      Plasma = i_blood$activity / i_bpr$bpr,
      ParentFraction = i_pf$parentFraction,
      AIF = i_aif$aif
    )

    class(input) <- c("interpblood", class(input))

    return(input)
  }
}


#' Update an old blooddata object to the new structure
#'
#' This is a convenience function for updating old blooddata to a new structure.
#' This is not the ideal way of doing this: better would be to re-create a
#' blooddata object from the data saved according to the new BIDS specification.
#' This way is faster, but keeps only the vital information.
#'
#' @param blooddata An old blooddata object in the old structure based on the
#' old PET BIDS specification
#'
#' @return A blooddata object in the new structure based on the newer PET BIDS
#' specification
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' update_blooddata(blooddata)
#' }
update_blooddata <- function(blooddata) {


  # Discrete Blood


  if( !is.null(blooddata$Data$Blood$Discrete$Data) ) {

    blooddata$Data$Blood$Discrete <-
      blooddata$Data$Blood$Discrete$Data

    blooddata$Data$Blood$Discrete$Avail = TRUE

    blooddata$Data$Blood$Discrete$time <- list(
      Description = "Transferred from old blooddata",
      Units = blooddata$Data$Blood$Discrete$units[1])

    blooddata$Data$Blood$Discrete$activity <- list(
      Description = "Transferred from old blooddata",
      Units = blooddata$Data$Blood$Discrete$units[3])

    blooddata$Data$Blood$Discrete$Values <-
      dplyr::mutate(blooddata$Data$Blood$Discrete$Values,
                    time = sampleStartTime + 0.5*sampleDuration) %>%
      dplyr::select(time, activity)
  } else {

    if(is.null(blooddata$Data$Blood$Discrete$Avail)) {
      blooddata$Data$Blood$Discrete$Avail <- FALSE
    }

  }

  # Plasma

  if( !is.null(blooddata$Data$Plasma$Data) ) {

    blooddata$Data$Plasma <-
      blooddata$Data$Plasma$Data

    blooddata$Data$Plasma$Avail = TRUE

    blooddata$Data$Plasma$time <- list(
      Description = "Transferred from old blooddata",
      Units = blooddata$Data$Plasma$units[1])

    blooddata$Data$Plasma$activity <- list(
      Description = "Transferred from old blooddata",
      Units = blooddata$Data$Plasma$units[3])

    blooddata$Data$Plasma$Values <-
      dplyr::mutate(blooddata$Data$Plasma$Values,
                    time = sampleStartTime + 0.5*sampleDuration) %>%
      dplyr::select(time, activity)
  } else {

    if(is.null(blooddata$Data$Plasma$Avail)) {
      blooddata$Data$Plasma$Avail <- FALSE
    }

  }

  # Parent Fraction

  if( !is.null(blooddata$Data$Metabolite$Data)) {

    blooddata$Data$Metabolite <-
      blooddata$Data$Metabolite$Data

    blooddata$Data$Metabolite$Avail = TRUE

    blooddata$Data$Metabolite$time <- list(
      Description = "Transferred from old blooddata",
      Units = blooddata$Data$Metabolite$units[1])

    blooddata$Data$Metabolite$activity <- list(
      Description = "Transferred from old blooddata",
      parentFraction = "unitless")

    blooddata$Data$Metabolite$Values <-
      dplyr::mutate(blooddata$Data$Metabolite$Values,
                    time = sampleStartTime + 0.5*sampleDuration) %>%
      dplyr::select(time, parentFraction)
  } else {

    if(is.null(blooddata$Data$Metabolite$Avail)) {
      blooddata$Data$Metabolite$Avail <- FALSE
    }

  }

  # Continuous Blood

  if( !is.null(blooddata$Data$Blood$Continuous$Data)) {

    blooddata$Data$Blood$Continuous <- c(
      blooddata$Data$Blood$Continuous$Data,
      blooddata$Data$Blood$Continuous[-which(
        names(blooddata$Data$Blood$Continuous)=="Data")])

    blooddata$Data$Blood$Continuous$Avail = TRUE

    blooddata$Data$Blood$Continuous$time <- list(
      Description = "Transferred from old blooddata",
      Units = blooddata$Data$Blood$Continuous$units[1])

    blooddata$Data$Blood$Continuous$activity <- list(
      Description = "Transferred from old blooddata",
      Units = blooddata$Data$Blood$Continuous$units[2])

  } else {

    if(is.null(blooddata$Data$Blood$Continuous$Avail)) {
      blooddata$Data$Blood$Continuous$Avail <- FALSE
    }

  }

  return(blooddata)


}

#' Perform Dispersion Correction on a blooddata object
#'
#' Dispersion corrected is performed on the continuous blood data if it has not
#' already been performed.
#'
#' @param blooddata A blooddata object created using one of the
#'   create_blooddata_* functions.
#' @param timedelta The time difference between each measured sample. Defaults
#'   to the most common time difference between the first 20 measurements.
#' @param keep_interpolated Defaults to TRUE: should interpolated samples which
#'   were not included in the original input be included in the output.
#' @param smooth_iterations The number of times that the smoothing of each pair
#'   of observations should be performed using blood_smooth. Defaults to 0.
#' @param tau The dispersion constant in unit seconds. If NULL (default),
#'   this value will be extracted from the constant in the blooddata object.
#'   This parameter is included as a manual override.
#'
#' @return A blooddata object after dispersion correction.
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' blooddata <- bd_blood_dispcor(blooddata)
#' }
#'
bd_blood_dispcor <- function(blooddata, timedelta = NULL,
                              keep_interpolated = T,
                              smooth_iterations = 0,
                              tau = NULL) {

  if (is.null(blooddata$Data$Blood$Continuous)) {

    stop("There is no continuous blood data to
         perform dispersion correction on.")
  }

  if (!blooddata$Data$Blood$Continuous$DispersionCorrected) { # not corrected

    if( is.null(blooddata$Data$Blood$Continuous$Avail) ) { # old blooddata
      blooddata <- update_blooddata(blooddata)
      message("
         It appears that you are using a blooddata object created using an
         older version of kinfitr, based on an older version of the PET BIDS
         standard. The blooddata object has been updated automatically. Ideally,
         please create your blooddata objects afresh from the raw data.
         Otherwise, you can use the update_blooddata_bids() command to keep only
         the bare essentials.")
    }

    blood_continuous <- blooddata$Data$Blood$Continuous$Values

    if(is.null(tau)) {
      tau = blooddata$Data$Blood$Continuous$DispersionConstant
    }

    blood_dc <- blood_dispcor(blood_continuous$time,
      blood_continuous$activity,
      tau = tau,
      timedelta, keep_interpolated = T,
      smooth_iterations = smooth_iterations
    )

    blooddata$Data$Blood$Continuous$Values <- blood_dc
    blooddata$Data$Blood$Continuous$DispersionCorrected <- TRUE
  }

  return(blooddata)
}

#' Plot a blooddata object
#'
#' @param blooddata A blooddata object, to which all the desired fits have been applied and added.
#' @param startTime The starting time for the interpolation. Defaults to zero. If, after application of the TimeShift value in the blooddata object, the startTime is still after zero, it will be set to zero.
#' @param stopTime The end time for the interpolation. Defaults to the maximum measured time.
#' @param interpPoints The number of points to interpolate over between the start and stop times.
#' @param colours The colours for the plot.
#' @param select_line Plot only one of the measurement lines without the others.
#'
#' @return A ggplot2 object
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' plot_blooddata(blooddata)
#' }
plot_blooddata <- function(blooddata,
                           startTime = 0,
                           stopTime = NULL,
                           interpPoints = 6000,
                           colours = c(
                             "#ff321b",
                             "#1bff7e",
                             "#501bff",
                             "#ffa81b"
                           ),
                           select_line = NULL) {


  if( is.null(blooddata$Data$Plasma$Avail) ) { # old blooddata
    blooddata <- update_blooddata_bids(blooddata)
    message("
         It appears that you are using a blooddata object created using an
         older version of kinfitr, based on an older version of the PET BIDS
         standard. The blooddata object has been updated for this function.
         Ideally, please create your blooddata objects afresh from the raw data.
         Otherwise, please use the update_blooddata_bids() command to keep only
         the bare essentials.")
  }

  # Predicted

  input <- bd_create_input(blooddata, startTime, stopTime, interpPoints)
  input$`Blood-Plasma Ratio` <- input$Blood / input$Plasma
  input <- dplyr::rename(input, "Parent Fraction" = ParentFraction)

  bloodmax <- max(c(input$Blood, input$AIF))
  input$Blood <- input$Blood / bloodmax
  input$AIF <- input$AIF / bloodmax
  input <- dplyr::select(input, -Plasma)

  pred <- suppressWarnings(
    tidyr::gather(input, Outcome, Value, -Time)
  )
  pred <- dplyr::arrange(pred, Outcome, Time)
  # pred <- dplyr::rename(pred, "Parent Fraction" = ParentFraction)

  # Measured data

  pf <- bd_extract(blooddata, output = "parentFraction", startTime, stopTime, interpPoints)
  pf <- dplyr::select(pf, Time = time, Value = parentFraction)
  pf <- dplyr::mutate(pf, Outcome = "Parent Fraction", Measurement = "Discrete")

  bpr <- bd_extract(blooddata, output = "BPR", startTime, stopTime, interpPoints)
  bpr <- dplyr::select(bpr, Time = time, Value = bpr)
  bpr <- dplyr::mutate(bpr, Outcome = "Blood-Plasma Ratio", Measurement = "Discrete")

  blood <- bd_extract(blooddata, output = "Blood", startTime, stopTime, interpPoints)
  blood <- dplyr::select(blood, Time = time, Value = activity, Measurement = Method)
  blood <- dplyr::mutate(blood, Outcome = "Blood", Value = Value / bloodmax)

  aif <- bd_extract(blooddata, output = "AIF", startTime, stopTime, interpPoints)
  aif <- dplyr::select(aif, Time = time, Value = aif, Measurement = Method)
  aif <- dplyr::mutate(aif, Outcome = "AIF", Value = Value / bloodmax)

  measured <- dplyr::bind_rows(list(pf, bpr, blood, aif))
  measured <- dplyr::arrange(measured, Outcome, Time)
  measured <- dplyr::mutate(measured,
                            Time = Time / 60,
                            dotsize = ifelse(Measurement == "Continuous", 1, 2)
  )


  # Check for BPR outside the plot

  if( mean(input$`Blood-Plasma Ratio` > 1.2, na.rm = T) > 0.5 ) {

    measured <- measured %>%
      dplyr::mutate( Value = ifelse(Outcome=="Blood-Plasma Ratio",
                                    yes = 1 / Value,
                                    no = Value),
                     Outcome = ifelse(Outcome=="Blood-Plasma Ratio",
                                      yes = "Plasma-Blood Ratio\n( = 1/BPR )",
                                      no  = Outcome))

    pred <- pred %>%
      dplyr::mutate( Value = ifelse(Outcome=="Blood-Plasma Ratio",
                                    yes = 1 / Value,
                                    no = Value),
                     Outcome = ifelse(Outcome=="Blood-Plasma Ratio",
                                      yes = "Plasma-Blood Ratio\n( = 1/BPR )",
                                      no  = Outcome))

    bpr$Value <- 1 / bpr$Value
  }


  # Plot

  plotmax <- min(c(1.2, max(bpr$Value, na.rm = T)))
  if (plotmax < 1) {
    plotmax <- 1
  }

  if (!( "Continuous" %in% measured$Measurement )) {
    shapes = 16
  } else {
    shapes = c(1, 16)
  }

  if(is.null(select_line)) {

    ggplot(data = measured, aes(x = Time, y = Value, colour = Outcome)) +
      geom_point(aes(shape = Measurement, size = dotsize)) +
      scale_shape_manual(values = shapes) +
      scale_size(range = c(1, 2.5)) +
      geom_line(
        data = pred, linewidth = 0.8, colour = "grey46",
        aes(group = Outcome), alpha = 0.5
      ) +
      geom_line(data = pred) +
      guides(size = "none") +
      coord_cartesian(ylim = c(0, plotmax))

  } else{

    # Plot only one line

    if(select_line == "Blood") {
      measured <- measured %>%
        dplyr::filter(Outcome=="Blood")

      pred <- pred %>%
        dplyr::filter(Outcome=="Blood")
    }

    if(select_line == "BPR") {
      measured <- measured %>%
        dplyr::filter(stringr::str_detect(Outcome, " Ratio"))

      pred <- pred %>%
        dplyr::filter(stringr::str_detect(Outcome, " Ratio"))
    }

    if(select_line == "parentFraction") {
      measured <- measured %>%
        dplyr::filter(Outcome == "Parent Fraction")

      pred <- pred %>%
        dplyr::filter(Outcome == "Parent Fraction")
    }

    if(select_line == "AIF") {
      measured <- measured %>%
        dplyr::filter(Outcome == "AIF")

      pred <- pred %>%
        dplyr::filter(Outcome == "AIF")
    }

    ggplot(data = measured, aes(x = Time, y = Value)) +
      geom_point(aes(shape = Measurement, size = dotsize)) +
      scale_shape_manual(values = shapes) +
      scale_size(range = c(1, 2.5)) +
      geom_line(
        data = pred, linewidth = 0.8, colour = "grey46",
        aes(group = Outcome), alpha = 0.5
      ) +
      geom_line(data = pred) +
      guides(size = "none") +
      coord_cartesian(ylim = c(0, plotmax))

  }


}


create_interpinput <- function(blood_input, startTime = 0,
                               stopTime = NULL,
                               interpPoints = 6000) {
  if ("interpblood" %in% class(blood_input)) {
    return(blood_input)
  }

  if ("blooddata" %in% class(blood_input)) {
    blood_input <- bd_create_input(
      blood_input,
      startTime,
      stopTime,
      interpPoints
    )
    return(blood_input)
  }
}
