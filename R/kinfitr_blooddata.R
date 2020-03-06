#' Create a blooddata object from data vectors
#'
#' This function creates a blooddata object from data vectors. This function
#' creates data objects which contain the minimum amount of data for modelling.
#' Ideally, I recommend storing your data according to the BIDS specification
#' and using the create_blooddata_bids command instead. Please refer to the PET
#' BIDS standard for further details about the inputs.
#'
#' @param Blood.Discrete.Data.Values.sampleStartTime In seconds. Sample start
#'   times (if durations known), or simply sample times (if durations unkown).
#' @param Blood.Discrete.Data.Values.sampleDuration In seconds. Sample durations
#'   (if known). Defaults to zero if durations unknown.
#' @param Blood.Discrete.Data.Values.activity In kBq/ml. Measured radioactivity.
#' @param Plasma.Data.Values.sampleStartTime In seconds. Sample start times (if
#'   durations known), or simply sample times (if durations unkown).
#' @param Plasma.Data.Values.sampleDuration In seconds. Sample durations (if
#'   known). Defaults to zero if durations unknown.
#' @param Plasma.Data.Values.activity In kBq/ml. Measured radioactivity.
#' @param Metabolite.Data.Values.sampleStartTime In seconds. Sample start times
#'   (if durations known), or simply sample times (if durations unkown).
#' @param Metabolite.Data.Values.sampleDuration In seconds. Sample durations (if
#'   known). Defaults to zero if durations unknown.
#' @param Metabolite.Data.Values.parentFraction Measured fraction.
#' @param Blood.Continuous.Data.Values.time In seconds.
#' @param Blood.Continuous.Data.Values.activity in kBq/ml.
#' @param Blood.Continuous.WithdrawalRate The rate at which the blood was
#'   withdrawn from the subject.
#' @param Blood.Continuous.DispersionConstant External dispersion time constant
#'   resulting from tubing.
#' @param Blood.Continuous.DispersionConstantUnits The units of the dispersion
#'   constant.
#' @param Blood.Continuous.DispersionCorrected Boolean flag specifying whether
#'   the continuous blood data have been dispersion-corrected.
#' @param TimeShift The extent to which all the times in the data should be
#'   shifted (in seconds). Defaults to 0.
#'
#' @return a blooddata object
#' @export
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @examples
#' \dontrun{
#' blooddata2 <- create_blooddata_components(
#'    Blood.Discrete.Data.Values.sampleStartTime =
#'      blooddata$Data$Blood$Discrete$Data$Values$sampleStartTime,
#'    Blood.Discrete.Data.Values.sampleDuration =
#'      blooddata$Data$Blood$Discrete$Data$Values$sampleDuration,
#'    Blood.Discrete.Data.Values.activity =
#'      blooddata$Data$Blood$Discrete$Data$Values$activity,
#'    Plasma.Data.Values.sampleStartTime =
#'      blooddata$Data$Plasma$Data$Values$sampleStartTime,
#'    Plasma.Data.Values.sampleDuration =
#'      blooddata$Data$Plasma$Data$Values$sampleDuration,
#'    Plasma.Data.Values.activity =
#'      blooddata$Data$Plasma$Data$Values$activity,
#'    Metabolite.Data.Values.sampleStartTime =
#'      blooddata$Data$Metabolite$Data$Values$sampleStartTime,
#'    Metabolite.Data.Values.sampleDuration =
#'      blooddata$Data$Metabolite$Data$Values$sampleDuration,
#'    Metabolite.Data.Values.parentFraction =
#'      blooddata$Data$Metabolite$Data$Values$parentFraction,
#'    Blood.Continuous.Data.Values.time =
#'      blooddata$Data$Blood$Continuous$Data$Values$time,
#'    Blood.Continuous.Data.Values.activity =
#'      blooddata$Data$Blood$Continuous$Data$Values$activity,
#'    Blood.Continuous.WithdrawalRate =
#'      blooddata$Data$Blood$Continuous$WithdrawalRate,
#'    Blood.Continuous.DispersionConstant =
#'      blooddata$Data$Blood$Continuous$DispersionConstant,
#'    Blood.Continuous.DispersionConstantUnits =
#'      blooddata$Data$Blood$Continuous$DispersionConstantUnits,
#'    Blood.Continuous.DispersionCorrected = FALSE,
#'    TimeShift = 0)
#'    }
create_blooddata_components <- function(
                                        Blood.Discrete.Data.Values.sampleStartTime,
                                        Blood.Discrete.Data.Values.sampleDuration = 0,
                                        Blood.Discrete.Data.Values.activity,
                                        Plasma.Data.Values.sampleStartTime,
                                        Plasma.Data.Values.sampleDuration = 0,
                                        Plasma.Data.Values.activity,
                                        Metabolite.Data.Values.sampleStartTime,
                                        Metabolite.Data.Values.sampleDuration = 0,
                                        Metabolite.Data.Values.parentFraction,
                                        Blood.Continuous.Data.Values.time = NULL,
                                        Blood.Continuous.Data.Values.activity = NULL,
                                        Blood.Continuous.WithdrawalRate = NULL,
                                        Blood.Continuous.DispersionConstant = NULL,
                                        Blood.Continuous.DispersionConstantUnits = NULL,
                                        Blood.Continuous.DispersionCorrected = TRUE,
                                        TimeShift = 0) {




  # Blood

  Blood <- list()

  # Continuous blood

  Blood$Continuous <- list(
    Data = list(
      Values = tibble::tibble(
        time = Blood.Continuous.Data.Values.time,
        activity = Blood.Continuous.Data.Values.activity
      )
    ),
    WithdrawalRate = Blood.Continuous.WithdrawalRate,
    DispersionConstant = Blood.Continuous.DispersionConstant,
    DispersionConstantUnits = Blood.Continuous.DispersionConstantUnits,
    DispersionCorrected = Blood.Continuous.DispersionCorrected
  )

  # Discrete blood
  Blood$Discrete <- list(
    Data = list(
      Values = tibble::tibble(
        sampleStartTime =
          Blood.Discrete.Data.Values.sampleStartTime,
        sampleDuration =
          Blood.Discrete.Data.Values.sampleDuration,
        activity = Blood.Discrete.Data.Values.activity
      )
    )
  )

  # Plasma
  Plasma <- list(
    Data = list(
      Values = tibble::tibble(
        sampleStartTime =
          Plasma.Data.Values.sampleStartTime,
        sampleDuration =
          Plasma.Data.Values.sampleDuration,
        activity = Plasma.Data.Values.activity
      )
    )
  )

  # Metabolite
  Metabolite <- list(
    Data = list(
      Values = tibble::tibble(
        sampleStartTime =
          Metabolite.Data.Values.sampleStartTime,
        sampleDuration =
          Metabolite.Data.Values.sampleDuration,
        parentFraction =
          Metabolite.Data.Values.parentFraction
      )
    )
  )

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
#' This function creates a blooddata object from data structured according to
#' the PET BIDS standard.
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
  if (!(length(predict(fit)) > 1)) {
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

  fitted <- tibble::tibble(time = time, predicted = predicted)

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









#' Create an input object from a blooddata object.
#'
#' Predicts and interpolates all of the functions to create an input object, or to a set of values for modelling of blood-related functions.
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
      activity = predict(blooddata$Models$Blood$Data,
        newdata = list(time = interptime)
      )
    )

    blood$activity <- predict(blooddata$Models$Blood$Data,
      newdata = list(time = blood$time)
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

  # Blood-to-Plasma Ratio

  plasma <- blooddata$Data$Plasma$Data$Values
  plasma <- dplyr::filter(plasma, !is.na(activity))

  plasma$time <- plasma$sampleStartTime +
    0.5 * plasma$sampleDuration
  plasma <- dplyr::arrange(plasma, time)

  commonvalues <- intersect(plasma$time, blood_discrete$time)

  bprvec <- blood_discrete$activity[blood_discrete$time %in% commonvalues] /
    plasma$activity[plasma$time %in% commonvalues]


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
      bpr = predict(blooddata$Models$BPR$Data,
        newdata = list(time = interptime)
      )
    )

    blood$bpr <- predict(blooddata$Models$BPR$Data,
      newdata = list(time = blood$time)
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
      parentFraction = predict(blooddata$Models$parentFraction$Data,
        newdata = list(time = interptime)
      )
    )

    blood$parentFraction <- predict(blooddata$Models$parentFraction$Data,
      newdata = list(time = blood$time)
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
    plasma = blood * (1 / bpr),
    aif = plasma * parentFraction
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
      aif = predict(blooddata$Models$AIF$Data,
        newdata = list(time = interptime)
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

interpends <- function(x, y, xi, method = "linear", yzero = NULL) {
  if (is.null(yzero)) {
    yzero <- head(y, 1)
  }

  if (head(xi, 1) < head(x, 1)) {
    x <- c(head(xi, 1), x)
    y <- c(yzero, y)
  }

  if (tail(xi, 1) > tail(x, 1)) {
    x <- c(x, tail(xi, 1))
    y <- c(y, tail(y, 1))
  }

  pracma::interp1(x, y, xi, method)
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
                             smooth_iterations = 0) {

  if (is.null(blooddata$Data$Blood$Continuous)) {

    stop("There is no continuous blood data to
         perform dispersion correction on.")
  }

  if (!blooddata$Data$Blood$Continuous$DispersionCorrected) { # not corrected

    blood_continuous <- blooddata$Data$Blood$Continuous$Data$Values

    blood_dc <- blood_dispcor(blood_continuous$time,
      blood_continuous$activity,
      tau = blooddata$Data$Blood$Continuous$DispersionConstant,
      timedelta, keep_interpolated = T,
      smooth_iterations = smooth_iterations
    )

    blooddata$Data$Blood$Continuous$Data$Values <- blood_dc
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
                           )) {

  # Predicted

  input <- bd_getdata(blooddata, startTime, stopTime, interpPoints)
  input$`Blood-Plasma Ratio` <- input$Blood / input$Plasma
  input <- dplyr::rename(input, "Parent Fraction" = ParentFraction)

  bloodmax <- max(c(input$Blood, input$AIF))
  input$Blood <- input$Blood / bloodmax
  input$AIF <- input$AIF / bloodmax
  input <- dplyr::select(input, -Plasma)

  pred <- tidyr::gather(input, Outcome, Value, -Time)
  pred <- dplyr::arrange(pred, Outcome, Time)
  # pred <- dplyr::rename(pred, "Parent Fraction" = ParentFraction)

  # Measured data

  pf <- bd_getdata(blooddata, startTime, stopTime, interpPoints, output = "parentFraction")
  pf <- dplyr::select(pf, Time = time, Value = parentFraction)
  pf <- dplyr::mutate(pf, Outcome = "Parent Fraction", Measurement = "Discrete")

  bpr <- bd_getdata(blooddata, startTime, stopTime, interpPoints, output = "BPR")
  bpr <- dplyr::select(bpr, Time = time, Value = bpr)
  bpr <- dplyr::mutate(bpr, Outcome = "Blood-Plasma Ratio", Measurement = "Discrete")

  blood <- bd_getdata(blooddata, startTime, stopTime, interpPoints, output = "Blood")
  blood <- dplyr::select(blood, Time = time, Value = activity, Measurement = Method)
  blood <- dplyr::mutate(blood, Outcome = "Blood", Value = Value / bloodmax)

  aif <- bd_getdata(blooddata, startTime, stopTime, interpPoints, output = "AIF")
  aif <- dplyr::select(aif, Time = time, Value = aif, Measurement = Method)
  aif <- dplyr::mutate(aif, Outcome = "AIF", Value = Value / bloodmax)

  measured <- dplyr::bind_rows(list(pf, bpr, blood, aif))
  measured <- dplyr::arrange(measured, Outcome, Time)
  measured <- dplyr::mutate(measured,
    Time = Time / 60,
    dotsize = ifelse(Measurement == "Continuous", 1, 2)
  )


  # Plot

  plotmax <- min(c(1.2, max(bpr$Value, na.rm = T)))
  if (plotmax < 1) {
    plotmax <- 1
  }

  ggplot(data = measured, aes(x = Time, y = Value, colour = Outcome)) +
    geom_point(aes(shape = Measurement, size = dotsize)) +
    scale_shape_manual(values = c(1, 16)) +
    scale_size(range = c(1, 2.5)) +
    geom_line(
      data = pred, size = 0.8, colour = "grey46",
      aes(group = Outcome), alpha = 0.5
    ) +
    geom_line(data = pred) +
    guides(size = FALSE) +
    coord_cartesian(ylim = c(0, plotmax))
}


create_interpinput <- function(blood_input, startTime = 0,
                               stopTime = NULL,
                               interpPoints = 6000) {
  if ("interpblood" %in% class(blood_input)) {
    return(blood_input)
  }

  if ("blooddata" %in% class(blood_input)) {
    blood_input <- bd_getdata(
      blood_input,
      startTime,
      stopTime,
      interpPoints
    )
    return(blood_input)
  }
}
