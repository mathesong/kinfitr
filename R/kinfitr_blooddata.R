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
  Blood.Continuous.Data.Values.time=NULL,
  Blood.Continuous.Data.Values.activity=NULL,
  Blood.Continuous.WithdrawalRate=NULL,
  Blood.Continuous.DispersionConstant=NULL,
  Blood.Continuous.DispersionConstantUnits=NULL,
  Blood.Continuous.DispersionCorrected = TRUE,
  TimeShift = 0) {


  # Continuous blood
  Blood.Continuous <- list(
    Data.Values = tibble::tibble(time = Blood.Continuous.Data.Values.time,
                                 activity = Blood.Continuous.Data.Values.activity),
    WithdrawalRate = Blood.Continuous.WithdrawalRate,
    DispersionConstant = Blood.Continuous.DispersionConstant,
    DispersionConstantUnits = Blood.Continuous.DispersionConstantUnits,
    DispersionCorrected = Blood.Continuous.DispersionCorrected
  )

  # Discrete blood
  Blood.Discrete <- list(
    Data.Values = tibble::tibble(sampleStartTime =
                                   Blood.Discrete.Data.Values.sampleStartTime,
                                 sampleDuration =
                                   Blood.Discrete.Data.Values.sampleDuration,
                                 activity = Blood.Discrete.Data.Values.activity))

  # Plasma
  Plasma <- list(
    Data.Values = tibble::tibble(sampleStartTime =
                                   Plasma.Data.Values.sampleStartTime,
                                 sampleDuration =
                                   Plasma.Data.Values.sampleDuration,
                                 activity = Plasma.Data.Values.activity))

  # Metabolite
  Metabolite <- list(
    Data.Values = tibble::tibble(sampleStartTime =
                                   Metabolite.Data.Values.sampleStartTime,
                                 sampleDuration =
                                   Metabolite.Data.Values.sampleDuration,
                                 parentFraction = Metabolite.Data.Values.activity))

  blooddata <- list(
    Data = list(
      Blood.Continuous = Blood.Continuous,
      Blood.Discrete = Blood.Discrete,
      Plasma = Plasma,
      Metabolite = Metabolite
    ),
    Models = list(
      Blood = list(Method = "interp", Data = NULL),
      BPR = list(Method = "interp", Data = NULL),
      parentFraction = list(Method = "interp", Data = NULL),
      AIF = list(Method = "interp", Data = NULL)
    ),
    TimeShift = TimeShift)

  blooddata <- rapply(blooddata, function(x) ifelse(x=="true",TRUE,x), how = "replace")
  blooddata <- rapply(blooddata, function(x) ifelse(x=="false",FALSE,x), how = "replace")

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
#' @examples
#' \dontrun{
#' a <- create_blooddata_bids('bids_sidecar.json')
#' }
create_blooddata_bids <- function(bids_data, TimeShift = 0) {

  if( !is.list(bids_data) ) {
    if( file.exists(bids_data) &
        tools::file_ext(bids_data) == "json" ) {
      bids_data <- jsonlite::fromJSON(bids_data)
    }
  }

  tibblify_bidsjson <- function(list) {
    list$Data.Values <- tibble::as_tibble(list$Data.Values)
    colnames(list$Data.Values) <- list$Data.Labels
    return(list)
  }


  # Continuous blood
  Blood.Continuous <- tibblify_bidsjson(bids_data$Blood.Continuous)

  # Discrete blood
  Blood.Discrete <- tibblify_bidsjson(bids_data$Blood.Discrete)

  # Plasma
  Plasma <- tibblify_bidsjson(bids_data$Plasma)

  # Metabolite
  Metabolite <- tibblify_bidsjson(bids_data$Metabolite)

  blooddata <- list(
    Data = list(
      Blood.Continuous = Blood.Continuous,
      Blood.Discrete = Blood.Discrete,
      Plasma = Plasma,
      Metabolite = Metabolite
    ),
    Models = list(
      Blood = list(Method = "interp", Data = NULL),
      BPR = list(Method = "interp", Data = NULL),
      parentFraction = list(Method = "interp", Data = NULL),
      AIF = list(Method = "interp", Data = NULL)
    ),
    TimeShift = 0)

  blooddata <- rapply(blooddata, function(x) ifelse(x=="true",TRUE,x), how = "replace")
  blooddata <- rapply(blooddata, function(x) ifelse(x=="false",FALSE,x), how = "replace")

  class(blooddata) <- "blooddata"

  return(blooddata)
}


blood_addfit <- function(blooddata, fit, modeltype = c("Blood",
                                                  "BPR",
                                                  "parentFraction",
                                                  "AIF") ) {

  # Verify fit object
  if( !( length(predict(fit)) > 1 ) ) {
    stop("The fit object should be capable of producing
         predicted values with the predict command.")
  }

  if(length(modeltype) > 1) {
    stop("modeltype should be defined for the model as one of
         the following: Blood, BPR, parentFraction or AIF")
  }
  modeltype <- match.arg(modeltype, c("Blood",
                                      "BPR",
                                      "parentFraction",
                                      "AIF"))

  blooddata$Models[[modeltype]] <- list(Method = "fit", Data=fit)

  return(blooddata)

}

blood_addfitpars <- function(blooddata, modelname, fitpars,
                             modeltype = c("Blood", "BPR", "parentFraction", "AIF")) {


  # Verify fit type
  if( !exists(modelname,
              mode = "function") ) {
    stop("The modelname should match the model name of
         the corresponding model function.")
  }

  if(length(modeltype) > 1) {
    stop("modeltype should be defined for the model as one of
         the following: Blood, BPR, parentFraction or AIF")
  }
  modeltype <- match.arg(modeltype, c("Blood",
                                      "BPR",
                                      "parentFraction",
                                      "AIF"))



  blooddata$Models[[modeltype]] <- list(Method = "fitpars",
                                        Data=list(Model = modelname,
                                                  Pars = fitpars))

  return(blooddata)

}


blood_addfitted <- function(blooddata, time, predicted,
                             modeltype = c("Blood", "BPR", "parentFraction", "AIF")) {


  if(length(modeltype) > 1) {
    stop("modeltype should be defined for the model as one of
         the following: Blood, BPR, parentFraction or AIF")
  }
  modeltype <- match.arg(modeltype, c("Blood",
                                      "BPR",
                                      "parentFraction",
                                      "AIF"))

  fitted <- tibble::tibble(time = time, predicted = predicted)

  blooddata$Models[[modeltype]] <- list(Method = "fitted", Data=fitted)

  return(blooddata)

}

blooddata2input <- function(blooddata,
                            startTime = 0,
                            stopTime = NULL,
                            interpPoints = 6000,
                            output = c("input",
                                       "Blood",
                                       "BPR",
                                       "parentFraction",
                                       "AIF")) {

  if(is.null(stopTime)) {

    stopTime = max( c(
      with(blooddata$Data$Blood.Discrete$Data.Values,
           sampleStartTime + sampleDuration),
      with(blooddata$Data$Plasma$Data.Values,
           sampleStartTime + sampleDuration),
      with(blooddata$Data$Metabolite$Data.Values,
           sampleStartTime + sampleDuration),
      blooddata$Data$Blood.Continuous$Data.Values$time) )

  }

  if(startTime > blooddata$TimeShift*(-1)) {
    startTime = 0
  }

  interptime <- seq(startTime, stopTime, length.out = interpPoints)

  output <- match.arg(output, c("input",
                                "Blood",
                                "BPR",
                                "parentFraction",
                                "AIF"))

  # Blood

  blood_discrete <- blooddata$Data$Blood.Discrete$Data.Values
  blood_discrete$time <- blood_discrete$sampleStartTime +
    0.5*blood_discrete$sampleDuration
  blood_discrete <- dplyr::arrange(blood_discrete, time)

  blood_continuous <- blooddata$Data$Blood.Continuous$Data.Values

  blood <- dplyr::bind_rows(blood_discrete, blood_continuous)
  blood <- dplyr::arrange(blood, time)
  blood$Method <- ifelse(is.na(blood$sampleDuration),
                         yes = "Continuous",
                         no = "Discrete")

  if(output == "Blood") {
    return(blood)
  }

    ## Interp
    if(blooddata$Models$Blood$Method == "interp") {

      suppressWarnings(
        i_blood <- tibble::tibble(time = interptime,
                                  activity = interpends(blood$time,
                                                        blood$activity,
                                                        interptime,
                                                        method = "linear",
                                                        yzero=0))
      )
    }

    ## Fit
    if(blooddata$Models$Blood$Method == "fit") {

      i_blood <- tibble::tibble(time = interptime,
                                activity = predict(blooddata$Models$Blood$Method$Data,
                                                   newdata = list(time = interptime) ))
    }

    ## Fit pars
    if(blooddata$Models$Blood$Method == "fitpars") {

      modelname <- blooddata$Models$Blood$Data$Model
      Pars <- append(list(time=interptime),
                     as.list(blooddata$Models$Blood$Data$Pars))

      i_blood <- tibble::tibble(time = interptime,
                                activity = do.call(what = modelname,
                                                   args = Pars) )
    }

    ## Fitted
    if(blooddata$Models$Blood$Method == "fitted") {

      i_blood <- tibble::tibble(time = interptime,
                                activity = interpends(blooddata$Models$Blood$Data$time,
                                                      blooddata$Models$Blood$Data$activity,
                                                      interptime,
                                                      method = "linear",
                                                      yzero=0))
    }

  # Blood-to-Plasma Ratio

  plasma <- blooddata$Data$Plasma$Data.Values
  plasma$time <- plasma$sampleStartTime +
    0.5*plasma$sampleDuration
  plasma <- dplyr::arrange(plasma, time)

  commonvalues <- intersect(plasma$time, blood_discrete$time)

  bprvec <- blood_discrete$activity[blood_discrete$time %in% commonvalues] /
    plasma$activity[plasma$time %in% commonvalues]


  bpr <- tibble::tibble(time = commonvalues, bpr = bprvec)

  if(output == "BPR") {
    return(bpr)
  }

    ## Interp
    if(blooddata$Models$BPR$Method == "interp") {

      i_bpr <- tibble::tibble(time = interptime,
                                bpr = interpends(bpr$time,
                                                 bpr$bpr,
                                                 interptime,
                                                 method = "linear"))
      blood$bpr <- interpends(bpr$time, bpr$bpr, blood$time,
                              method = "linear")
    }

    ## Fit
    if(blooddata$Models$BPR$Method == "fit") {

      i_bpr <- tibble::tibble(time = interptime,
                                bpr = predict(blooddata$Models$BPR$Method$Data,
                                                   newdata = list(time = interptime) ))

      blood$bpr <- predict(blooddata$Models$BPR$Method$Data,
                           newdata = list(time = blood$time) )
    }

    ## Fit pars
    if(blooddata$Models$BPR$Method == "fitpars") {

      modelname <- blooddata$Models$BPR$Data$Model
      Pars <- append(list(time=interptime),
                     as.list(blooddata$Models$BPR$Data$Pars))

      i_bpr <- tibble::tibble(time = interptime,
                                activity = do.call(what = modelname,
                                                   args = Pars) )
      Pars <- append(list(time=blood$time),
                     as.list(blooddata$Models$BPR$Data$Pars))

      blood$bpr <- do.call(what = modelname,
                           args = Pars)
    }

    ## Fitted
    if(blooddata$Models$BPR$Method == "fitted") {

      i_bpr <- tibble::tibble(time = interptime,
                                activity = interpends(blooddata$Models$BPR$Data$time,
                                                      blooddata$Models$BPR$Data$activity,
                                                      interptime,
                                                      method = "linear"))
      blood$bpr <- interpends(blooddata$Models$BPR$Data$time,
                              blooddata$Models$BPR$Data$activity,
                              blood$time,
                              method = "linear")

    }


  # Parent Fraction

  pf <- blooddata$Data$Metabolite$Data.Values
  pf$time <- pf$sampleStartTime + 0.5*pf$sampleDuration

  if(output == "parentFraction") {
    return(pf)
  }

    ## Interp
    if(blooddata$Models$parentFraction$Method == "interp") {

      i_pf <- tibble::tibble(time = interptime,
                              parentFraction = interpends(pf$time,
                                               pf$parentFraction,
                                               interptime,
                                               method = "linear",
                                               yzero = 1))
      blood$parentFraction <- interpends(pf$time, pf$parentFraction,
                                               blood$time,
                                               method = "linear")

    }

    ## Fit
    if(blooddata$Models$parentFraction$Method == "fit") {

      i_pf <- tibble::tibble(time = interptime,
                              parentFraction = predict(blooddata$Models$parentFraction$Method$Data,
                                            newdata = list(time = interptime) ))

      blood$parentFraction <- predict(blooddata$Models$parentFraction$Method$Data,
                           newdata = list(time = blood$time) )
    }

    ## Fit pars
    if(blooddata$Models$parentFraction$Method == "fitpars") {

      modelname <- blooddata$Models$parentFraction$Data$Model
      Pars <- append(list(time=interptime),
                     as.list(blooddata$Models$parentFraction$Data$Pars))

      i_pf <- tibble::tibble(time = interptime,
                              activity = do.call(what = modelname,
                                                 args = Pars) )

      Pars <- append(list(time=blood$time),
                     as.list(blooddata$Models$parentFraction$Data$Pars))

      blood$parentFraction <- do.call(what = modelname,
                           args = Pars)
    }

    ## Fitted
    if(blooddata$Models$parentFraction$Method == "fitted") {

      i_pf <- tibble::tibble(time = interptime,
                              parentFraction = interpends(blooddata$Models$parentFraction$Data$time,
                                                    blooddata$Models$parentFraction$Data$parentFraction,
                                                    interptime,
                                                    method = "linear",
                                                    yzero = 1))

      blood$parentFraction <- interpends(blooddata$Models$parentFraction$Data$time,
                              blooddata$Models$parentFraction$Data$parentFraction,
                              blood$time,
                              method = "linear",
                              yzero=1)
    }

  # AIF

  aif <- blood
  aif <- dplyr::rename(aif, blood = activity)
  aif <- dplyr::mutate(aif,
                       plasma = blood * (1/bpr),
                       aif = plasma * parentFraction)

  if(output == "AIF") {
    return(aif)
  }


    ## Interp
    if(blooddata$Models$AIF$Method == "interp") {

      suppressWarnings(
        i_aif <- tibble::tibble(time = interptime,
                              aif = interpends(aif$time,
                                               aif$aif,
                                               interptime,
                                               method = "linear",
                                               yzero = 0))
      )

    }

    ## Fit
    if(blooddata$Models$AIF$Method == "fit") {

      i_aif <- tibble::tibble(time = interptime,
                              aif = predict(blooddata$Models$AIF$Method$Data,
                                           newdata = list(time = interptime) ))
    }

    ## Fit pars
    if(blooddata$Models$AIF$Method == "fitpars") {

      modelname <- blooddata$Models$AIF$Data$Model
      Pars <- append(list(time=interptime),
                     as.list(blooddata$Models$AIF$Data$Pars))

      i_aif <- tibble::tibble(time = interptime,
                             activity = do.call(what = modelname,
                                                args = Pars) )
    }

    ## Fitted
    if(blooddata$Models$AIF$Method == "fitted") {

      i_aif <- tibble::tibble(time = interptime,
                              aif = interpends(blooddata$Models$AIF$Data$time,
                                               blooddata$Models$AIF$Data$aif,
                                               interptime,
                                               method = "linear",
                                               yzero = 0))
    }

  if(output == "input") {

    input <- tibble::tibble(
      Time = (interptime + TimeShift)/60,
      Blood = i_blood$activity,
      Plasma = i_blood$activity * i_bpr$bpr,
      ParentFraction = i_pf$parentFraction,
      AIF = i_aif$aif
    )

    return(input)

  }

}

interpends <- function(x, y, xi, method="linear", yzero = NULL) {

  if(is.null(yzero)) {
    yzero = head(y , 1)
  }

  if( head(xi, 1) < head(x, 1)) {
    x <- c( head(xi, 1), x)
    y <- c( yzero, y)
  }

  if( tail(xi, 1) < tail(x, 1)) {
    x <- c( x, tail(xi, 1) )
    y <- c( y, tail(y , 1) )
  }

  pracma::interp1(x, y, xi, method)

}

