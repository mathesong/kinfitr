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
      parentFrac = list(Method = "interp", Data = NULL),
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

  # Shouldn't be here
  bids_data <- "Z:\\kipet\\atlas_images\\human\\vff\\Granville\\TRT_ShareDat\\PBR28\\PBR28_mini\\sub-03\\pet\\sub-03_ses-1_task-rest_pet.json"
  bids_data <- jsonlite::fromJSON(bids_data)
  bids_data$Blood.Continuous$Data.Labels <- c("time", "activity")


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
      parentFrac = list(Method = "interp", Data = NULL),
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
                                                  "parentFrac",
                                                  "AIF") ) {

  # Verify fit object
  if( !( length(predict(fit)) > 1 ) ) {
    stop("The fit object should be capable of producing
         predicted values with the predict command.")
  }

  if(length(modeltype) > 1) {
    stop("modeltype should be defined for the model as one of
         the following: Blood, BPR, parentFrac or AIF")
  }
  modeltype <- match.arg(modeltype)

  blooddata$Models[[modeltype]] <- list(Method = "fit", Data=fit)

  return(blooddata)

}

blood_addfitpars <- function(blooddata, modelname, fitpars,
                             modeltype = c("Blood", "BPR", "parentFrac", "AIF")) {


  # Verify fit type
  if( !exists(modelname,
              where = "package:kinfitr",
              mode = "function") ) {
    stop("The modelname should match the model name of
         the corresponding model within the kinfitr package.")
  }

  if(length(modeltype) > 1) {
    stop("modeltype should be defined for the model as one of
         the following: Blood, BPR, parentFrac or AIF")
  }
  modeltype <- match.arg(modeltype)

  blooddata$Models[[modeltype]] <- list(Method = "fitpars", Data=fitpars)

  return(blooddata)

}


blood_addfitted <- function(blooddata, time, predicted,
                             modeltype = c("Blood", "BPR", "parentFrac", "AIF")) {


  if(length(modeltype) > 1) {
    stop("modeltype should be defined for the model as one of
         the following: Blood, BPR, parentFrac or AIF")
  }
  modeltype <- match.arg(modeltype)

  fitted <- tibble::tibble(time = time, predicted = predicted)

  blooddata$Models[[modeltype]] <- list(Method = "fitted", Data=fitted)

  return(blooddata)

}

blooddata2input <- function(blooddata, startTime, stopTime, interpPoints = 6000) {



}
