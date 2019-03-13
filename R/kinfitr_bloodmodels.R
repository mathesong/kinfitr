#' Blood Model: Splines
#'
#' Fits two or three (if both discrete and continuous) splines to blood or AIF data to smooth it
#'
#' @param time The time of each measurement in seconds
#' @param activity The radioactivity of each measurement
#' @param Method Optional. The method of collection, i.e. "Discrete" or "Continuous"
#' @param weights Optional. Weights of each measurement.
#'
#' @return A model fit including all of the individual models.
#' @export
#'
#' @examples
#' blooddata <- create_blooddata_bids(pbr28$jsondata[[1]])
#' blood <- bd_getdata(blooddata, output = "Blood")
#' blood_fit <- blmod_splines(blood$time,
#'                            blood$activity,
#'                            Method = blood$Method)
blmod_splines <- function(time, activity, Method = NULL, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(time))
  }

  blood <- tibble::tibble(time = time, activity = activity, weights = weights)

  blood <- dplyr::filter(blood, !is.na(activity))
  blood <- dplyr::arrange(blood, time)


  if (is.null(Method)) {
    blood$Method <- "Discrete"
  } else {
    blood$Method <- Method
  }

  if (!is.null(Method)) {
    if (!all(Method %in% c("Continuous", "Discrete"))) {
      stop("Unrecognised Method input - it should
           either be Continuous or Discrete")
    }
  }


  peaktime <- blood$time[blood$activity == max(blood$activity)]

  before_peak <- dplyr::filter(blood, time <= peaktime)
  before_peak$time [ nrow(before_peak) ] <-
    before_peak$time [ nrow(before_peak) ] - 0.001 # predict until nearly there

  after_peak <- dplyr::filter(blood, time >= peaktime)

  if ("Continuous" %in% Method) {
    before_peak <- dplyr::filter(before_peak, Method == "Continuous")
    after_peak_d <- dplyr::filter(after_peak, Method == "Discrete")
    after_peak_c <- dplyr::filter(after_peak, Method == "Continuous")

    before <- pspline::sm.spline(before_peak$time, before_peak$activity, w = before_peak$weights)
    after_d <- pspline::sm.spline(after_peak_d$time, after_peak_d$activity, w = after_peak_d$weights)
    after_c <- pspline::sm.spline(after_peak_c$time, after_peak_c$activity, w = after_peak_c$weights)

    start_overlap <- min(after_peak_d$time)
    stop_overlap <- max(after_peak_c$time)
  } else { # i.e. no Continuous

    before <- pspline::sm.spline(before_peak$time, before_peak$activity, w = before_peak$weights)
    after_d <- pspline::sm.spline(after_peak$time, after_peak$activity, w = after_peak$weights)
    after_c <- pspline::sm.spline(after_peak$time, after_peak$activity, w = after_peak$weights)

    start_overlap <- peaktime
    stop_overlap <- max(after_peak$time)
  }

  out <- list(
    before = before,
    after_d = after_d,
    after_c = after_c,
    peaktime = peaktime,
    start_overlap = start_overlap,
    stop_overlap = stop_overlap
  )

  class(out) <- c("blood_splines", class(out))

  return(out)
}


# #' Predict method for blood splines
# #'
# #' Predicts values for new times for blood splines fits.
# #'
# #' @param object Blood splines fit.
# #' @param newdata A new data list, including times to be predicted for.
# #'
# #' @return Model predictions
# #' @export
# #'
# #' @examples
# #' blooddata <- create_blooddata_bids(pbr28$jsondata[[1]])
# #' blood <- bd_getdata(blooddata, output = "Blood")
# #' blood_fit <- blmod_splines(blood$time,
# #'                            blood$activity,
# #'                            Method = blood$Method)
# #'
# #' predict(blood_fit, newdata=list(time=1:10))
predict.blood_splines <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    pred_before <- predict(object$before)
    pred_x_before <- pred_before$x

    pred_after_d <- predict(object$after_d)
    pred_x_after_d <- pred_after_d$x

    pred_after_c <- predict(object$after_c)
    pred_x_after_c <- pred_after_c$x

    pred_x_after <- unique(c(pred_x_after_d, pred_x_after_c))
    pred_x_after <- pred_x_after[order(pred_x_after)]

    newdata <- list(time = c(pred_x_before, pred_x_after))
  }

  pred_before <- predict(object$before, x = newdata$time)[, 1]
  pred_after_d <- predict(object$after_d, x = newdata$time)[, 1]
  pred_after_c <- predict(object$after_c, x = newdata$time)[, 1]

pred_before <- tibble::tibble(time = newdata$time, activity = pred_before)
pred_before <- dplyr::mutate_all(pred_before, dplyr::funs(replace(., is.nan(.), 0)))
  pred_after <- tibble::tibble(
    time = newdata$time,
    activity_c = pred_after_c,
    activity_d = pred_after_d,
    overlapfrac =
      (newdata$time -
        object$start_overlap) /
        (object$stop_overlap -
          object$start_overlap)
  )

  pred_after <- dplyr::mutate_all(pred_after, dplyr::funs(replace(., is.nan(.), 0)))

  pred_after$cweights <- dplyr::case_when(
    pred_after$overlapfrac < 0 ~ 1,
    pred_after$overlapfrac > 1 ~ 0,
    TRUE ~ 1 - pred_after$overlapfrac
  )

  pred_after$dweights <- 1 - pred_after$cweights



  pred_after$activity <- apply(pred_after,
    MARGIN = 1,
    function(x)
      stats::weighted.mean(
        x = c(x["activity_c"], x["activity_d"]),
        w = c(x["cweights"], x["dweights"])
      )
  )

  pred_before <- dplyr::filter(pred_before, time < object$peaktime)
  pred_after <- dplyr::filter(pred_after, time >= object$peaktime)

  pred <- dplyr::bind_rows(pred_before, pred_after)

  preds <- dplyr::pull(pred, activity)

  return(preds)
}
