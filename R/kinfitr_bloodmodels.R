#' Blood Model: Tidy inputs
#'
#' Tidies the inputs to blood models ready for modelling
#'
#' @param time The time of each measurement in seconds
#' @param activity The radioactivity of each measurement
#' @param Method Optional. The method of collection, i.e. "Discrete" or "Continuous"
#' @param weights Optional. Weights of each measurement.
#'
#' @return A blood data frame ready for modelling
#' @export
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' blood <- bd_extract(blooddata, output = "Blood")
#' blood <- blmod_tidyinput(blood$time,
#'                            blood$activity,
#'                            Method = blood$Method)
blmod_tidyinput <- function(time, activity, Method = NULL, weights = NULL) {

  if (is.null(Method)) {
    Method <- rep("Discrete", length(time))
  }

  # if (is.null(weights)) {
  #
  #   if(length(unique(Method))==2) {
  #     wc <- as.data.frame(table(Method))
  #     wc$weights <- 1/wc$Freq
  #     wc$weights <- wc$weights / max(wc$weights)
  #     wc <- wc[,-2]
  #   }
  #
  #   weights <- rep(1, length(time))
  # }

  blood <- tibble::tibble(time = time, activity = activity,
                          Method = Method)

  if (is.null(weights)) {
    methodsummary <- as.data.frame(table(Method=blood$Method))
    methodsummary$weights <- 1/methodsummary$Freq
    methodsummary$weights <- methodsummary$weights / max(methodsummary$weights)
    methodsummary <- methodsummary[,-2]

    blood <- merge(blood, methodsummary, by="Method")
  } else {
    blood$weights <- weights
  }

  blood <- dplyr::filter(blood, !is.na(activity))
  blood <- dplyr::arrange(blood, time)
  blood$activity <- ifelse(blood$activity < 0, yes = 0, no = blood$activity)

  # Peak fraction - for weighing continuos and discrete samples
  peaktime <- blood$time[which.max(blood$activity)]

  if("Continuous" %in% blood$Method) {
    cont_samples <- dplyr::filter(blood, Method=="Continuous")
    end_cont <- cont_samples$time[which.max(cont_samples$time)]

    if(end_cont > peaktime) {
      peakfrac <- (blood$time - peaktime) / (end_cont - peaktime)
      peakfrac[peakfrac < 0] <- 0
      peakfrac[peakfrac > 1] <- 1
      blood$peakfrac <- peakfrac
    } else {
      blood$peakfrac <- ifelse(blood$time < peaktime,
                               yes = 0, no=1)
    }
  } else {
    blood$peakfrac <- 1
  }

  if (!all(blood$Method %in% c("Continuous", "Discrete"))) {
    stop("Unrecognised Method input - it should
         either be Continuous or Discrete")
  }

  return(blood)

}


#' Blood Model: Splines
#'
#' Fits two or three (if both discrete and continuous) splines to blood or AIF data to smooth it
#'
#' @param time The time of each measurement in seconds
#' @param activity The radioactivity of each measurement
#' @param Method Optional. The method of collection, i.e. "Discrete" or "Continuous"
#' @param weights Optional. Weights of each measurement.
#' @param bs_before Optional. Defines the basis function for the points before the peak.
#' @param bs_after_c Optional. Defines the basis function for the continuous points after the peak.
#' @param bs_after_d Optional. Defines the basis function for the discrete points after the peak.
#' @param k_before Optional. Defines the dimension of the basis for the points before the peak.
#' @param k_after_c Optional. Defines the dimension of the basis for the continuous points after the peak.
#' @param k_after_d Optional. Defines the dimension of the basis for the discrete points after the peak.
#'
#' @return A model fit including all of the individual models of class blood_splines.
#' @export
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' blood <- bd_extract(blooddata, output = "Blood")
#' blood_fit <- blmod_splines(blood$time,
#'                            blood$activity,
#'                            Method = blood$Method)
blmod_splines <- function(time, activity, Method = NULL, weights = NULL,
                          bs_before="cr", bs_after_c="cr", bs_after_d="cr",
                          k_before=-1, k_after_c=-1, k_after_d=-1) {

  blood <- blmod_tidyinput(time, activity, Method, weights)

  peaktime <- blood$time[blood$activity == max(blood$activity)]

  before_peak <- dplyr::filter(blood, time <= peaktime)
  before_peak$time [ nrow(before_peak) ] <-
    before_peak$time [ nrow(before_peak) ] - 0.001 # predict until nearly there

  # Check knots
  if( nrow(before_peak) < 12 ) {
    if(k_before == -1 || k_before > (nrow(before_peak) - 2)) {
      warning("k_before reduced on account of few pre-peak samples")
      k_before <- nrow(before_peak) - 2
    }
  }

  after_peak <- dplyr::filter(blood, time >= peaktime)

  if ("Continuous" %in% Method) {
    before_peak <- dplyr::filter(before_peak, Method == "Continuous")
    after_peak_d <- dplyr::filter(after_peak, Method == "Discrete")
    after_peak_c <- dplyr::filter(after_peak, Method == "Continuous")


    # Check knots
    if( nrow(after_peak_d) < 12 ) {
      if(k_after_d == -1 || k_after_d > (nrow(after_peak_d) - 2)) {
        warning("k_after_d reduced on account of few post-peak discrete samples")
        k_after_d <- nrow(after_peak_d) - 2
      }
    }

    if( nrow(after_peak_c) < 12 ) {
      if(k_after_c == -1 || k_after_c > (nrow(after_peak_c) - 2)) {
        warning("k_after_c reduced on account of few post-peak continuous samples")
        k_after_c <- nrow(after_peak_c) - 2
      }
    }

    # Fit
    before <- mgcv::gam(activity ~ s(time, bs = bs_before, k=k_before),
                        weights = weights,
                        data = before_peak)

    after_d <- mgcv::gam(activity ~ s(log(time), bs = bs_after_d, k=k_after_d),
                         weights = weights,
                         data = after_peak_d)

    after_c <- mgcv::gam(activity ~ s(log(time), bs = bs_after_c, k=k_after_c),
                         weights = weights,
                         data = after_peak_c)

    start_overlap <- min(after_peak_d$time)
    stop_overlap <- max(after_peak_c$time)


  } else { # i.e. no Continuous

    # Check knots
    if( nrow(after_peak) < 12 ) {
      if(k_after_d == -1 || k_after_d > (nrow(after_peak) - 2)) {
        warning("k_after_d reduced on account of few post-peak samples")
        k_after_d <- nrow(after_peak) - 2
      }
    }

    # Fit
    before <- mgcv::gam(activity ~ s(time, bs = bs_before, k=k_before),
                        weights = weights,
                        data = before_peak)

    after_d <- mgcv::gam(activity ~ s(log(time), bs = bs_after_d, k=k_after_d),
                         weights = weights,
                         data = after_peak)

    after_c <- mgcv::gam(activity ~ s(log(time), bs = bs_after_d, k=k_after_d),
                         weights = weights,
                         data = after_peak)

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

  class(out) <- c("blood_splines", "blmod", class(out))

  return(out)
}


#' Predicted values for blood splines
#'
#' Predicts values for new times for blood splines fits.
#'
#' @param object Blood splines fit.
#' @param newdata A new data list, including times to be predicted for.
#'
#' @return Model predictions
#' @export
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blood <- bd_extract(blooddata, output = "Blood")
#' blood_fit <- blmod_splines(blood$time,
#'                            blood$activity,
#'                            Method = blood$Method)
#'
#' predict_blood_splines(blood_fit, newdata=list(time=1:10))
#' @export
predict_blood_splines <- function(object, newdata = NULL) {
  if (is.null(newdata)) {
    pred_before <- as.numeric(predict(object$before))
    pred_x_before <- object$before$model$time

    # Remove our extra point before the peak
    pred_before <- pred_before[-length(pred_before)]
    pred_x_before <- pred_x_before[-length(pred_x_before)]

    pred_after_d <- as.numeric(predict(object$after_d))
    pred_x_after_d <- exp(object$after_d$model$`log(time)`)

    pred_after_c <- as.numeric(predict(object$after_c))
    pred_x_after_c <- exp(object$after_c$model$`log(time)`)

    pred_x_after <- unique(c(pred_x_after_d, pred_x_after_c))
    pred_x_after <- pred_x_after[order(pred_x_after)]

    newdata <- list(time = c(pred_x_before, pred_x_after))

#     warning("Without the newdata argument, this function predicts only for all
#             unique values of time. If there are duplicated time values (e.g.
#             manual and continuous), there will be fewer predicted values than times.")
  }

  newdata_nozero <- newdata
  newdata_nozero$time <- ifelse(newdata_nozero$time < 0.00001,
                                yes = 0.00001, no = newdata_nozero$time)

  pred_before <-  as.numeric(predict(object$before, newdata = newdata))
  pred_after_d <- as.numeric(predict(object$after_d, newdata = newdata_nozero))
  pred_after_c <- as.numeric(predict(object$after_c, newdata = newdata_nozero))

  pred_before <- tibble::tibble(time = newdata$time,
                                activity = pred_before)
  pred_before <- dplyr::mutate_all(pred_before,
                                   ~replace(., is.nan(.), 0))
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

  pred_after <- dplyr::mutate_all(pred_after,
                                  ~replace(., is.nan(.), 0))

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


#' Create starting parameters for an exponential blood model
#'
#' This function guesses reasonable starting parameters for exponential blood
#' models. I use the suggestions from page 515 of Pinheiro & Bates Mixed-Effects
#' Models in S and S-PLUS for creating the exponential starting parameters.
#'
#' @param time The time of each measurement.
#' @param activity The radioactivity of each measurement.
#' @param fit_exp3 Should a third exponential be fitted, or is a bi-exponential
#'   fit desired? Default is TRUE.
#' @param expdecay_props What proportions of the decay should be used for
#'   choosing starting parameters for the exponential decay. Defaults to 1/60
#'   and 1/10, i.e. start to 1/60, 1/60 to 1/10 and 1/10 to end. If fitting only
#'   two exponentials, the second value will be used.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @return A list of the starting parameters
#' @export
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' aif <- bd_extract(blooddata, output = "AIF")
#' start <- blmod_exp_startpars(aif$time,
#'                            aif$aif)
blmod_exp_startpars <- function(time, activity, fit_exp3=TRUE,
                                expdecay_props = c(1/60, 0.1)) {

  startpars <- list()

  if(!fit_exp3) {
    expdecay_props[1] <- expdecay_props[2]
  }

  blood <- blmod_tidyinput(time, activity, Method=NULL, weights = NULL)

  # Peaktime
  startpars$peaktime <- blood$time[which.max(blood$activity)][1]

  # Peakval
  startpars$peakval <- blood$activity[blood$time==startpars$peaktime][1]

  # t0
  bloodrise <- blood[blood$time <= startpars$peaktime, ]

  if( nrow(bloodrise) < 2 ) {
    stop("There are an insufficient number of data points in the AIF curve prior to the peak to fit this model.")
  }

  rise_lm <- lm(activity ~ time, data=bloodrise)

  rise_seg <- try(segmented::segmented(rise_lm, npsi=1), silent = T)

  if( !any(class(rise_seg) == "try-error") ) {  # If success
    if(is.numeric(rise_seg$psi[2])) {     # And if numeric

      suppressWarnings(
        interppred <- pracma::interp1(x = bloodrise$time,
                                      y = bloodrise$activity,
                                      xi = rise_seg$psi[2])
      )

      if(interppred < 0.2 * startpars$peakval) { # And if < 20% of peakval
        startpars$t0 <- rise_seg$psi[2]
      } else {
        startpars$t0 <- NA
      }
    } else {
      startpars$t0 <- NA
    }
  } else {
    startpars$t0 <- NA
  }

  #### Old Method
  if( is.na(startpars$t0) ) {

    rise_lm <- lm(activity ~ time, data=bloodrise, weights=activity)
    rise_coef <- coef(rise_lm)
    startpars$t0 <- as.numeric(-1*(rise_coef[1] / rise_coef[2]))

    if(startpars$t0 < 0) { # This can happen if the peak is very dispersed
      bloodrise <- bloodrise[blood$activity >= 0.1*startpars$peakval, ]
      rise_lm <- lm(activity ~ time, data=bloodrise)
      rise_coef <- coef(rise_lm)
      startpars$t0 <- as.numeric(-1*(rise_coef[1]/rise_coef[2]))
    }

    if(startpars$t0 < 0) { # If still less than 0
      startpars$t0 <- 0
    }
  }



  # Decay
  blood_decay <- blood[blood$time > startpars$peaktime,]
  blood_decay$time <- blood_decay$time - startpars$peaktime


  # Exponentials

  ## Third Exponential

  blood_exp_part3 <- dplyr::filter(blood_decay,
                                   dplyr::between(time,
                                                  expdecay_props[2]*max(time),
                                                  max(time)) &
                                     activity>0)

  exp3_mod <- lm(log(abs(activity)) ~ time,
                 data=blood_exp_part3)

  exp3_coef <-  as.numeric(coef(exp3_mod))

  C <- exp(exp3_coef[1])
  gamma <- exp(log(abs(exp3_coef[2])))

  if(fit_exp3) {
    startpars$C     <- C
    startpars$gamma <- gamma
  } else {
    startpars$C     <- 0
    startpars$gamma <- 0

    startpars$B     <- C
    startpars$beta  <-  gamma
  }

  blood_decay$activity_2ex <- blood_decay$activity -
    C*exp( -gamma * (blood_decay$time) )

  ## Second Exponential
  blood_exp_part2 <- dplyr::filter(blood_decay,
                                   dplyr::between(time,
                                                  expdecay_props[1]*max(time),
                                                  expdecay_props[2]*max(time)) &
                                     activity_2ex>0)

  if(fit_exp3) {

    exp2_mod <- lm(log(abs(activity_2ex)) ~ time,
                   data=blood_exp_part2)

    exp2_coef <-  as.numeric(coef(exp2_mod))

    B <- exp(exp2_coef[1])
    beta <- exp(log(abs(exp2_coef[2])))

    startpars$B     <- B
    startpars$beta  <- beta

    blood_decay$activity_1ex <- blood_decay$activity_2ex -
      B*exp( -beta * (blood_decay$time) )

  } else {
    blood_decay$activity_1ex <- blood_decay$activity_2ex
  }

  ## First Exponential
  blood_exp_part1 <- dplyr::filter(blood_decay,
                                   dplyr::between(time,
                                                  min(time),
                                                  expdecay_props[1]*max(time)) &
                                     activity_1ex>0)

  if( nrow(blood_exp_part1) > 3 ) {

  exp1_mod <- lm(log(abs(activity_1ex)) ~ time,
                 data=blood_exp_part1)

  exp1_coef <-  as.numeric(coef(exp1_mod))

  A <- exp(exp1_coef[1])
  alpha <- exp(log(abs(exp1_coef[2])))

  startpars$A <- A
  startpars$alpha <- alpha

  } else {

    startpars$A <- startpars$B
    startpars$alpha <- startpars$beta
  }

  # return the list in correct order
  out <- list(
    t0 = startpars$t0,
    peaktime = startpars$peaktime,
    peakval = startpars$peakval,
    A = startpars$A,
    alpha = startpars$alpha,
    B = startpars$B,
    beta = startpars$beta,
    C = startpars$C,
    gamma = startpars$gamma
  )

  return(out)

}


#' Blood Model: Exponential
#'
#' Fit an exponential model to AIF data with the ability to model the peak. In
#' other words, this model fits a single NLS model which describes the rise and
#' the fall of the AIF. This approach is more flexible, but also more sensitive
#' to fitting failures.
#'
#' This model fits a bi- or tri-exponential model to AIF data. This model can be
#' specified in a fairly large number of ways, depending primarily on the
#' quality of the input data. Reasonably conservative defaults are provided.
#'
#' The model can fit the time zero point (\code{fit_t0}), otherwise the rise
#' starts at time point zero. It can fit two or three exponentials to the curve
#' after the peak (\code{fit_exp0}). When it comes to fitting the peak itself,
#' the time of the peak can be fit or set (\code{fit_peaktime}). When the
#' peaktime is set, it is either set to the time point of the maximal measured
#' value, or it can be set to another value using \code{peaktime_val}. For
#' fitting the peak value, it can be fit as a unique parameter or set using
#' \code{fit_peakval}. If it is set, then it can either be set to the largest
#' measured value, or it can be set to the fitted combination of \code{A+B+C},
#' i.e. the initial part of the decay after the peak.
#'
#'
#' @param time The time of each measurement in seconds
#' @param activity The radioactivity of each measurement
#' @param Method Optional. The method of collection, i.e. "Discrete" or
#'   "Continuous"
#' @param weights Optional. Weights of each measurement.
#' @param fit_t0 Should time point zero be fitted? If not, it is set to 0.
#'   Default is TRUE.
#' @param fit_exp3 Should the third exponential be fitted, or should a
#'   bi-exponential model be used? Default is TRUE for a tri-exponential model.
#' @param fit_peaktime Should the time of the peak be fitted? Default is FALSE.
#'   This is potentially useful for data where sampling frequency was low around
#'   the peak.
#' @param fit_peakval Should the value of the peak be fitted? Default is FALSE.
#'   This is potentially useful for data where sampling frequency was low around
#'   the peak.
#' @param peaktime_val Optional. If \code{fit_peaktime} is FALSE, then this
#'   parameter will define the peaktime. If \code{fit_peaktime} is FALSE and
#'   this parameter is left as NULL, then the empirical measured peaktime will
#'   be used.
#' @param peakval_set Optional. If \code{fit_peakval} is FALSE, then this
#'   parameter defines the peakval as either the empirical measured peakval
#'   (TRUE), or the combination of A+B+c (FALSE).
#' @param lower Optional. The lower limits of the fit. If left as NULL, they
#'   will be given reasonable defaults (mostly 50\% of the starting parameters).
#' @param upper Optional. The upper limits of the fit. If left as NULL, they
#'   will be given reasonable defaults (mostly 150\% of the starting
#'   parameters).
#' @param start Optional. The starting parameters for the fit. If left as NULL,
#'   they will be selected using \code{blmod_exp_startpars}.
#' @param multstart_lower Optional. The lower limits of the starting parameters.
#' @param multstart_upper Optional. The upper limits of the starting parameters.
#' @param multstart_iter The number of fits to perform with different starting
#'   parameters. If set to 1, then the starting parameters will be used for a
#'   single fit.
#' @param taper_weights Should the weights be tapered to gradually trade off
#'   between the continuous and discrete samples after the peak?
#' @param check_startpars Optional. Return only the starting parameters. Useful
#'   for debugging fits which do not work.
#' @param expdecay_props What proportions of the decay should be used for
#'   choosing starting parameters for the exponential decay. Defaults to 1/60
#'   and 1/10, i.e. start to 1/60, 1/60 to 1/10 and 1/10 to end. If fitting only
#'   two exponentials, the second value will be used.
#'
#' @return A model fit including all of the individual parameters, fit details,
#'   and model fit object of class blood_exp.
#' @export
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' aif <- bd_extract(blooddata, output = "AIF")
#' blood_fit <- blmod_exp(aif$time,
#'                            aif$aif,
#'                            Method = aif$Method, multstart_iter = 1)
blmod_exp <- function(time, activity, Method = NULL,
                         weights = NULL,
                         fit_t0=TRUE, fit_exp3=TRUE,
                         fit_peaktime=FALSE, fit_peakval=FALSE,
                         peaktime_val = NULL, peakval_set = TRUE,
                         lower = NULL,
                         upper = NULL,
                         start = NULL,
                         multstart_lower = NULL,
                         multstart_upper = NULL,
                         multstart_iter = 100,
                         taper_weights = TRUE,
                         check_startpars = FALSE,
                         expdecay_props = c(1/60, 0.1)) {


  # Tidy up
  blood <- blmod_tidyinput(time, activity, Method, weights)


  # Create starting parameters
  ## Note: start contains actual starting parameters. Startvals useful for other
  ## things even if starting parameters are defined.
  startvals <- blmod_exp_startpars(time, activity,
                                   fit_exp3,
                                   expdecay_props)
  if(is.null(start)) {
    start <- startvals
  }

  # Fix weights
  if( length(unique(blood$Method)) == 2 ) { # i.e. both cont and discrete
    discrete_before_peak <- which(blood$Method=="Discrete" &
                                    blood$time < startvals$peaktime)

    blood$weights[discrete_before_peak] <- 0
  }

  if(taper_weights) {
    blood$weights <- ifelse(blood$Method=="Continuous",
                            yes = blood$weights * (1-blood$peakfrac),
                            no = blood$weights * blood$peakfrac)
  }


  # Set up start, upper and lower parameter lists
  if(is.null(lower)) {
    lower <- purrr::map(startvals, ~(.x - abs(.x*0.8)))
    lower$t0 <- 0
    lower$peaktime <- startvals$t0

    lower$A <- 0.5*startvals$peakval
    lower$B <- 0
    lower$C <- 0
  }

  if(is.null(upper)) {
    upper <- purrr::map(startvals, ~(.x + abs(.x*2)))
    upper$t0 <- startvals$peaktime
  }


  if(!fit_peaktime) {
    if(is.null(peaktime_val)) {
      peaktime_val <- startvals$peaktime
    }

    lower$peaktime <- NULL
    upper$peaktime <- NULL
    start$peaktime <- NULL
  }

  if (!fit_t0) {
    lower$t0 <- NULL
    upper$t0 <- NULL
    start$t0 <- NULL
  }

  if (!fit_exp3) {
    lower$C <- NULL
    upper$C <- NULL
    start$C <- NULL

    lower$gamma <- NULL
    upper$gamma <- NULL
    start$gamma <- NULL
  }

  if (!fit_peakval) {
    if(peakval_set) {
      formula_peakval <- startvals$peakval
    } else {
      formula_peakval <- "A+B+C"
    }

    lower$peakval <- NULL
    upper$peakval <- NULL
    start$peakval <- NULL
  }

  if(is.null(multstart_lower)) {
    multstart_lower <- lower
  }

  if(is.null(multstart_upper)) {
    multstart_upper <- upper
  }

  lower <- as.numeric(as.data.frame(lower))
  upper <- as.numeric(as.data.frame(upper))

  multstart_lower <- as.numeric(as.data.frame(multstart_lower))
  multstart_upper <- as.numeric(as.data.frame(multstart_upper))


  if(check_startpars) {
    out <- list(start = start, startvals = startvals)
    return(out)
  }


  # Set up the formula
  formula <- paste("activity ~ blmod_triexp_model(time, ",
                   "t0", ifelse(fit_t0, "", "=0"),", ",
                   "peaktime", ifelse(fit_peaktime, "", paste0("=", peaktime_val)),
                   ", ",
                   "peakval", ifelse(fit_peakval, "", paste0("=", formula_peakval)),
                   ", ",
                   "A, alpha, B, beta, ",
                   "C", ifelse(fit_exp3, "", "=0"), ", ",
                   "gamma", ifelse(fit_exp3, "", "=0"), ")",
                   sep = "")


  # Fit the model, with normal NLS or with multstart
  if( multstart_iter == 1 ) {
    modelout <- minpack.lm::nlsLM(as.formula(formula),
                      data = blood,
                      lower = lower,
                      upper = upper,
                      start = start,
                      weights = weights)
  } else {

    modelout_single <- try(minpack.lm::nlsLM(as.formula(formula),
                                  data = blood,
                                  lower = lower,
                                  upper = upper,
                                  start = start,
                                  weights = weights))

    if( !any(class(modelout_single) == "try-error") ) {  # If success
      AIC_single <- AIC(modelout_single)
    } else {
      AIC_single <- Inf
    }

    multmodelout_multi <- nls.multstart::nls_multstart(
      formula = as.formula(formula), modelweights = weights,
      data = blood, iter = multstart_iter,
      start_lower = multstart_lower, start_upper = multstart_upper,
      supp_errors = "Y", lower = lower, upper = upper)

    if( !is.null(multmodelout_multi) ) {  # If success
      AIC_multi <- AIC(multmodelout_multi)
    } else {
        AIC_multi <- Inf
    }

    if( min( c(AIC_single, AIC_multi) ) == Inf ) {
      stop("None of the model fits were able to converge")
    } else {
      if(AIC_multi < AIC_single) {
        modelout <- multmodelout_multi
      } else {
        modelout <- modelout_single
      }
    }


  }

  # Prepare output
  coefficients <- as.data.frame(as.list(coef(modelout)))

  if(!fit_peaktime) {
    coefficients$peaktime <- ifelse(is.null(peaktime_val),
                                    startvals$peaktime,
                                    peaktime_val)
  }

  if (!fit_t0) {
    coefficients$t0 <- 0
  }

  if (!fit_exp3) {
    coefficients$C <- 0
    coefficients$gamma <- 0
  }

  if (!fit_peakval) {

    abc <- with(coefficients, A+B+C)

    coefficients$peakval <- ifelse(peakval_set,
                                   startvals$peakval,
                                   abc)
  }

  coefficients <- coefficients[order(names(coefficients))]

  start <- as.data.frame(start)
  upper <- as.data.frame(as.list(upper), col.names = names(start))
  lower <- as.data.frame(as.list(lower), col.names = names(start))

  fit_details <- as.data.frame(
    list(
      fit_t0 = fit_t0,
      fit_exp3 = fit_exp3,
      fit_peaktime = fit_peaktime,
      fit_peakval = fit_peakval,
      peakval_set = peakval_set
    ))

  out <- list(
    par = coefficients,
    fit = modelout,
    start = start,
    lower = lower,
    upper = upper,
    fit_details = fit_details,
    blood = blood
  )

  class(out) <- c("blood_exp", "blmod", class(out))

  return(out)

}


# #' Fit an exponential model to AIF data using the measured peak
# #'
# #' This model fits a bi- or tri-exponential model to AIF data. This model, in
# #' contrast to the \code{blmod_expsep()} model, uses the measured peak as
# #' the peak, and describes the rise by interpolating between t0 and the peak,
# #' and fits an exponential model only to the fall of the curve. This model is
# #' more robust, but less flexible than the fitpeak model.
# #'
# #' The model can fit two or three exponentials to the curve
# #' after the peak (\code{fit_exp0}). For \code{t0}, \code{peaktime} and
# #' \code{peakval}, the values can be specified by the user, otherwise they are
# #' determined from the data. \code{t0} is determined by fitting a regression
# #' line through the rise, while \code{peaktime} and \code{peakval} are just
# #' the x and y values of the maximal point in the measured data.
# #'
# #'
# #'
# #' @param time The time of each measurement in seconds
# #' @param activity The radioactivity of each measurement
# #' @param Method Optional. The method of collection, i.e. "Discrete" or
# #'   "Continuous"
# #' @param weights Optional. Weights of each measurement.
# #' @param fit_exp3 Should the third exponential be fitted, or should a
# #'   bi-exponential model be used? Default is TRUE for a tri-exponential model.
# #' @param peaktime_val Optional. This allows the user to specify the peaktime
# #'   value which will be used. Otherwise it will be specified from the data.
# #' @param peakval Optional. This allows the user to specify the peak value
# #'   which will be used. Otherwise it will be specified from the data.
# #' @param t0 Optional. This allows the user to specify the t0 value
# #'   which will be used. Otherwise it will be estimated from the data using a
# #'   linear regression.
# #' @param lower Optional. The lower limits of the fit. If left as NULL, they
# #'   will be given reasonable defaults (mostly 50\% of the starting parameters).
# #' @param upper Optional. The upper limits of the fit. If left as NULL, they
# #'   will be given reasonable defaults (mostly 150\% of the starting
# #'   parameters).
# #' @param start Optional. The starting parameters for the fit. If left as NULL,
# #'   they will be selected using \code{blmod_exp_startpars}.
# #' @param multstart_lower Optional. The lower limits of the starting parameters.
# #' @param multstart_upper Optional. The upper limits of the starting parameters.
# #' @param multstart_iter The number of fits to perform with different starting
# #'   parameters. If set to 1, then the starting parameters will be used for a
# #'   single fit.
# #' @param taper_weights Should the weights be tapered to gradually trade off
# #'   between the continuous and discrete samples after the peak?
# #' @param check_startpars Optional. Return only the starting parameters. Useful
# #'   for debugging fits which do not work.
# #' @param expdecay_props What proportions of the decay should be used for
# #'   choosing starting parameters for the exponential decay. Defaults to 1/60
# #'   and 1/10, i.e. start to 1/60, 1/60 to 1/10 and 1/10 to end. If fitting only
# #'   two exponentials, the second value will be used.
# #'
# #' @return A model fit including all of the individual parameters, fit details,
# #'   and model fit object of class blood_exp.
# #' @export
# #'
# #' @examples
# #' blooddata <- create_blooddata_bids(pbr28$jsondata[[3]])
# #' blooddata <- bd_blood_dispcor(blooddata)
# #' aif <- bd_extract(blooddata, output = "AIF")
# #' blood_fit <- blmod_exp_sep(aif$time,
# #'                            aif$aif,
# #'                            Method = aif$Method, multstart_iter = 1)
# blmod_exp_sep <- function(time, activity, Method = NULL,
#                              weights = NULL,
#                              fit_exp3=TRUE,
#                              peaktime_val = NULL, peakval = NULL, t0 = NULL,
#                              lower = NULL,
#                              upper = NULL,
#                              start = NULL,
#                              multstart_lower = NULL,
#                              multstart_upper = NULL,
#                              multstart_iter = 100,
#                              taper_weights = TRUE,
#                              check_startpars = FALSE,
#                              expdecay_props = c(1/60, 0.1)) {
#
#
#   # Tidy up
#   blood <- blmod_tidyinput(time, activity, Method, weights)
#
#
#   # Create starting parameters
#   ## Note: start contains actual starting parameters. Startvals useful for other
#   ## things even if starting parameters are defined.
#   startvals <- blmod_exp_startpars(time, activity,
#                                    fit_exp3,
#                                    expdecay_props)
#   if(is.null(start)) {
#     start <- startvals
#   }
#
#   # Fix weights
#   if( length(unique(blood$Method)) == 2 ) { # i.e. both cont and discrete
#     discrete_before_peak <- which(blood$Method=="Discrete" &
#                                     blood$time < startvals$peaktime)
#
#     blood$weights[discrete_before_peak] <- 0
#   }
#
#   if(taper_weights) {
#     blood$weights <- ifelse(blood$Method=="Continuous",
#                             yes = blood$weights * (1-blood$peakfrac),
#                             no = blood$weights * blood$peakfrac)
#   }
#
#
#
#   # Set up parameter bounds
#
#   ## lower
#   if(is.null(lower)) {
#     lower <- purrr::map(startvals, ~(.x - abs(.x*0.5)))
#
#     lower$A <- 0.5*startvals$peakval
#     lower$B <- 0
#     lower$C <- 0
#   }
#
#   ## upper
#   if(is.null(upper)) {
#     upper <- purrr::map(startvals, ~(.x + abs(.x*0.5)))
#   }
#
#   ## Fitting exp3
#   if (!fit_exp3) {
#     lower$C <- NULL
#     upper$C <- NULL
#     start$C <- NULL
#
#     lower$gamma <- NULL
#     upper$gamma <- NULL
#     start$gamma <- NULL
#   }
#
#   ## Fix starting parameters for sep
#
#   ### Peaktime
#   peaktime_val <- startvals$peaktime
#   lower$peaktime <- NULL
#   upper$peaktime <- NULL
#   start$peaktime <- NULL
#
#   ### t0
#   if(is.null(t0)) {
#     t0 <- startvals$t0
#   }
#   lower$t0 <- NULL
#   upper$t0 <- NULL
#   start$t0 <- NULL
#
#   ### Peakval
#   if(is.null(peakval)) {
#     peakval <- startvals$peakval
#   }
#   lower$peakval <- NULL
#   upper$peakval <- NULL
#   start$peakval <- NULL
#
#   ## multstart
#   if(is.null(multstart_lower)) {
#     multstart_lower <- lower
#   }
#
#   if(is.null(multstart_upper)) {
#     multstart_upper <- upper
#   }
#
#   ## Defining them
#   lower <- as.numeric(as.data.frame(lower))
#   upper <- as.numeric(as.data.frame(upper))
#
#   multstart_lower <- as.numeric(as.data.frame(multstart_lower))
#   multstart_upper <- as.numeric(as.data.frame(multstart_upper))
#
#
#   if(check_startpars) {
#     out <- list(start = start, startvals = startvals)
#     return(out)
#   }
#
#
#   # Set up the formula for the exponential fit
#   expfit_formula <- paste0("activity ~ A * exp(-alpha * (time-",
#                            peaktime_val,"))",
#                            " + B * exp(-beta * (time-",
#                            peaktime_val,"))",
#                            ifelse(fit_exp3, paste0(" + C * exp(-gamma * (time-",
#                                                    peaktime_val,"))"), ""))
#
#
#   # Fit the model, with normal NLS or with multstart
#   if( multstart_iter == 1 ) {
#     modelout <- minpack.lm::nlsLM(as.formula(expfit_formula),
#                                   data = blood,
#                                   lower = lower,
#                                   upper = upper,
#                                   start = start,
#                                   weights = weights,
#                                   control = minpack.lm::nls.lm.control(
#                                     maxiter = 100))
#   } else {
#     modelout <- nls.multstart::nls_multstart(
#       formula = as.formula(expfit_formula), modelweights = weights,
#       data = blood, iter = multstart_iter,
#       start_lower = multstart_lower, start_upper = multstart_upper,
#       supp_errors = "Y", lower = lower, upper = upper)
#
#     if(is.null(modelout)) {
#       stop("None of the model fits were able to converge")
#     }
#   }
#
#   # Prepare output
#   coefficients <- as.data.frame(as.list(coef(modelout)))
#
#   coefficients$peaktime <- peaktime_val
#   coefficients$peakval <- peakval
#
#   coefficients$t0 <- t0
#
#   if (!fit_exp3) {
#     coefficients$C <- 0
#     coefficients$gamma <- 0
#   }
#
#   coefficients <- coefficients[order(names(coefficients))]
#
#   start <- as.data.frame(start)
#   upper <- as.data.frame(as.list(upper), col.names = names(start))
#   lower <- as.data.frame(as.list(lower), col.names = names(start))
#
#   fit_details <- as.data.frame(
#     list(
#       fit_exp3 = fit_exp3
#     ))
#
#   out <- list(
#     par = coefficients,
#     fit = modelout,
#     start = start,
#     lower = lower,
#     upper = upper,
#     fit_details = fit_details,
#     blood = blood
#   )
#
#   class(out) <- c("blood_exp_sep", class(out))
#
#   return(out)
#
# }



#' Tri-exponential model for fitting of arterial input functions
#'
#' This is the model itself for the tri-exponential model of AIF decay, in which
#' both the rise and fall are modelled together in the same nonlinear model.
#'
#' @param time Time of each sample.
#' @param t0 The delay time. This is the point at which the linear rise begins.
#' @param peaktime The time of the peak.
#' @param peakval The radioactivity value of the peak.
#' @param A The multiplier of the first exponential.
#' @param alpha The rate of the first exponential.
#' @param B The multiplier of the second exponential.
#' @param beta The rate of the second exponential.
#' @param C The multiplier of the third exponential.
#' @param gamma The rate of the third exponential.
#'
#' @return Model predictions
#' @export
#'
#' @examples
#' blmod_triexp_model(1:100, 2, 5, 100, 80, 0.5, 15, 0.1, 5, 0.01)
blmod_triexp_model <- function(time, t0, peaktime, peakval, A, alpha,
                               B, beta, C, gamma) {

  tcorr <- time - t0
  peaktime <- peaktime - t0

  t_before <- tcorr[ which(tcorr < 0) ]
  t_after <- tcorr[ which(tcorr >= 0) ]

  t_beforepeak <- t_after[ which(t_after <= peaktime) ]
  t_afterpeak <- t_after[ which(t_after >  peaktime) ]


  # Before t0
  Cpl_0_t0 <- rep(0, length.out = length(t_before))

  # Rise to the peak
  Cpl_t0_peaktime <- (peakval / peaktime) * t_beforepeak

  # Descent
  Cpl_peaktime_end <-
    A * exp(-alpha * (t_afterpeak-peaktime)) +
    B * exp(-beta * (t_afterpeak-peaktime)) +
    C * exp(-gamma * (t_afterpeak-peaktime))

  out <- c(Cpl_0_t0, Cpl_t0_peaktime, Cpl_peaktime_end)

  return(out)

}

#' Predicted values for tri-exponential fits
#'
#' Predicts values for new times for tri-exponential blood fits.
#'
#' @param object Tri-exponential fit.
#' @param newdata A new data list, including times to be predicted for.
#'
#' @return Model predictions
#' @export
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blood <- bd_extract(blooddata, output = "AIF")
#' blood_fit <- blmod_exp(blood$time,
#'                            blood$aif,
#'                            Method = blood$Method)
#'
#' predict_blood_exp(blood_fit, newdata=list(time=1:10))
predict_blood_exp <- function(object, newdata = NULL) {

  if(is.null(newdata)) {
    newdata <- list()
    newdata$time = object$blood$time
  }

  pars <- object$par

  blmod_triexp_model(time = newdata$time,
                     t0 = pars$t0,
                     peaktime = pars$peaktime,
                     peakval = pars$peakval,
                     A = pars$A,
                     alpha = pars$alpha,
                     B = pars$B,
                     beta = pars$beta,
                     C = pars$C,
                     gamma = pars$gamma)

}

# predict.blood_exp_sep <- function(object, newdata = NULL) {
#
#   if(is.null(newdata)) {
#     newdata <- list(
#       time = object$blood$time
#     )
#   }
#
#   # Solve for linear portion
#   gradient <- with(object$par, peakval / ( peaktime - t0))
#   intercept <- -1*(object$par$t0 * gradient)
#
#   # Let's tibble
#   dat <- tibble::tibble(time = newdata$time,
#                 curve = dplyr::case_when(
#                   time < object$par$t0 ~ "Before t0",
#                   time > object$par$peaktime ~ "Descent",
#                   TRUE ~ "Rise"))
#
#   pred <- dplyr::case_when(
#     dat$curve == "Before t0" ~ 0,
#     dat$curve == "Descent"   ~ predict(object$fit, newdata),
#     dat$curve == "Rise"      ~ newdata$time*gradient + intercept
#   )
#
#   return(pred)
#
# }



#' Create starting parameters for an Feng blood model
#'
#' This function estimates reasonable starting parameters for the Feng blood
#' models. I use the suggestions from Tonietto (2012).
#'
#' @param time The time of each measurement.
#' @param activity The radioactivity of each measurement.
#' @param expdecay_props What proportions of the decay should be used for
#'   choosing starting parameters for the exponential decay. Defaults to 0
#'   and 0.5, i.e. the latter two exponentials, B, C, beta, gamma are estimated
#'   using halfway to the end of the decay, and the beginning to halfway through
#'   the decay. The first parameters, A and alpha, are estimated from the ascent.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @return A list of the starting parameters
#' @export
#'
#' @references Tonietto, M. (2012). Methods and models for the characterization
#' of arterial input function in dynamic PET studies. Master thesis.
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' aif <- bd_extract(blooddata, output = "AIF")
#' start <- blmod_feng_startpars(aif$time,
#'                            aif$aif)
blmod_feng_startpars <- function(time, activity,
                                 expdecay_props = c(0, 0.5)) {

  startpars <- list()

  blood <- blmod_tidyinput(time, activity, Method=NULL, weights = NULL)

  # Peaktime
  startpars$peaktime <- blood$time[which.max(blood$activity)][1]

  # Peakval
  startpars$peakval <- blood$activity[blood$time==startpars$peaktime][1]

  # t0
  bloodrise <- blood[blood$time <= startpars$peaktime, ]
  rise_lm <- lm(activity ~ time, data=bloodrise)

  rise_seg <- try(segmented::segmented(rise_lm, npsi=1), silent = T)

  if( !any(class(rise_seg) == "try-error") ) {  # If success
    if(is.numeric(rise_seg$psi[2])) {     # And if numeric

      suppressWarnings(
          interppred <- pracma::interp1(x = bloodrise$time,
                                    y = bloodrise$activity,
                                    xi = rise_seg$psi[2])
      )

      if(interppred < 0.2 * startpars$peakval) { # And if < 20% of peakval
        startpars$t0 <- rise_seg$psi[2]
      } else {
        startpars$t0 <- NA
      }
    } else {
      startpars$t0 <- NA
    }
  } else {
    startpars$t0 <- NA
  }

  #### Old Method
  if( is.na(startpars$t0) ) {

    rise_lm <- lm(activity ~ time, data=bloodrise, weights=activity)
    rise_coef <- coef(rise_lm)
    startpars$t0 <- as.numeric(-1*(rise_coef[1] / rise_coef[2]))

    if(startpars$t0 < 0) { # This can happen if the peak is very dispersed
      bloodrise <- bloodrise[blood$activity >= 0.1*startpars$peakval, ]
      rise_lm <- lm(activity ~ time, data=bloodrise)
      rise_coef <- coef(rise_lm)
      startpars$t0 <- as.numeric(-1*(rise_coef[1]/rise_coef[2]))
    }

    if(startpars$t0 < 0) { # If still less than 0
      startpars$t0 <- 0
    }
  }

    if(startpars$t0 < 0) { # If still less than 0
      startpars$t0 <- 0
    }



  # Decay
  blood_decay <- blood[blood$time > startpars$peaktime,]
  blood_decay$time <- blood_decay$time - startpars$peaktime

  blood_risedecay <- blood
  blood_risedecay$time <- blood_risedecay$time - startpars$peaktime


  # Exponentials

  ## Third Exponential

  blood_exp_part3 <- dplyr::filter(blood_decay,
                                   dplyr::between(time,
                                                  expdecay_props[2]*max(time),
                                                  max(time)) &
                                     activity>0)

  exp3_mod <- lm(log(abs(activity)) ~ time,
                 data=blood_exp_part3)

  exp3_coef <-  as.numeric(coef(exp3_mod))

  C <- exp(exp3_coef[1])
  gamma <- exp(log(abs(exp3_coef[2])))

  startpars$C     <- C
  startpars$gamma <- gamma

  blood_decay$activity_2ex <- blood_decay$activity -
    C*exp( -gamma * (blood_decay$time) )

  blood_risedecay$activity_2ex <- blood_risedecay$activity -
    C*exp( -gamma * (blood_risedecay$time) )

  ## Second Exponential
  blood_exp_part2 <- dplyr::filter(blood_decay,
                                   dplyr::between(time,
                                                  expdecay_props[1]*max(time),
                                                  expdecay_props[2]*max(time)) &
                                     activity_2ex>0)


  exp2_mod <- lm(log(abs(activity_2ex)) ~ time,
                 data=blood_exp_part2)

  exp2_coef <-  as.numeric(coef(exp2_mod))

  B <- exp(exp2_coef[1])
  beta <- exp(log(abs(exp2_coef[2])))

  startpars$B     <- B
  startpars$beta  <- beta

  blood_decay$activity_1ex <- blood_decay$activity_2ex -
    B*exp( -beta * (blood_decay$time) )

  blood_risedecay$activity_1ex <- blood_risedecay$activity_2ex -
    B*exp( -beta * (blood_risedecay$time) )


  ## First Exponential
  # blood_exp_part1 <- dplyr::filter(blood_risedecay,
  #                                  dplyr::between(time,
  #                                                 0,
  #                                                 startpars$peaktime) &
  #                                    activity_1ex>0)

  blood_exp_part1 <- bloodrise %>%
    dplyr::filter(time > 0)

  if( nrow(blood_exp_part1) > 3 ) {

    risetime = startpars$peaktime - startpars$t0


    tmp1 = - ( B + C ) / ( C*exp(risetime*beta) + exp(risetime*gamma) *
                             (B - startpars$peakval*exp(risetime*beta) ) )

    tmp2_exptop = C*exp(risetime*beta)*(-1+risetime*beta) + exp(risetime*gamma) *
      ( B*(-1+risetime*gamma) - exp(risetime*beta)*startpars$peakval*(-1+risetime*(beta + gamma)) )
    tmp2_expbot = C*exp(risetime*beta) + exp(risetime*gamma)*
      (B - startpars$peakval*exp(risetime*beta))
    tmp2 = exp( tmp2_exptop / tmp2_expbot )

    tmp = tmp1 * tmp2

    alpha1 = 1 / ( risetime*( C*exp(risetime*beta) + exp(risetime*gamma)*
                                (B - startpars$peakval*exp(risetime*beta) ) ) )

    alpha2.1 = (exp(risetime*gamma)) * (B - startpars$peakval*exp(risetime*beta) +
                                          B*risetime*beta)
    alpha2.2 = C*exp(risetime*beta)*(1+risetime*gamma)
    alpha2.3 = C*exp(risetime*beta) + ( exp(risetime*gamma) *
                                          (B - startpars$peakval*exp(risetime*beta)) )
    alpha2.4 = LambertW::W(tmp)

    alpha = alpha1 * ( alpha2.1 + alpha2.2 + (alpha2.3 * alpha2.4) )

    A1 = B*(alpha - beta*exp(risetime*(alpha-beta)) ) +
      C*( alpha - gamma*exp(risetime * (alpha-gamma) ) )
    A2 = -1 + risetime*alpha
    A = A1 / A2

    startpars$A <- A
    startpars$alpha <- alpha

  } else {

    startpars$A <- startpars$B
    startpars$alpha <- startpars$beta

    warning("not enough points left to define A and alpha: consider
            revising the expdecay_props")
  }

  # return the list in correct order
  out <- list(
    t0 = startpars$t0,
    peaktime = startpars$peaktime,
    peakval = startpars$peakval,
    A = startpars$A,
    alpha = startpars$alpha,
    B = startpars$B,
    beta = startpars$beta,
    C = startpars$C,
    gamma = startpars$gamma
  )

  return(out)

}



#' Blood Model: Feng
#'
#' Fit a Feng model to AIF data. This model fits the conventional Feng model to
#' AIF data. Note: this model is pretty terrible for just about all tracers
#' except for FDG. The convolved Feng model, fengconv is usually much more
#' effective.
#'
#'
#' @param time The time of each measurement in seconds
#' @param activity The radioactivity of each measurement
#' @param Method Optional. The method of collection, i.e. "Discrete" or
#'   "Continuous"
#' @param weights Optional. Weights of each measurement.
#' @param fit_t0 Should time point zero be fitted? If not, it is set to 0.
#'   Default is TRUE.
#' @param lower Optional. The lower limits of the fit. If left as NULL, they
#'   will be given reasonable defaults (mostly 10\% of the starting parameters).
#' @param upper Optional. The upper limits of the fit. If left as NULL, they
#'   will be given reasonable defaults (mostly 150\% of the starting
#'   parameters).
#' @param start Optional. The starting parameters for the fit. If left as NULL,
#'   they will be selected using \code{blmod_exp_startpars}.
#' @param multstart_lower Optional. The lower limits of the starting parameters.
#' @param multstart_upper Optional. The upper limits of the starting parameters.
#' @param multstart_iter The number of fits to perform with different starting
#'   parameters. If set to 1, then the starting parameters will be used for a
#'   single fit.
#' @param taper_weights Should the weights be tapered to gradually trade off
#'   between the continuous and discrete samples after the peak?
#' @param check_startpars Optional. Return only the starting parameters. Useful
#'   for debugging fits which do not work.
#' @param expdecay_props What proportions of the decay should be used for
#'   choosing starting parameters for the exponential decay. Defaults to 0 and
#'   0.5, i.e. the latter two exponentials, B, C, beta, gamma are estimated
#'   using halfway to the end of the decay, and the beginning to halfway through
#'   the decay. The first parameters, A and alpha, are estimated from the
#'   ascent.
#'
#' @return A model fit including all of the individual parameters, fit details,
#'   and model fit object of class blood_feng.
#' @export
#'
#' @references Wang, Xinmin, and Dagan Feng. "A study on physiological parameter
#'   estimation accuracy for tracer kinetic modeling with positron emission
#'   tomography (pet)." 1992 American Control Conference. IEEE, 1992.
#'   Feng, D., Z. Wang, and S. C. Huang. "Tracer plasma time-activity curves in
#'   circulatory system for positron emission tomography kinetic modeling
#'   studies." IFAC Proceedings Volumes 26.2 (1993): 175-178.
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' aif <- bd_extract(blooddata, output = "AIF")
#' blood_fit <- blmod_feng(aif$time,
#'                            aif$aif,
#'                            Method = aif$Method, multstart_iter = 1)
blmod_feng <- function(time, activity, Method = NULL,
                      weights = NULL,
                      fit_t0=TRUE,
                      lower = NULL,
                      upper = NULL,
                      start = NULL,
                      multstart_lower = NULL,
                      multstart_upper = NULL,
                      multstart_iter = 100,
                      taper_weights = TRUE,
                      check_startpars = FALSE,
                      expdecay_props = c(1/60, 0.1)) {


  # Tidy up
  blood <- blmod_tidyinput(time, activity, Method, weights)


  # Create starting parameters
  ## Note: start contains actual starting parameters. Startvals useful for other
  ## things even if starting parameters are defined.
  startvals <- blmod_feng_startpars(time, activity,
                                   expdecay_props)
  if(is.null(start)) {
    start <- startvals
  }

  # Fix weights
  if( length(unique(blood$Method)) == 2 ) { # i.e. both cont and discrete
    discrete_before_peak <- which(blood$Method=="Discrete" &
                                    blood$time < startvals$peaktime)

    blood$weights[discrete_before_peak] <- 0
  }

  if(taper_weights) {
    blood$weights <- ifelse(blood$Method=="Continuous",
                            yes = blood$weights * (1-blood$peakfrac),
                            no = blood$weights * blood$peakfrac)
  }


  # Set up start, upper and lower parameter lists
  if(is.null(lower)) {
    lower <- purrr::map(startvals, ~(.x - abs(.x*0.5)))
    lower$t0 <- 0
    lower$peaktime <- startvals$t0

    lower$A <- 0.1*startvals$A
    lower$B <- 0.1*startvals$B
    lower$C <- 0.1*startvals$C
  }

  if(is.null(upper)) {
    upper <- purrr::map(startvals, ~(20*abs(.x)))
    upper$t0 <- startvals$peaktime
  }


  lower$peaktime <- NULL
  upper$peaktime <- NULL
  start$peaktime <- NULL

  lower$peakval <- NULL
  upper$peakval <- NULL
  start$peakval <- NULL

  if (!fit_t0) {
    lower$t0 <- NULL
    upper$t0 <- NULL
    start$t0 <- NULL
  }

  if(is.null(multstart_lower)) {
    multstart_lower <- lower
  }

  if(is.null(multstart_upper)) {
    multstart_upper <- upper
  }

  lower <- as.numeric(as.data.frame(lower))
  upper <- as.numeric(as.data.frame(upper))

  multstart_lower <- as.numeric(as.data.frame(multstart_lower))
  multstart_upper <- as.numeric(as.data.frame(multstart_upper))


  if(check_startpars) {
    out <- list(start = start, startvals = startvals)
    return(out)
  }


  # Set up the formula
  formula <- paste("activity ~ blmod_feng_model(time, ",
                   "t0", ifelse(fit_t0, "", "=0"),", ",
                   "A, alpha, B, beta, C, gamma)",
                   sep = "")


  # Fit the model, with normal NLS or with multstart
  if( multstart_iter == 1 ) {
    modelout <- minpack.lm::nlsLM(as.formula(formula),
                                  data = blood,
                                  lower = lower,
                                  upper = upper,
                                  start = start,
                                  weights = weights)
  } else {

    modelout_single <- try(minpack.lm::nlsLM(as.formula(formula),
                                             data = blood,
                                             lower = lower,
                                             upper = upper,
                                             start = start,
                                             weights = weights))

    if( !any(class(modelout_single) == "try-error") ) {  # If success
      AIC_single <- AIC(modelout_single)
    } else {
      AIC_single <- Inf
    }

    multmodelout_multi <- nls.multstart::nls_multstart(
      formula = as.formula(formula), modelweights = weights,
      data = blood, iter = multstart_iter,
      start_lower = multstart_lower, start_upper = multstart_upper,
      supp_errors = "Y", lower = lower, upper = upper)

    if( !is.null(multmodelout_multi) ) {  # If success
      AIC_multi <- AIC(multmodelout_multi)
    } else {
      AIC_multi <- Inf
    }

    if( min( c(AIC_single, AIC_multi) ) == Inf ) {
      stop("None of the model fits were able to converge")
    } else {
      if(AIC_multi < AIC_single) {
        modelout <- multmodelout_multi
      } else {
        modelout <- modelout_single
      }
    }


  }

  # Prepare output
  coefficients <- as.data.frame(as.list(coef(modelout)))


  if (!fit_t0) {
    coefficients$t0 <- 0
  }

  coefficients <- coefficients[order(names(coefficients))]

  start <- as.data.frame(start)
  upper <- as.data.frame(as.list(upper), col.names = names(start))
  lower <- as.data.frame(as.list(lower), col.names = names(start))

  fit_details <- as.data.frame(
    list(
      fit_t0 = fit_t0
    ))

  out <- list(
    par = coefficients,
    fit = modelout,
    start = start,
    lower = lower,
    upper = upper,
    fit_details = fit_details,
    blood = blood
  )

  class(out) <- c("blood_feng", "blmod", class(out))

  return(out)

}


#' Feng model for fitting of arterial input functions
#'
#' This is the model itself for the Feng model of the AIF.
#'
#' @param time Time of each sample.
#' @param t0 The delay time. This is the point at which the linear rise begins.
#' @param A The multiplier of the first exponential.
#' @param alpha The rate of the first exponential.
#' @param B The multiplier of the second exponential.
#' @param beta The rate of the second exponential.
#' @param C The multiplier of the third exponential.
#' @param gamma The rate of the third exponential.
#'
#' @return Model predictions
#' @export
#'
#' @examples
#' blmod_feng_model(1:1000, 30, 7, 0.02, 2, 0.005, 4.5, 0.0005)
blmod_feng_model <- function(time, t0, A, alpha, B, beta, C, gamma) {

  tcorr <- time - t0

  t_before <- tcorr[ which(tcorr <= 0) ]
  t_after <- tcorr[ which(tcorr > 0) ]

  # Before t0
  Cpl_0_t0 <- rep(0, length.out = length(t_before))

  # After t0
  term1 = ( A*t_after - B - C )*exp(-alpha*t_after)
  term2 = B * exp(-beta*t_after)
  term3 = C * exp(-gamma*t_after)

  Cpl_t0_end = term1 + term2 + term3

  out = c( Cpl_0_t0 , Cpl_t0_end )

  return(out)
}

#' Predicted values for Feng fits
#'
#' Predicts values for new times for Feng blood fits.
#'
#' @param object Feng fit.
#' @param newdata A new data list, including times to be predicted for.
#'
#' @return Model predictions
#' @export
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blood <- bd_extract(blooddata, output = "AIF")
#' blood_fit <- blmod_feng(blood$time,
#'                            blood$aif,
#'                            Method = blood$Method)
#'
#' predict_blood_feng(blood_fit, newdata=list(time=1:10))
predict_blood_feng <- function(object, newdata = NULL) {

  if(is.null(newdata)) {
    newdata <- list()
    newdata$time = object$blood$time
  }

  pars <- object$par

  blmod_feng_model(time = newdata$time,
                     t0 = pars$t0,
                     A = pars$A,
                     alpha = pars$alpha,
                     B = pars$B,
                     beta = pars$beta,
                     C = pars$C,
                     gamma = pars$gamma)

}



#' Blood Model: FengConv
#'
#' Fits a Feng model, convolved with an infusion injection to AIF data. This
#' model is significantly more effective than the original Feng model in the
#' vast majority of applications. The infusion duration can be a little bit
#' tricky to fit. I recommend trying to find out the approximate duration of the
#' injection and including this as a fixed parameter for most stable
#' performance.
#'
#' @param time The time of each measurement in seconds
#' @param activity The radioactivity of each measurement
#' @param inftime The infusion time during which the tracer is administered. It
#'   is recommended to provide this value if known. Leaving this argument empty
#'   will result in the inftime being fitted. Alternatively, two values can be
#'   specified, e.g. inftime=c(25, 35) to fit the infusion time with specified
#'   upper and lower limits if the known infusion time is approximate.
#' @param Method Optional. The method of collection, i.e. "Discrete" or
#'   "Continuous"
#' @param weights Optional. Weights of each measurement.
#' @param fit_t0 Should time point zero be fitted? If not, it is set to 0.
#'   Default is TRUE.
#' @param lower Optional. The lower limits of the fit. If left as NULL, they
#'   will be given reasonable defaults (mostly 10\% of the starting parameters).
#' @param upper Optional. The upper limits of the fit. If left as NULL, they
#'   will be given reasonable defaults (mostly 20x the starting parameters).
#' @param start Optional. The starting parameters for the fit. If left as NULL,
#'   they will be selected using \code{blmod_feng_startpars}.
#' @param multstart_lower Optional. The lower limits of the starting parameters.
#' @param multstart_upper Optional. The upper limits of the starting parameters.
#' @param multstart_iter The number of fits to perform with different starting
#'   parameters. If set to 1, then the starting parameters will be used for a
#'   single fit.
#' @param taper_weights Should the weights be tapered to gradually trade off
#'   between the continuous and discrete samples after the peak?
#' @param check_startpars Optional. Return only the starting parameters. Useful
#'   for debugging fits which do not work.
#' @param expdecay_props What proportions of the decay should be used for
#'   choosing starting parameters for the exponential decay. Defaults to 0 and
#'   0.5, i.e. the latter two exponentials, B, C, beta, gamma are estimated
#'   using halfway to the end of the decay, and the beginning to halfway through
#'   the decay. The first parameters, A and alpha, are estimated from the
#'   ascent.
#'
#' @return A model fit including all of the individual parameters, fit details,
#'   and model fit object of class blood_fengconv.
#' @export
#'
#' @references Wong, Koon-Pong, Sung-Cheng Huang, and Michael J. Fulham.
#'   "Evaluation of an input function model that incorporates the injection
#'   schedule in FDG-PET studies." 2006 IEEE Nuclear Science Symposium
#'   Conference Record. Vol. 4. IEEE, 2006. Tonietto, M. (2012). Methods and
#'   models for the characterization of arterial input function in dynamic PET
#'   studies. Master thesis.
#'
#'
#' @examples
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' aif <- bd_extract(blooddata, output = "AIF")
#' blood_fit <- blmod_fengconv(aif$time,
#'                            aif$aif,
#'                            Method = aif$Method,
#'                            multstart_iter = 1,
#'                            inftime = 22)
blmod_fengconv <- function(time, activity, inftime = NULL,
                           Method = NULL,
                           weights = NULL,
                           fit_t0=TRUE,
                           lower = NULL,
                           upper = NULL,
                           start = NULL,
                           multstart_lower = NULL,
                           multstart_upper = NULL,
                           multstart_iter = 500,
                           taper_weights = TRUE,
                           check_startpars = FALSE,
                           expdecay_props = c(0, 0.5)) {


  # Tidy up
  blood <- blmod_tidyinput(time, activity, Method, weights)


  # Create starting parameters
  ## Note: start contains actual starting parameters. Startvals useful for other
  ## things even if starting parameters are defined.
  startvals <- blmod_feng_startpars(time, activity,
                                    expdecay_props)
  if(is.null(start)) {
    start <- startvals
  }

  # Fix weights
  if( length(unique(blood$Method)) == 2 ) { # i.e. both cont and discrete
    discrete_before_peak <- which(blood$Method=="Discrete" &
                                    blood$time < startvals$peaktime)

    blood$weights[discrete_before_peak] <- 0
  }

  if(taper_weights) {
    blood$weights <- ifelse(blood$Method=="Continuous",
                            yes = blood$weights * (1-blood$peakfrac),
                            no = blood$weights * blood$peakfrac)
  }


  # Set up start, upper and lower parameter lists
  if(is.null(lower)) {
    lower <- purrr::map(startvals, ~(.x - abs(.x*0.05)))
    lower$t0 <- 0
    lower$peaktime <- startvals$t0

    lower$A <- 0
    lower$B <- 0
    lower$C <- 0

    lower$alpha <- 0
    lower$beta <- 0
    lower$gamma <- 0
  }

  if(is.null(upper)) {
    upper <- purrr::map(startvals, ~(100*.x))
    upper$t0 <- startvals$peaktime
  }


  lower$peaktime <- NULL
  upper$peaktime <- NULL
  start$peaktime <- NULL

  lower$peakval <- NULL
  upper$peakval <- NULL
  start$peakval <- NULL

  if (!fit_t0) {
    lower$t0 <- NULL
    upper$t0 <- NULL
    start$t0 <- NULL
  }

  if (is.null(inftime)) {

    if(max(blood$time) > 300) {

      # Assuming seconds
      start$ti <- 30
      lower$ti <- 0
      upper$ti <- 300

    } else {

      # Assuming minutes
      start$ti <- 0.5
      lower$ti <- 0
      upper$ti <- 5

    }
  }

  if(length(inftime == 2)) {

    lower$ti <- min(inftime)
    upper$ti <- max(inftime)
    start$ti <- mean(inftime)

    inftime <- NULL

  }

  if(is.null(multstart_lower)) {
    multstart_lower <- lower
  }

  if(is.null(multstart_upper)) {
    multstart_upper <- upper
  }

  lower <- as.numeric(as.data.frame(lower))
  upper <- as.numeric(as.data.frame(upper))

  multstart_lower <- as.numeric(as.data.frame(multstart_lower))
  multstart_upper <- as.numeric(as.data.frame(multstart_upper))


  if(check_startpars) {
    out <- list(start = start, startvals = startvals)
    return(out)
  }


  # Set up the formula
  formula <- paste("activity ~ blmod_fengconv_model(time, ",
                   "t0", ifelse(fit_t0, "", "=0"),", ",
                   "A, alpha, B, beta, C, gamma, ",
                   "ti", ifelse(is.null(inftime), "", paste0("=", inftime)),")",
                   sep = "")


  # Fit the model, with normal NLS or with multstart
  if( multstart_iter == 1 ) {
    modelout <- minpack.lm::nlsLM(as.formula(formula),
                                  data = blood,
                                  lower = lower,
                                  upper = upper,
                                  start = start,
                                  weights = weights,
                                  control  = minpack.lm::nls.lm.control(maxiter = 200))
  } else {

    modelout_single <- try(minpack.lm::nlsLM(as.formula(formula),
                                             data = blood,
                                             lower = lower,
                                             upper = upper,
                                             start = start,
                                             weights = weights,
                                             control  = minpack.lm::nls.lm.control(
                                               maxiter = 200)))

    if( !any(class(modelout_single) == "try-error") ) {  # If success
      AIC_single <- AIC(modelout_single)
    } else {
      AIC_single <- Inf
    }

    multmodelout_multi <- nls.multstart::nls_multstart(
      formula = as.formula(formula), modelweights = weights,
      data = blood, iter = multstart_iter,
      start_lower = multstart_lower, start_upper = multstart_upper,
      supp_errors = "Y",
      lower = lower, upper = upper,
      convergence_count = FALSE,
      control  = minpack.lm::nls.lm.control(
        maxiter = 200))

    if( !is.null(multmodelout_multi) ) {  # If success
      AIC_multi <- AIC(multmodelout_multi)
    } else {
      AIC_multi <- Inf
    }

    if( min( c(AIC_single, AIC_multi) ) == Inf ) {
      stop("None of the model fits were able to converge")
    } else {
      if(AIC_multi < AIC_single) {
        modelout <- multmodelout_multi
      } else {
        modelout <- modelout_single
      }
    }


  }

  # Prepare output
  coefficients <- as.data.frame(as.list(coef(modelout)))


  if (!fit_t0) {
    coefficients$t0 <- 0
  }

  if (!is.null(inftime)) {
    coefficients$ti <- inftime
  }

  coefficients <- coefficients[order(names(coefficients))]

  start <- as.data.frame(start)
  upper <- as.data.frame(as.list(upper), col.names = names(start))
  lower <- as.data.frame(as.list(lower), col.names = names(start))

  fit_details <- as.data.frame(
    list(
      fit_t0 = fit_t0,
      fit_inftime = is.null(inftime)
    ))

  out <- list(
    par = coefficients,
    fit = modelout,
    start = start,
    lower = lower,
    upper = upper,
    fit_details = fit_details,
    blood = blood
  )

  class(out) <- c("blood_fengconv", "blmod", class(out))

  return(out)

}

#' Convolved Feng model for fitting of arterial input functions
#'
#' This is the model itself for the Fengconv model of the AIF with an extended
#' injection.
#'
#' @param time Time of each sample.
#' @param t0 The delay time. This is the point at which the linear rise begins.
#' @param A The multiplier of the first exponential.
#' @param alpha The rate of the first exponential.
#' @param B The multiplier of the second exponential.
#' @param beta The rate of the second exponential.
#' @param C The multiplier of the third exponential.
#' @param gamma The rate of the third exponential.
#' @param ti The infusion time of the bolus.
#'
#' @return Model predictions
#' @export
#'
#' @examples
#' blmod_fengconv_model(1:1000, 30, 220, 0.4, 100, 0.05, 22, 0.003, 30)
blmod_fengconv_model <- function(time, t0, A, alpha, B, beta, C, gamma, ti) {

  # Integral
  g <- function(time, A, alpha, B, beta, C, gamma) {
    term1 <- ( (A - alpha*(B + C)) / alpha^2 ) * (1-exp(-alpha*time))
    term2 <- (A / alpha) * time * exp(-alpha * time)
    term3 <- (B / beta) * (1 - exp(-beta * time))
    term4 <- (C / gamma) * (1 - exp(-gamma * time))

    out <- term1 - term2 + term3 + term4
    return(out)
  }

  tcorr <- time - t0

  # Before and after t0
  t_before <- tcorr[ which(tcorr <= 0) ]
  t_after <- tcorr[ which(tcorr > 0) ]


  # Before and after t0 + ti
  t_middle <- t_after[which(t_after <= ti)]
  t_end <- t_after[which(t_after > ti)]
  t_end_corr <- t_end - ti

  # Before t0
  Cpl_0_t0 <- rep(0, length.out = length(t_before))

  # After t0 but before t0+ti
  Cpl_t0_t0ti <- g(t_middle,
                   A, alpha, B, beta, C, gamma)


  # After t0
  Cpl_t0ti_end <- g(t_end,
                    A, alpha, B, beta, C, gamma) -
    g(t_end_corr,
      A, alpha, B, beta, C, gamma)


  out <- c(Cpl_0_t0, Cpl_t0_t0ti, Cpl_t0ti_end)

  out <- out / ti

  return(out)
}


#' Predicted values for Fengconv fits
#'
#' Predicts values for new times for Fengconv blood fits.
#'
#' @param object Fengconv fit.
#' @param newdata A new data list, including times to be predicted for.
#'
#' @return Model predictions
#' @export
#'
#' @examples
#' \dontrun{
#' blooddata <- pbr28$blooddata[[1]]
#' blood <- bd_extract(blooddata, output = "AIF")
#' blood_fit <- blmod_fengconv(blood$time,
#'                            blood$aif,
#'                            Method = blood$Method)
#'
#' predict_blood_fengconv(blood_fit, newdata=list(time=1:10))
#' }
predict_blood_fengconv <- function(object, newdata = NULL) {

  if(is.null(newdata)) {
    newdata <- list()
    newdata$time = object$blood$time
  }

  pars <- object$par

  blmod_fengconv_model(time = newdata$time,
                   t0 = pars$t0,
                   A = pars$A,
                   alpha = pars$alpha,
                   B = pars$B,
                   beta = pars$beta,
                   C = pars$C,
                   gamma = pars$gamma,
                   ti = pars$ti)

}


#' Blood Model: FengConv plus a rise term
#'
#' Fits a Feng model, convolved with an infusion injection, and then with an
#' additional rise to an asymptote term, to AIF data. This model is a slight
#' modification of the fengconv model which just allows the model to be a little
#' bit more capable of catching a rising tail: the model was designed for
#' bolus+infusion paradigms - however it should be stated that FengConv is
#' usually pretty capable in these cases too. By bolus+infusion, I mean a bolus,
#' injected as a short infusion, followed by a long infusion. The infusion
#' duration can be a little bit tricky to fit. I recommend trying to find out
#' the approximate duration of the injection and including this as a fixed
#' parameter for most stable performance.
#'
#'
#' @param time The time of each measurement in seconds
#' @param activity The radioactivity of each measurement
#' @param inftime The infusion time during which the tracer is administered. It
#'   is recommended to provide this value if known. Leaving this argument empty
#'   will result in the inftime being fitted. Alternatively, two values can be
#'   specified, e.g. inftime=c(25, 35) to fit the infusion time with specified
#'   upper and lower limits if the known infusion time is approximate.
#' @param Method Optional. The method of collection, i.e. "Discrete" or
#'   "Continuous"
#' @param weights Optional. Weights of each measurement.
#' @param fit_t0 Should time point zero be fitted? If not, it is set to 0.
#'   Default is TRUE.
#' @param lower Optional. The lower limits of the fit. If left as NULL, they
#'   will be given reasonable defaults (mostly 50\% of the starting parameters).
#' @param upper Optional. The upper limits of the fit. If left as NULL, they
#'   will be given reasonable defaults (mostly 150\% of the starting
#'   parameters).
#' @param start Optional. The starting parameters for the fit. If left as NULL,
#'   they will be selected using \code{blmod_exp_startpars}.
#' @param multstart_lower Optional. The lower limits of the starting parameters.
#' @param multstart_upper Optional. The upper limits of the starting parameters.
#' @param multstart_iter The number of fits to perform with different starting
#'   parameters. If set to 1, then the starting parameters will be used for a
#'   single fit.
#' @param taper_weights Should the weights be tapered to gradually trade off
#'   between the continuous and discrete samples after the peak?
#' @param check_startpars Optional. Return only the starting parameters. Useful
#'   for debugging fits which do not work.
#' @param expdecay_props What proportions of the decay should be used for
#'   choosing starting parameters for the exponential decay. Defaults to 0 and
#'   0.5, i.e. the latter two exponentials, B, C, beta, gamma are estimated
#'   using halfway to the end of the decay, and the beginning to halfway through
#'   the decay. The first parameters, A and alpha, are estimated from the
#'   ascent.
#'
#' @return A model fit including all of the individual parameters, fit details,
#'   and model fit object of class blood_fengconvplus.
#' @export
#'
#' @examples
#' \dontrun{
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' aif <- bd_extract(blooddata, output = "AIF")
#' blood_fit <- blmod_fengconvplus(aif$time,
#'                            aif$aif,
#'                            Method = aif$Method,
#'                            multstart_iter = 100)
#' }
blmod_fengconvplus <- function(time, activity, inftime = NULL,
                           Method = NULL,
                           weights = NULL,
                           fit_t0=TRUE,
                           lower = NULL,
                           upper = NULL,
                           start = NULL,
                           multstart_lower = NULL,
                           multstart_upper = NULL,
                           multstart_iter = 500,
                           taper_weights = TRUE,
                           check_startpars = FALSE,
                           expdecay_props = c(0, 0.5)) {


  # Tidy up
  blood <- blmod_tidyinput(time, activity, Method, weights)


  # Create starting parameters
  ## Note: start contains actual starting parameters. Startvals useful for other
  ## things even if starting parameters are defined.
  startvals <- blmod_feng_startpars(time, activity,
                                    expdecay_props)
  if(is.null(start)) {
    start <- startvals
  }

  # Fix weights
  if( length(unique(blood$Method)) == 2 ) { # i.e. both cont and discrete
    discrete_before_peak <- which(blood$Method=="Discrete" &
                                    blood$time < startvals$peaktime)

    blood$weights[discrete_before_peak] <- 0
  }

  if(taper_weights) {
    blood$weights <- ifelse(blood$Method=="Continuous",
                            yes = blood$weights * (1-blood$peakfrac),
                            no = blood$weights * blood$peakfrac)
  }


  # Set up start, upper and lower parameter lists
  if(is.null(lower)) {
    lower <- purrr::map(startvals, ~(.x - abs(.x*0.5)))
    lower$t0 <- 0
    lower$peaktime <- startvals$t0

    lower$A <- 0
    lower$B <- 0
    lower$C <- 0

    lower$alpha <- 0
    lower$beta <- 0
    lower$gamma <- 0

    lower$asymptote <- 0
    lower$slope <- 0
  }

  if(is.null(upper)) {
    upper <- purrr::map(startvals, ~(100*.x))
    upper$t0 <- startvals$peaktime

    upper$asymptote <- 20*median(tail(activity, 3))
    upper$slope <- 0.05
  }


  start$asymptote <- 0.5*median(tail(activity, 3))
  start$slope <- 0.1

  lower$peaktime <- NULL
  upper$peaktime <- NULL
  start$peaktime <- NULL

  lower$peakval <- NULL
  upper$peakval <- NULL
  start$peakval <- NULL

  if (!fit_t0) {
    lower$t0 <- NULL
    upper$t0 <- NULL
    start$t0 <- NULL
  }

  if (is.null(inftime)) {
    start$ti <- 30
    lower$ti <- 1
    upper$ti <- 300
  }

  if(length(inftime == 2)) {

    lower$ti <- min(inftime)
    upper$ti <- max(inftime)
    start$ti <- mean(inftime)

    inftime <- NULL

  }

  if(is.null(multstart_lower)) {
    multstart_lower <- lower
  }

  if(is.null(multstart_upper)) {
    multstart_upper <- upper
  }

  lower <- as.numeric(as.data.frame(lower))
  upper <- as.numeric(as.data.frame(upper))

  multstart_lower <- as.numeric(as.data.frame(multstart_lower))
  multstart_upper <- as.numeric(as.data.frame(multstart_upper))


  if(check_startpars) {
    out <- list(start = start, startvals = startvals)
    return(out)
  }


  # Set up the formula
  formula <- paste("activity ~ blmod_fengconvplus_model(time, ",
                   "t0", ifelse(fit_t0, "", "=0"),", ",
                   "A, alpha, B, beta, C, gamma, ",
                   "ti", ifelse(is.null(inftime), "", paste0("=", inftime)),", ",
                   "asymptote, slope",")",
                   sep = "")


  # Fit the model, with normal NLS or with multstart
  if( multstart_iter == 1 ) {
    modelout <- minpack.lm::nlsLM(as.formula(formula),
                                  data = blood,
                                  lower = lower,
                                  upper = upper,
                                  start = start,
                                  weights = weights,
                                  control  = minpack.lm::nls.lm.control(maxiter = 200))
  } else {

    modelout_single <- try(minpack.lm::nlsLM(as.formula(formula),
                                             data = blood,
                                             lower = lower,
                                             upper = upper,
                                             start = start,
                                             weights = weights,
                                             control  = minpack.lm::nls.lm.control(
                                               maxiter = 200)))

    if( !any(class(modelout_single) == "try-error") ) {  # If success
      AIC_single <- AIC(modelout_single)
    } else {
      AIC_single <- Inf
    }

    multmodelout_multi <- nls.multstart::nls_multstart(
      formula = as.formula(formula), modelweights = weights,
      data = blood, iter = multstart_iter,
      start_lower = multstart_lower, start_upper = multstart_upper,
      supp_errors = "Y", lower = lower, upper = upper,
      convergence_count = FALSE,
      control  = minpack.lm::nls.lm.control(
        maxiter = 200))

    if( !is.null(multmodelout_multi) ) {  # If success
      AIC_multi <- AIC(multmodelout_multi)
    } else {
      AIC_multi <- Inf
    }

    if( min( c(AIC_single, AIC_multi) ) == Inf ) {
      stop("None of the model fits were able to converge")
    } else {
      if(AIC_multi < AIC_single) {
        modelout <- multmodelout_multi
      } else {
        modelout <- modelout_single
      }
    }


  }

  # Prepare output
  coefficients <- as.data.frame(as.list(coef(modelout)))


  if (!fit_t0) {
    coefficients$t0 <- 0
  }

  if (!is.null(inftime)) {
    coefficients$ti <- inftime
  }

  coefficients <- coefficients[order(names(coefficients))]

  start <- as.data.frame(start)
  upper <- as.data.frame(as.list(upper), col.names = names(start))
  lower <- as.data.frame(as.list(lower), col.names = names(start))

  fit_details <- as.data.frame(
    list(
      fit_t0 = fit_t0,
      fit_inftime = is.null(inftime)
    ))

  out <- list(
    par = coefficients,
    fit = modelout,
    start = start,
    lower = lower,
    upper = upper,
    fit_details = fit_details,
    blood = blood
  )

  class(out) <- c("blood_fengconvplus", "blmod", class(out))

  return(out)

}

#' Convolved Feng model for fitting of arterial input functions with an additional rise term
#'
#' This is the model itself for the Fengconvplus model of the AIF with an extended
#' injection and additional rising term.
#'
#' @param time Time of each sample.
#' @param t0 The delay time. This is the point at which the linear rise begins.
#' @param A The multiplier of the first exponential.
#' @param alpha The rate of the first exponential.
#' @param B The multiplier of the second exponential.
#' @param beta The rate of the second exponential.
#' @param C The multiplier of the third exponential.
#' @param gamma The rate of the third exponential.
#' @param ti The infusion time of the bolus.
#' @param asymptote The value of the asymptote to which the rise rises.
#' @param slope The slope of the rise term, i.e. the rate at which it rises.
#'
#' @return Model predictions
#' @export
#'
#' @examples
#' blmod_fengconvplus_model(1:1000, 30, 220, 0.4, 100, 0.05, 22, 0.003, 30, 40, 0.001)
blmod_fengconvplus_model <- function(time, t0, A, alpha, B, beta, C, gamma, ti,
                                     asymptote, slope) {

  fengconv_out <- blmod_fengconv_model(time, t0, A, alpha, B, beta, C, gamma, ti)

  risefunc <- ifelse(time > t0,
                     yes=-asymptote + (2*asymptote)/(1 + exp(-slope*(time-t0))),
                     no = 0)

  out <- fengconv_out + risefunc

  return(out)
}

#' Predicted values for Fengconvplus fits
#'
#' Predicts values for new times for Fengconvplus blood fits.
#'
#' @param object Fengconvplus fit.
#' @param newdata A new data list, including times to be predicted for.
#'
#' @return Model predictions
#' @export
#'
#' @examples
#' \dontrun{
#' blooddata <- pbr28$blooddata[[1]]
#' blood <- bd_extract(blooddata, output = "AIF")
#' blood_fit <- blmod_fengconvplus(blood$time,
#'                            blood$aif,
#'                            Method = blood$Method)
#'
#' predict_blood_fengconvplus(blood_fit, newdata=list(time=1:10))
#' }
predict_blood_fengconvplus <- function(object, newdata = NULL) {

  if(is.null(newdata)) {
    newdata <- list()
    newdata$time = object$blood$time
  }

  pars <- object$par

  blmod_fengconvplus_model(time = newdata$time,
                       t0 = pars$t0,
                       A = pars$A,
                       alpha = pars$alpha,
                       B = pars$B,
                       beta = pars$beta,
                       C = pars$C,
                       gamma = pars$gamma,
                       ti = pars$ti,
                       asymptote = pars$asymptote,
                       slope = pars$slope)

}
