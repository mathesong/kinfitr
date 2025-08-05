#' Convolution Function
#'
#' Function to convolve two vectors of equal length for kinetic modelling.
#'
#' @param a A numeric vector of TAC function at equally spaced intervals
#' @param b A numeric vector of TAC function at equally spaced intervals of equal length to \code{a}
#' @param stepsize The time interval between each value of \code{a} and \code{b}
#'
#' @return A numeric vector of the convolved function of the same length as \code{a} and \code{b}
#'
#' @examples
#' a <- rnorm(20)
#' b <- rnorm(20)
#' stepsize <- 1
#'
#' kinfit_convolve(a, b, stepsize)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

kinfit_convolve <- function(a, b, stepsize) {
  if (length(a) != length(b)) {
    stop("The vectors are not of the same length")
  }

  out <- pracma::conv(a, b)[1:length(a)] * stepsize
  return(out)
}

kinfit_uninterpbydur <- function(interptime, interppred, times, durtimes) {
  starttimes <- times - durtimes / 2
  interpstart <- sapply(starttimes, function(x) which.min(abs(interptime - x)))[-1]
  interpend <- c((interpstart - 1)[-1], length(interptime))

  startstop <- data.frame(start = interpstart, stop = interpend)

  predvals <- c(0, mapply(function(start, stop) mean(interppred[start:stop]), startstop$start, startstop$stop))

  return(predvals)
}

zerotrim_input <- function(t_inp, inp) {

  # Create a zero point by drawing a regression between the peak and the peak second differential

  input <- data.frame(Time = t_inp, Value = inp)

  input$Time <- input$Time

  input_start <- subset(input, input$Time < which.max(input$Value))
  # qplot(input_start$Time, input_start$Value)
  input_start$dif1 <- c(0, diff(input_start$Value))
  input_start$dif2 <- c(0, diff(input_start$dif1))
  lineseg <- input_start[which.max(input_start$dif2):which.max(input_start$Value), ]

  abc <- lm(lineseg$Value ~ lineseg$Time)

  # qplot(input_start$Time, input_start$Value) + geom_point() + geom_abline(slope = abc$coefficients[2], intercept=abc$coefficients[1])

  starttime <- -1 * abc$coefficients[1] / abc$coefficients[2]

  input <- subset(input, input$Time > starttime)
  input <- rbind(c(starttime, 0), input)

  input$Time <- input$Time - input$Time[1]

  return(input)
}


plot_fitpoints <- function(ymeasured, yfitted, xmeasured, xfitted = NULL) {
  if (is.null(xfitted)) {
    xfitted <- xmeasured
  }

  # plot the values
  plot(xmeasured, ymeasured)

  # draw the curve
  lines(xfitted, yfitted, col = "blue")
}

plot_fitfunc <- function(ymeasured, xmeasured, fitfunction, parameters) {

  # plot the values
  plot(xmeasured, ymeasured)

  # generate a range of values for xmeasured in small increments to create a smooth line
  xRange <- seq(min(xmeasured), max(xmeasured), length.out = 1024)

  # generate the predicted y values
  yValues <- predict(fitfunction, par = parameters, x = xRange)

  # draw the curve
  lines(xRange, yValues, col = "blue")
}


#' Plot Residuals of a Model
#'
#' Function to visualise the residuals of a model fit.  This function only works for model outputs for
#' which a fit object is available called \code{outputobject$fit}.
#'
#' @param outputobject The output of a kinetic model, including a fit object \code{outputobject$fit}
#'
#' @return A ggplot2 object of the residual plot
#'
#' @examples
#' # Note: Reference region models, and irreversible binding models, should not
#' # be used for PBR28 - this is just to demonstrate function
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' reftac <- pbr28$tacs[[2]]$CBL
#' roitac <- pbr28$tacs[[2]]$STR
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' srtmout <- srtm(t_tac, reftac, roitac)
#'
#' plot_residuals(srtmout)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

plot_residuals <- function(outputobject) {
  Time <- outputobject$tacs$Time
  Residuals <- residuals(outputobject$fit)
  Weights <- weights(outputobject$fit)

  if (length(Time) > length(Residuals)) {
    Time <- tail(Time, length(Residuals))
  }

  plotdf <- data.frame(Time = Time, Residuals = Residuals, Weights = Weights)

  ylimits <- c(-max(abs(plotdf$Residuals[plotdf$Weights > 0])), max(abs(plotdf$Residuals[plotdf$Weights > 0])))

  outplot <- ggplot(plotdf, aes(x = Time, y = Residuals)) +
    geom_point(aes(size = Weights)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") + ylim(ylimits) +
    scale_size(range = c(1, 3)) + geom_line()

  return(outplot)
}


#' Tidy Up for Reference Region Methods
#'
#' Function to tidy up the input argument vectors for reference region models.
#'
#' @param t_tac Numeric vector of times for each frame in minutes.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param roitac Numeric vector of radioactivity concentrations in the target tissue for each frame.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting.
#' @param frameStartEnd Optional. This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#'
#' @return A dataframe containing the tidied up t_tac, reftac, roitac and weights ready for modelling.
#'
#' @details This function i) adds uniform weights if weights are not specified, ii) checks that the
#' lengths of t_tac, reftac, roitac and weights are of the same length, iii) shortens the vectors
#' if a frameStartEnd is specified, iv) adds a zero frame if there is not one, and v) checks that
#' times are in minutes and not in seconds.
#'
#' @examples
#'
#' #' # Note: Reference region models, and irreversible binding models, should not
#' # be used for PBR28 - this is just to demonstrate function
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' reftac <- pbr28$tacs[[2]]$CBL
#' roitac <- pbr28$tacs[[2]]$STR
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd = c(1, 10))
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

tidyinput_ref <- function(t_tac, reftac, roitac, weights, frameStartEnd) {
  if (is.null(weights) == T) {
    weights <- rep(1, length(reftac))
  }

  lengths <- c(length(t_tac), length(reftac), length(roitac), length(weights))
  if (!all(lengths == lengths[1])) {
    stop("The lengths of t_tac, reftac, roitac and/or weights are not equal")
  }

  if (!is.null(frameStartEnd)) {
    t_tac <- t_tac[ frameStartEnd[1]:frameStartEnd[2] ]
    reftac <- reftac[ frameStartEnd[1]:frameStartEnd[2] ]
    roitac <- roitac[ frameStartEnd[1]:frameStartEnd[2] ]
    weights <- weights[ frameStartEnd[1]:frameStartEnd[2] ]
  }

  if (min(t_tac) > 0) {
    t_tac <- c(0, t_tac)
    roitac <- c(0, roitac)
    reftac <- c(0, reftac)
    weights <- c(0, weights)
  }

  if (max(t_tac) > 300) {
    warning("\n      ******************************************************************************\n          It looks like you have included seconds instead of minutes for time:\n          this can cause wrong/weird results, and should be avoided. Ignore this\n          warning if you just have really long measurements (over 300 minutes).\n      ******************************************************************************")
  }

  out <- data.frame(
    t_tac = t_tac, reftac = reftac, roitac = roitac,
    weights = weights
  )

  return(out)
}


#' Tidy Up for Models with Arterial Input
#'
#' Function to tidy up the input argument vectors for arterial input models.
#'
#' @param t_tac Numeric vector of times for each frame in minutes.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting.
#' @param frameStartEnd Optional. This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#'
#' @return A list containing the tidied up t_tac, tac and weights ready for modelling.
#'
#' @details This function i) adds uniform weights if weights are not specified, ii) checks that the
#' lengths of t_tac, tac and weights are of the same length, iii) shortens the vectors
#' if a frameStartEnd is specified, iv) Checks whether there are negative values in t_tac due to
#' adjusting for delay iv) adds a zero frame if there is not one, and v) checks that
#' times are in minutes and not in seconds.
#'
#' @examples
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' tac <- pbr28$tacs[[2]]$STR
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' tidyinput_art(t_tac, tac, weights, frameStartEnd = c(1, 10))
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

tidyinput_art <- function(t_tac, tac, weights, frameStartEnd) {
  if (is.null(weights)) {
    weights <- rep(1, length(tac))
  }

  lengths <- c(length(t_tac), length(tac), length(weights))
  if (!all(lengths == lengths[1])) {
    stop("The lengths of t_tac, tac and/or weights are not equal")
  }

  if (!is.null(frameStartEnd)) {
    t_tac <- t_tac[ frameStartEnd[1]:frameStartEnd[2] ]
    tac <- tac[ frameStartEnd[1]:frameStartEnd[2] ]
    weights <- weights[ frameStartEnd[1]:frameStartEnd[2] ]
  }

  if (min(t_tac) < 0) {
    stop("There are negative times in the TAC")
  }

  if (min(t_tac) > 0) {
    t_tac <- c(0, t_tac)
    tac <- c(0, tac)
    weights <- c(0, weights)
  }

  if (max(t_tac) > 300) {
    warning("\n      ******************************************************************************\n          It looks like you have included seconds instead of minutes for time:\n          this can cause wrong/weird results, and should be avoided. Ignore this\n          warning if you just have really long measurements (over 300 minutes).\n      ******************************************************************************")
  }

  out <- data.frame(
    t_tac = t_tac, tac = tac,
    weights = weights
  )

  return(out)
}

#' Maximum Percentage Residual
#'
#' Function to determine the maximum percentage residual from the fitted value, i.e. the maximum
#' value of a residual as a percentage of its fitted value.  This is useful for obtaining estimates
#' of t*.
#'
#' @param outputobject The output of a kinetic model, including a fit object \code{outputobject$fit}
#'
#' @return The maximum percentage of the fitted value of a residual in the fit.
#'
#' @examples
#'
#' # Note: Reference region models, and irreversible binding models, should not
#' # be used for PBR28 - this is just to demonstrate function
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' reftac <- pbr28$tacs[[2]]$CBL
#' roitac <- pbr28$tacs[[2]]$STR
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' refloganout <- refLogan(t_tac, reftac, roitac, 0.1, tstarIncludedFrames = 9)
#'
#' maxpercres(refloganout)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

maxpercres <- function(outputobject) {
  Residuals <- abs(residuals(outputobject$fit))
  Fitted <- abs(fitted(outputobject$fit))
  PercRes <- 100 * Residuals / Fitted
  MaxPercRes <- max(PercRes)
  return(MaxPercRes)
}

#' Plot Kinetic Model Fit: Generic Function
#'
#' Function to plot the output of a kinetic model. This function calls the specific plotting functions
#' for each model based on the output of the model.
#'
#' @param modelout The output object of the model fitting procedure.
#' @param ... Additional optional arguments.
#'
#' @return A ggplot2 object of the plot.
#'
#' @details This function uses the \code{out$model} name to call the correct function to plot the model fit.
#'
#' @examples
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' tac <- pbr28$tacs[[2]]$FC
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1
#' )
#'
#' fit <- ma1(t_tac, tac, input, 10, weights)
#'
#' plot_kinfit(fit, roiname = "FC")
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

plot_kinfit <- function(modelout, ...) {
  modelname <- modelout$model
  input_list <- as.list(substitute(list(...)))
  do.call(what = paste0("plot_", modelname, "fit"), args = list(modelout, ...))
  # eval(expr = paste0('plot_', modelname, 'fit'))
}

#' Check and fix multstart parameters
#'
#' Checks and fixes for inputted multstart parameters. Checks missing multstart boundaries,
#' checks if some parameters aren't to be tried with different starting parameters, and checks
#' the lengths of the multstart bounds
#'
#' @param start Original starting values for without multstart
#' @param lower Lower fitting bounds
#' @param upper Upper fitting bounds
#' @param multstart_iter Iterations
#' @param multstart_lower Lower starting bounds
#' @param multstart_upper Upper starting bounds
#'
#' @return A list containing the parameters
#' @export
#'
#' @examples
#' start <- list(a = 1, b = 2, c = 3)
#' upper <- lapply(start, function(x) x + 1)
#' lower <- lapply(start, function(x) x - 1)
#' multstart_iter <- 5
#'
#' multstart_lower <- NULL
#' multstart_upper <- NULL
#'
#' fix_multstartpars(start, lower, upper, multstart_iter, multstart_lower, multstart_upper)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
fix_multstartpars <- function(start, lower, upper, multstart_iter, multstart_lower, multstart_upper) {
  if (length(multstart_iter) == length(start) || length(multstart_iter) == 1) {
    multstart_l <- as.list(lower)
    multstart_u <- as.list(upper)

    parnames <- names(start)

    ### Adding multstart boundaries
    if (!is.null(multstart_lower)) {
      if (class(multstart_lower) != "list" ||
        sum(!(names(multstart_lower) %in% parnames)) > 0) { # Are there any names not in start?
        stop("multstart_lower should be a named list whose names match the parameters to be fitted")
      }

      multstart_lower <- modifyList(multstart_l, multstart_lower)
    } else {
      multstart_lower <- multstart_l
    }

    if (!is.null(multstart_upper)) {
      if (class(multstart_upper) != "list" ||
        sum(!(names(multstart_upper) %in% parnames)) > 0) { # Are there any names not in start?
        stop("multstart_upper should be a named list whose names match the parameters to be fitted")
      }

      multstart_upper <- modifyList(multstart_u, multstart_upper)
    } else {
      multstart_upper <- multstart_u
    }

    ### No multstart for some variables ###
    if (any(multstart_iter == 1)) {
      non_iterable <- which(multstart_iter == 1)
      multstart_lower[non_iterable] <- start[non_iterable]
      multstart_upper[non_iterable] <- start[non_iterable]
    }

    # Check order
    multstart_lower <- multstart_lower[names(start)]
    multstart_upper <- multstart_upper[names(start)]


    # Turn into named vectors again

    multstart_lower <- as.numeric(multstart_lower)
    multstart_upper <- as.numeric(multstart_upper)

    names(multstart_lower) <- names(start)
    names(multstart_upper) <- names(start)
  } else {
    stop("multstart_iter should be of length 1 or of the same length as the number of parameters")
  }

  out <- list(
    multstart_lower = multstart_lower,
    multstart_upper = multstart_upper
  )

  return(out)
}

#' Estimate the standard error of parameters using the delta method
#'
#' This function returns the SE of parameters, as a proportion of their
#' parameter estimates, calculated using the delta method as implemented in the
#' car package. If the parameter does not exist, this function will return NA.
#'
#' @param fit A fit object for a particular model - in this case the actual lm
#'   or nls object
#' @param expression The expression for which the SE is being calculated
#'
#' @return The SE as a proportion of the parameter estimate.
#' @export
#'
#' @examples
#' data(mtcars)
#' a <- lm(mpg ~ cyl + disp, data=mtcars)
#'
#' get_se(a, "cyl")
#' get_se(a, "cyl/disp")
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
get_se <- function(fit, expression) {

  secov <- function(fit, expression) {
    carout <- car::deltaMethod(fit, expression)
    abs( carout$SE / carout$Estimate )
  }

  possibly_se <- purrr::possibly(secov, otherwise = NA)

  possibly_se(fit, expression)

}

#' Convert between different units of radioactivity
#'
#' @param values Radioactivity values which should be converted. This can be a single value or a vector of values.
#' @param from_units The units from which the radioactivity should be converted.
#' @param to_units The units to which the radioactivity should be converted.
#'
#' @return The converted radioactivity values in the new units
#' @export
#'
#' @examples
#' unit_convert(1, "nCi", "Bq")
#' unit_convert(1:5, "nCi", "Bq")
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
unit_convert <- function(values,
                         from_units = c("GBq", "MBq", "kBq","Bq",
                                       "mBq", "uBq", "nBq","pBq",
                                       "GCi", "MCi", "kCi","Ci",
                                       "mCi", "uCi", "nCi", "pCi"),
                         to_units = c("GBq", "MBq", "kBq","Bq",
                                     "mBq", "uBq", "nBq","pBq",
                                     "GCi", "MCi", "kCi","Ci",
                                     "mCi", "uCi", "nCi", "pCi")) {

  if(length(from_units) > 1 | length(to_units) > 1) {
    stop("Units must be selected to convert between")
  }

  # Match the units

  from <- match.arg(from_units, c("GBq", "MBq", "kBq", "Bq",
                                    "mBq", "uBq", "nBq", "pBq",
                                    "GCi", "MCi", "kCi", "Ci",
                                    "mCi", "uCi", "nCi", "pCi"))

  to <- match.arg(to_units, c("GBq", "MBq", "kBq","Bq",
                                  "mBq", "uBq", "nBq","pBq",
                                  "GCi", "MCi", "kCi","Ci",
                                  "mCi", "uCi", "nCi", "pCi"))

  # Extract the meaning from the units

  from_bqci <- ifelse( grepl("Bq", from), "Bq", "Ci")
  from_mult <- strsplit(x = from, split = from_bqci)[[1]]

  to_bqci <- ifelse( grepl("Bq", to), "Bq", "Ci")
  to_mult <- strsplit(x = to, split = to_bqci)[[1]]

  from_multno <- dplyr::case_when(
    from_mult == "G" ~ 1e9,
    from_mult == "M" ~ 1e6,
    from_mult == "k" ~ 1e3,
    from_mult == ""  ~ 1,
    from_mult == "m" ~ 1e-3,
    from_mult == "u" ~ 1e-6,
    from_mult == "n" ~ 1e-9,
    from_mult == "n" ~ 1e-12,
    TRUE ~ 1
  )

  to_multno <- dplyr::case_when(
    to_mult == "G" ~ 1e9,
    to_mult == "M" ~ 1e6,
    to_mult == "k" ~ 1e3,
    to_mult == ""  ~ 1,
    to_mult == "m" ~ 1e-3,
    to_mult == "u" ~ 1e-6,
    to_mult == "n" ~ 1e-9,
    to_mult == "p" ~ 1e-12,
  )

  # Removing the from multiplier

  values <- values * from_multno

  # Converting between Bq and Ci

  values <- dplyr::case_when(
    from_bqci == to_bqci ~ values,
    from_bqci == "Ci" & to_bqci == "Bq" ~ values * 37000000000,
    from_bqci == "Bq" & to_bqci == "Ci" ~ values / 37000000000,
  )

  # Adding the to multiplier

  values <- values / to_multno

  # Return

  return(values)

}

#' Reverse Decay Correction of TAC values
#'
#' This function uses the equation described here to perform the calculations: http://www.turkupetcentre.net/petanalysis/decay.html
#'
#' @param t_start The starting times of the frames in minutes.
#' @param t_end The end times of the frames in minutes.
#' @param tac The uncorrected radioactivity values.
#' @param radioisotope The radioisotope.
#' @param radionuclide_halflife The half-life of the radionuclide in minutes. If specified, this overrides the radioisotope parameter. Defaults to NULL.
#'
#' @return The tac before decay correction
#' @export
#'
#' @examples
#' data(pbr28)
#' s1 <- pbr28$tacs[[1]]
#'
#' # Assuming the data were not decay-corrected (they are)
#' decay_uncorrect(
#'   s1$StartTime/60,
#'  (s1$StartTime + s1$Duration)/60,
#'  tac = s1$FC)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
decay_uncorrect <- function(t_start, t_end, tac, radioisotope = c("C11", "O15", "F18"), radionuclide_halflife = NULL) {

  if (!is.null(radionuclide_halflife)) {
    hl <- radionuclide_halflife
  } else {
    radioisotope <- match.arg(radioisotope, c("C11", "O15", "F18"))
    hl <- dplyr::case_when(
      radioisotope == "C11" ~ 20.4,
      radioisotope == "O15" ~ 2.05,
      radioisotope == "F18" ~ 109.8
    )
  }

  lambda <- log(2) / hl

  term1 <- exp(lambda*t_start) * lambda * (t_end - t_start)
  term2 <- 1 - exp( -1 * lambda * (t_end - t_start) )

  tac_uncor <- tac * ( term2 / term1 )

  tac_uncor[is.nan(tac_uncor)] <- 0

  return(tac_uncor)

}

#' Perform Decay Correction of TAC values
#'
#' This function uses the equation described here to perform the calculations: http://www.turkupetcentre.net/petanalysis/decay.html
#'
#' @param t_start The starting times of the frames in minutes.
#' @param t_end The end times of the frames in minutes.
#' @param tac_uncor The uncorrected radioactivity values.
#' @param radioisotope The radioisotope.
#' @param radionuclide_halflife The half-life of the radionuclide in minutes. If specified, this overrides the radioisotope parameter. Defaults to NULL.
#'
#' @return The tac after decay correction
#' @export
#'
#' @examples
#' data(pbr28)
#' s1 <- pbr28$tacs[[1]]
#'
#' # Assuming the data were not decay-corrected (they are)
#' decay_correct(
#'   s1$StartTime/60,
#'  (s1$StartTime + s1$Duration)/60,
#'  tac_uncor = s1$FC)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
decay_correct <- function(t_start, t_end, tac_uncor,
                          radioisotope = c("C11", "O15", "F18"),
                          radionuclide_halflife = NULL) {

  if (!is.null(radionuclide_halflife)) {
    hl <- radionuclide_halflife
  } else {
    radioisotope <- match.arg(radioisotope, c("C11", "O15", "F18"))
    hl <- dplyr::case_when(
      radioisotope == "C11" ~ 20.4,
      radioisotope == "O15" ~ 2.05,
      radioisotope == "F18" ~ 109.8
    )
  }

  lambda <- log(2) / hl

  term1 <- exp(lambda*t_start) * lambda * (t_end - t_start)
  term2 <- 1 - exp( -1 * lambda * (t_end - t_start) )

  tac <- tac_uncor * ( term1 / term2 )

  tac[is.nan(tac)] <- 0

  return(tac)

}

#' Cumulative integral over frames
#'
#' This function uses the frame radioactivity values and durations to create a
#' cumulative integral. It sums the multiple of radioactivity and duration
#' for all frames preceding the current frame, and halves that of the current
#' frame. This is an alternative to using the trapezoidal method.
#'
#' @param dur Frame durations in minutes.
#' @param tac Radioactivity concentrations for each frame.
#'
#' @return The cumulative integral of radioactivity.
#' @export
#'
#' @examples
#' data(pbr28)
#'
#' tac <- pbr28$tacs[[2]]$WB
#' dur <- pbr28$tacs[[2]]$Duration / 60
#'
#' frame_integ <- frame_cumsum(dur, tac)
#'
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
frame_cumsum <- function(dur, tac) {

  if( max(dur) > 60 ) {
    warning("The frame durations are very long. Make sure they are in minutes and not seconds")
  }

  durtac <- dur*tac
  cumdurtac <- cumsum(durtac)

  before_cumdurtac <- c(0, cumdurtac[-(length(cumdurtac))])
  halfway_durtac <- 0.5 * durtac

  out <- before_cumdurtac + halfway_durtac

  return(out)

}

camel <- function(x){ #function for camel case
  capit <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2, nchar(x)))
  first <- strsplit(x, "\\_")[[1]][1]
  rest <- sapply(strsplit(x, "\\_")[[1]][-1],
                 function(x) paste(capit(x), sep=""))

  paste(first, rest, collapse="")
}

get_units_radioactivity <- function(x) {

  rad <- gsub("[\\.\\/].*", "", x)
  vol <- gsub("\\w*[\\.\\/]", "", x)
  vol <- gsub("\\-1$", "", vol)

  list(rad=rad, vol=vol)

}


#' Predict Blood Model Values: Generic Function
#'
#' Function to predict values from blood models
#'
#' @param modelout The output object of the model fitting procedure.
#' @param ... Additional optional arguments.
#'
#' @return The predicted values
#'
#' @details This function uses the \code{class(model)} name to call the correct function to plot the model fit.
#'
#' @examples
#'
#' blooddata <- pbr28$blooddata[[1]]
#' blooddata <- bd_blood_dispcor(blooddata)
#' blood <- bd_extract(blooddata, output = "Blood")
#' blood_fit <- blmod_splines(blood$time,
#'                            blood$activity,
#'                            Method = blood$Method)
#'
#' predict_blmod(blood_fit)
#'
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
predict_blmod <- function(modelout, ...) {
  modelname <- class(modelout)[1]
  input_list <- as.list(substitute(list(...)))

  out <- do.call(what = paste0("predict_", modelname), args = list(modelout, ...))

  as.numeric(out)

}



#' Extract coefficients
#'
#' Function to extract the coefficients from kinfit objects
#'
#' @param x The output object of the model fitting procedure.
#' @param ... Additional optional arguments.
#'
#' @examples
#' \dontrun{
#' loganout <- Loganplot(t_tac, tac, input, 10, weights)
#' coef(loganout)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
coef_kinfit <- function(x, ...) {
  out <- x$par
  return(out)
}

interpends <- function(x, y, xi, method = "linear", yzero = NULL) {
  if (is.null(yzero)) { # start point
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

