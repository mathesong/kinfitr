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
#' stepsize = 1
#' 
#' kinfit_convolve(a,b,stepsize)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

kinfit_convolve <- function(a, b, stepsize) {
  
  if(length(a) != length(b)) {
    stop('The vectors are not of the same length')
  }
  
  out <- pracma::conv(a, b)[1:length(a)]*stepsize
  return(out)
  
}

kinfit_uninterpbydur <- function(interptime, interppred, times, durtimes) {
  
  starttimes = times - durtimes/2
  interpstart = sapply(starttimes, function(x) which.min(abs(interptime-x)))[-1]
  interpend = c((interpstart-1)[-1], length(interptime))
  
  startstop <- data.frame(start = interpstart, stop = interpend)
  
  predvals = c(0, mapply(function(start, stop) mean(interppred[start:stop]), startstop$start, startstop$stop))
  
  return(predvals)
  
}

zerotrim_input <- function(t_inp, inp) {
  
  # Create a zero point by drawing a regression between the peak and the peak second differential

  input <- data.frame(Time = t_inp, Value = inp)
  
  input$Time <- input$Time
  
  input_start <- subset(input, input$Time < which.max(input$Value))
  #qplot(input_start$Time, input_start$Value)
  input_start$dif1 <- c(0, diff(input_start$Value))
  input_start$dif2 <- c(0, diff(input_start$dif1))
  lineseg = input_start[which.max(input_start$dif2):which.max(input_start$Value),]
  
  abc <- lm(lineseg$Value ~ lineseg$Time)
  
  #qplot(input_start$Time, input_start$Value) + geom_point() + geom_abline(slope = abc$coefficients[2], intercept=abc$coefficients[1])
  
  starttime <- -1*abc$coefficients[1]/abc$coefficients[2]
  
  input <- subset(input, input$Time > starttime)
  input <- rbind(c(starttime, 0), input)
  
  input$Time <- input$Time - input$Time[1]
  
  return(input)
  
}


plot_fitpoints <- function(ymeasured, yfitted, xmeasured, xfitted = NULL) {
  
  if(is.null(xfitted)) { xfitted = xmeasured }
  
  # plot the values
  plot(xmeasured, ymeasured)
  
  #draw the curve
  lines(xfitted, yfitted, col="blue")
  
}

plot_fitfunc <- function(ymeasured, xmeasured, fitfunction, parameters) {
  
  # plot the values
  plot(xmeasured, ymeasured)
  
  # generate a range of values for xmeasured in small increments to create a smooth line
  xRange <- seq(min(xmeasured), max(xmeasured), length.out = 1024)
  
  # generate the predicted y values
  yValues <- predict(fitfunction, par = parameters, x = xRange)
  
  #draw the curve
  lines(xRange, yValues, col="blue")
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
#' srtmout <- srtm(t_tac, reftac, roitac)
#' plot_residuals(srtmout)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

plot_residuals <- function(outputobject) {
  
  Time = outputobject$tacs$Time
  Residuals = residuals(outputobject$fit)
  Weights = weights(outputobject$fit)
  
  if(length(Time) > length(Residuals)) {
    Time = tail(Time, length(Residuals))
  }
  
  plotdf <- data.frame(Time = Time, Residuals = Residuals, Weights = Weights)
  
  ylimits = c( -max(abs(plotdf$Residuals[plotdf$Weights > 0])) , max(abs(plotdf$Residuals[plotdf$Weights > 0])) )
  
  outplot <- ggplot(plotdf, aes(x = Time, y = Residuals)) + 
    geom_point(aes(size = Weights)) + 
    geom_hline(aes(yintercept = 0), linetype="dashed") + ylim(ylimits) +
    scale_size(range = c(1,3)) + geom_line()
    
  return(outplot)
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
#' refloganout <- refLogan(times, reftac, roitac, k2prime, tstarIncludedFrames = 9, weights=weights)
#' maxpercres(refloganout)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

maxpercres <- function(outputobject) {
  
  Residuals = abs(residuals(outputobject$fit))
  Fitted = abs(fitted(outputobject$fit))
  PercRes = 100 * Residuals / Fitted
  MaxPercRes = max(PercRes)
  return(MaxPercRes)
  
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
#' tidyinput_ref(t_tac, reftac, roitac, weights)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

tidyinput_ref <- function(t_tac, reftac, roitac, weights, frameStartEnd) {
  
  if(missing(weights) == T) {
    weights = rep(1, length(reftac))
  }
  
  lengths <- c(length(t_tac), length(reftac),  length(roitac), length(weights))
  if(!all(lengths == lengths[1])) {
    stop('The lengths of t_tac, reftac, roitac and/or weights are not equal')
  }
  
  if(!missing(frameStartEnd)) {
    t_tac <- t_tac[ frameStartEnd[1] : frameStartEnd[2] ]
    reftac <- reftac[ frameStartEnd[1] : frameStartEnd[2] ]
    roitac <- roitac[ frameStartEnd[1] : frameStartEnd[2] ]
    weights <- weights[ frameStartEnd[1] : frameStartEnd[2] ]
  }
  
  if(min(t_tac) > 0) {
    t_tac = c(0, t_tac)
    roitac = c(0, roitac)
    reftac = c(0, reftac)
    weights = c(0, weights)
  }
  
  if(max(t_tac) > 300) {
    warning('
      ******************************************************************************
          It looks like you have included seconds instead of minutes for time:
          this can cause wrong/weird results, and should be avoided. Ignore this
          warning if you just have really long measurements (over 300 minutes).
      ******************************************************************************')
  }
  
  out <- data.frame(t_tac = t_tac, reftac = reftac, roitac = roitac, 
              weights = weights)
  
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
#' tidyinput_art(t_tac, tac, weights)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

tidyinput_art <- function(t_tac, tac, weights, frameStartEnd) {
  
  if(missing(weights)) {
    weights = rep(1, length(tac))
  }
  
  lengths <- c(length(t_tac), length(tac), length(weights))
  if(!all(lengths == lengths[1])) {
    stop('The lengths of t_tac, tac and/or weights are not equal')
  }
  
  if(!missing(frameStartEnd)) {
    t_tac <- t_tac[ frameStartEnd[1] : frameStartEnd[2] ]
    tac <- tac[ frameStartEnd[1] : frameStartEnd[2] ]
    weights <- weights[ frameStartEnd[1] : frameStartEnd[2] ]
  }
  
  if(min(t_tac) < 0) {
    stop('There are negative times in the TAC')
  }
  
  if(min(t_tac) > 0) {
    t_tac = c(0, t_tac)
    tac = c(0, tac)
    weights = c(0, weights)
  }
  
  if(max(t_tac) > 300) {
    warning('
      ******************************************************************************
          It looks like you have included seconds instead of minutes for time:
          this can cause wrong/weird results, and should be avoided. Ignore this
          warning if you just have really long measurements (over 300 minutes).
      ******************************************************************************')
  }
  
  out <- data.frame(t_tac = t_tac, tac = tac, 
              weights = weights)
  
  return(out)
  
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
#' loganout <- Loganplot(t_tac, tac, input, 10, weights)
#' plot_kinfit(loganout)
#' 
#' srtmout <- srtm(t_tac, reftac, roitac)
#' plot_kinfit(srtmout, roiname = 'Putamen', refname = 'Cerebellum')
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

plot_kinfit <- function(modelout, ...) {
  modelname <- modelout$model
  input_list <- as.list(substitute(list(...)))
  do.call(what = paste0('plot_', modelname, 'fit'), args = list(modelout, ...))
  #eval(expr = paste0('plot_', modelname, 'fit'))
}