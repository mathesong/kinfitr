#' Non-Invasive Logan Plot
#'
#' Function to fit the non-invasive Logan plot model of Logan et al. (1996) to
#' data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param roitac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param k2prime Value of k2prime to be used for the fitting, i.e. the average
#'   tissue-to-plasma clearance rate. This can be obtained from another model,
#'   or set at a specified value. If using SRTM to estimate this value, it is
#'   equal to k2 / R1.
#' @param tstar The t* specification for regression. If tstar_type="frames",
#'   this is the number of frames from the end to include (e.g., 10 means last
#'   10 frames). If tstar_type="time", this is the time point (in minutes) after
#'   which all frames with midpoints later than this time are included. This
#'   value can be estimated using \code{refLogan_tstar}.
#' @param tstar_type Either "frames" (default) or "time", specifying how to
#'   interpret tstar.
#' @param tstarIncludedFrames Deprecated. Use 'tstar' with 'tstar_type="frames"'
#'   instead.
#' @param weights Optional. Numeric vector of the conventional frame-wise weights
#'   assigned to each frame. If not specified, uniform weights will be used.
#'   Specified weights are internally transformed to account for the dependent
#'   variable transformation in the Logan plot. If \code{dur} is not provided,
#'   weights cannot be transformed and uniform weights will be used. We include
#'   zero at time zero: if not included, it is added.
#' @param dur Optional. Numeric vector of the time durations of the frames. If
#'   not included, the integrals will be calculated using trapezoidal
#'   integration.
#' @param frameStartEnd Optional: This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This can be used to assess
#'   time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and
#'   end time point instead of defining the frame numbers using frameStartEnd.
#'   This function will restrict the model to all time frames whose t_tac is
#'   between the values, i.e. c(0,5) will select all frames with midtimes during
#'   the first 5 minutes.

#'
#' @return A list with a data frame of the fitted parameters \code{out$par},
#'   their percentage standard errors (scaled so that 1 represents 100\%)
#'   \code{out$par.se}, the
#'   model fit object \code{out$fit}, a dataframe containing the TACs of the
#'   data \code{out$tacs}, a dataframe containing the TACs of the fitted values
#'   \code{out$fitvals}, a vector of the weights \code{out$weights}, the
#'   specified k2prime value \code{out$k2prime}, and the specified
#'   tstarIncludedFrames value \code{out$tstarIncludedFrames}
#'
#' @examples
#'
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' reftac <- simref$tacs[[2]]$Reference
#' roitac <- simref$tacs[[2]]$ROI1
#' weights <- simref$tacs[[2]]$Weights
#'
#' fit <- refLogan(t_tac, reftac, roitac, k2prime = 0.1, tstar = 15, weights = weights)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Logan J, Fowler JS, Volkow ND, Wang GJ, Ding YS, Alexoff DL.
#'   Distribution volume ratios without blood sampling from graphical analysis
#'   of PET data. Journal of Cerebral Blood Flow & Metabolism. 1996 Sep
#'   1;16(5):834-40.
#'
#' @export

refLogan <- function(t_tac, reftac, roitac, k2prime, tstar, weights = NULL,
                     dur = NULL, frameStartEnd = NULL, timeStartEnd = NULL,
                     tstar_type = "frames", tstarIncludedFrames = NULL) {

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  # Handle deprecated parameter
  if (!is.null(tstarIncludedFrames)) {
    warning("tstarIncludedFrames is deprecated and will be removed in a future version. Use 'tstar' with 'tstar_type=\"frames\"' instead", call. = FALSE)
    if (!missing(tstar)) {
      stop("Cannot specify both 'tstar' and 'tstarIncludedFrames'")
    }
    tstar <- tstarIncludedFrames
    tstar_type <- "frames"
  }

  # Validate tstar_type
  tstar_type <- match.arg(tstar_type, c("frames", "time"))

  # Handle missing tstar
  if (missing(tstar)) {
    tstar <- length(t_tac)
    tstar_type <- "frames"
    warning("No value specified for tstar: defaulting to including all frames. This may produce biased outcomes.", call. = FALSE)
  }

  # Tidying

  tidyinput <- tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd)

  if (!is.null(dur)) {
    tidyinput_dur <- tidyinput_ref(dur, reftac, roitac, weights, frameStartEnd)
    dur <- tidyinput_dur$t_tac
  }

  t_tac <- tidyinput$t_tac
  reftac <- tidyinput$reftac
  roitac <- tidyinput$roitac
  weights <- tidyinput$weights

  # Convert tstar based on type
  if (tstar_type == "time") {
    frames_after_tstar <- which(t_tac >= tstar)
    tstarIncludedFrames <- length(frames_after_tstar)
  } else {
    tstarIncludedFrames <- tstar
  }

  # Transform weights for graphical analysis (if provided)
  # Check if real weights were provided (more than just 0s and 1s from tidyinput)
  unique_weights <- unique(weights[is.finite(weights)])
  real_weights_provided <- length(setdiff(unique_weights, c(0, 1))) > 0
  if (!is.null(weights) && real_weights_provided) {
    if (!is.null(dur)) {
      weights <- weights_Logan_transform(t_tac, dur, roitac, weights)
      # Center weights so mean of equilibrium weights equals 1
      equil_weights <- tail(weights, tstarIncludedFrames)
      equil_mean <- mean(equil_weights[is.finite(equil_weights)])
      if (is.finite(equil_mean) && equil_mean > 0) {
        weights <- weights / equil_mean
      }
      # Set pre-equilibrium weights to 1
      pre_equil_idx <- seq_len(length(weights) - tstarIncludedFrames)
      weights[pre_equil_idx] <- 1
    } else {
      message("Weights provided but dur is NULL. Frame durations are required for ",
              "weight transformation in the Logan plot. Using uniform weights.")
      weights <- rep(1, length(t_tac))
    }
  }

  # Parameters

  if (!is.null(dur)) {

    logan_roi <- frame_cumsum(dur, roitac) / roitac
    logan_ref <- (frame_cumsum(dur, reftac) + reftac / k2prime) / roitac

  } else {

    logan_roi <- pracma::cumtrapz(t_tac, roitac) / roitac
    logan_ref <- (pracma::cumtrapz(t_tac, reftac) + reftac / k2prime) / roitac

  }

  logan_equil_roi <- tail(logan_roi, tstarIncludedFrames)
  logan_equil_ref <- tail(logan_ref, tstarIncludedFrames)
  weights_equil <- tail(weights, tstarIncludedFrames)


  # Solution

  logan_model <- lm(logan_equil_roi ~ logan_equil_ref, weights = weights_equil)

  # Output

  par <- as.data.frame(list(bp = as.numeric(logan_model$coefficients[2]) - 1))

  par.se <- par
  names(par.se) <- paste0(names(par.se), ".se")
  par.se$bp.se <- get_se(logan_model, "logan_equil_ref - 1")

  fit <- logan_model

  tacs <- data.frame(Time = t_tac, Reference = reftac, Target = roitac)

  if (!is.null(dur)) { tacs$Duration <- dur }

  fitvals <- data.frame(Logan_ROI = logan_roi, Logan_Ref = logan_ref)

  out <- list(
    par = par, par.se = par.se, fit = fit, tacs = tacs,
    fitvals = fitvals, weights = weights, k2prime = k2prime,
    tstarIncludedFrames = tstarIncludedFrames, model = "refLogan"
  )

  class(out) <- c("refLogan", "kinfit")

  return(out)
}

#' Plot: Non-Invasive Logan Plot
#'
#' Function to visualise the fit of the refLogan model to data.
#'
#' @param refloganout The output object of the refLogan fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' reftac <- simref$tacs[[2]]$Reference
#' roitac <- simref$tacs[[2]]$ROI1
#' weights <- simref$tacs[[2]]$Weights
#'
#' fit <- refLogan(t_tac, reftac, roitac, k2prime = 0.1, tstar = 10, weights = weights)
#'
#' plot_refLoganfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_refLoganfit <- function(refloganout, roiname = NULL) {
  plotdf <- data.frame(
    Weights = refloganout$weights,
    Logan_ref = refloganout$fitvals$Logan_Ref,
    Logan_roi = refloganout$fitvals$Logan_ROI,
    Equilibrium = as.character("Before")
  )

  plotdf$Equilibrium <- as.character(plotdf$Equilibrium)
  plotdf$Equilibrium [ (nrow(plotdf) - (refloganout$tstarIncludedFrames - 1)):nrow(plotdf)  ] <- "After"

  # Set pre-tstar weights to 1 for display (so points are visible but don't affect scale)
  plotdf$Weights[plotdf$Equilibrium == "Before"] <- 1

  plotdf$Equilibrium <- forcats::fct_inorder(factor(plotdf$Equilibrium))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Equilibrium)
  colScale <- scale_colour_manual(name = paste0(roiname, "\nEquilibrium"), values = myColors)

  xlabel <- expression(paste("(", "",
                                    paste(integral(, paste("0"), paste("", "t")),
                                          "C", phantom()[{ paste("R") }],"(",tau,")d",tau, " + ",
                                          frac(paste("C", phantom()[{ paste("R") }],"(t)"),
                                               paste("k",phantom()[{ paste("2") }],"\'",))),")",
                                      " / ",
                                 "C", phantom()[{ paste("T") }],"(t)"))

  ylabel <- expression(paste("", "", integral(, paste("0"), paste("", "t")),
                                 "C", phantom()[{ paste("T") }],"(",tau,")d",tau, " / ",
                                 "C", phantom()[{ paste("T") }],"(t)"))


  # Fix limits

  after_equil <- plotdf %>%
    dplyr::filter(Equilibrium == "After")

  finite_x <- after_equil$Logan_ref[is.finite(after_equil$Logan_ref)]
  finite_y <- after_equil$Logan_roi[is.finite(after_equil$Logan_roi)]

  if (length(finite_x) > 0 && length(finite_y) > 0) {
    meanval_x <- mean(finite_x)
    meanval_y <- mean(finite_y)

    minval_x <- min(finite_x)
    xmin <- min(0, minval_x - (0.1 * abs(meanval_x)))
    maxval_x <- max(finite_x)
    xmax <- maxval_x + (0.1 * abs(meanval_x))

    minval_y <- min(finite_y)
    ymin <- min(0, minval_y - (0.1 * abs(meanval_y)))
    maxval_y <- max(finite_y)
    ymax <- maxval_y + (0.1 * abs(meanval_y))
  } else {
    xmin <- NULL; xmax <- NULL
    ymin <- NULL; ymax <- NULL
  }

  # Plot

  outplot <- ggplot(data = plotdf, aes(x = Logan_ref, y = Logan_roi, colour = Equilibrium)) +
    geom_point(aes(shape = "a", size = Weights), na.rm = TRUE) +
    geom_abline(
      slope = as.numeric(refloganout$fit$coefficients[2]),
      intercept = as.numeric(refloganout$fit$coefficients[1])
    ) +
    xlab(xlabel) + ylab(ylabel) + colScale +
    guides(shape = "none", color = guide_legend(order = 1)) +
    scale_size(range = c(1, 3)) +
    coord_cartesian(xlim=c(xmin, xmax),
                    ylim=c(ymin, ymax))

  return(outplot)
}


#' Tstar Finder: Non-Invasive Logan Plot
#'
#' Function to identify where t* is for the non-invasive Logan plot.
#'
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param k2prime Value of k2prime to be used for the fitting, i.e. the average tissue-to-plasma clearance rate. This can be
#' obtained from another model, or set at a specified value. If using SRTM to estimate this value, it is equal to k2 / R1.
#' @param filename The name of the output image: filename_refLogan.jpeg
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_refLogan.jpeg
#'
#' @examples
#' \dontrun{
#' refLogan_tstar(t_tac, reftac, taclow, tacmed, tachigh, k2prime, "demonstration")
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

refLogan_tstar <- function(t_tac, reftac, lowroi, medroi, highroi, k2prime, filename = NULL, frameStartEnd = NULL, timeStartEnd = NULL, gridbreaks = 2) {
  frameStartEnd <- tstar_frameStartEnd(t_tac, frameStartEnd, timeStartEnd)

  frames <- length(reftac)
  lowroi_fit <- refLogan(t_tac, reftac, lowroi, k2prime, length(reftac), frameStartEnd = frameStartEnd)
  medroi_fit <- refLogan(t_tac, reftac, medroi, k2prime, length(reftac), frameStartEnd = frameStartEnd)
  highroi_fit <- refLogan(t_tac, reftac, highroi, k2prime, length(reftac), frameStartEnd = frameStartEnd)

  xlabel <- expression(paste("(", "",
                             paste(integral(, paste("0"), paste("", "t")),
                                   "C", phantom()[{ paste("R") }],"(",tau,")d",tau, " + ",
                                   frac(paste("C", phantom()[{ paste("R") }],"(t)"),
                                        paste("k",phantom()[{ paste("2") }],"\'",))),")",
                             " / ",
                             "C", phantom()[{ paste("T") }],"(t)"))

  ylabel <- expression(paste("", "", integral(, paste("0"), paste("", "t")),
                             "C", phantom()[{ paste("T") }],"(",tau,")d",tau, " / ",
                             "C", phantom()[{ paste("T") }],"(t)"))

  # The first frame's linearised point is 0/0 and so NA by construction:
  # both Logan and Patlak divide by a tissue or plasma concentration that
  # is zero at t = 0. na.rm declares that, rather than warning once per
  # panel about a value that is never going to be there.
  low_linplot <- ggplot(lowroi_fit$fitvals, aes(x = Logan_Ref, y = Logan_ROI)) + geom_point(na.rm = TRUE) + ggtitle("Low") + xlab(xlabel) + ylab(ylabel)
  med_linplot <- ggplot(medroi_fit$fitvals, aes(x = Logan_Ref, y = Logan_ROI)) + geom_point(na.rm = TRUE) + ggtitle("Medium") + xlab(xlabel) + ylab(ylabel)
  high_linplot <- ggplot(highroi_fit$fitvals, aes(x = Logan_Ref, y = Logan_ROI)) + geom_point(na.rm = TRUE) + ggtitle("High") + xlab(xlabel) + ylab(ylabel)

  tstarInclFrames <- 3:frames

  fitfunc <- function(roitac, tstar) {
    refLogan(t_tac, reftac, roitac, k2prime = k2prime, tstar = tstar, frameStartEnd = frameStartEnd)
  }

  comp <- tstar_compute(t_tac, lowroi, medroi, highroi, fitfunc,
                        function(f) f$par$bp, tstarInclFrames)

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow = 1)

  totalplot <- tstar_finalise_plot(
    linrow, lowroi_fit, medroi_fit, highroi_fit,
    comp$r2_df, comp$maxperc_df, comp$outcome_df,
    outcome_ylab = expression(BP[ND]), t_tac = t_tac,
    tstarInclFrames = tstarInclFrames, gridbreaks = gridbreaks,
    outcome_ylim = NULL, filename = filename, modelname = "refLogan"
  )

  return(totalplot)
}
