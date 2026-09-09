#' Patlak Plot
#'
#' Function to fit the Patlak Plot of of Patlak et al (1983) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param input Data frame containing the blood, plasma, and parent fraction
#'   concentrations over time.  This can be generated using the
#'   \code{blood_interp} function.
#' @param tstar The t* specification for regression. If tstar_type="frames",
#'   this is the number of frames from the end to include (e.g., 10 means last
#'   10 frames). If tstar_type="time", this is the time point (in minutes) after
#'   which all frames with midpoints later than this time are included. This
#'   value can be estimated using \code{Patlak_tstar}.
#' @param tstar_type Either "frames" (default) or "time", specifying how to
#'   interpret tstar.
#' @param tstarIncludedFrames Deprecated. Use 'tstar' with 'tstar_type="frames"'
#'   instead.
#' @param weights Optional. Numeric vector of the conventional frame-wise weights
#'   assigned to each frame. If not specified, uniform weights will be used.
#'   Specified weights are internally transformed to account for the dependent
#'   variable transformation in the Patlak plot. We include zero at time zero:
#'   if not included, it is added.
#' @param inpshift Optional. The number of minutes by which to shift the timing
#'   of the input data frame forwards or backwards. If not specified, this will
#'   be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will
#'   be ignored and assumed to be 0%. If specified, it will be corrected for
#'   prior to parameter estimation using the following equation: \deqn{C_{T}(t)
#'   = \frac{C_{Measured}(t) - vB\times C_{B}(t)}{1-vB}}
#' @param frameStartEnd Optional: This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This can be used to assess
#'   time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and
#'   end time point instead of defining the frame numbers using frameStartEnd.
#'   This function will restrict the model to all time frames whose t_tac is
#'   between the values, i.e. c(0,5) will select all frames with midtimes during
#'   the first 5 minutes.
#'
#'
#' @return A list with a data frame of the fitted parameters \code{out$par},,
#'   their percentage standard errors (scaled so that 1 represents 100\%)
#'   \code{out$par.se}, the model fit object \code{out$fit}, a dataframe
#'   containing the TACs of the data \code{out$tacs}, a dataframe containing the
#'   fitted values \code{out$fitvals}, the blood input data frame after time
#'   shifting \code{input}, a vector of the weights \code{out$weights}, the
#'   inpshift value used \code{inpshift}, the specified vB value \code{out$vB},
#'   and the specified tstarIncludedFrames value \code{out$tstarIncludedFrames}.
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
#' fit <- Patlakplot(t_tac, tac, input, tstar=10, weights, inpshift = 0.1)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Patlak CS, Blasberg RG, Fenstermacher JD. Graphical evaluation of
#'   blood-to-brain transfer constants from multiple-time uptake data. Journal
#'   of Cerebral Blood Flow & Metabolism. 1983 Mar 1;3(1):1-7.
#'
#' @export

Patlakplot <- function(t_tac, tac, input, tstar, weights = NULL,
                       inpshift = 0, vB = 0, frameStartEnd = NULL, timeStartEnd = NULL,
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

  tidyinput <- tidyinput_art(t_tac, tac, weights, frameStartEnd)

  t_tac <- tidyinput$t_tac
  tac <- tidyinput$tac
  weights <- tidyinput$weights


  # Convert tstar based on type
  if (tstar_type == "time") {
    frames_after_tstar <- which(t_tac >= tstar)
    tstarIncludedFrames <- length(frames_after_tstar)
  } else {
    tstarIncludedFrames <- tstar
  }


  newvals <- shift_timings(
    t_tac = t_tac,
    tac = tac,
    input = input,
    inpshift = inpshift
  )


  t_tac <- newvals$t_tac
  tac <- newvals$tac

  t_inp <- newvals$input$Time
  blood <- newvals$input$Blood
  aif <- newvals$input$AIF

  # Parameters

  interptime <- newvals$input$Time
  i_tac <- pracma::interp1(t_tac, tac, interptime, method = "linear")

  # Blood Volume Correction (nothing happens if vB = 0)
  i_tac <- (i_tac - vB * blood) / (1 - vB)
  tac_uncor <- tac
  tac <- pracma::interp1(interptime, i_tac, t_tac, method = "linear")

  # Transform weights for graphical analysis (if provided)
  # Check if real weights were provided (more than just 0s and 1s from tidyinput)
  unique_weights <- unique(weights[is.finite(weights)])
  real_weights_provided <- length(setdiff(unique_weights, c(0, 1))) > 0
  if (!is.null(weights) && real_weights_provided) {
    weights <- weights_Patlak_transform(t_tac, tac, newvals$input, weights)
    # Center weights so mean of equilibrium weights equals 1
    equil_weights <- tail(weights, tstarIncludedFrames)
    equil_mean <- mean(equil_weights[is.finite(equil_weights)])
    if (is.finite(equil_mean) && equil_mean > 0) {
      weights <- weights / equil_mean
    }
    # Set pre-equilibrium weights to 1
    pre_equil_idx <- seq_len(length(weights) - tstarIncludedFrames)
    weights[pre_equil_idx] <- 1
  }

  patlak_roi <- i_tac / aif
  patlak_plasma <- as.numeric((pracma::cumtrapz(interptime, aif)) / aif)

  patlak_roi <- pracma::interp1(interptime, patlak_roi, t_tac, method = "linear")
  patlak_plasma <- pracma::interp1(interptime, patlak_plasma, t_tac, method = "linear")

  patlak_equil_roi <- tail(patlak_roi, tstarIncludedFrames)
  patlak_equil_plasma <- tail(patlak_plasma, tstarIncludedFrames)
  weights_equil <- tail(weights, tstarIncludedFrames)


  # Solution

  patlak_model <- lm(patlak_equil_roi ~ patlak_equil_plasma, weights = weights_equil)

  # Output

  par <- as.data.frame(list(Ki = as.numeric(patlak_model$coefficients[2])))
  fit <- patlak_model

  tacs <- data.frame(Time = t_tac, Target = tac, Target_uncor = tac_uncor) # uncorrected for blood volume

  fitvals <- data.frame(Patlak_Plasma = patlak_plasma, Patlak_ROI = patlak_roi)

  input <- newvals$input

  par.se <- par
  par.se[1,] <- purrr::map_dbl(names(coef(patlak_model)), ~ get_se(patlak_model, .x))[2]
  names(par.se) <- paste0(names(par.se), ".se")

  out <- list(
    par = par, par.se = par.se, fit = patlak_model, tacs = tacs,
    fitvals = fitvals, input = input, weights = weights,
    inpshift = inpshift, vB = vB, tstarIncludedFrames = tstarIncludedFrames,
    model = "Patlak"
  )

  class(out) <- c("Patlak", "kinfit")

  return(out)
}

#' Plot: Patlak Plot
#'
#' Function to visualise the fit of the Patlak Plot model to data.
#'
#' @param patlakout The output object of the Patlak Plot fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
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
#' fit <- Patlakplot(t_tac, tac, input, 10, weights, inpshift = 0.1)
#'
#' plot_Patlakfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_Patlakfit <- function(patlakout, roiname = NULL) {
  plotdf <- data.frame(
    Weights = patlakout$weights,
    Patlak_Plasma = patlakout$fitvals$Patlak_Plasma,
    Patlak_ROI = patlakout$fitvals$Patlak_ROI,
    Equilibrium = as.character("Before")
  )

  plotdf$Equilibrium <- as.character(plotdf$Equilibrium)
  plotdf$Equilibrium [ (nrow(plotdf) - (patlakout$tstarIncludedFrames - 1)):nrow(plotdf)  ] <- "After"

  # Set pre-tstar weights to 1 for display (so points are visible but don't affect scale)
  plotdf$Weights[plotdf$Equilibrium == "Before"] <- 1

  plotdf$Equilibrium <- forcats::fct_inorder(factor(plotdf$Equilibrium))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Equilibrium)
  colScale <- scale_colour_manual(name = paste0(roiname, "\nLinear"),
                                  values = myColors)

  xlimits <- c(0, tail(plotdf$Patlak_Plasma, 1))

  xlabel <- expression(paste("", "", integral(, paste("0"), paste("", "t")),
                             "C", phantom()[{ paste("P") }],"(",tau,")d",tau, " / ",
                             "C", phantom()[{ paste("P") }],"(t)"))

  ylabel <- expression(paste("C", phantom()[{ paste("T") }],"(t)", " / ",
                             "C", phantom()[{ paste("P") }],"(t)"))

  outplot <- ggplot(data = plotdf, aes(x = Patlak_Plasma, y = Patlak_ROI, colour = Equilibrium)) +
    geom_point(aes(shape = "a", size = Weights), na.rm = TRUE) +
    geom_abline(
      slope = as.numeric(patlakout$fit$coefficients[2]),
      intercept = as.numeric(patlakout$fit$coefficients[1])
    ) +
    xlab(xlabel) + ylab(ylabel) + coord_cartesian(xlim = xlimits) + colScale +
    guides(shape = "none", color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

  return(outplot)
}

#' Tstar Finder: Patlak Plot
#'
#' Function to identify where t* is for the Patlak Plot.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param filename The name of the output image: filename_Patlak.jpeg
#' @param inpshift Optional. The number of minutes by which to shift the timing of the input data frame forwards or backwards.
#' If not specified, this will be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will
#'   be ignored and assumed to be 0%. If specified, it will be corrected for
#'   prior to parameter estimation using the following equation: \deqn{C_{T}(t)
#'   = \frac{C_{Measured}(t) - vB\times C_{B}(t)}{1-vB}}
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_Patlakplot.jpeg
#'
#' @examples
#' \dontrun{
#' Patlak_tstar(t_tac, lowroi, medroi, highroi, input,
#'   filename = "demonstration",
#'   inpshift = onetcmout$par$inpshift, vB = 0.05, gridbreaks = 4
#' )
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

Patlak_tstar <- function(t_tac, lowroi, medroi, highroi, input, filename = NULL, inpshift = 0, vB = 0, frameStartEnd = NULL, timeStartEnd = NULL, gridbreaks = 2) {
  frameStartEnd <- tstar_frameStartEnd(t_tac, frameStartEnd, timeStartEnd)

  frames <- length(t_tac)
  lowroi_fit <- Patlakplot(t_tac, lowroi, input, tstar = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  medroi_fit <- Patlakplot(t_tac, medroi, input, tstar = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  highroi_fit <- Patlakplot(t_tac, highroi, input, tstar = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)

  low_xlimits <- c(0, tail(lowroi_fit$fitvals$Patlak_Plasma, 1))
  med_xlimits <- c(0, tail(medroi_fit$fitvals$Patlak_Plasma, 1))
  high_xlimits <- c(0, tail(highroi_fit$fitvals$Patlak_Plasma, 1))

  xlabel <- expression(paste("", "", integral(, paste("0"), paste("", "t")),
                             "C", phantom()[{ paste("P") }],"(",tau,")d",tau, " / ",
                             "C", phantom()[{ paste("P") }],"(t)"))

  ylabel <- expression(paste("C", phantom()[{ paste("T") }],"(t)", " / ",
                             "C", phantom()[{ paste("P") }],"(t)"))

  # The first frame's linearised point is 0/0 and so NA by construction:
  # both Logan and Patlak divide by a tissue or plasma concentration that
  # is zero at t = 0. na.rm declares that, rather than warning once per
  # panel about a value that is never going to be there.
  low_linplot <- ggplot(lowroi_fit$fitvals, aes(x = Patlak_Plasma, y = Patlak_ROI)) + geom_point(na.rm = TRUE) + ggtitle("Low") + xlab(xlabel) + ylab(ylabel) + coord_cartesian(xlim = low_xlimits)
  med_linplot <- ggplot(medroi_fit$fitvals, aes(x = Patlak_Plasma, y = Patlak_ROI)) + geom_point(na.rm = TRUE) + ggtitle("Medium") + xlab(xlabel) + ylab(ylabel) + coord_cartesian(xlim = med_xlimits)
  high_linplot <- ggplot(highroi_fit$fitvals, aes(x = Patlak_Plasma, y = Patlak_ROI)) + geom_point(na.rm = TRUE) + ggtitle("High") + xlab(xlabel) + ylab(ylabel) + coord_cartesian(xlim = high_xlimits)

  tstarInclFrames <- 3:frames

  fitfunc <- function(roitac, tstar) {
    Patlakplot(t_tac, roitac, input, tstar = tstar, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  }

  comp <- tstar_compute(t_tac, lowroi, medroi, highroi, fitfunc,
                        function(f) f$par$Ki, tstarInclFrames)

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow = 1)

  totalplot <- tstar_finalise_plot(
    linrow, lowroi_fit, medroi_fit, highroi_fit,
    comp$r2_df, comp$maxperc_df, comp$outcome_df,
    outcome_ylab = expression(K[i]), t_tac = t_tac,
    tstarInclFrames = tstarInclFrames, gridbreaks = gridbreaks,
    outcome_ylim = NULL, filename = filename, modelname = "Patlakplot"
  )

  return(totalplot)
}
