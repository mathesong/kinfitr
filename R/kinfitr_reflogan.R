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

  outplot <- ggplot(data = plotdf, aes(x = Logan_ref, y = Logan_roi, colour = Equilibrium)) +
    geom_point(aes(shape = "a", size = Weights)) +
    geom_abline(
      slope = as.numeric(refloganout$fit$coefficients[2]),
      intercept = as.numeric(refloganout$fit$coefficients[1])
    ) +
    xlab(xlabel) + ylab(ylabel) + colScale +
    guides(shape = "none", color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

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
  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

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

  low_linplot <- qplot(lowroi_fit$fitvals$Logan_Ref, lowroi_fit$fitvals$Logan_ROI) + ggtitle("Low") + xlab(xlabel) + ylab(ylabel)
  med_linplot <- qplot(medroi_fit$fitvals$Logan_Ref, medroi_fit$fitvals$Logan_ROI) + ggtitle("Medium") + xlab(xlabel) + ylab(ylabel)
  high_linplot <- qplot(highroi_fit$fitvals$Logan_Ref, highroi_fit$fitvals$Logan_ROI) + ggtitle("High") + xlab(xlabel) + ylab(ylabel)

  tstarInclFrames <- 3:frames
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  bp_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)

  for (i in 1:length(tstarInclFrames)) {
    lowfit <- refLogan(t_tac, reftac, lowroi, k2prime = k2prime, tstar = tstarInclFrames[i], frameStartEnd = frameStartEnd)
    medfit <- refLogan(t_tac, reftac, medroi, k2prime = k2prime, tstar = tstarInclFrames[i], frameStartEnd = frameStartEnd)
    highfit <- refLogan(t_tac, reftac, highroi, k2prime = k2prime, tstar = tstarInclFrames[i], frameStartEnd = frameStartEnd)

    r2_df$Low[i] <- summary(lowfit$fit)$r.squared
    r2_df$Medium[i] <- summary(medfit$fit)$r.squared
    r2_df$High[i] <- summary(highfit$fit)$r.squared

    maxperc_df$Low[i] <- maxpercres(lowfit)
    maxperc_df$Medium[i] <- maxpercres(medfit)
    maxperc_df$High[i] <- maxpercres(highfit)

    bp_df$Low[i] <- lowfit$par$bp
    bp_df$Medium[i] <- medfit$par$bp
    bp_df$High[i] <- highfit$par$bp
  }

  xlabel <- "Number of Included Frames"
  ylab_r2 <- expression(R^2)
  ylab_mp <- "Maximum Percentage Deviation"


  # R Squared plots

  low_r2plot <- ggplot(r2_df, aes(x = Frames, y = Low)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim = c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)
  med_r2plot <- ggplot(r2_df, aes(x = Frames, y = Medium)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim = c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)
  high_r2plot <- ggplot(r2_df, aes(x = Frames, y = High)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim = c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)

  # Max Percentage Variation Plots

  maxperc_df$inclmins <- rev(max(t_tac) - t_tac)[-c(1, 2)]
  maxperc_df$tstar <- rev(t_tac)[-c(1, 2)]

  low_mpplot <- ggplot(maxperc_df, aes(x = Frames, y = Low)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim = c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text", x = 3, y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Low + 1.4, label = round(maxperc_df$tstar, 1), size = 3, colour = "red") + annotate("text", x = 3, y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Low + 0.7, label = round(maxperc_df$inclmins, 1), size = 3, colour = "blue")
  med_mpplot <- ggplot(maxperc_df, aes(x = Frames, y = Medium)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim = c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text", x = 3, y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Medium + 1.4, label = round(maxperc_df$tstar, 1), size = 3, colour = "red") + annotate("text", x = 3, y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Medium + 0.7, label = round(maxperc_df$inclmins, 1), size = 3, colour = "blue")
  high_mpplot <- ggplot(maxperc_df, aes(x = Frames, y = High)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim = c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text", x = 3, y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$High + 1.4, label = round(maxperc_df$tstar, 1), size = 3, colour = "red") + annotate("text", x = 3, y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$High + 0.7, label = round(maxperc_df$inclmins, 1), size = 3, colour = "blue")


  # TAC Plot

  tacplotdf <- data.frame(cbind(Time = lowroi_fit$tacs$Time, Reference = lowroi_fit$tacs$Reference, Low = lowroi_fit$tacs$Target, Medium = medroi_fit$tacs$Target, High = highroi_fit$tacs$Target))
  tacplotdf <- tidyr::gather(tacplotdf, key = Region, value = Radioactivity, -Time)

  tacplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(tacplotdf$Region)))

  myColors <- RColorBrewer::brewer.pal(4, "Set1")
  names(myColors) <- levels(tacplotdf$Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  tacplot <- ggplot(tacplotdf, aes(x = Time, y = Radioactivity, colour = Region)) + geom_point() + geom_line() + colScale


  # BP Plot

  bpplotdf <- tidyr::gather(bp_df, key = Region, value = BP, -Frames, -Time)
  bpplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(bpplotdf$Region)))

  bpplot <- ggplot(bpplotdf, aes(x = Frames, y = BP, colour = Region)) +
    geom_point() + geom_line() +
    scale_x_continuous(breaks = seq(min(tstarInclFrames),
                                    max(tstarInclFrames), by = gridbreaks)) +
    ylab(expression(BP[ND])) +
    colScale


  # Output

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow = 1)
  r2row <- cowplot::plot_grid(low_r2plot, med_r2plot, high_r2plot, nrow = 1)
  mprow <- cowplot::plot_grid(low_mpplot, med_mpplot, high_mpplot, nrow = 1)
  outrow <- cowplot::plot_grid(tacplot, bpplot, rel_widths = c(2, 1))

  totalplot <- cowplot::plot_grid(linrow, r2row, mprow, outrow, nrow = 4)

  if (!is.null(filename)) {
    jpeg(filename = paste0(filename, "_refLogan.jpeg"), width = 300, height = 400, units = "mm", res = 600)
    totalplot
    dev.off()
  }

  return(totalplot)
}
