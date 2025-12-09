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
    geom_point(aes(shape = "a", size = Weights)) +
    geom_abline(
      slope = as.numeric(patlakout$fit$coefficients[2]),
      intercept = as.numeric(patlakout$fit$coefficients[1])
    ) +
    xlab(xlabel) + ylab(ylabel) + xlim(xlimits) + colScale +
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
  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

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

  low_linplot <- qplot(lowroi_fit$fitvals$Patlak_Plasma, lowroi_fit$fitvals$Patlak_ROI) + ggtitle("Low") + xlab(xlabel) + ylab(ylabel) + xlim(low_xlimits)
  med_linplot <- qplot(medroi_fit$fitvals$Patlak_Plasma, medroi_fit$fitvals$Patlak_ROI) + ggtitle("Medium") + xlab(xlabel) + ylab(ylabel) + xlim(med_xlimits)
  high_linplot <- qplot(highroi_fit$fitvals$Patlak_Plasma, highroi_fit$fitvals$Patlak_ROI) + ggtitle("High") + xlab(xlabel) + ylab(ylabel) + xlim(high_xlimits)

  tstarInclFrames <- 3:frames
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  Ki_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)

  for (i in 1:length(tstarInclFrames)) {
    lowfit <- Patlakplot(t_tac, lowroi, input, tstar = tstarInclFrames[i], inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
    medfit <- Patlakplot(t_tac, medroi, input, tstar = tstarInclFrames[i], inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
    highfit <- Patlakplot(t_tac, highroi, input, tstar = tstarInclFrames[i], inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)

    r2_df$Low[i] <- summary(lowfit$fit)$r.squared
    r2_df$Medium[i] <- summary(medfit$fit)$r.squared
    r2_df$High[i] <- summary(highfit$fit)$r.squared

    maxperc_df$Low[i] <- maxpercres(lowfit)
    maxperc_df$Medium[i] <- maxpercres(medfit)
    maxperc_df$High[i] <- maxpercres(highfit)

    Ki_df$Low[i] <- lowfit$par$Ki
    Ki_df$Medium[i] <- medfit$par$Ki
    Ki_df$High[i] <- highfit$par$Ki
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

  tacplotdf <- data.frame(cbind(Time = lowroi_fit$tacs$Time, Low = lowroi_fit$tacs$Target, Medium = medroi_fit$tacs$Target, High = highroi_fit$tacs$Target))
  tacplotdf <- tidyr::gather(tacplotdf, key = Region, value = Radioactivity, -Time)

  tacplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(tacplotdf$Region)))

  myColors <- RColorBrewer::brewer.pal(4, "Set1")
  names(myColors) <- levels(tacplotdf$Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  tacplot <- ggplot(tacplotdf, aes(x = Time, y = Radioactivity, colour = Region)) + geom_point() + geom_line() + colScale


  # Ki Plot

  Kiplotdf <- tidyr::gather(Ki_df, key = Region, value = Ki, -Frames, -Time)
  Kiplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(Kiplotdf$Region)))

  Kiplot <- ggplot(Kiplotdf, aes(x = Frames, y = Ki, colour = Region)) + geom_point() + geom_line() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylab(expression(K[i])) + colScale


  # Output

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow = 1)
  r2row <- cowplot::plot_grid(low_r2plot, med_r2plot, high_r2plot, nrow = 1)
  mprow <- cowplot::plot_grid(low_mpplot, med_mpplot, high_mpplot, nrow = 1)
  outrow <- cowplot::plot_grid(tacplot, Kiplot, rel_widths = c(2, 1))

  totalplot <- cowplot::plot_grid(linrow, r2row, mprow, outrow, nrow = 4)

  if (!is.null(filename)) {
    jpeg(filename = paste0(filename, "_Patlakplot.jpeg"), width = 300, height = 400, units = "mm", res = 600)
    totalplot
    dev.off()
  }

  return(totalplot)
}
