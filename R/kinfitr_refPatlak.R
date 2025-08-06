#' Patlak Reference Tissue Model
#'
#' Function to fit the Patlak Reference Tissue Model of Patlak & Blasbert (1985) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame. We include zero at
#' time zero: if not included, it is added.
#' @param roitac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param tstarIncludedFrames The number of frames to be used in the regression model, i.e. the number of frames for which
#' the function is linear after pseudo-equilibrium is reached. This is a count from the end of the measurement, so a value of
#' 10 means that last 10 frames will be used. This value can be estimated using \code{refPatlak_tstar}.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param dur Optional. Numeric vector of the time durations of the frames. If
#' not included, the integrals will be calculated using trapezoidal integration.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit}, a dataframe
#' containing the TACs of the data \code{out$tacs}, a dataframe containing the TACs of the fitted values \code{out$fitvals},
#' a vector of the weights \code{out$weights}, and the specified tstarIncludedFrames value \code{out$tstarIncludedFrames}
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
#' fit <- refPatlak(t_tac, reftac, roitac, tstarIncludedFrames = 10, weights = weights)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Patlak CS, Blasberg RG. Graphical evaluation of blood-to-brain transfer constants from multiple-time uptake data. Generalizations. Journal of Cerebral Blood Flow & Metabolism. 1985 Dec 1;5(4):584-90.
#'
#' @export

refPatlak <- function(t_tac, reftac, roitac, tstarIncludedFrames,
                      weights = NULL, dur = NULL, frameStartEnd = NULL, timeStartEnd = NULL) {

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1], 
                       tail(which(t_tac <= timeStartEnd[2]), 1))
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


  # Parameters

  if (!is.null(dur)) {

    patlak_roi <- roitac / reftac
    patlak_ref <- frame_cumsum(dur, reftac) / reftac

  } else {

    patlak_roi <- roitac / reftac
    patlak_ref <- pracma::cumtrapz(t_tac, reftac) / reftac

  }

  patlak_equil_roi <- tail(patlak_roi, tstarIncludedFrames)
  patlak_equil_ref <- tail(patlak_ref, tstarIncludedFrames)
  weights_equil <- tail(weights, tstarIncludedFrames)


  # Solution

  patlak_model <- lm(patlak_equil_roi ~ patlak_equil_ref, weights = weights_equil)

  # Output

  par <- as.data.frame(list(K = as.numeric(patlak_model$coefficients[2])))
  fit <- patlak_model

  tacs <- data.frame(Time = t_tac, Reference = reftac, Target = roitac )

  if (!is.null(dur)) {
    tacs$Duration <- dur
  }

  fitvals <- data.frame(Patlak_ROI = patlak_roi, Patlak_Ref = patlak_ref)

  out <- list(
    par = par, fit = fit, tacs = tacs, fitvals = fitvals, weights = weights,
    tstarIncludedFrames = tstarIncludedFrames, model = "refPatlak"
  )

  class(out) <- c("refPatlak", "kinfit")

  return(out)
}

#' Plot: Patlak Reference Tissue Model
#'
#' Function to visualise the fit of the refPatlak model to data.
#'
#' @param refpatlakout The output object of the refPatlak fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
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
#' fit <- refPatlak(t_tac, reftac, roitac, tstarIncludedFrames = 10, weights = weights)
#'
#' plot_refPatlakfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_refPatlakfit <- function(refpatlakout, roiname = NULL) {
  plotdf <- data.frame(
    Weights = refpatlakout$weights,
    Patlak_ref = refpatlakout$fitvals$Patlak_Ref,
    Patlak_roi = refpatlakout$fitvals$Patlak_ROI,
    Equilibrium = as.character("Before")
  )


  plotdf$Equilibrium <- as.character(plotdf$Equilibrium)
  plotdf$Equilibrium [ (nrow(plotdf) - (refpatlakout$tstarIncludedFrames - 1)):nrow(plotdf)  ] <- "After"

  plotdf$Equilibrium <- forcats::fct_inorder(factor(plotdf$Equilibrium))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Equilibrium)
  colScale <- scale_colour_manual(name = paste0(roiname, "\nEquilibrium"), values = myColors)

  ylabel <- expression(paste("C", phantom()[{ paste("T") }],"(t)", " / ",
                                  "C", phantom()[{ paste("R") }],"(t)"))

  xlabel <- expression(paste("", "", integral(, paste("0"), paste("", "t")),
                                  "C", phantom()[{ paste("R") }],"(",tau,")d",tau, " / ",
                                  "C", phantom()[{ paste("R") }],"(t)"))

  outplot <- ggplot(data = plotdf, aes(x = Patlak_ref, y = Patlak_roi, colour = Equilibrium)) +
    geom_point(aes(shape = "a", size = Weights)) +
    geom_abline(
      slope = as.numeric(refpatlakout$fit$coefficients[2]),
      intercept = as.numeric(refpatlakout$fit$coefficients[1])
    ) +
    xlab(xlabel) + ylab(ylabel) + colScale +
    guides(shape = "none", color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

  return(outplot)
}


#' Tstar Finder: Patlak Reference Tissue Model
#'
#' Function to identify where t* is for the Patlak Reference Tissue Model.
#'
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param filename The name of the output image: filename_refPatlak.jpeg
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_refPatlak.jpeg
#'
#' @examples
#' \dontrun{
#' refPatlak_tstar(t_tac, reftac, taclow, tacmed, tachigh, "demonstration")
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

refPatlak_tstar <- function(t_tac, reftac, lowroi, medroi, highroi, filename = NULL, frameStartEnd = NULL, timeStartEnd = NULL, gridbreaks = 2) {
  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1], 
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  frames <- length(reftac)
  lowroi_fit <- refPatlak(t_tac, reftac, lowroi, length(reftac), frameStartEnd = frameStartEnd)
  medroi_fit <- refPatlak(t_tac, reftac, medroi, length(reftac), frameStartEnd = frameStartEnd)
  highroi_fit <- refPatlak(t_tac, reftac, highroi, length(reftac), frameStartEnd = frameStartEnd)

  ylabel <- expression(paste("C", phantom()[{ paste("T") }],"(t)", " / ",
                             "C", phantom()[{ paste("R") }],"(t)"))

  xlabel <- expression(paste("", "", integral(, paste("0"), paste("", "t")),
                             "C", phantom()[{ paste("R") }],"(",tau,")d",tau, " / ",
                             "C", phantom()[{ paste("R") }],"(t)"))

  low_linplot <- qplot(lowroi_fit$fitvals$Patlak_Ref, lowroi_fit$fitvals$Patlak_ROI) + ggtitle("Low") + xlab(xlabel) + ylab(ylabel)
  med_linplot <- qplot(medroi_fit$fitvals$Patlak_Ref, medroi_fit$fitvals$Patlak_ROI) + ggtitle("Medium") + xlab(xlabel) + ylab(ylabel)
  high_linplot <- qplot(highroi_fit$fitvals$Patlak_Ref, highroi_fit$fitvals$Patlak_ROI) + ggtitle("High") + xlab(xlabel) + ylab(ylabel)

  tstarInclFrames <- 3:frames
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  k_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)

  for (i in 1:length(tstarInclFrames)) {
    lowfit <- refPatlak(t_tac, reftac, lowroi, tstarIncludedFrames = tstarInclFrames[i], frameStartEnd = frameStartEnd)
    medfit <- refPatlak(t_tac, reftac, medroi, tstarIncludedFrames = tstarInclFrames[i], frameStartEnd = frameStartEnd)
    highfit <- refPatlak(t_tac, reftac, highroi, tstarIncludedFrames = tstarInclFrames[i], frameStartEnd = frameStartEnd)

    r2_df$Low[i] <- summary(lowfit$fit)$r.squared
    r2_df$Medium[i] <- summary(medfit$fit)$r.squared
    r2_df$High[i] <- summary(highfit$fit)$r.squared

    maxperc_df$Low[i] <- maxpercres(lowfit)
    maxperc_df$Medium[i] <- maxpercres(medfit)
    maxperc_df$High[i] <- maxpercres(highfit)

    k_df$Low[i] <- lowfit$par$K
    k_df$Medium[i] <- medfit$par$K
    k_df$High[i] <- highfit$par$K
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


  # K Plot

  kplotdf <- tidyr::gather(k_df, key = Region, value = K, -Frames, -Time)
  kplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(kplotdf$Region)))

  kplot <- ggplot(kplotdf, aes(x = Frames, y = K, colour = Region)) + geom_point() + geom_line() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylab(expression(K[i])) + colScale


  # Output

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow = 1)
  r2row <- cowplot::plot_grid(low_r2plot, med_r2plot, high_r2plot, nrow = 1)
  mprow <- cowplot::plot_grid(low_mpplot, med_mpplot, high_mpplot, nrow = 1)
  outrow <- cowplot::plot_grid(tacplot, kplot, rel_widths = c(2, 1))

  totalplot <- cowplot::plot_grid(linrow, r2row, mprow, outrow, nrow = 4)

  if (!is.null(filename)) {
    jpeg(filename = paste0(filename, "_refPatlak.jpeg"), width = 300, height = 400, units = "mm", res = 600)
    totalplot
    dev.off()
  }

  return(totalplot)
}
