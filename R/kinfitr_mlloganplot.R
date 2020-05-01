#' Multilinear Logan Plot
#'
#' Function to fit the multilinear Logan Plot of Turkheimer et al (2003) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param tstarIncludedFrames The number of frames to be used in the regression model, i.e. the number of frames for which
#' the function is linear after pseudo-equilibrium is reached. This is a count from the end of the measurement, so a value of
#' 10 means that last 10 frames will be used. This value can be estimated using \code{Logan_tstar}.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param inpshift Optional. The number of minutes by which to shift the timing of the input data frame forwards or backwards.
#' If not specified, this will be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will be ignored and assumed to be 0%. If specified, it
#' will be corrected for prior to parameter estimation using the following equation:
#' \deqn{C_{T}(t) = \frac{C_{Measured}(t) - vB\times C_{B}(t)}{1-vB}}
#' @param dur Optional. Numeric vector of the time durations of the frames. If
#' not included, the integrals will be calculated using trapezoidal integration.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#'
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit},
#' a dataframe containing the TACs of the data \code{out$tacs}, a dataframe containing the fitted values \code{out$fitvals},
#' the blood input data frame after time shifting \code{input}, a vector of the weights \code{out$weights},
#' the inpshift value used \code{inpshift}, the specified vB value \code{out$vB}, and the specified tstarIncludedFrames
#' value \code{out$tstarIncludedFrames}.
#'
#' @examples
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' tac <- pbr28$tacs[[2]]$FC
#' weights <- pbr28$tacs[[2]]$Weights
#' dur <- pbr28$tacs[[2]]$Duration/60
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1
#' )
#'
#' fit1 <- mlLoganplot(t_tac, tac, input, 10, weights)
#' fit2 <- mlLoganplot(t_tac, tac, input, 10, weights, inpshift = 0.1, vB = 0.05)
#' fit3 <- mlLoganplot(t_tac, tac, input, 10, weights, inpshift = 0.1, vB = 0.05, dur=dur)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Turkheimer FE, Aston JA, Banati RB, Riddell C, Cunningham VJ. A linear wavelet filter for parametric imaging with dynamic PET. IEEE transactions on medical imaging. 2003 Mar;22(3):289-301.
#'
#' @export


mlLoganplot <- function(t_tac, tac, input, tstarIncludedFrames, weights = NULL,
                        inpshift = 0, vB = 0, dur = NULL, frameStartEnd = NULL) {


  # Tidying

  tidyinput <- tidyinput_art(t_tac, tac, weights, frameStartEnd)

  if (!is.null(dur)) {
    tidyinput_dur <- tidyinput_art(dur, tac, weights, frameStartEnd)
    dur <- tidyinput_dur$t_tac
  }

  t_tac <- tidyinput$t_tac
  tac <- tidyinput$tac
  weights <- tidyinput$weights


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

  if (!is.null(dur)) {

    term1_dv <- frame_cumsum(dur, tac)
    term2 <- as.numeric(pracma::cumtrapz(interptime, aif))
    term3 <- -tac

    term2 <- pracma::interp1(interptime, term2, t_tac, method = "linear")

  } else {

    term1_dv <- as.numeric(pracma::cumtrapz(interptime, i_tac))
    term2 <- as.numeric(pracma::cumtrapz(interptime, aif))
    term3 <- -i_tac

    term1_dv <- pracma::interp1(interptime, term1_dv, t_tac, method = "linear")
    term2 <- pracma::interp1(interptime, term2, t_tac, method = "linear")
    term3 <- pracma::interp1(interptime, term3, t_tac, method = "linear")

  }

  fitvals <- data.frame(
    Time = t_tac, Weights = weights,
    Term1_DV = term1_dv, Term2 = term2, Term3 = term3
  )

  equil <- rep("Before", length(tac))
  equil[(length(equil) - tstarIncludedFrames + 1):length(equil)] <- "After"

  fitvals$Equilibrium <- equil

  fitvals_equil <- subset(fitvals, fitvals$Equilibrium == "After")


  # Solution

  mllogan_model <- lm(Term1_DV ~ Term2 + Term3 - 1, weights = Weights, data = fitvals_equil)

  # Output

  Vt <- as.numeric(mllogan_model$coefficients[1])
  k2 <- as.numeric(1 / mllogan_model$coefficients[2])
  par <- as.data.frame(list(Vt = Vt, k2 = k2))

  tacs <- data.frame(Time = t_tac, Target = tac, Target_uncor = tac_uncor) # uncorrected for blood volume

  fittedvals <- (as.numeric(mllogan_model$coefficients[1]) * term2) + (as.numeric(mllogan_model$coefficients[2] * term3))
  fitvals <- data.frame(Term1_DV = term1_dv, Term2 = term2, Term3 = term3, Fitted = fittedvals)

  input <- newvals$input

  out <- list(
    par = par, fit = mllogan_model, tacs = tacs,
    fitvals = fitvals, input = input, weights = weights,
    inpshift = inpshift, vB = vB, tstarIncludedFrames = tstarIncludedFrames,
    model = "mlLogan"
  )

  class(out) <- c("mlLogan", "kinfit")

  return(out)
}


#' Plot: Multilinear Logan Plot
#'
#' Function to visualise the fit of the multilinear Logan Plot model to data.
#'
#' @param mlloganout The output object of the multilinear Logan Plot fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
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
#' fit <- mlLoganplot(t_tac, tac, input, 10, weights)
#' plot_mlLoganfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_mlLoganfit <- function(mlloganout, roiname = NULL) {
  plotdf <- data.frame(
    Weights = mlloganout$weights,
    Term1_DV = mlloganout$fitvals$Term1_DV,
    Fitted = mlloganout$fitvals$Fitted,
    Equilibrium = as.character("Before")
  )


  plotdf$Equilibrium <- as.character(plotdf$Equilibrium)
  plotdf$Equilibrium [ (nrow(plotdf) - (mlloganout$tstarIncludedFrames - 1)):nrow(plotdf)  ] <- "After"

  plotdf$Equilibrium <- forcats::fct_inorder(factor(plotdf$Equilibrium))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Equilibrium)
  colScale <- scale_colour_manual(name = paste0(roiname, "\nEquilibrium"), values = myColors)

  outplot <- ggplot(data = plotdf, aes(x = Fitted, y = Term1_DV, colour = Equilibrium)) +
    geom_point(aes(shape = "a", size = Weights)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("Fitted Values") + ylab("Integ(C_Tissue)") + colScale +
    guides(shape = FALSE, color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

  return(outplot)
}


#' Tstar Finder: Multilinear Logan Plot
#'
#' Function to identify where t* is for the multilinear Logan Plot.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param filename The name of the output image: filename_Logan.jpeg
#' @param inpshift Optional. The number of minutes by which to shift the timing of the input data frame forwards or backwards.
#' If not specified, this will be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will
#'   be ignored and assumed to be 0%. If specified, it will be corrected for
#'   prior to parameter estimation using the following equation: \deqn{C_{T}(t)
#'   = \frac{C_{Measured}(t) - vB\times C_{B}(t)}{1-vB}}
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_mlLogan.jpeg
#'
#' @examples
#' \dontrun{
#' mlLogan_tstar(t_tac, lowroi, medroi, highroi, input,
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

mlLogan_tstar <- function(t_tac, lowroi, medroi, highroi, input, filename = NULL, inpshift = 0, vB = 0, frameStartEnd = NULL, gridbreaks = 2) {
  frames <- length(t_tac)
  lowroi_fit <- mlLoganplot(t_tac, lowroi, input, tstarIncludedFrames = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  medroi_fit <- mlLoganplot(t_tac, medroi, input, tstarIncludedFrames = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  highroi_fit <- mlLoganplot(t_tac, highroi, input, tstarIncludedFrames = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)

  xlabel <- "Fitted Values"
  ylabel <- "Integ(C_Tissue)"

  low_linplot <- qplot(lowroi_fit$fitvals$Fitted, lowroi_fit$fitvals$Term1_DV) + ggtitle("Low") + xlab(xlabel) + ylab(ylabel)
  med_linplot <- qplot(medroi_fit$fitvals$Fitted, medroi_fit$fitvals$Term1_DV) + ggtitle("Medium") + xlab(xlabel) + ylab(ylabel)
  high_linplot <- qplot(highroi_fit$fitvals$Fitted, highroi_fit$fitvals$Term1_DV) + ggtitle("High") + xlab(xlabel) + ylab(ylabel)

  tstarInclFrames <- 3:frames
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  vt_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)

  for (i in 1:length(tstarInclFrames)) {
    lowfit <- mlLoganplot(t_tac, lowroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
    medfit <- mlLoganplot(t_tac, medroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
    highfit <- mlLoganplot(t_tac, highroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)

    r2_df$Low[i] <- summary(lowfit$fit)$r.squared
    r2_df$Medium[i] <- summary(medfit$fit)$r.squared
    r2_df$High[i] <- summary(highfit$fit)$r.squared

    maxperc_df$Low[i] <- maxpercres(lowfit)
    maxperc_df$Medium[i] <- maxpercres(medfit)
    maxperc_df$High[i] <- maxpercres(highfit)

    vt_df$Low[i] <- lowfit$par$Vt
    vt_df$Medium[i] <- medfit$par$Vt
    vt_df$High[i] <- highfit$par$Vt
  }

  xlabel <- "Number of Included Frames"
  ylab_r2 <- expression(R^2)
  ylab_mp <- "Maximum Percentage Variance"


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


  # Vt Plot

  vtplotdf <- tidyr::gather(vt_df, key = Region, value = Vt, -Frames, -Time)
  vtplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(vtplotdf$Region)))

  vtplot <- ggplot(vtplotdf, aes(x = Frames, y = Vt, colour = Region)) + geom_point() + geom_line() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylab(expression(V[T])) + colScale


  # Output

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow = 1)
  r2row <- cowplot::plot_grid(low_r2plot, med_r2plot, high_r2plot, nrow = 1)
  mprow <- cowplot::plot_grid(low_mpplot, med_mpplot, high_mpplot, nrow = 1)
  outrow <- cowplot::plot_grid(tacplot, vtplot, rel_widths = c(2, 1))

  totalplot <- cowplot::plot_grid(linrow, r2row, mprow, outrow, nrow = 4)

  if (!is.null(filename)) {
    jpeg(filename = paste0(filename, "_mlLoganplot.jpeg"), width = 300, height = 400, units = "mm", res = 600)
    totalplot
    dev.off()
  }

  return(totalplot)
}
