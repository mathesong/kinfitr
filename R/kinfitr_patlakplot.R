#' Patlak Plot
#'
#' Function to fit the Patlak Plot of of Patlak et al (1983) to data.
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
#' Patlakplot(t_tac, tac, input, 10, weights, inpshift = onetcmout$par$inpshift)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Patlak CS, Blasberg RG, Fenstermacher JD. Graphical evaluation of blood-to-brain transfer constants from multiple-time uptake data. Journal of Cerebral Blood Flow & Metabolism. 1983 Mar 1;3(1):1-7.
#'
#' @export

Patlakplot <- function(t_tac, tac, input, tstarIncludedFrames, weights, inpshift = 0, vB = 0, frameStartEnd) {


  # Tidying

  tidyinput <- tidyinput_art(t_tac, tac, weights, frameStartEnd)

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
  blood <- newvals$input$blood
  plasma <- newvals$input$plasma
  parentfrac <- newvals$input$parentfrac

  # Parameters

  corrplasma <- newvals$input$plasma * newvals$input$parentfrac

  interptime <- newvals$input$Time
  i_tac <- pracma::interp1(t_tac, tac, interptime, method = "linear")

  # Blood Volume Correction (nothing happens if vB = 0)
  i_tac <- (i_tac - vB * blood) / (1 - vB)

  patlak_roi <- i_tac / corrplasma
  patlak_plasma <- as.numeric((pracma::cumtrapz(interptime, corrplasma)) / corrplasma)

  patlak_roi <- pracma::interp1(interptime, patlak_roi, t_tac, method = "linear")
  patlak_plasma <- pracma::interp1(interptime, patlak_plasma, t_tac, method = "linear")

  patlak_equil_roi <- tail(patlak_roi, tstarIncludedFrames)
  patlak_equil_plasma <- tail(patlak_plasma, tstarIncludedFrames)
  weights_equil <- tail(weights, tstarIncludedFrames)


  # Solution

  patlak_model <- lm(patlak_equil_roi ~ patlak_equil_plasma, weights = weights_equil)

  # Output

  par <- as.data.frame(list(K = as.numeric(patlak_model$coefficients[2])))
  fit <- patlak_model

  tacs <- data.frame(Time = t_tac, Target = tac)

  fitvals <- data.frame(Patlak_Plasma = patlak_plasma, Patlak_ROI = patlak_roi)

  input <- newvals$input

  par.se <- as.data.frame(as.list(sqrt(abs(vcov(patlak_model)[, 1]))))[-1]
  names(par.se) <- paste0("K", ".se")

  out <- list(
    par = par, par.se = par.se, fit = patlak_model, tacs = tacs,
    fitvals = fitvals, input = input, weights = weights,
    inpshift = inpshift, vB = vB, tstarIncludedFrames = tstarIncludedFrames
  )
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
#' plot_Patlakfit(patlakout)
#'
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

  if (is.null(roiname)) {
    roiname <- "ROI"
  }

  if (roiname != "ROI") {
    plotdf <- plyr::rename(plotdf, c("ROI_measured" = roiname))
  }

  plotdf$Equilibrium <- as.character(plotdf$Equilibrium)
  plotdf$Equilibrium [ (nrow(plotdf) - (patlakout$tstarIncludedFrames - 1)):nrow(plotdf)  ] <- "After"

  plotdf$Equilibrium <- forcats::fct_inorder(factor(plotdf$Equilibrium))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Equilibrium)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  xlimits <- c(0, tail(plotdf$Patlak_Plasma, 1))

  xlabel <- "Integ(C_Plasma)/C_Tissue"
  ylabel <- "Integ(C_Tissue)/C_Tissue"

  outplot <- ggplot(data = plotdf, aes(x = Patlak_Plasma, y = Patlak_ROI, colour = Equilibrium)) +
    geom_point(aes(shape = "a", size = Weights)) +
    geom_abline(
      slope = as.numeric(patlakout$fit$coefficients[2]),
      intercept = as.numeric(patlakout$fit$coefficients[1])
    ) +
    xlab(xlabel) + ylab(ylabel) + xlim(xlimits) + colScale +
    guides(shape = FALSE, color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

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
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_Patlakplot.jpeg
#'
#' @examples
#' Patlak_tstar(t_tac, lowroi, medroi, highroi, input, filename='demonstration', inpshift = onetcmout$par$inpshift, vB = 0.05, gridbreaks=4)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

Patlak_tstar <- function(t_tac, lowroi, medroi, highroi, input, filename, inpshift = 0, frameStartEnd, gridbreaks=2) {
  frames <- length(t_tac)
  lowroi_fit <- Patlakplot(t_tac, lowroi, input, tstarIncludedFrames = frames, inpshift = inpshift, frameStartEnd = frameStartEnd)
  medroi_fit <- Patlakplot(t_tac, medroi, input, tstarIncludedFrames = frames, inpshift = inpshift, frameStartEnd = frameStartEnd)
  highroi_fit <- Patlakplot(t_tac, highroi, input, tstarIncludedFrames = frames, inpshift = inpshift, frameStartEnd = frameStartEnd)

  xlabel <- "Integ(C_Plasma) / C_Plasma"
  ylabel <- "C_Tissue / C_Plasma"

  low_xlimits <- c(0, tail(lowroi_fit$fitvals$Patlak_Plasma, 1))
  med_xlimits <- c(0, tail(medroi_fit$fitvals$Patlak_Plasma, 1))
  high_xlimits <- c(0, tail(highroi_fit$fitvals$Patlak_Plasma, 1))

  low_linplot <- qplot(lowroi_fit$fitvals$Patlak_Plasma, lowroi_fit$fitvals$Patlak_ROI) + ggtitle("Low") + xlab(xlabel) + ylab(ylabel) + xlim(low_xlimits)
  med_linplot <- qplot(medroi_fit$fitvals$Patlak_Plasma, medroi_fit$fitvals$Patlak_ROI) + ggtitle("Medium") + xlab(xlabel) + ylab(ylabel) + xlim(med_xlimits)
  high_linplot <- qplot(highroi_fit$fitvals$Patlak_Plasma, highroi_fit$fitvals$Patlak_ROI) + ggtitle("High") + xlab(xlabel) + ylab(ylabel) + xlim(high_xlimits)

  tstarInclFrames <- 3:frames
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  K_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)

  for (i in 1:length(tstarInclFrames)) {
    lowfit <- Patlakplot(t_tac, lowroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, frameStartEnd = frameStartEnd)
    medfit <- Patlakplot(t_tac, medroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, frameStartEnd = frameStartEnd)
    highfit <- Patlakplot(t_tac, highroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, frameStartEnd = frameStartEnd)

    r2_df$Low[i] <- summary(lowfit$fit)$r.squared
    r2_df$Medium[i] <- summary(medfit$fit)$r.squared
    r2_df$High[i] <- summary(highfit$fit)$r.squared

    maxperc_df$Low[i] <- maxpercres(lowfit)
    maxperc_df$Medium[i] <- maxpercres(medfit)
    maxperc_df$High[i] <- maxpercres(highfit)

    K_df$Low[i] <- lowfit$par$K
    K_df$Medium[i] <- medfit$par$K
    K_df$High[i] <- highfit$par$K
  }

  xlabel <- "Number of Included Frames"
  ylab_r2 <- expression(R ^ 2)
  ylab_mp <- "Maximum Percentage Variance"


  # R Squared plots

  low_r2plot <- ggplot(r2_df, aes(x = Frames, y = Low)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)
  med_r2plot <- ggplot(r2_df, aes(x = Frames, y = Medium)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)
  high_r2plot <- ggplot(r2_df, aes(x = Frames, y = High)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)


  # Max Percentage Variation Plots

  maxperc_df$inclmins <- rev(max(t_tac) - t_tac)[-c(1, 2)]
  maxperc_df$tstar <- rev(t_tac)[-c(1, 2)]

  low_mpplot <- ggplot(maxperc_df, aes(x = Frames, y = Low)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text", x = 3, y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Low + 1.4, label = round(maxperc_df$tstar, 1), size = 3, colour = "red") + annotate("text", x = 3, y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Low + 0.7, label = round(maxperc_df$inclmins, 1), size = 3, colour = "blue")
  med_mpplot <- ggplot(maxperc_df, aes(x = Frames, y = Medium)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text", x = 3, y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Medium + 1.4, label = round(maxperc_df$tstar, 1), size = 3, colour = "red") + annotate("text", x = 3, y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Medium + 0.7, label = round(maxperc_df$inclmins, 1), size = 3, colour = "blue")
  high_mpplot <- ggplot(maxperc_df, aes(x = Frames, y = High)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text", x = 3, y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$High + 1.4, label = round(maxperc_df$tstar, 1), size = 3, colour = "red") + annotate("text", x = 3, y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$High + 0.7, label = round(maxperc_df$inclmins, 1), size = 3, colour = "blue")



  # TAC Plot

  tacplotdf <- data.frame(cbind(Time = lowroi_fit$tacs$Time, Low = lowroi_fit$tacs$Target, Medium = medroi_fit$tacs$Target, High = highroi_fit$tacs$Target))
  tacplotdf <- tidyr::gather(tacplotdf, key = Region, value = Radioactivity, -Time)

  tacplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(tacplotdf$Region)))

  myColors <- RColorBrewer::brewer.pal(4, "Set1")
  names(myColors) <- levels(tacplotdf$Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  tacplot <- ggplot(tacplotdf, aes(x = Time, y = Radioactivity, colour = Region)) + geom_point() + geom_line() + colScale


  # K Plot

  Kplotdf <- tidyr::gather(K_df, key = Region, value = K, -Frames, -Time)
  Kplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(Kplotdf$Region)))

  Kplot <- ggplot(Kplotdf, aes(x = Frames, y = K, colour = Region)) + geom_point() + geom_line() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylab(expression(K[i])) + colScale


  # File Output

  jpeg(filename = paste0(filename, "_Patlakplot.jpeg"), width = 300, height = 400, units = "mm", res = 600)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(4, 3, heights = unit(c(1, 1), "null"))))
  print(low_linplot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(med_linplot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(high_linplot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  print(low_r2plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(med_r2plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(high_r2plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 3))
  print(low_mpplot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  print(med_mpplot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
  print(high_mpplot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 3))
  print(tacplot, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2))
  print(Kplot, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 3))
  dev.off()
}
