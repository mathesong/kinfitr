#' Logan Plot
#'
#' Function to fit the Logan Plot to data.
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
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit},
#' a dataframe containing the TACs of the data \code{out$tacs}, a dataframe containing the fitted values \code{out$fitvals},
#' the blood input data frame after time shifting \code{input}, a vector of the weights \code{out$weights},
#' the inpshift value used \code{inpshift}, the specified vB value \code{out$vB} and the specified
#' tstarIncludedFrames value \code{out$tstarIncludedFrames}.
#'
#' @examples
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times/60
#' tac <- pbr28$tacs[[2]]$FC
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time/60 , pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time /60 , pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1 )
#'
#' fit1 <- Loganplot(t_tac, tac, input, 10, weights)
#' fit2 <- Loganplot(t_tac, tac, input, 10, weights, inpshift = 0.1, vB=0.05)
#'
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Logan J, Fowler JS, Volkow ND, Wolf AP, Dewey SL, Schlyer DJ, MacGregor RR, Hitzemann R, Bendriem B, Gatley SJ, Christman DR. Graphical Analysis of Reversible Radioligand Binding from Time-Activity Measurements Applied to [N-11C-Methyl]-(-)-Cocaine PET Studies in Human Subjects. Journal of Cerebral Blood Flow & Metabolism. 1990 Sep 1;10(5):740-7.
#'
#' @export

Loganplot <- function(t_tac, tac, input, tstarIncludedFrames, weights=NULL,
                      inpshift = 0, vB = 0, frameStartEnd=NULL) {

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
  blood <- newvals$input$Blood
  aif <- newvals$input$AIF

  # Parameters

  interptime <- newvals$input$Time
  i_tac <- pracma::interp1(t_tac, tac, interptime, method = "linear")

  # Blood Volume Correction (nothing happens if vB = 0)
  i_tac <- (i_tac - vB * blood) / (1 - vB)

  logan_roi <- as.numeric(pracma::cumtrapz(interptime, i_tac) / i_tac)
  logan_plasma <- as.numeric((pracma::cumtrapz(interptime, aif)) / i_tac)

  logan_roi <- pracma::interp1(interptime, logan_roi, t_tac, method = "linear")
  logan_plasma <- pracma::interp1(interptime, logan_plasma, t_tac, method = "linear")

  logan_equil_roi <- tail(logan_roi, tstarIncludedFrames)
  logan_equil_plasma <- tail(logan_plasma, tstarIncludedFrames)
  weights_equil <- tail(weights, tstarIncludedFrames)


  # Solution

  logan_model <- lm(logan_equil_roi ~ logan_equil_plasma, weights = weights_equil)

  # Output

  par <- as.data.frame(list(Vt = as.numeric(logan_model$coefficients[2])))
  fit <- logan_model

  tacs <- data.frame(Time = t_tac, Target = tac)

  fitvals <- data.frame(Logan_Plasma = logan_plasma, Logan_ROI = logan_roi)

  input <- newvals$input

  par.se <- as.data.frame(as.list(sqrt(abs(vcov(logan_model)[, 1]))))[-1]
  names(par.se) <- paste0("Vt", ".se")

  out <- list(
    par = par, par.se = par.se, fit = logan_model, tacs = tacs,
    fitvals = fitvals, input = input, weights = weights,
    inpshift = inpshift, vB = vB, tstarIncludedFrames = tstarIncludedFrames,
    model = "Logan"
  )

  class(out) <- c("Logan", "kinfit")

  return(out)
}


#' Plot: Logan Plot
#'
#' Function to visualise the fit of the Logan Plot model to data.
#'
#' @param loganout The output object of the Logan Plot fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times/60
#' tac <- pbr28$tacs[[2]]$FC
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time/60 , pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time /60 , pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1 )
#'
#' fit <- Loganplot(t_tac, tac, input, 10, weights)
#' plot_Loganfit(fit)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export


plot_Loganfit <- function(loganout, roiname = NULL) {
  plotdf <- data.frame(
    Weights = loganout$weights,
    Logan_Plasma = loganout$fitvals$Logan_Plasma,
    Logan_ROI = loganout$fitvals$Logan_ROI,
    Equilibrium = as.character("Before")
  )

  plotdf$Equilibrium <- as.character(plotdf$Equilibrium)
  plotdf$Equilibrium [ (nrow(plotdf) - (loganout$tstarIncludedFrames - 1)):nrow(plotdf)  ] <- "After"

  plotdf$Equilibrium <- forcats::fct_inorder(factor(plotdf$Equilibrium))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Equilibrium)
  colScale <- scale_colour_manual(name = paste0(roiname, "\nEquilibrium"), values = myColors)

  xlimits <- c(0, tail(plotdf$Logan_Plasma, 1))

  xlabel <- "Integ(C_Plasma)/C_Tissue"
  ylabel <- "Integ(C_Tissue)/C_Tissue"

  outplot <- ggplot(data = plotdf, aes(x = Logan_Plasma, y = Logan_ROI, colour = Equilibrium)) +
    geom_point(aes(shape = "a", size = Weights)) +
    geom_abline(
      slope = as.numeric(loganout$fit$coefficients[2]),
      intercept = as.numeric(loganout$fit$coefficients[1])
    ) +
    xlab(xlabel) + ylab(ylabel) + xlim(xlimits) + colScale +
    guides(shape = FALSE, color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

  return(outplot)
}


#' Tstar Finder: Logan Plot
#'
#' Function to identify where t* is for the Logan Plot.
#'
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param lowroi Numeric vector of radioactivity concentrations in a target
#'   tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target
#'   tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target
#'   tissue for each frame. This should be from a ROI with high binding.
#' @param input Data frame containing the blood, plasma, and parent fraction
#'   concentrations over time.  This can be generated using the
#'   \code{blood_interp} function.
#' @param filename The name of the output image: filename_Logan.jpeg
#' @param inpshift Optional. The number of minutes by which to shift the timing
#'   of the input data frame forwards or backwards. If not specified, this will
#'   be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will
#'   be ignored and assumed to be 0%. If specified, it will be corrected for
#'   prior to parameter estimation using the following equation: \deqn{C_{T}(t)
#'   = \frac{C_{Measured}(t) - vB\times C_{B}(t)}{1-vB}}
#' @param frameStartEnd Optional: This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This is to assess time
#'   stability.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_Logan.jpeg
#'
#' @examples
#' \dontrun{
#' Logan_tstar(t_tac, lowroi, medroi, highroi, input, filename='demonstration',
#'             inpshift = onetcmout$par$inpshift, vB = 0.05, gridbreaks=4)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

Logan_tstar <- function(t_tac, lowroi, medroi, highroi, input, filename = NULL, inpshift = 0, vB = 0.05, frameStartEnd = NULL, gridbreaks=2) {

  frames <- length(t_tac)
  lowroi_fit <- Loganplot(t_tac, lowroi, input, tstarIncludedFrames = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  medroi_fit <- Loganplot(t_tac, medroi, input, tstarIncludedFrames = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  highroi_fit <- Loganplot(t_tac, highroi, input, tstarIncludedFrames = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)

  xlabel <- "Integ(C_Plasma)/C_Tissue"
  ylabel <- "Integ(C_Tissue)/C_Tissue"

  low_xlimits <- c(0, tail(lowroi_fit$fitvals$Logan_Plasma, 1))
  med_xlimits <- c(0, tail(medroi_fit$fitvals$Logan_Plasma, 1))
  high_xlimits <- c(0, tail(highroi_fit$fitvals$Logan_Plasma, 1))

  low_linplot <- qplot(lowroi_fit$fitvals$Logan_Plasma, lowroi_fit$fitvals$Logan_ROI) + ggtitle("Low") + xlab(xlabel) + ylab(ylabel) + xlim(low_xlimits)
  med_linplot <- qplot(medroi_fit$fitvals$Logan_Plasma, medroi_fit$fitvals$Logan_ROI) + ggtitle("Medium") + xlab(xlabel) + ylab(ylabel) + xlim(med_xlimits)
  high_linplot <- qplot(highroi_fit$fitvals$Logan_Plasma, highroi_fit$fitvals$Logan_ROI) + ggtitle("High") + xlab(xlabel) + ylab(ylabel) + xlim(high_xlimits)

  tstarInclFrames <- 3:frames
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  vt_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)

  for (i in 1:length(tstarInclFrames)) {
    lowfit <- Loganplot(t_tac, lowroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
    medfit <- Loganplot(t_tac, medroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
    highfit <- Loganplot(t_tac, highroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)

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
  ylab_r2 <- expression(R ^ 2)
  ylab_mp <- "Maximum Percentage Variance"


  # R Squared plots

  low_r2plot <- ggplot(r2_df, aes(x = Frames, y = Low)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim=c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)
  med_r2plot <- ggplot(r2_df, aes(x = Frames, y = Medium)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim=c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)
  high_r2plot <- ggplot(r2_df, aes(x = Frames, y = High)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim=c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)


  # Max Percentage Variation Plots

  maxperc_df$inclmins <- rev(max(t_tac) - t_tac)[-c(1, 2)]
  maxperc_df$tstar <- rev(t_tac)[-c(1, 2)]

  low_mpplot <- ggplot(maxperc_df, aes(x = Frames, y = Low)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim=c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text", x = 3, y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Low + 1.4, label = round(maxperc_df$tstar, 1), size = 3, colour = "red") + annotate("text", x = 3, y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Low + 0.7, label = round(maxperc_df$inclmins, 1), size = 3, colour = "blue")
  med_mpplot <- ggplot(maxperc_df, aes(x = Frames, y = Medium)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim=c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text", x = 3, y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Medium + 1.4, label = round(maxperc_df$tstar, 1), size = 3, colour = "red") + annotate("text", x = 3, y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$Medium + 0.7, label = round(maxperc_df$inclmins, 1), size = 3, colour = "blue")
  high_mpplot <- ggplot(maxperc_df, aes(x = Frames, y = High)) + geom_point() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + coord_cartesian(ylim=c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text", x = 3, y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$High + 1.4, label = round(maxperc_df$tstar, 1), size = 3, colour = "red") + annotate("text", x = 3, y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0) + annotate("text", x = maxperc_df$Frames, y = maxperc_df$High + 0.7, label = round(maxperc_df$inclmins, 1), size = 3, colour = "blue")



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

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow=1)
  r2row <- cowplot::plot_grid(low_r2plot, med_r2plot, high_r2plot, nrow=1)
  mprow <- cowplot::plot_grid(low_mpplot, med_mpplot, high_mpplot, nrow=1)
  outrow <- cowplot::plot_grid(tacplot, vtplot, rel_widths = c(2,1))

  totalplot <- cowplot::plot_grid(linrow, r2row, mprow, outrow, nrow = 4)

  if(!is.null(filename)) {
    jpeg(filename = paste0(filename, "_Loganplot.jpeg"), width = 300, height = 400, units = "mm", res = 600)
    totalplot
    dev.off()
  }

  return(totalplot)
}
