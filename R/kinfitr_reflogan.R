#' Non-Invasive Logan Plot
#'
#' Function to fit the non-invasive Logan plot model of Logan et al. (1996) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame. We include zero at
#' time zero: if not included, it is added.
#' @param roitac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param k2prime Value of k2prime to be used for the fitting, i.e. the average tissue-to-plasma clearance rate. This can be
#' obtained from another model, or set at a specified value. If using SRTM to estimate this value, it is equal to k2 / R1.
#' @param tstarIncludedFrames The number of frames to be used in the regression model, i.e. the number of frames for which
#' the function is linear after pseudo-equilibrium is reached. This is a count from the end of the measurement, so a value of
#' 10 means that last 10 frames will be used. This value can be estimated using \code{refLogan_tstar}.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.

#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit}, a dataframe
#' containing the TACs of the data \code{out$tacs}, a dataframe containing the TACs of the fitted values \code{out$fitvals},
#' a vector of the weights \code{out$weights}, the specified k2prime value \code{out$k2prime}, and the specified
#' tstarIncludedFrames value \code{out$tstarIncludedFrames}
#'
#' @examples
#' refLogan(t_tac, reftac, roitac, k2prime=0.1, tstarIncludedFrames=10, weights=weights)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Logan J, Fowler JS, Volkow ND, Wang GJ, Ding YS, Alexoff DL. Distribution volume ratios without blood sampling from graphical analysis of PET data. Journal of Cerebral Blood Flow & Metabolism. 1996 Sep 1;16(5):834-40.
#'
#' @export

refLogan <- function(t_tac, reftac, roitac, k2prime, tstarIncludedFrames, weights, frameStartEnd) {


  # Tidying

  tidyinput <- tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd)

  t_tac <- tidyinput$t_tac
  reftac <- tidyinput$reftac
  roitac <- tidyinput$roitac
  weights <- tidyinput$weights


  # Parameters

  logan_roi <- pracma::cumtrapz(t_tac, roitac) / roitac
  logan_ref <- (pracma::cumtrapz(t_tac, reftac) + reftac / k2prime) / roitac

  logan_equil_roi <- tail(logan_roi, tstarIncludedFrames)
  logan_equil_ref <- tail(logan_ref, tstarIncludedFrames)
  weights_equil <- tail(weights, tstarIncludedFrames)


  # Solution

  logan_model <- lm(logan_equil_roi ~ logan_equil_ref, weights = weights_equil)

  # Output

  par <- as.data.frame(list(bp = as.numeric(logan_model$coefficients[2]) - 1))
  fit <- logan_model

  tacs <- data.frame(Time = t_tac, Reference = reftac, Target = roitac)

  fitvals <- data.frame(Logan_ROI = logan_roi, Logan_Ref = logan_ref)

  out <- list(
    par = par, fit = fit, tacs = tacs, fitvals = fitvals,
    weights = weights, k2prime = k2prime,
    tstarIncludedFrames = tstarIncludedFrames, model = "refLogan"
  )

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
#' plot_refLoganfit(refloganoutout)
#'
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

  plotdf$Equilibrium <- forcats::fct_inorder(factor(plotdf$Equilibrium))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Equilibrium)
  colScale <- scale_colour_manual(name = paste0(roiname, "\nEquilibrium"), values = myColors)

  outplot <- ggplot(data = plotdf, aes(x = Logan_ref, y = Logan_roi, colour = Equilibrium)) +
    geom_point(aes(shape = "a", size = Weights)) +
    geom_abline(
      slope = as.numeric(refloganout$fit$coefficients[2]),
      intercept = as.numeric(refloganout$fit$coefficients[1])
    ) +
    xlab("[Integ(C_Ref)+C_Ref/k2prime] / C_Tissue") + ylab("Integ(C_Tissue) / C_Tissue") + colScale +
    guides(shape = FALSE, color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

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
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_refLogan.jpeg
#'
#' @examples
#' refLogan_tstar(t_tac, reftac, taclow, tacmed, tachigh, k2prime = k2prime, 'demonstration')
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

refLogan_tstar <- function(t_tac, reftac, lowroi, medroi, highroi, k2prime, filename, gridbreaks = 2) {
  frames <- length(reftac)
  lowroi_fit <- refLogan(t_tac, reftac, lowroi, k2prime, length(reftac))
  medroi_fit <- refLogan(t_tac, reftac, medroi, k2prime, length(reftac))
  highroi_fit <- refLogan(t_tac, reftac, highroi, k2prime, length(reftac))

  logan_xlab <- "[Integ(C_Ref)+C_Ref/k2prime] / C_Tissue"
  logan_ylab <- "Integ(C_Tissue) / C_Tissue"

  low_linplot <- qplot(lowroi_fit$fitvals$Logan_Ref, lowroi_fit$fitvals$Logan_ROI) + ggtitle("Low") + xlab(logan_xlab) + ylab(logan_ylab)
  med_linplot <- qplot(medroi_fit$fitvals$Logan_Ref, medroi_fit$fitvals$Logan_ROI) + ggtitle("Medium") + xlab(logan_xlab) + ylab(logan_ylab)
  high_linplot <- qplot(highroi_fit$fitvals$Logan_Ref, highroi_fit$fitvals$Logan_ROI) + ggtitle("High") + xlab(logan_xlab) + ylab(logan_ylab)

  tstarInclFrames <- 3:frames
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  bp_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)

  for (i in 1:length(tstarInclFrames)) {
    lowfit <- refLogan(t_tac, reftac, lowroi, k2prime = k2prime, tstarIncludedFrames = tstarInclFrames[i])
    medfit <- refLogan(t_tac, reftac, medroi, k2prime = k2prime, tstarIncludedFrames = tstarInclFrames[i])
    highfit <- refLogan(t_tac, reftac, highroi, k2prime = k2prime, tstarIncludedFrames = tstarInclFrames[i])

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

  bpplot <- ggplot(bpplotdf, aes(x = Frames, y = BP, colour = Region)) + geom_point() + geom_line() + scale_x_continuous(breaks = seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylab(expression(BP[ND])) + colScale


  # File Output

  jpeg(filename = paste0(filename, "_refLogan.jpeg"), width = 300, height = 400, units = "mm", res = 600)
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
  print(bpplot, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 3))
  dev.off()
}
