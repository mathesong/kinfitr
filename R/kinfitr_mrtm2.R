#' Ichise's Multilinear Reference Tissue Model 2
#'
#' Function to fit MRTM2 of Ichise et al. (2003) to data.  This model is often used alongside MRTM1: MRTM1 is run on a
#' high-binding region to obtain an estimate of k2prime, and this k2prime value is used for MRTM2.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame. We include zero at
#' time zero: if not included, it is added.
#' @param roitac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param k2prime Value of k2prime to be used for the fitting, i.e. the average tissue-to-plasma clearance rate. This can be
#' obtained from MRTM1 of SRTM, or set at a specified value. If using SRTM to estimate this value, it is equal to k2 / R1.
#' @param tstarIncludedFrames Optional. The number of frames to be used in the multiple regression. This is a count from the end of the
#' measurement, so a value of 10 means that last 10 frames will be used. This value can be estimated using \code{mrtm2_tstar}. Note
#' that this tstar differs from that of the non-invasive Logan plot (which is the point at which pseudo-equilibrium is reached).
#' Rather, with MRTM1 and MRTM2, all frames can be used provided that the kinetics in the target tissue can be described by a one
#' tissue compartment model (like SRTM).  If this is not the case, a t* is required.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit}, a dataframe
#' containing the TACs of the data \code{out$tacs}, a dataframe containing the TACs of the fitted values \code{out$fitvals},
#' a vector of the weights \code{out$weights}, the specified k2prime value \code{out$k2prime}, and the specified
#' tstarIncludedFrames value \code{out$tstarIncludedFrames}.
#'
#' @examples
#' # Note: Reference region models, and irreversible binding models, should not
#' # be used for PBR28 - this is just to demonstrate function
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times/60
#' reftac <- pbr28$tacs[[2]]$CBL
#' roitac <- pbr28$tacs[[2]]$STR
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' fit <- mrtm2(t_tac, reftac, roitac, 0.001, weights=weights)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references @references Ichise M, Liow JS, Lu JQ, Takano A, Model K, Toyama H, Suhara T, Suzuki K, Innis RB, Carson RE. Linearized Reference Tissue Parametric Imaging Methods: Application to [11C]DASB Positron Emission Tomography Studies of the Serotonin Transporter in Human Brain. Journal of Cerebral Blood Flow & Metabolism. 2003 Sep 1;23(9):1096-112.
#'
#' @export

mrtm2 <- function(t_tac, reftac, roitac, k2prime, tstarIncludedFrames=NULL, weights=NULL, frameStartEnd=NULL) {


  # Tidying

  tidyinput <- tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd)

  t_tac <- tidyinput$t_tac
  reftac <- tidyinput$reftac
  roitac <- tidyinput$roitac
  weights <- tidyinput$weights

  if (is.null(tstarIncludedFrames)) {
    tstarIncludedFrames <- length(reftac)
  }


  # Parameters

  term1 <- pracma::cumtrapz(t_tac, reftac) + (1 / k2prime) * reftac
  term2 <- pracma::cumtrapz(t_tac, roitac)

  fitvals <- data.frame(
    Time = t_tac, Reference = reftac, Target = roitac, Weights = weights,
    Term1 = term1, Term2 = term2
  )

  if (is.null(tstarIncludedFrames) != T) {
    equil <- rep("Before", length(roitac))
    equil[(length(equil) - tstarIncludedFrames + 1):length(equil)] <- "After"
  } else {
    equil <- rep("After", length(roitac))
  }

  fitvals$Equilibrium <- equil

  fitvals_equil <- subset(fitvals, fitvals$Equilibrium == "After")

  # Solution

  mrtm2_model <- lm(Target ~ Term1 + Term2 - 1, weights = Weights, data = fitvals_equil)

  # Output

  ic <- 0 # is there an intercept?

  bp <- as.numeric(-((coef(mrtm2_model)[1 + ic] / coef(mrtm2_model)[2 + ic]) + 1))

  par <- as.data.frame(list(bp = bp))
  fit <- mrtm2_model

  tacs <- data.frame(Time = t_tac, Reference = reftac, Target = roitac)

  fitvals <- fitvals_equil
  fitvals$Target_fitted <- as.numeric(predict(mrtm2_model))

  out <- list(
    par = par, fit = fit, tacs = tacs, fitvals = fitvals, weights = weights,
    k2prime = k2prime, tstarIncludedFrames = tstarIncludedFrames, model = "mrtm2"
  )

  class(out) <- c("mrtm2", "kinfit")

  return(out)
}

#' Plot: MRTM2
#'
#' Function to visualise the fit of the MRTM2 model to data.
#'
#' @param mrtm2out The output object of the mrtm2 fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#' @param refname Optional. The name of the Reference Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times/60
#' reftac <- pbr28$tacs[[2]]$CBL
#' roitac <- pbr28$tacs[[2]]$STR
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' fit <- mrtm2(t_tac, reftac, roitac, 0.001, weights=weights)
#'
#' plot_mrtm2fit(fit)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_mrtm2fit <- function(mrtm2out, roiname = NULL, refname = NULL) {
  measured <- data.frame(
    Time = mrtm2out$tacs$Time,
    Reference = mrtm2out$tacs$Reference,
    ROI.measured = mrtm2out$tacs$Target,
    Weights = mrtm2out$weights
  )

  fitted <- data.frame(
    Time = mrtm2out$fitvals$Time,
    ROI.fitted = mrtm2out$fitvals$Target_fitted,
    Weights = mrtm2out$fitvals$Weights
  )

  if (is.null(roiname)) {
    roiname <- "ROI"
  }
  if (is.null(refname)) {
    refname <- "Reference"
  }

  measured <- plyr::rename(measured, c(
    "ROI.measured" = paste0(roiname, ".measured"),
    "Reference" = refname
  ))

  fitted <- plyr::rename(fitted, c("ROI.fitted" = paste0(roiname, ".fitted")))

  tidymeasured <- tidyr::gather(
    measured, key = Region, value = Radioactivity,
    -Time, -Weights, factor_key = F
  )

  tidyfitted <- tidyr::gather(
    fitted, key = Region, value = Radioactivity,
    -Time, -Weights, factor_key = F
  )

  Region <- forcats::fct_inorder(factor(c(tidymeasured$Region, tidyfitted$Region)))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  output <- ggplot(tidymeasured, aes(x = Time, y = Radioactivity, colour = Region)) +
    geom_point(data = tidymeasured, aes(shape = "a", size = Weights)) +
    geom_line(data = tidyfitted) + colScale +
    guides(shape = FALSE, color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

  return(output)
}

#' Tstar Finder: MRTM2
#'
#' Function to identify where t* is for MRTM2. Remember that this is unnecessary if the kinetics in the target tissue can be described by a one
#' tissue compartment model (like SRTM).
#'
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param k2prime Value of k2prime to be used for the fitting, i.e. the average tissue-to-plasma clearance rate. This can be
#' obtained from another model, such as MRTM1, SRTM or set at a specified value. If using SRTM to estimate this value, it is equal to k2 / R1.
#' @param filename The name of the output image: filename_mrtm1.jpeg
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_mrtm2.jpeg
#'
#' @examples
#' \dontrun{
#' mrtm2_tstar(t_tac, reftac, taclow, tacmed, tachigh, 'demonstration')
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

mrtm2_tstar <- function(t_tac, reftac, lowroi, medroi, highroi, k2prime, filename, frameStartEnd = NULL, gridbreaks = 2) {
  frames <- length(reftac)

  lowroi_fit <- mrtm2(t_tac, reftac, lowroi, k2prime = k2prime, frameStartEnd = frameStartEnd)
  medroi_fit <- mrtm2(t_tac, reftac, medroi, k2prime = k2prime, frameStartEnd = frameStartEnd)
  highroi_fit <- mrtm2(t_tac, reftac, highroi, k2prime = k2prime, frameStartEnd = frameStartEnd)

  low_linplot <- plot_mrtm2fit(lowroi_fit) + ggtitle("Low") + ylim(0, max(c(lowroi_fit$tacs$Reference, lowroi_fit$tacs$Target)) * 1.1) + theme(legend.position = "none")
  med_linplot <- plot_mrtm2fit(medroi_fit) + ggtitle("Medium") + ylim(0, max(c(medroi_fit$tacs$Reference, medroi_fit$tacs$Target)) * 1.1) + theme(legend.position = "none")
  high_linplot <- plot_mrtm2fit(highroi_fit) + ggtitle("High") + ylim(0, max(c(highroi_fit$tacs$Reference, highroi_fit$tacs$Target)) * 1.1) + theme(legend.position = "none")

  tstarInclFrames <- 3:frames
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  bp_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)

  for (i in 1:length(tstarInclFrames)) {
    lowfit <- mrtm2(t_tac, reftac, lowroi, k2prime = k2prime, tstarIncludedFrames = tstarInclFrames[i], frameStartEnd = frameStartEnd)
    medfit <- mrtm2(t_tac, reftac, medroi, k2prime = k2prime, tstarIncludedFrames = tstarInclFrames[i], frameStartEnd = frameStartEnd)
    highfit <- mrtm2(t_tac, reftac, highroi, k2prime = k2prime, tstarIncludedFrames = tstarInclFrames[i], frameStartEnd = frameStartEnd)

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

  jpeg(filename = paste0(filename, "_mrtm2.jpeg"), width = 300, height = 400, units = "mm", res = 600)
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
