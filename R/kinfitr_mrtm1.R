#' Ichise's Multilinear Reference Tissue Model
#'
#' Function to fit MRTM1 of Ichise et al. (2003) to data.  This model is often used alongside MRTM2: MRTM1 is run on a
#' high-binding region to obtain an estimate of k2prime, and this k2prime value is used for MRTM2.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame. We include zero at
#' time zero: if not included, it is added.
#' @param roitac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param tstarIncludedFrames Optional. The number of frames to be used in the multiple regression. This is a count from the end of the
#' measurement, so a value of 10 means that last 10 frames will be used. This value can be estimated using \code{mrtm1_tstar}. Note
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
#' a vector of the weights \code{out$weights}, and the specified tstarIncludedFrames value \code{out$tstarIncludedFrames}
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
#' fit <- mrtm1(t_tac, reftac, roitac, weights=weights)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Ichise M, Liow JS, Lu JQ, Takano A, Model K, Toyama H, Suhara T, Suzuki K, Innis RB, Carson RE. Linearized Reference Tissue Parametric Imaging Methods: Application to [11C]DASB Positron Emission Tomography Studies of the Serotonin Transporter in Human Brain. Journal of Cerebral Blood Flow & Metabolism. 2003 Sep 1;23(9):1096-112.
#'
#' @export

mrtm1 <- function(t_tac, reftac, roitac, weights=NULL, tstarIncludedFrames=NULL, frameStartEnd=NULL) {

  # Tidying

  tidyinput <- tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd)

  t_tac <- tidyinput$t_tac
  reftac <- tidyinput$reftac
  roitac <- tidyinput$roitac
  weights <- tidyinput$weights

  if (is.null(tstarIncludedFrames)) {
    tstarIncludedFrames <- length(reftac)
  }

  if(tstarIncludedFrames > length(reftac)) {
    warning("tstarIncludedFrames is greater than the number of frames.
            Setting it to include all the frames.")
    tstarIncludedFrames <- length(reftac)
  }


  # Parameters

  term1 <- pracma::cumtrapz(t_tac, reftac)
  term2 <- pracma::cumtrapz(t_tac, roitac)
  term3 <- reftac

  fitvals <- data.frame(
    Time = t_tac, Reference = reftac, Target = roitac, Weights = weights,
    Term1 = term1, Term2 = term2, Term3 = term3
  )

  if (!is.null(tstarIncludedFrames)) {
    equil <- rep("Before", length(roitac))
    equil[(length(equil) - tstarIncludedFrames + 1):length(equil)] <- "After"
  } else {
    equil <- rep("After", length(roitac))
  }

  fitvals$Equilibrium <- equil

  fitvals_equil <- subset(fitvals, fitvals$Equilibrium == "After")

  # Solution

  mrtm1_model <- lm(Target ~ Term1 + Term2 + Term3 - 1, weights = Weights, data = fitvals_equil)

  # Output

  ic <- 0 # is there an intercept?

  bp <- as.numeric(-((coef(mrtm1_model)[1 + ic] / coef(mrtm1_model)[2 + ic]) + 1))
  k2prime <- as.numeric(coef(mrtm1_model)[1 + ic] / coef(mrtm1_model)[3 + ic])

  par <- as.data.frame(list(bp = bp, k2prime = k2prime))
  fit <- mrtm1_model

  tacs <- data.frame(Time = t_tac, Reference = reftac, Target = roitac)

  fitvals <- fitvals_equil
  fitvals$Target_fitted <- as.numeric(predict(mrtm1_model))

  out <- list(
    par = par, fit = fit, tacs = tacs, fitvals = fitvals, weights = weights,
    tstarIncludedFrames = tstarIncludedFrames, model = "mrtm1"
  )

  class(out) <- c("mrtm1", "kinfit")

  return(out)
}

#' Plot: MRTM1
#'
#' Function to visualise the fit of the MRTM1 model to data.
#'
#' @param mrtm1out The output object of the mrtm1 fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#' @param refname Optional. The name of the Reference Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
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
#' fit <- mrtm1(t_tac, reftac, roitac, weights=weights)
#'
#' plot_mrtm1fit(fit)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_mrtm1fit <- function(mrtm1out, roiname = NULL, refname = NULL) {
  measured <- data.frame(
    Time = mrtm1out$tacs$Time,
    Reference = mrtm1out$tacs$Reference,
    ROI.measured = mrtm1out$tacs$Target,
    Weights = mrtm1out$weights
  )

  fitted <- data.frame(
    Time = mrtm1out$fitvals$Time,
    ROI.fitted = mrtm1out$fitvals$Target_fitted,
    Weights = mrtm1out$fitvals$Weights
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

  tidymeasured$Region[tidymeasured$Region == "ROI.measured"] <- paste0(roiname, ".measured")
  tidymeasured$Region[tidymeasured$Region == "Reference"] <- paste0(refname)
  tidyfitted$Region[tidyfitted$Region == "ROI.fitted"] <- paste0(roiname, ".fitted")

  Region <- forcats::fct_inorder(factor(c(tidymeasured$Region, tidyfitted$Region)))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  outplot <- ggplot(tidymeasured, aes(x = Time, y = Radioactivity, colour = Region)) +
    geom_point(data = tidymeasured, aes(shape = "a", size = Weights)) +
    geom_line(data = tidyfitted) +
    guides(shape = FALSE, color = guide_legend(order = 1)) + colScale +
    scale_size(range = c(1, 3))

  # print(outplot)
  return(outplot)
}

#' Tstar Finder: MRTM1
#'
#' Function to identify where t* is for MRTM1. Remember that this is unnecessary if the kinetics in the target tissue can be described by a one
#' tissue compartment model (like SRTM).
#'
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param filename The name of the output image: filename_mrtm1.jpeg
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_mrtm1.jpeg
#'
#' @examples
#' \dontrun{
#' mrtm1_tstar(t_tac, reftac, taclow, tacmed, tachigh, 'demonstration')
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

mrtm1_tstar <- function(t_tac, reftac, lowroi, medroi, highroi, filename=NULL, frameStartEnd=NULL, gridbreaks = 2) {
  frames <- length(reftac)

  lowroi_fit <- mrtm1(t_tac, reftac, lowroi, frameStartEnd = frameStartEnd)
  medroi_fit <- mrtm1(t_tac, reftac, medroi, frameStartEnd = frameStartEnd)
  highroi_fit <- mrtm1(t_tac, reftac, highroi, frameStartEnd = frameStartEnd)

  low_linplot <- plot_mrtm1fit(lowroi_fit) + ggtitle("Low") + ylim(0, max(c(lowroi_fit$tacs$Reference, lowroi_fit$tacs$Target)) * 1.1) + theme(legend.position = "none")
  med_linplot <- plot_mrtm1fit(medroi_fit) + ggtitle("Medium") + ylim(0, max(c(medroi_fit$tacs$Reference, medroi_fit$tacs$Target)) * 1.1) + theme(legend.position = "none")
  high_linplot <- plot_mrtm1fit(highroi_fit) + ggtitle("High") + ylim(0, max(c(highroi_fit$tacs$Reference, highroi_fit$tacs$Target)) * 1.1) + theme(legend.position = "none")

  tstarInclFrames <- 3:frames
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  bp_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)

  for (i in 1:length(tstarInclFrames)) {
    lowfit <- mrtm1(t_tac, reftac, lowroi, tstarIncludedFrames = tstarInclFrames[i], frameStartEnd = frameStartEnd)
    medfit <- mrtm1(t_tac, reftac, medroi, tstarIncludedFrames = tstarInclFrames[i], frameStartEnd = frameStartEnd)
    highfit <- mrtm1(t_tac, reftac, highroi, tstarIncludedFrames = tstarInclFrames[i], frameStartEnd = frameStartEnd)

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


  # Output

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow=1)
  r2row <- cowplot::plot_grid(low_r2plot, med_r2plot, high_r2plot, nrow=1)
  mprow <- cowplot::plot_grid(low_mpplot, med_mpplot, high_mpplot, nrow=1)
  outrow <- cowplot::plot_grid(tacplot, bpplot, rel_widths = c(2,1))

  totalplot <- cowplot::plot_grid(linrow, r2row, mprow, outrow, nrow = 4)

  if(!is.null(filename)) {
    jpeg(filename = paste0(filename, "_mrtm1.jpeg"), width = 300, height = 400, units = "mm", res = 600)
    totalplot
    dev.off()
  }

  return(totalplot)
}
