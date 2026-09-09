# Internal helpers shared by all of the `*_tstar` t* finding plot functions
# (mrtm1, mrtm2, refLogan, refmlLogan, refPatlak, Logan, mlLogan, ma1, Patlak).
# These collapse the large amount of plotting and looping boilerplate that used
# to be copy-pasted into every `*_tstar` function.


#' Resolve frameStartEnd from a timeStartEnd specification
#'
#' @param t_tac Time vector.
#' @param frameStartEnd Frame-based subset (returned unchanged if supplied).
#' @param timeStartEnd Time-based subset, converted to frames if frameStartEnd
#'   is not supplied.
#'
#' @return A frameStartEnd vector (or NULL).
#'
#' @noRd
tstar_frameStartEnd <- function(t_tac, frameStartEnd, timeStartEnd) {
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(
      which(t_tac >= timeStartEnd[1])[1],
      tail(which(t_tac <= timeStartEnd[2]), 1)
    )
  }
  frameStartEnd
}


#' Fit each ROI across all candidate t* values
#'
#' Runs the model for the low, medium and high binding ROIs across every value
#' of \code{tstarInclFrames}, collecting the R-squared, maximum percentage
#' residual and outcome measure for each.
#'
#' @param t_tac Time vector (used only to label the data frames).
#' @param lowroi,medroi,highroi The three ROI TACs.
#' @param fitfunc A function of \code{(roitac, tstar)} returning a model fit.
#' @param get_outcome A function of a fit returning the outcome measure (e.g.
#'   \code{function(f) f$par$Vt}).
#' @param tstarInclFrames The candidate numbers of included frames.
#'
#' @return A list with \code{r2_df}, \code{maxperc_df} and \code{outcome_df}.
#'
#' @noRd
tstar_compute <- function(t_tac, lowroi, medroi, highroi, fitfunc, get_outcome,
                          tstarInclFrames) {
  zeros <- rep(0, length(tstarInclFrames))

  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[tstarInclFrames], Low = zeros, Medium = zeros, High = zeros)
  outcome_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[tstarInclFrames], Low = zeros, Medium = zeros, High = zeros)

  for (i in seq_along(tstarInclFrames)) {
    lowfit <- fitfunc(lowroi, tstarInclFrames[i])
    medfit <- fitfunc(medroi, tstarInclFrames[i])
    highfit <- fitfunc(highroi, tstarInclFrames[i])

    r2_df$Low[i] <- summary(lowfit$fit)$r.squared
    r2_df$Medium[i] <- summary(medfit$fit)$r.squared
    r2_df$High[i] <- summary(highfit$fit)$r.squared

    maxperc_df$Low[i] <- maxpercres(lowfit)
    maxperc_df$Medium[i] <- maxpercres(medfit)
    maxperc_df$High[i] <- maxpercres(highfit)

    outcome_df$Low[i] <- get_outcome(lowfit)
    outcome_df$Medium[i] <- get_outcome(medfit)
    outcome_df$High[i] <- get_outcome(highfit)
  }

  list(r2_df = r2_df, maxperc_df = maxperc_df, outcome_df = outcome_df)
}


#' x-axis breaks for the frame-based t* plots
#'
#' Breaks start at 4 (the first even number) rather than at the first frame (3),
#' so that the gridlines fall on even numbers of included frames: 4, 6, 8 ...
#'
#' @noRd
tstar_framebreaks <- function(tstarInclFrames, gridbreaks) {
  scale_x_continuous(breaks = seq(4, max(tstarInclFrames), by = gridbreaks))
}


#' Build the three R-squared plots (low, medium, high binding)
#'
#' @noRd
tstar_r2plots <- function(r2_df, tstarInclFrames, gridbreaks) {
  xbreaks <- tstar_framebreaks(tstarInclFrames, gridbreaks)
  ybounds <- coord_cartesian(ylim = c(0.99, 1))
  xlabel <- "Number of Included Frames"
  ylab_r2 <- expression(R^2)

  one <- function(region) {
    ggplot(r2_df, aes(x = Frames, y = .data[[region]])) + geom_point() +
      xbreaks + ybounds + xlab(xlabel) + ylab(ylab_r2)
  }

  list(low = one("Low"), med = one("Medium"), high = one("High"))
}


#' Build the three maximum-percentage-deviation plots (low, medium, high)
#'
#' The t* (red) and included-minutes (blue) labels are drawn vertically in a
#' fixed band rather than directly above each point, so that they remain legible
#' and paired even when there are many frames. To keep each panel sparse, only
#' every \code{label_step}-th frame is labelled, and the three panels are
#' staggered (e.g. Low: 3, 6, 9 ...; Medium: 4, 7, 10 ...; High: 5, 8, 11 ...)
#' so that across the panels the time for every frame is shown exactly once.
#'
#' @param label_step Spacing between labelled frames within a panel. Defaults to
#'   3 (one panel per offset), so the three panels jointly label every frame.
#'
#' @noRd
tstar_maxpercplots <- function(maxperc_df, t_tac, tstarInclFrames, gridbreaks, label_step = 3) {
  maxperc_df$inclmins <- rev(max(t_tac) - t_tac)[-c(1, 2)]
  maxperc_df$tstar <- rev(t_tac)[-c(1, 2)]

  xbreaks <- tstar_framebreaks(tstarInclFrames, gridbreaks)
  ybounds <- coord_cartesian(ylim = c(0, 20))
  xlabel <- "Number of Included Frames"
  ylab_mp <- "Maximum Percentage Deviation"

  # Keys (top left) explaining the two coloured number rows
  keys <- list(
    annotate("text", x = min(tstarInclFrames), y = 20, label = "t* Minutes", colour = "red", size = 3, hjust = 0),
    annotate("text", x = min(tstarInclFrames), y = 20 - 0.7, label = "Included Minutes", colour = "blue", size = 3, hjust = 0)
  )

  # Per-frame numbers, vertical and in a fixed band, drawn only for the frames
  # requested for this panel so they never overlap
  make_labels <- function(label_frames) {
    df <- maxperc_df[maxperc_df$Frames %in% label_frames, , drop = FALSE]
    list(
      annotate("text", x = df$Frames, y = 14, label = round(df$tstar, 1), colour = "red", size = 2.5, angle = 90),
      annotate("text", x = df$Frames, y = 11, label = round(df$inclmins, 1), colour = "blue", size = 2.5, angle = 90)
    )
  }

  one <- function(region, label_frames) {
    ggplot(maxperc_df, aes(x = Frames, y = .data[[region]])) + geom_point() +
      xbreaks + ybounds + xlab(xlabel) + ylab(ylab_mp) + keys + make_labels(label_frames)
  }

  n <- length(tstarInclFrames)
  low_frames <- tstarInclFrames[seq(1, n, by = label_step)]
  med_frames <- tstarInclFrames[seq(2, n, by = label_step)]
  high_frames <- tstarInclFrames[seq(3, n, by = label_step)]

  list(low = one("Low", low_frames),
       med = one("Medium", med_frames),
       high = one("High", high_frames))
}


#' Build the bottom output row: TAC plot beside the outcome-measure plot
#'
#' The two plots share a single colour legend, which is shown only on the TAC
#' plot (left); it is removed from the outcome plot (right) to avoid showing the
#' same legend twice.
#'
#' @noRd
tstar_outrow <- function(lowroi_fit, medroi_fit, highroi_fit, outcome_df,
                         outcome_ylab, tstarInclFrames, gridbreaks,
                         outcome_ylim = NULL) {

  # TAC plot, including the reference TAC when the model has one
  tac_cols <- list(
    Time = lowroi_fit$tacs$Time,
    Low = lowroi_fit$tacs$Target,
    Medium = medroi_fit$tacs$Target,
    High = highroi_fit$tacs$Target
  )
  if ("Reference" %in% names(lowroi_fit$tacs)) {
    tac_cols <- c(list(Time = lowroi_fit$tacs$Time, Reference = lowroi_fit$tacs$Reference),
                  tac_cols[c("Low", "Medium", "High")])
  }
  tacplotdf <- data.frame(do.call(cbind, tac_cols))
  tacplotdf <- tidyr::gather(tacplotdf, key = Region, value = Radioactivity, -Time)
  tacplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(tacplotdf$Region)))

  myColors <- RColorBrewer::brewer.pal(4, "Set1")
  names(myColors) <- levels(tacplotdf$Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  tacplot <- ggplot(tacplotdf, aes(x = Time, y = Radioactivity, colour = Region)) +
    geom_point() + geom_line() + colScale

  # Outcome measure plot (shares the colour scale, but drops the legend)
  outcomeplotdf <- tidyr::gather(outcome_df, key = Region, value = Outcome, -Frames, -Time)
  outcomeplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(outcomeplotdf$Region)))

  outcomeplot <- ggplot(outcomeplotdf, aes(x = Frames, y = Outcome, colour = Region)) +
    geom_point() + geom_line() +
    tstar_framebreaks(tstarInclFrames, gridbreaks) +
    ylab(outcome_ylab) + colScale + theme(legend.position = "none")

  if (!is.null(outcome_ylim)) {
    outcomeplot <- outcomeplot + coord_cartesian(ylim = outcome_ylim)
  }

  cowplot::plot_grid(tacplot, outcomeplot, rel_widths = c(2, 1))
}


#' Assemble the full t* plot grid from the linearised-fit row downwards
#'
#' Builds the R-squared and maximum-percentage rows, the TAC/outcome row, and
#' stacks them under the supplied linearised-fit row. Optionally saves a jpeg.
#'
#' @noRd
tstar_finalise_plot <- function(linrow, lowroi_fit, medroi_fit, highroi_fit,
                                r2_df, maxperc_df, outcome_df, outcome_ylab,
                                t_tac, tstarInclFrames, gridbreaks,
                                outcome_ylim = NULL, filename = NULL,
                                modelname = NULL) {

  r2 <- tstar_r2plots(r2_df, tstarInclFrames, gridbreaks)
  mp <- tstar_maxpercplots(maxperc_df, t_tac, tstarInclFrames, gridbreaks)

  r2row <- cowplot::plot_grid(r2$low, r2$med, r2$high, nrow = 1)
  mprow <- cowplot::plot_grid(mp$low, mp$med, mp$high, nrow = 1)
  outrow <- tstar_outrow(lowroi_fit, medroi_fit, highroi_fit, outcome_df,
                         outcome_ylab, tstarInclFrames, gridbreaks, outcome_ylim)

  totalplot <- cowplot::plot_grid(linrow, r2row, mprow, outrow, nrow = 4)

  if (!is.null(filename)) {
    jpeg(filename = paste0(filename, "_", modelname, ".jpeg"), width = 300, height = 400, units = "mm", res = 600)
    print(totalplot)
    dev.off()
  }

  totalplot
}
