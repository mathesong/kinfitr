#' Standardised Uptake Value
#'
#' Function to calculate the standardised uptake value and its integral.
#'
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param dur Numeric vector of durations for each frame in minutes. This will be used instead of the middle points of the frames if provided.
#' @param injRad The injected radioactivity.  If not included, this will be set to 1 in case one is using SUV ratios.
#' @param bodymass The body mass of the participant. If not included, this will be set to 1 in case one is using SUV ratios.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes. Note that this requires t_tac.
#'
#' @return A list with a data frame of the calculated parameters \code{out$par} and a dataframe containing the TACs both of the
#' original data and the transformed values \code{out$tacs}. \code{out$tacs} contains the whole measured TAC, with a logical
#' \code{Included} column marking the frames which the integral was taken over, i.e. those selected by
#' \code{frameStartEnd} or \code{timeStartEnd}.
#'
#' \code{out$par} contains \code{SUV}, the mean SUV over the frames which were included; \code{SUV_AUC}, its integral over
#' those frames; \code{SUV_denominator}, the injected radioactivity divided by the body mass which was applied; and
#' \code{window_duration} and \code{n_frames} describing what was integrated. With the default \code{injRad} and
#' \code{bodymass} of 1, the denominator is 1 and the outcomes are in units of radioactivity concentration.
#'
#' @examples
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' dur <- pbr28$tacs[[2]]$Duration / 60
#' tac <- pbr28$tacs[[2]]$FC
#'
#'
#' fit1 <- suv(tac, t_tac, injRad = 150, bodymass = 85)
#' fit2 <- suv(tac, dur = dur, injRad = 150, bodymass = 85)
#' fit3 <- suv(tac, t_tac = t_tac, dur = dur, injRad = 150, bodymass = 85)
#' fit4 <- suv(tac, t_tac = t_tac, dur = dur, injRad = 150, bodymass = 85,
#'    frameStartEnd = c(30,33))
#' fit5 <- suv(tac, t_tac = t_tac, dur = dur, injRad = 150, bodymass = 85,
#'    timeStartEnd = c(40,61))
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

suv <- function(tac, t_tac = NULL, dur = NULL, injRad = 1, bodymass = 1,
                frameStartEnd = NULL, timeStartEnd = NULL) {

  # Tidying

  if( is.null(frameStartEnd) && !is.null(timeStartEnd) && is.null(t_tac) ) {
    stop("timeStartEnd can only be used if t_tac is provided.")
  }

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  # Which of the measured frames are included, recorded against the input before
  # tidyinput_art() subsets it and potentially adds a frame at time zero, so that
  # the whole TAC can be returned marked with what was used.
  if (is.null(frameStartEnd)) {
    included <- rep(TRUE, length(tac))
  } else {
    included <- seq_along(tac) %in% (frameStartEnd[1]:frameStartEnd[2])
  }

  n_frames <- sum(included)

  # Kept for the output: below, tac/t_tac/dur become the tidied subset which
  # the integral is taken over.
  tac_measured <- tac
  t_tac_measured <- t_tac
  dur_measured <- dur

  if (!is.null(dur) & is.null(t_tac)) {

    tidyinput <- tidyinput_art(dur, tac, tac,
                               frameStartEnd) # Don't need weights, thus just set to same as tac
    dur <- tidyinput$t_tac

    tac <- tidyinput$tac
  }

  if (is.null(dur) & !is.null(t_tac)) {
    tidyinput <- tidyinput_art(t_tac, tac, tac,
                               frameStartEnd) # Don't need weights, thus just set to same as tac
    t_tac <- tidyinput$t_tac

    tac <- tidyinput$tac
  }

  if (is.null(dur) & is.null(t_tac)) {
    stop("Either t_tac or dur must be provided")
  }

  if (!is.null(dur) & !is.null(t_tac)) {
    tidyinput1 <- tidyinput_art(dur, tac, tac,
                               frameStartEnd) # Don't need weights, thus just set to same as tac
    dur <- tidyinput1$t_tac

    tidyinput2 <- tidyinput_art(t_tac, tac, tac,
                               frameStartEnd) # Don't need weights, thus just set to same as tac
    t_tac <- tidyinput2$t_tac

    tac <- tidyinput1$tac
  }


  # 'Model'

  denominator <- injRad / bodymass

  # Integrate over the measured frames which were included, rather than over the
  # vectors tidied above. tidyinput_art() subsets to the window and then prepends
  # a frame at time zero, so for a window which does not begin at the first frame
  # it would otherwise place a point at the origin ahead of the window: the
  # trapezoidal integral would then include the whole area between injection and
  # the start of the window, and the duration would be measured from zero. The
  # tidying is kept above for the input checking it performs.
  inc <- which(included)
  tac_inc <- tac_measured[inc]

  if (is.null(dur_measured)) {
    t_inc <- t_tac_measured[inc]

    # A frame at time zero belongs only when integrating from the start of the
    # acquisition, which is what an unwindowed TAC has always done.
    if (inc[1] == 1L && min(t_inc) > 0) {
      t_inc <- c(0, t_inc)
      tac_inc <- c(0, tac_inc)
    }

    intSUV <- pracma::trapz(t_inc, tac_inc) / denominator
    window_duration <- max(t_inc) - min(t_inc)
  } else {
    dur_inc <- dur_measured[inc]

    intSUV <- sum((tac_inc * dur_inc) / denominator)
    window_duration <- sum(dur_inc)
  }


  # Output

  par <- as.data.frame(list(
    SUV = intSUV / window_duration,
    SUV_AUC = intSUV,
    SUV_denominator = denominator,
    window_duration = window_duration,
    n_frames = n_frames
  ))

  # The whole measured TAC, marked with the frames which the integral was taken
  # over, so that a plot can show both the curve and the window within it.
  tacs <- data.frame(
    Target = tac_measured,
    SUV = tac_measured / denominator
  )

  if (is.null(dur_measured)) {
    tacs$Duration <- NA
  } else {
    tacs$Duration <- dur_measured
  }

  if (is.null(t_tac_measured)) {
    tacs$Time <- NA
  } else {
    tacs$Time <- t_tac_measured
  }

  tacs$Included <- included

  out <- list(par = par, tacs = tacs, model = "suv")

  class(out) <- c("suv", "kinfit")

  return(out)
}


#' Standardised Uptake Value Ratio
#'
#' Function to calculate the standardised uptake value ratio (SUVR): the target
#' region integral over a frame or time window divided by the reference region
#' integral over the same window. This follows the regional TAC ratio approach
#' described by the Turku PET Centre's \code{dftratio}.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param roitac Numeric vector of radioactivity concentrations in the target tissue for each frame.
#' @param dur Numeric vector of durations for each frame in minutes. If provided, the integral is a frame-duration
#' weighted sum rather than a trapezoidal approximation, matching \code{\link{suv}}.
#' @param injRad The injected radioactivity. If not included, this will be set to 1, as in \code{\link{suv}}, so that
#' the SUV outcomes are in units of radioactivity concentration. The SUVR is unaffected either way.
#' @param bodymass The body mass of the participant. If not included, this will be set to 1, as in \code{\link{suv}}.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' Either element may be \code{NA} for a one-sided window.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame
#' numbers using frameStartEnd. As in \code{\link{suv}}, this restricts the calculation to all time frames whose t_tac is
#' between the values: whole frames are selected by their midpoint, and never a part of a frame. Either element may be
#' \code{NA} for a one-sided window.
#'
#' @return A list with a data frame of the calculated parameters \code{out$par}, a data frame of the whole TACs with a
#' logical \code{Included} column marking the frames which were integrated \code{out$tacs}, and the resolved frame window
#' \code{out$window}.
#'
#' \code{out$par} contains \code{SUVR}, the ratio itself; \code{SUV} and \code{SUV_ref}, the mean SUV of each region over
#' the window; \code{SUV_AUC} and \code{SUV_ref_AUC}, their integrals; \code{SUV_denominator}, the injected radioactivity
#' divided by the body mass which was applied; and \code{window_duration} and \code{n_frames} describing what was
#' integrated.
#'
#' Because SUVR is a ratio of two integrals over the same frames, the injected radioactivity and body mass cancel: the
#' SUVR is identical whether it is calculated from radioactivity concentrations or from SUV, so it is correct whether or
#' not a dose was provided. The SUV outcomes do not cancel, and are in units of radioactivity concentration without one,
#' which \code{SUV_denominator} records.
#'
#' @examples
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' reftac <- simref$tacs[[2]]$Reference
#' roitac <- simref$tacs[[2]]$ROI1
#' dur <- simref$tacs[[2]]$Duration
#'
#' fit1 <- suvr(t_tac, reftac, roitac, dur = dur)
#' fit2 <- suvr(t_tac, reftac, roitac, dur = dur, timeStartEnd = c(20, 60))
#' fit3 <- suvr(t_tac, reftac, roitac, dur = dur, timeStartEnd = c(20, 60),
#'    injRad = 150, bodymass = 85)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

suvr <- function(t_tac, reftac, roitac, dur = NULL,
                 injRad = 1, bodymass = 1,
                 frameStartEnd = NULL, timeStartEnd = NULL) {

  n <- length(t_tac)

  if (n == 0) {
    stop("t_tac is empty.")
  }

  if (length(roitac) != n || length(reftac) != n) {
    stop("The lengths of t_tac, reftac and/or roitac are not equal")
  }

  if (!is.null(dur) && length(dur) != n) {
    stop("The lengths of t_tac and dur are not equal")
  }

  window <- suvr_window(t_tac, n, frameStartEnd, timeStartEnd)
  frames <- window$start:window$end
  included <- seq_len(n) %in% frames


  # 'Model'

  # Both regions go through suv() with the resolved window passed as frame
  # numbers, so that the two integrals are taken over exactly the same frames
  # and the reported window cannot disagree with what was integrated.
  roi_suv <- suv(
    tac = roitac, t_tac = t_tac, dur = dur,
    injRad = injRad, bodymass = bodymass,
    frameStartEnd = c(window$start, window$end)
  )

  ref_suv <- suv(
    tac = reftac, t_tac = t_tac, dur = dur,
    injRad = injRad, bodymass = bodymass,
    frameStartEnd = c(window$start, window$end)
  )


  # Output

  par <- as.data.frame(list(
    SUVR = roi_suv$par$SUV_AUC / ref_suv$par$SUV_AUC,
    SUV = roi_suv$par$SUV,
    SUV_ref = ref_suv$par$SUV,
    SUV_AUC = roi_suv$par$SUV_AUC,
    SUV_ref_AUC = ref_suv$par$SUV_AUC,
    SUV_denominator = roi_suv$par$SUV_denominator,
    window_duration = roi_suv$par$window_duration,
    n_frames = roi_suv$par$n_frames
  ))

  tacs <- data.frame(
    Time = t_tac,
    Reference = reftac,
    Target = roitac,
    Included = included
  )

  if (!is.null(dur)) {
    tacs$Duration <- dur
  }

  out <- list(
    par = par, tacs = tacs, window = window, model = "suvr"
  )

  class(out) <- c("suvr", "kinfit")

  return(out)
}

# Resolve a frame or time window down to a first and last frame number, matching
# the convention of suv() that a time window selects whole frames by their
# midpoints. Either end of the window may be missing, giving a one-sided window.
# This is normalised here, rather than left to suv(), so that a one-sided window
# does not become an NA subscript, and so that the window which was integrated
# can be reported.
suvr_window <- function(t_tac, n, frameStartEnd, timeStartEnd) {
  scalar <- function(x, i) {
    if (is.null(x) || length(x) < i) {
      return(NA_real_)
    }
    v <- suppressWarnings(as.numeric(x[[i]]))
    if (length(v) != 1 || is.na(v)) NA_real_ else v
  }

  start <- 1L
  end <- n

  if (!is.null(frameStartEnd)) {
    lo <- scalar(frameStartEnd, 1)
    hi <- scalar(frameStartEnd, 2)
    # A window of c(0, 0) means no window, i.e. all frames.
    if (!(is.na(lo) && is.na(hi)) && !isTRUE(all.equal(c(lo, hi), c(0, 0)))) {
      if (!is.na(lo) && lo >= 1) start <- as.integer(min(lo, n))
      if (!is.na(hi) && hi >= 1) end <- as.integer(min(hi, n))
    }
  } else if (!is.null(timeStartEnd)) {
    lo <- scalar(timeStartEnd, 1)
    hi <- scalar(timeStartEnd, 2)
    if (!(is.na(lo) && is.na(hi)) && !isTRUE(all.equal(c(lo, hi), c(0, 0)))) {
      if (is.na(lo)) lo <- -Inf
      if (is.na(hi)) hi <- Inf
      sel <- which(t_tac >= lo & t_tac <= hi)
      if (length(sel) == 0) {
        stop(
          "No frames fall within the requested time window of ", lo, " to ", hi,
          " minutes. The frame times run from ", round(min(t_tac), 2), " to ",
          round(max(t_tac), 2), " minutes."
        )
      }
      start <- min(sel)
      end <- max(sel)
    }
  }

  if (start > end) {
    stop(
      "The requested window is empty: it starts at frame ", start,
      " and ends at frame ", end, "."
    )
  }

  list(start = as.integer(start), end = as.integer(end))
}

#' Plot: Standardised Uptake Value Ratio
#'
#' Function to visualise the SUVR calculation. The target and reference TACs are
#' shown, with the frames which were integrated shaded underneath: the two shaded
#' areas are the integrals whose ratio is the SUVR.
#'
#' The shading follows the integration. Where durations were provided, the
#' integrals are frame-duration weighted sums, so each included frame is drawn as
#' a rectangle spanning its own start and end times, whose area is that frame's
#' contribution. Where only times were provided, the integrals are trapezoidal,
#' and the areas beneath the curves are shaded instead.
#'
#' @param suvrout The output object of the suvr fitting procedure.
#' @param roiname Optional. The name of the Target region to see it on the plot.
#' @param refname Optional. The name of the Reference region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' data(simref)
#'
#' t_tac <- simref$tacs[[2]]$Times
#' reftac <- simref$tacs[[2]]$Reference
#' roitac <- simref$tacs[[2]]$ROI1
#' dur <- simref$tacs[[2]]$Duration
#'
#' fit <- suvr(t_tac, reftac, roitac, dur = dur, timeStartEnd = c(20, 60))
#'
#' plot_suvrfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_suvrfit <- function(suvrout, roiname = NULL, refname = NULL) {
  if (is.null(roiname)) {
    roiname <- "ROI"
  }
  if (is.null(refname)) {
    refname <- "Reference"
  }

  n <- nrow(suvrout$tacs)

  duration <- if ("Duration" %in% colnames(suvrout$tacs)) {
    suvrout$tacs$Duration
  } else {
    rep(NA_real_, n)
  }

  measured <- data.frame(
    Time = rep(suvrout$tacs$Time, 2),
    Duration = rep(duration, 2),
    Radioactivity = c(suvrout$tacs$Target, suvrout$tacs$Reference),
    Included = rep(suvrout$tacs$Included, 2),
    Region = c(rep(roiname, n), rep(refname, n))
  )

  measured$Region <- forcats::fct_inorder(factor(measured$Region))

  shaded <- measured[measured$Included, , drop = FALSE]

  myColors <- RColorBrewer::brewer.pal(3, "Set1")[1:2]
  names(myColors) <- levels(measured$Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)
  fillScale <- scale_fill_manual(name = "Region", values = myColors)

  outplot <- ggplot(measured, aes(x = Time, y = Radioactivity, colour = Region))

  if (!all(is.na(duration))) {
    # One rectangle per included frame per region, spanning that frame in time:
    # its area is the frame's term in sum(tac * dur), and the ratio of the two
    # sets of areas is the SUVR.
    shaded$xmin <- shaded$Time - shaded$Duration / 2
    shaded$xmax <- shaded$Time + shaded$Duration / 2

    outplot <- outplot +
      geom_rect(
        data = shaded,
        aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Radioactivity, fill = Region),
        inherit.aes = FALSE,
        alpha = 0.3, colour = NA
      )
  } else {
    outplot <- outplot +
      geom_ribbon(
        data = shaded,
        aes(ymin = 0, ymax = Radioactivity, fill = Region),
        alpha = 0.3, colour = NA
      )
  }

  outplot <- outplot +
    geom_line() +
    geom_point(size = 1.5) +
    expand_limits(y = 0) +
    colScale + fillScale

  return(outplot)
}


#' Plot: Standardised Uptake Value
#'
#' Function to visualise the SUV calculation. The whole SUV time activity curve
#' is shown, with the frames which were included shaded beneath it. The shaded
#' area is \code{SUV_AUC}, and the dashed line is the mean SUV over those frames,
#' i.e. that area divided by the duration which was integrated. Frames outside
#' the window are drawn but not shaded.
#'
#' The shading follows the integration. Where durations were provided, the
#' integral is a frame-duration weighted sum, so each included frame is drawn as
#' a rectangle spanning its own start and end times: the area of each rectangle
#' is that frame's contribution to the total, making the relative contribution of
#' long and short frames apparent. Where only times were provided, the integral
#' is trapezoidal, and the area beneath the curve itself is shaded instead.
#'
#' @param suvout The output object of the SUV fitting procedure.
#' @param roiname Optional. The name of the Target region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' dur <- pbr28$tacs[[2]]$Duration / 60
#' tac <- pbr28$tacs[[2]]$FC
#'
#' fit <- suv(tac, t_tac = t_tac, dur = dur, injRad = 150, bodymass = 85)
#'
#' plot_suvfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_suvfit <- function(suvout, roiname = NULL) {
  if (is.null(roiname)) {
    roiname <- "ROI"
  }

  plotdata <- data.frame(
    Time = suvout$tacs$Time,
    Duration = suvout$tacs$Duration,
    SUV = suvout$tacs$SUV,
    Included = suvout$tacs$Included,
    Region = roiname
  )

  has_dur <- !all(is.na(plotdata$Duration))

  # suv() can be called with durations alone, in which case it records no times.
  # Place the frames end to end so that the curve can still be drawn: this is the
  # same ordering the integral assumes.
  if (all(is.na(plotdata$Time))) {
    plotdata$Time <- cumsum(plotdata$Duration) - plotdata$Duration / 2
    xlabel <- "Cumulative frame duration (min)"
  } else {
    xlabel <- "Time (min)"
  }

  # Without a dose the denominator is 1, and the values are radioactivity
  # concentrations rather than SUV. Say which is being shown.
  ylabel <- if (isTRUE(all.equal(suvout$par$SUV_denominator, 1))) {
    "Radioactivity"
  } else {
    "SUV"
  }

  shaded <- plotdata[plotdata$Included, , drop = FALSE]

  myColor <- RColorBrewer::brewer.pal(3, "Set1")[1]

  outplot <- ggplot(plotdata, aes(x = Time, y = SUV))

  if (has_dur) {
    # One rectangle per included frame, spanning that frame in time and running
    # from zero to its value: its area is the frame's term in sum(tac * dur).
    shaded$xmin <- shaded$Time - shaded$Duration / 2
    shaded$xmax <- shaded$Time + shaded$Duration / 2

    outplot <- outplot +
      geom_rect(
        data = shaded,
        aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = SUV),
        inherit.aes = FALSE,
        fill = myColor, alpha = 0.3, colour = myColor, linewidth = 0.2
      )
  } else {
    # Trapezoidal integration: the area under the curve itself is what was taken.
    outplot <- outplot +
      geom_ribbon(
        data = shaded,
        aes(ymin = 0, ymax = SUV),
        fill = myColor, alpha = 0.3, colour = NA
      )
  }

  outplot <- outplot +
    geom_line(colour = myColor) +
    geom_point(aes(shape = Included), colour = myColor, size = 1.5) +
    geom_hline(
      yintercept = suvout$par$SUV,
      linetype = "dashed", colour = myColor
    ) +
    scale_shape_manual(
      name = "Included",
      values = c(`TRUE` = 16, `FALSE` = 1),
      breaks = c("TRUE", "FALSE"),
      limits = c("TRUE", "FALSE")
    ) +
    expand_limits(y = 0) +
    labs(x = xlabel, y = ylabel)

  return(outplot)
}
