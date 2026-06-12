#' Ichise Multilinear Analysis 1
#'
#' Function to fit the MA1 of Ichise et al (2002) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it  will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param tstar The t* specification for regression. If tstar_type="frames",
#' this is the number of frames from the end to include (e.g., 10 means last 10 frames).
#' If tstar_type="time", this is the time point (in minutes) after which all frames
#' with midpoints later than this time are included. This value can be estimated using \code{ma1_tstar}.
#' @param tstar_type Either "frames" (default) or "time", specifying how to interpret tstar.
#' @param tstarIncludedFrames Deprecated. Use 'tstar' with 'tstar_type="frames"' instead.
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
#' This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit},
#' a dataframe containing the TACs of the data \code{out$tacs}, a dataframe containing the fitted values \code{out$fitvals},
#' the blood input data frame after time shifting \code{input}, a vector of the weights \code{out$weights},
#' the inpshift value used \code{inpshift}, the specified vB value \code{out$vB} and the specified tstarIncludedFrames
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
#' fit1 <- ma1(t_tac, tac, input, 10, weights)
#' fit2 <- ma1(t_tac, tac, input, 10, weights, inpshift = 0.1, vB = 0.05)
#' fit3 <- ma1(t_tac, tac, input, 10, weights, inpshift = 0.1, vB = 0.05, dur = dur)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Ichise M, Toyama H, Innis RB, Carson RE. Strategies to improve neuroreceptor parameter estimation by linear regression analysis. Journal of Cerebral Blood Flow & Metabolism. 2002 Oct 1;22(10):1271-81.
#'
#' @export

ma1 <- function(t_tac, tac, input, tstar, weights = NULL,
                inpshift = 0, vB = 0, dur = NULL, frameStartEnd = NULL, timeStartEnd = NULL,
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

  if (!is.null(dur)) {
    tidyinput_dur <- tidyinput_art(dur, tac, weights, frameStartEnd)
    dur <- tidyinput_dur$t_tac
  }

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

  if (!is.null(dur)) {

    term1 <- as.numeric(pracma::cumtrapz(interptime, aif))
    term1 <- pracma::interp1(interptime, term1, t_tac, method = "linear")

    term2 <- frame_cumsum(dur, tac)

  } else {

    term1 <- as.numeric(pracma::cumtrapz(interptime, aif))
    term2 <- as.numeric(pracma::cumtrapz(interptime, i_tac))

    term1 <- pracma::interp1(interptime, term1, t_tac, method = "linear")
    term2 <- pracma::interp1(interptime, term2, t_tac, method = "linear")

  }

  tac_equil <- tail(tac, tstarIncludedFrames)
  t_tac_equil <- tail(t_tac, tstarIncludedFrames)
  term1_equil <- tail(term1, tstarIncludedFrames)
  term2_equil <- tail(term2, tstarIncludedFrames)
  weights_equil <- tail(weights, tstarIncludedFrames)


  # Solution

  ma1_model <- lm(tac_equil ~ term1_equil + term2_equil - 1, weights = weights_equil)

  # Output

  Vt <- as.numeric(-ma1_model$coefficients[1] / ma1_model$coefficients[2])

  par <- as.data.frame(list(Vt = Vt))
  fit <- ma1_model

  tacs <- data.frame(Time = t_tac, Target = tac, Target_uncor = tac_uncor) # uncorrected for blood volume

  if(!is.null(dur)) { tacs$Duration = dur }

  fitvals <- data.frame(
    Time = t_tac_equil, Target = tac_equil, Term1 = term1_equil, Term2 = term2_equil,
    Target_fitted = as.numeric(predict(ma1_model)), Weights = weights_equil
  )

  input <- newvals$input

  par.se <- par
  names(par.se) <- paste0(names(par.se), ".se")
  par.se$Vt.se <- get_se(ma1_model, "-term1_equil / term2_equil")


  out <- list(
    par = par, par.se = par.se, fit = ma1_model, tacs = tacs, fitvals = fitvals,
    input = input, weights = weights, inpshift = inpshift, vB = vB,
    tstarIncludedFrames = tstarIncludedFrames, model = "ma1"
  )

  class(out) <- c("ma1", "kinfit")

  return(out)
}

#' Plot: Ichise Multilinear Analysis 1
#'
#' Function to visualise the fit of the MA1 model to data.
#'
#' @param ma1out The output object of the MA1 fitting procedure.
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
#' fit <- ma1(t_tac, tac, input, 10, weights)
#' plot_ma1fit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_ma1fit <- function(ma1out, roiname = NULL) {
  measured <- data.frame(
    Time = ma1out$tacs$Time,
    Target.measured = ma1out$tacs$Target,
    Weights = ma1out$weights
  )

  fitted <- data.frame(
    Time = ma1out$fitvals$Time,
    Target.fitted = ma1out$fitvals$Target_fitted,
    Weights = ma1out$fitvals$Weights
  )

  if (is.null(roiname)) {
    roiname <- "ROI"
  }

  measured <- dplyr::rename(measured, !!paste0(roiname, ".measured") := Target.measured)

  fitted <- dplyr::rename(fitted, !!paste0(roiname, ".fitted") := Target.fitted)

  tidymeasured <- tidyr::gather(
    measured,
    key = Region, value = Radioactivity,
    -Time, -Weights, factor_key = F
  )

  tidyfitted <- tidyr::gather(
    fitted,
    key = Region, value = Radioactivity,
    -Time, -Weights, factor_key = F
  )

  Region <- forcats::fct_inorder(factor(c(tidymeasured$Region, tidyfitted$Region)))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  outplot <- ggplot(tidymeasured, aes(x = Time, y = Radioactivity, colour = Region)) +
    geom_point(data = tidymeasured, aes(shape = "a", size = Weights)) +
    geom_line(data = tidyfitted) + colScale +
    guides(shape = "none", color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

  return(outplot)
}


#' Tstar Finder: Ichise Multilinear Analysis 1
#'
#' Function to identify where t* is for MA1.
#'
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param filename The name of the output image: filename_ma1.jpeg
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
#' @return Saves a jpeg of the plots as filename_ma1.jpeg
#'
#' @examples
#' \dontrun{
#' ma1_tstar(t_tac, lowroi, medroi, highroi, input,
#'   filename = "demonstration",
#'   inpshift = onetcmout$par$inpshift, frameStartEnd, gridbreaks = 4
#' )
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export


ma1_tstar <- function(t_tac, lowroi, medroi, highroi, input, filename = NULL, inpshift = 0, vB = 0, frameStartEnd = NULL, timeStartEnd = NULL, gridbreaks = 2) {
  frameStartEnd <- tstar_frameStartEnd(t_tac, frameStartEnd, timeStartEnd)

  frames <- length(t_tac)
  lowroi_fit <- ma1(t_tac, lowroi, input, tstar = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  medroi_fit <- ma1(t_tac, medroi, input, tstar = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  highroi_fit <- ma1(t_tac, highroi, input, tstar = frames, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)

  low_linplot <- plot_ma1fit(lowroi_fit) + ggtitle("Low") + ylim(0, max(lowroi_fit$tacs$Target * 1.1)) + theme(legend.position = "none")
  med_linplot <- plot_ma1fit(medroi_fit) + ggtitle("Medium") + ylim(0, max(medroi_fit$tacs$Target * 1.1)) + theme(legend.position = "none")
  high_linplot <- plot_ma1fit(highroi_fit) + ggtitle("High") + ylim(0, max(highroi_fit$tacs$Target * 1.1)) + theme(legend.position = "none")

  tstarInclFrames <- 3:frames

  fitfunc <- function(roitac, tstar) {
    ma1(t_tac, roitac, input, tstar = tstar, inpshift = inpshift, vB = vB, frameStartEnd = frameStartEnd)
  }

  comp <- tstar_compute(t_tac, lowroi, medroi, highroi, fitfunc,
                        function(f) f$par$Vt, tstarInclFrames)

  outcome_vals <- c(comp$outcome_df$Low, comp$outcome_df$Medium, comp$outcome_df$High)
  outcome_ylim <- c(min(outcome_vals), max(outcome_vals))
  if (outcome_ylim[2] > 20 || outcome_ylim[1] < 0) {
    outcome_ylim <- c(0, 20)
  }

  linrow <- cowplot::plot_grid(low_linplot, med_linplot, high_linplot, nrow = 1)

  totalplot <- tstar_finalise_plot(
    linrow, lowroi_fit, medroi_fit, highroi_fit,
    comp$r2_df, comp$maxperc_df, comp$outcome_df,
    outcome_ylab = expression(V[T]), t_tac = t_tac,
    tstarInclFrames = tstarInclFrames, gridbreaks = gridbreaks,
    outcome_ylim = outcome_ylim, filename = filename, modelname = "ma1"
  )

  return(totalplot)
}
