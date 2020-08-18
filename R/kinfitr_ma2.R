#' Ichise Multilinear Analysis 2
#'
#' Function to fit the MA2 of Ichise et al (2002) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
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
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit},
#' a dataframe containing the TACs of the data \code{out$tacs}, a dataframe containing the fitted values \code{out$fitvals},
#' the blood input data frame after time shifting \code{input}, a vector of the weights \code{out$weights},
#' the inpshift value used \code{inpshift} and the specified vB value \code{out$vB}.
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
#' fit1 <- ma2(t_tac, tac, input, weights)
#' fit2 <- ma2(t_tac, tac, input, weights, inpshift = 0.1, vB = 0.05)
#' fit3 <- ma2(t_tac, tac, input, weights, inpshift = 0.1, vB = 0.05, dur = dur)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Ichise M, Toyama H, Innis RB, Carson RE. Strategies to improve neuroreceptor parameter estimation by linear regression analysis. Journal of Cerebral Blood Flow & Metabolism. 2002 Oct 1;22(10):1271-81.
#'
#' @export


ma2 <- function(t_tac, tac, input, weights = NULL, inpshift = 0, vB = 0,
                dur = NULL, frameStartEnd = NULL) {


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

    term1 <- as.numeric(pracma::cumtrapz(interptime,
                                         pracma::cumtrapz(interptime, aif)))
    term2 <- frame_cumsum(dur, frame_cumsum(dur, tac))
    term3 <- frame_cumsum(dur, tac)
    term4 <- as.numeric(pracma::cumtrapz(interptime, aif))

    term1 <- pracma::interp1(interptime, term1, t_tac, method = "linear")
    term4 <- pracma::interp1(interptime, term4, t_tac, method = "linear")

  } else {

    term1 <- as.numeric(pracma::cumtrapz(interptime,
                                         pracma::cumtrapz(interptime, aif)))
    term2 <- as.numeric(pracma::cumtrapz(interptime,
                                         pracma::cumtrapz(interptime, i_tac)))
    term3 <- as.numeric(pracma::cumtrapz(interptime, i_tac))
    term4 <- as.numeric(pracma::cumtrapz(interptime, aif))

    term1 <- pracma::interp1(interptime, term1, t_tac, method = "linear")
    term2 <- pracma::interp1(interptime, term2, t_tac, method = "linear")
    term3 <- pracma::interp1(interptime, term3, t_tac, method = "linear")
    term4 <- pracma::interp1(interptime, term4, t_tac, method = "linear")

  }


  # Solution

  ma2_model <- lm(tac ~ term1 + term2 + term3 + term4 - 1, weights = weights)

  # Output

  Vt <- as.numeric(-ma2_model$coefficients[1] / ma2_model$coefficients[2])

  coefs <- as.numeric(ma2_model$coefficients)

  Vt <- as.numeric(-coefs[1] / coefs[2])

  Vs <- (-coefs[1] * (coefs[1] + coefs[3] * coefs[4]) + coefs[2] * coefs[4]^2) / (coefs[2] * (coefs[1] + coefs[3] * coefs[4]))

  Vnd <- Vt - Vs


  # Attempting to get rate constants

  K1 <- coefs[4]
  k2 <- -((coefs[1] / coefs[4]) + coefs[3])
  k4 <- -coefs[2] / k2
  k3 <- -coefs[3] - k4 - k2


  par <- as.data.frame(list(Vt = Vt, Vs = Vs, Vnd = Vnd, K1 = K1, k2 = k2, k3 = k3, k4 = k4))
  fit <- ma2_model

  tacs <- data.frame(Time = t_tac, Target = tac, Target_uncor = tac_uncor) # uncorrected for blood volume

  if (!is.null(dur)) {
    tacs$Duration <- dur
  }

  fitvals <- data.frame(
    Time = t_tac, Target = tac, Term1 = term1, Term2 = term2,
    Term3 = term3, Term4 = term4,
    Target_fitted = as.numeric(predict(ma2_model)), Weights = weights
  )

  input <- newvals$input

  out <- list(
    par = par, fit = ma2_model, tacs = tacs,
    fitvals = fitvals, input = input, weights = weights,
    inpshift = inpshift, vB = vB, model = "ma2"
  )

  class(out) <- c("ma2", "kinfit")

  return(out)
}


#' Plot: Ichise Multilinear Analysis 2
#'
#' Function to visualise the fit of the MA2 model to data.
#'
#' @param ma2out The output object of the MA2 fitting procedure.
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
#' fit <- ma2(t_tac, tac, input, weights)
#' plot_ma2fit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_ma2fit <- function(ma2out, roiname = NULL) {
  measured <- data.frame(
    Time = ma2out$tacs$Time,
    Target.measured = ma2out$tacs$Target,
    Weights = ma2out$weights
  )

  fitted <- data.frame(
    Time = ma2out$fitvals$Time,
    Target.fitted = ma2out$fitvals$Target_fitted,
    Weights = ma2out$fitvals$Weights
  )

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
    guides(shape = FALSE, color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

  return(outplot)
}
