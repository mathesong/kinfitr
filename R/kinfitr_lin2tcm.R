#' Linear Two Tissue Compartment Model
#'
#' Function to fit the linearised 2TCM function.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param input Data frame containing the blood, plasma, and parent fraction
#'   concentrations over time.  This can be generated using the
#'   \code{blood_interp} function.
#' @param weights Optional. Numeric vector of the weights assigned to each frame
#'   in the fitting. We include zero at time zero: if not included, it is added.
#'   If not specified, uniform weights will be used.
#' @param inpshift Optional. The number of minutes by which to shift the timing
#'   of the input data frame forwards or backwards. If not specified, this will
#'   be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will
#'   be fitted. If specified as a number (e.g. 0.05 for 5%), then that value
#'   will be used.
#' @param dur Optional. Numeric vector of the time durations of the frames. If
#'   not included, the integrals will be calculated using trapezoidal
#'   integration.
#' @param frameStartEnd Optional: This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This can be used to assess
#'   time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and
#'   end time point instead of defining the frame numbers using frameStartEnd.
#'   This function will restrict the model to all time frames whose t_tac is
#'   between the values, i.e. c(0,5) will select all frames with midtimes during
#'   the first 5 minutes.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the
#'   model fit object \code{out$fit}, a dataframe containing the TACs of the
#'   data \code{out$tacs}, a dataframe containing the fitted values
#'   \code{out$fitvals}, the blood input data frame after time shifting
#'   \code{input}, a vector of the weights \code{out$weights}, the inpshift
#'   value used \code{inpshift} and the specified vB value \code{out$vB}.
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
#' fit1 <- lin2tcm(t_tac, tac, input, weights, inpshift=0.1)
#' fit2 <- lin2tcm(t_tac, tac, input, weights, inpshift = 0.1, vB=0.05)
#' fit3 <- lin2tcm(t_tac, tac, input, weights, inpshift = 0.1, vB=0.05, dur = dur)
#' fit4 <- lin2tcm(t_tac, tac, input, weights, inpshift = 0.1, dur = dur)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Oikonen, V (2003). Multilinear solution for 4-compartment model:
#'   I. Tissue compartments in series. Gjedde A, Wong DF 1990. Modeling
#'   neuroreceptor binding of radioligands in vivo. In: Quantitative imaging:
#'   neuroreceptors, neurotransmitters, and enzymes. (Eds. Frost JJ, Wagner HM
#'   Jr). Raven Press, 51-79.
#'
#' @export
lin2tcm <- function(t_tac, tac, input, weights = NULL, inpshift = 0,
                    vB = NULL, dur = NULL, frameStartEnd = NULL, timeStartEnd = NULL) {


  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
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

  # Blood Volume Correction

  if(!is.null(vB)) {

    vB_fitted = FALSE
    vB = vB
    i_tac <- (i_tac - vB * blood) / (1 - vB)

  } else {

    vB_fitted = TRUE

  }

  if(vB_fitted) {

    if (!is.null(dur)) {

      term1 <- as.numeric(pracma::cumtrapz(interptime,
                                           pracma::cumtrapz(interptime, aif)))
      term2 <- as.numeric(pracma::cumtrapz(interptime, aif))
      term3 <- aif
      term4 <- -1*frame_cumsum(dur, frame_cumsum(dur, tac))
      term5 <- -1*frame_cumsum(dur, tac)

      term1 <- pracma::interp1(interptime, term1, t_tac, method = "linear")
      term2 <- pracma::interp1(interptime, term2, t_tac, method = "linear")
      term3 <- pracma::interp1(interptime, term3, t_tac, method = "linear")

    } else {

      term1 <- as.numeric(pracma::cumtrapz(interptime,
                                           pracma::cumtrapz(interptime, aif)))
      term2 <- as.numeric(pracma::cumtrapz(interptime, aif))
      term3 <- aif
      term4 <- -1*as.numeric(pracma::cumtrapz(interptime,
                                              pracma::cumtrapz(interptime, i_tac)))
      term5 <- -1*as.numeric(pracma::cumtrapz(interptime, i_tac))

      term1 <- pracma::interp1(interptime, term1, t_tac, method = "linear")
      term2 <- pracma::interp1(interptime, term2, t_tac, method = "linear")
      term3 <- pracma::interp1(interptime, term3, t_tac, method = "linear")
      term4 <- pracma::interp1(interptime, term4, t_tac, method = "linear")
      term5 <- pracma::interp1(interptime, term5, t_tac, method = "linear")

    }

    terms <- data.frame(Term1 = term1, Term2 = term2, Term3 = term3,
                        Term4 = term4, Term5 = term5)

    # Solution

    lin2tcm_model <- lm(tac ~ term1 + term2 + term3 + term4 + term5 - 1,
                        weights = weights)

    # Output

    coefs <- as.numeric(lin2tcm_model$coefficients)

    K1 <- coefs[2] - (coefs[3] * coefs[5])
    k2 <- coefs[5] - (coefs[1] - (coefs[3] * coefs[4]))/K1
    k4 <- coefs[4] / k2
    k3 <- coefs[5] - k2 - k4
    vB <- coefs[3]

  } else {

    if (!is.null(dur)) {

      term1 <- as.numeric(pracma::cumtrapz(interptime,
                                           pracma::cumtrapz(interptime, aif)))
      term2 <- as.numeric(pracma::cumtrapz(interptime, aif))
      term3 <- -1*frame_cumsum(dur, frame_cumsum(dur, tac))
      term4 <- -1*frame_cumsum(dur, tac)

      term1 <- pracma::interp1(interptime, term1, t_tac, method = "linear")
      term2 <- pracma::interp1(interptime, term2, t_tac, method = "linear")

    } else {

      term1 <- as.numeric(pracma::cumtrapz(interptime,
                                           pracma::cumtrapz(interptime, aif)))
      term2 <- as.numeric(pracma::cumtrapz(interptime, aif))
      term3 <- -1*as.numeric(pracma::cumtrapz(interptime,
                                              pracma::cumtrapz(interptime, i_tac)))
      term4 <- -1*as.numeric(pracma::cumtrapz(interptime, i_tac))

      term1 <- pracma::interp1(interptime, term1, t_tac, method = "linear")
      term2 <- pracma::interp1(interptime, term2, t_tac, method = "linear")
      term3 <- pracma::interp1(interptime, term3, t_tac, method = "linear")
      term4 <- pracma::interp1(interptime, term4, t_tac, method = "linear")

    }

    terms <- data.frame(Term1 = term1, Term2 = term2, Term3 = term3,
                        Term4 = term4)

    # Solution

    lin2tcm_model <- lm(tac ~ term1 + term2 + term3 + term4 - 1,
                        weights = weights)

    # Output

    coefs <- as.numeric(lin2tcm_model$coefficients)

    K1 <- coefs[2]
    k2 <- coefs[4] - coefs[1] / K1
    k4 <- coefs[3] / k2
    k3 <- coefs[4] - k2 - k4

  }


  # Output

  Vt <- as.numeric( (K1/k2)*(1 + (k3/k4)))

  par <- as.data.frame(list(K1 = K1, k2 = k2, k3 = k3, k4 = k4, vB = vB, Vt = Vt))
  fit <- lin2tcm_model

  tacs <- data.frame(Time = t_tac, Target = tac)

  fitvals <- data.frame(
    Time = t_tac, Target = tac,
    Target_fitted = as.numeric(predict(lin2tcm_model)), Weights = weights
  )

  if(!is.null(dur)) {
    fitvals$Duration = dur
  }

  fitvals <- dplyr::bind_cols(fitvals, terms)

  input <- newvals$input

  out <- list(
    par = par, fit = lin2tcm_model, tacs = tacs,
    fitvals = fitvals, input = input, weights = weights,
    inpshift = inpshift, vB = vB, vB_fitted = vB_fitted, model = "lin2tcm"
  )

  class(out) <- c("lin2tcm", "kinfit")

  return(out)
}


#' Plot: Linear 2TCM
#'
#' Function to visualise the fit of the linearised 2TCM model to data.
#'
#' @param lin2tcmout The output object of the linearised 2TCM fitting procedure.
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
#' fit <- lin2tcm(t_tac, tac, input, weights, inpshift=0.1)
#' plot_lin2tcmfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_lin2tcmfit <- function(lin2tcmout, roiname = NULL) {
  measured <- data.frame(
    Time = lin2tcmout$tacs$Time,
    Target.measured = lin2tcmout$tacs$Target,
    Weights = lin2tcmout$weights
  )

  fitted <- data.frame(
    Time = lin2tcmout$fitvals$Time,
    Target.fitted = lin2tcmout$fitvals$Target_fitted,
    Weights = lin2tcmout$fitvals$Weights
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
    guides(shape = "none", color = guide_legend(order = 1)) + scale_size(range = c(1, 3))

  return(outplot)
}


#' Profile the inpshift using the linearised 2TCM
#'
#' Function to fit the linearised 2TCM function with several different delay
#' values to find the optimal delay.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param input Data frame containing the blood, plasma, and parent fraction
#'   concentrations over time.  This can be generated using the
#'   \code{blood_interp} function.
#' @param weights Optional. Numeric vector of the weights assigned to each frame
#'   in the fitting. We include zero at time zero: if not included, it is added.
#'   If not specified, uniform weights will be used.
#' @param vB Optional. The blood volume fraction.  If not specified, this will
#'   be fitted. If specified as a number (e.g. 0.05 for 5%), then that value
#'   will be used.
#' @param dur Optional. Numeric vector of the time durations of the frames. If
#'   not included, the integrals will be calculated using trapezoidal
#'   integration.
#' @param frameStartEnd Optional: This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This can be used to assess
#'   time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and
#'   end time point instead of defining the frame numbers using frameStartEnd.
#'   This function will restrict the model to all time frames whose t_tac is
#'   between the values, i.e. c(0,5) will select all frames with midtimes during
#'   the first 5 minutes.
#' @param inpshift_vals Optional. The values of the inpshift to assess with the
#'   grid. By default, a grid between -1 and 1 with spacing of 0.01 will be
#'   used.
#'
#' @return A plot with the residual weighted sums of squares for each value of
#'   the input shift
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
#' lin2tcm_inpshiftProfile(t_tac, tac, input, weights,
#'    inpshift_vals = seq(-0.5, 0.5, length.out=101))
#' lin2tcm_inpshiftProfile(t_tac, tac, input, dur = dur)
#' lin2tcm_inpshiftProfile(t_tac, tac, input, vB=0.05,
#'   frameStartEnd = c(1,15))
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Oikonen, V (2003). Multilinear solution for 4-compartment model:
#'   I. Tissue compartments in series. Gjedde A, Wong DF 1990. Modeling
#'   neuroreceptor binding of radioligands in vivo. In: Quantitative imaging:
#'   neuroreceptors, neurotransmitters, and enzymes. (Eds. Frost JJ, Wagner HM
#'   Jr). Raven Press, 51-79.
#'
#' @export
lin2tcm_inpshiftProfile <- function(t_tac, tac, input, weights = NULL, vB = NULL,
                                    dur = NULL, frameStartEnd = NULL,
                                    timeStartEnd = NULL, inpshift_vals = NULL) {

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    frameStartEnd <- c(which(t_tac >= timeStartEnd[1])[1],
                       tail(which(t_tac <= timeStartEnd[2]), 1))
  }

  if ( is.null(inpshift_vals) ) {
    inpshift_vals <- seq(from = -1, to = 1, by = 0.01)
  }


  inpshift_RSS <- purrr::map_dbl(inpshift_vals, ~lin2tcm_RSS(
    t_tac, tac, input, weights, inpshift = .x,
    vB, dur, frameStartEnd))


  inpshift_profile <- tibble::tibble(
    inpshift = inpshift_vals,
    log_RSS = log(inpshift_RSS),
    lag = dplyr::lag(log_RSS, 1),
    lead = dplyr::lead(log_RSS, 1)
  )

  inpshift_labels <- inpshift_profile %>%
    dplyr::mutate(Label = ifelse(log_RSS < lag & log_RSS < lead,
                                 yes = inpshift, no = NA)) %>%
    dplyr::filter(!is.na(Label))  %>%
    dplyr::arrange(log_RSS) %>%
    dplyr::mutate(Label = round(Label, 2)) %>%
    dplyr::slice(1:3)

  inpshift_profile <- dplyr::left_join(inpshift_profile, inpshift_labels)


  ggplot(inpshift_profile, aes(x=inpshift, y=log_RSS)) +
    geom_line() +
    geom_text(aes(label=Label), nudge_y = -0.4) +
    geom_point(aes(x=Label), size=4, shape=1) +
    labs(y = "log(RSS)", x="Input Shift (min)")

}

lin2tcm_RSS <- function(t_tac, tac, input, weights = NULL, inpshift = 0,
                        vB = NULL, dur = NULL, frameStartEnd = NULL) {

  fit <- lin2tcm(t_tac, tac, input, weights, inpshift,
                 vB, dur, frameStartEnd)

  sum(weights(fit$fit) * residuals(fit$fit)^2)



}
