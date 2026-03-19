#' Simultaneous Estimation of Non-Displaceable Binding (SIME)
#'
#' Function to fit the SIME Model of Ogden et al (2015) to data to estimate Vnd
#' for a set of TACs using long-format input data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes (repeated for
#'   each region). We use the time halfway through the frame as well as a zero.
#'   If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations for each frame,
#'   stacked across regions.
#' @param region Character vector identifying the region for each observation.
#' @param input Data frame containing the blood, plasma, and parent fraction
#'   concentrations over time.  This can be generated using the
#'   \code{blood_interp} function.
#' @param Vndgrid The grid of Vnd values which will be tested to see which one
#'   has the best fit.
#' @param weights Optional. Numeric vector of the weights assigned to each frame
#'   in the fitting. Can be the same length as \code{tac} (one weight per
#'   observation) or the same length as the number of frames per region (weights
#'   are recycled for each region). If not specified, uniform weights will be
#'   used.
#' @param roiweights Optional. Named numeric vector of the weights assigned to
#'   each ROI in the fitting, keyed by region name (e.g.,
#'   \code{c(FC=1, STR=0.8)}). If unnamed, matched by position against unique
#'   regions. If not specified, uniform weights will be used.
#' @param inpshift Optional. The number of minutes by which to shift the timing
#'   of the input data frame forwards or backwards. If not specified, this will
#'   be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will
#'   be set to 0.05. This can be fitted using 1TCM or 2TCM.
#' @param twotcmstart Optional. Region name string specifying which region to
#'   use for fitting a 2TCM model to obtain starting parameters for k2, k3 and
#'   k4 (e.g., \code{"FC"}). Best to use the largest ROI.
#' @param frameStartEnd Optional: This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). Applied per region.
#' @param timeStartEnd Optional. This allows one to specify the beginning and
#'   end time point instead of defining the frame numbers using frameStartEnd.
#'   This function will restrict the model to all time frames whose t_tac is
#'   between the values, i.e. c(0,5) will select all frames with midtimes
#'   during the first 5 minutes.
#' @param k2.start Optional. Starting parameter for fitting of k2. Default is
#'   0.1.
#' @param k2.lower Optional. Lower bound for the fitting of k2. Default is 0.
#' @param k2.upper Optional. Upper bound for the fitting of k2. Default is 0.5.
#' @param k3.start Optional. Starting parameter for fitting of k3. Default is
#'   0.1.
#' @param k3.lower Optional. Lower bound for the fitting of k3. Default is 0.
#' @param k3.upper Optional. Upper bound for the fitting of k3. Default is 0.5.
#' @param k4.start Optional. Starting parameter for fitting of k4. Default is
#'   0.1.
#' @param k4.lower Optional. Lower bound for the fitting of k4. Default is 0.
#' @param k4.upper Optional. Upper bound for the fitting of k4. Default is 0.5.
#' @param success_cutoff Optional. Should values of Vnd for which a certain
#'   proportion of ROIs failed be included? Default is 0.5, i.e. 50% should have
#'   successfully fitted.
#'
#' @return A list with a data frame of the fitted parameter \code{out$par}, a
#'   long-format data frame containing the times, regions and TACs
#'   \code{out$tacs} (with columns Time, Region, Radioactivity), the mean cost
#'   values after fitting (after ROI weighting) \code{out$fitvals}, the ROI cost
#'   values after fitting (before ROI weighting) \code{out$roifits}, the blood
#'   input data frame after time shifting \code{input}, a vector of the weights
#'   \code{out$weights}, a vector of the ROI weights \code{out$roiweights}, the
#'   inpshift value used \code{inpshift} and the vB value used \code{out$vB},
#'   and the success cutoff \code{success_cutoff}.
#'
#'
#' @examples
#' \dontrun{
#' data(pbr28)
#'
#' sime_long <- pbr28$tacs[[2]] %>%
#'   dplyr::select(Times, Weights, FC:CBL) %>%
#'   dplyr::mutate(Times = Times / 60) %>%
#'   tidyr::pivot_longer(cols = -c(Times, Weights),
#'                       names_to = "region", values_to = "tac")
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1
#' )
#'
#' Vndgrid <- seq(from = 0, to = 3, by = 0.5)
#' SIMEout <- SIME(sime_long$Times, sime_long$tac, sime_long$region,
#'   input, Vndgrid,
#'   weights = sime_long$Weights,
#'   inpshift = 0.1, vB = 0.05
#' )
#' }
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
#'
#' @references Ogden RT, Zanderigo F, Parsey RV. Estimation of in vivo
#'   nonspecific binding in positron emission tomography studies without
#'   requiring a reference region. NeuroImage. 2015 Mar 31;108:234-42.
#'
#' @importFrom dplyr "%>%"

SIME <- function(t_tac, tac, region, input, Vndgrid,
                 weights = NULL, roiweights = NULL,
                 inpshift = 0, vB = NULL, twotcmstart = NULL,
                 frameStartEnd = NULL, timeStartEnd = NULL,
                 k2.start = 0.1, k2.lower = 0, k2.upper = 0.5,
                 k3.start = 0.1, k3.lower = 0, k3.upper = 0.5,
                 k4.start = 0.1, k4.lower = 0, k4.upper = 0.5,
                 success_cutoff = 0.5) {

  region <- as.character(region)
  regions <- unique(region)

  # Convert timeStartEnd to frameStartEnd if needed
  if (is.null(frameStartEnd) && !is.null(timeStartEnd)) {
    # Use time vector from first region to determine frame indices
    t_first <- t_tac[region == regions[1]]
    frameStartEnd <- c(which(t_first >= timeStartEnd[1])[1],
                       tail(which(t_first <= timeStartEnd[2]), 1))
  }

  # Tidying via tidyinput_long
  tidied <- tidyinput_long(t_tac, tac, region, weights, frameStartEnd)
  t_tac <- tidied$t_tac
  tac <- tidied$tac
  region <- tidied$region
  weights <- tidied$weights
  regions <- unique(region)

  # ROI weights
  if (is.null(roiweights)) {
    roiweights <- rep(1, length(regions))
    names(roiweights) <- regions
  } else if (!is.null(names(roiweights))) {
    roiweights <- roiweights[regions]
  } else {
    names(roiweights) <- regions
  }
  roiweights <- roiweights / max(roiweights)

  Regions <- data.frame(
    Region = regions,
    roiweights = as.numeric(roiweights[regions]),
    stringsAsFactors = FALSE
  )

  if (length(roiweights) != length(regions)) {
    stop("The number of ROIs and roiweights do not match")
  }

  # Shift timings
  newvals <- shift_timings_long(t_tac, tac, region, input, inpshift)

  t_tac <- newvals$t_tac
  tac <- newvals$tac
  region <- newvals$region
  input <- newvals$input
  t_tac_unique <- newvals$t_tac_unique

  # Extract per-region weights (from tidied data, matching shifted times)
  # After shift_timings_long, the time vector may have changed, so we use

  # the weights from the tidied data which has the same structure
  weights_per_region <- weights[tidied$region == regions[1]]


  # Parameters

  start <- c(k2 = k2.start, k3 = k3.start, k4 = k4.start)
  lower <- c(k2 = k2.lower, k3 = k3.lower, k4 = k4.lower)
  upper <- c(k2 = k2.upper, k3 = k3.upper, k4 = k4.upper)


  # 2tcm Starting Parameters

  if (!is.null(twotcmstart)) {
    twotcm_tac <- tac[region == twotcmstart]
    twotcmout <- twotcm(
      t_tac_unique,
      tac = twotcm_tac, input,
      weights = weights_per_region, inpshift = inpshift,
      vB = vB
    )
    start[1] <- twotcmout$par$k2
    start[2] <- twotcmout$par$k3
    start[3] <- twotcmout$par$k4
    if (is.null(vB)) {
      vB <- twotcmout$par$vB
    }
  } else {
    if (is.null(vB)) {
      vB <- 0.05
    }
  }

  # Solution — build tidytacs directly in long format

  tidytacs <- data.frame(
    Time = t_tac,
    Region = region,
    Radioactivity = tac,
    weights = rep(weights_per_region, times = length(regions)),
    stringsAsFactors = FALSE
  )

  frames <- nrow(tidytacs)

  tidytacs <- tidytacs[rep(1:nrow(tidytacs), times = length(Vndgrid)), ]
  tidytacs$Vnd <- rep(Vndgrid, each = frames)

  tidytacs_nested <- tidyr::nest(tidytacs, tacs = c(-Region, -Vnd))

  fit_SIMEroi <- function(tacs, input, Vnd, vB, start, upper, lower) {
    tac <- tacs$Radioactivity
    t_tac <- tacs$Time
    weights <- tacs$weights

    modeldata <- list(
      tac = tacs$Radioactivity,
      t_tac = tacs$Time,
      weights = tacs$weights,
      input = input
    )

    fit <- minpack.lm::nlsLM(
      tac ~ SIME_model(t_tac, input, Vnd, k2, k3, k4, vB = vB),
      data = modeldata,
      start = start, lower = lower, upper = upper,
      weights = weights, control = minpack.lm::nls.lm.control(maxiter = 200)
    )

    output <- broom::glance(fit)
    output$RSS <- sum(weights(fit) * residuals(fit)^2)

    return(output)
  }

  fit_SIMEroi_possibly <- purrr::possibly(fit_SIMEroi, otherwise = NA, quiet = T)


  tidytacs_nested <- dplyr::mutate(tidytacs_nested, fit = purrr::pmap(
    list(tacs, Vnd),
    fit_SIMEroi_possibly,
    input = input,
    vB = vB,
    start = start,
    lower = lower,
    upper = upper
  ))

  tidypars <- dplyr::mutate(
    tidytacs_nested,
    success = purrr::map_lgl(fit, ~ is.data.frame(.x))
  )

  fit_success = tidypars
  fit_success = dplyr::group_by(fit_success, Vnd)
  fit_success = dplyr::summarise(fit_success, prop_success=mean(success))

  tidypars <- dplyr::filter(tidypars, success == T)
  tidypars <- dplyr::select(tidypars, -tacs, -success)
  tidypars <- tidyr::unnest(tidypars, cols = c(fit))
  tidypars <- dplyr::left_join(Regions, tidypars, by = "Region")
  tidypars <- dplyr::mutate(tidypars, RSSw = RSS * roiweights)
  tidypars <- dplyr::left_join(tidypars, fit_success, by = "Vnd")

  # Calculating SSmean

  SSmean <- dplyr::select(tidypars, Region, roiweights, Vnd, RSSw, prop_success)
  SSmean <- dplyr::filter(SSmean, prop_success >= success_cutoff)
  SSmean <- dplyr::group_by(SSmean, Vnd)
  SSmean <- dplyr::summarise(SSmean, RSSw = mean(RSSw))
  SSmean <- dplyr::ungroup(SSmean)


  # Output

  Vnd <- SSmean$Vnd[which.min(SSmean$RSSw)]

  par <- as.data.frame(list(Vnd = Vnd))

  tacs <- data.frame(
    Time = t_tac,
    Region = region,
    Radioactivity = tac,
    stringsAsFactors = FALSE
  )

  fitvals <- SSmean

  out <- list(
    par = par, tacs = tacs, fitvals = fitvals, roifits = tidypars, input = input,
    weights = weights_per_region, roiweights = roiweights, inpshift = inpshift,
    vB = vB, success_cutoff = success_cutoff, model = "SIME"
  )

  class(out) <- c("SIME", "kinfit")

  return(out)
}


#' Model: SIME
#'
#' This is the SIME Model model itself by which predicted values are generated.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function. It should already be shifted if a shift is desired, using \code{shift_timings}.
#' @param Vnd The specified Vnd value for the given ROI.
#' @param k2 Parameter value for k2
#' @param k3 Parameter value for k3
#' @param k4 Parameter value for k4
#' @param vB Parameter value for vB
#'
#' @return A numeric vector of the predicted values of the TAC in the target region with the given Vnd.
#'
#'
#' @examples
#' \dontrun{
#' SIME_model(t_tac, input, Vnd = 5, k2 = 0.1, k3 = 0.05, k4 = 0.04, vB = 0.05)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Ogden RT, Zanderigo F, Parsey RV. Estimation of in vivo nonspecific binding in positron emission tomography studies without requiring a reference region. NeuroImage. 2015 Mar 31;108:234-42.
#'
#' @export

SIME_model <- function(t_tac, input, Vnd, k2, k3, k4, vB) {
  blood <- input$Blood
  aif <- input$AIF

  interptime <- input$Time
  step <- interptime[2] - interptime[1]

  K1 <- k2 * Vnd

  delta <- sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4)
  th1 <- (k2 + k3 + k4 + delta) / 2
  th2 <- (k2 + k3 + k4 - delta) / 2
  ph1 <- K1 * (th1 - k3 - k4) / delta
  ph2 <- K1 * (th2 - k3 - k4) / (-delta)

  a <- ph1 * exp(-th1 * interptime) + ph2 * exp(-th2 * interptime)
  b <- aif

  i_outtac <- kinfit_convolve(a, b, step)

  # Correction for vB
  i_outtac <- i_outtac * (1 - vB) + vB * aif

  outtac <- pracma::interp1(interptime, i_outtac, t_tac)

  return(outtac)
}

#' Plot: SIME
#'
#' Function to visualise the fit of the SIME Model to data.
#'
#' @param SIMEout The output object of the SIME fitting procedure.
#'
#' @return A ggplot2 object of the plot.
#'
#' @description The black line of this plot is the total cost function for the simultaneous estimation and is the true
#' outcome for the SIME model. I have also included, in coloured points, the costs for the individual ROIs such that
#' one can visualise whether the overall estimated Vnd value is overly influenced by one ROI in particular, or whether
#' one ROI shows a completely different pattern from the rest of the ROIs.
#'
#' @examples
#' \dontrun{
#' plot_SIMEfit(SIMEout)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_SIMEfit <- function(SIMEout) {
  fitvals <- SIMEout$fitvals

  roifits <- SIMEout$roifits

  minmax <- list(
    min = min(roifits$RSSw),
    max = max(roifits$RSSw)
  )

  roifits$Exclude = roifits$prop_success < SIMEout$success_cutoff

  outplot <- ggplot(roifits, aes(x = Vnd, y = RSSw)) +
    geom_point(aes(colour = Region, shape=Exclude)) +
    geom_line(data = fitvals, aes(x = Vnd, y = RSSw), linewidth = 1) +
    geom_vline(xintercept = SIMEout$par$Vnd, linetype = "dashed") +
    xlab(expression(V[ND])) +
    scale_y_log10(
      "Cost",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_shape_manual(values=c(19, 1)) +
    guides(shape="none")

  return(outplot)
}
