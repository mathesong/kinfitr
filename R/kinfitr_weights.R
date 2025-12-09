#' Create weights using durations and counts
#'
#' This set of methods allows for creating a weighting scheme based on TAC data
#' alone. It is advised to calculate weights once for each PET measurement,
#' using a large region (e.g. whole brain).
#'
#' @param t_start The starting times of the frames in minutes.
#' @param t_end The end times of the frames in minutes.
#' @param tac The time activity curve.
#' @param radioisotope The radioisotope.
#' @param radionuclide_halflife The half-life of the radionuclide in minutes. If specified, this overrides the radioisotope parameter. Defaults to NULL.
#' @param method Which method should be used? 1 represents duration /
#'   (tac_uncorrected). 2 represents sqrt(durations*tac_uncorrected). 3
#'   represents duration / tac. 4 represents sqrt(durations). 5 represents
#'   durations * exp((-ln(2)) / halflife ). 6 represents durations /
#'   tac. 7 represents durations. 8 represents duration^2 / tac_uncor.
#'   9 represents((durations^2 / (durations*tac))*corrections^2) (courtesy of
#'   Claus Svarer).
#'   Uncorrected refers to decay correction.
#' @param minweight The minimum weight. Weights will be calculated as a fraction
#'   between this value and 1. A zero frame with duration=0 will be set to 0
#'   though.
#' @param minweight_risetopeak Should there be a linear rise of the minimum weight to
#'   the peak weight value? This means that values before the maximum weight can
#'   be below the minweight specified. This is helpful for downweighting noisy
#'   values towards the start of the TAC. Defaults to FALSE.
#' @param weight_checkn The number of values of the weights to check to make
#'   sure that things haven't gone terribly wrong. It will check that this
#'   number of weights are all above half of the maximum weight.
#'
#' @return The vector of weights.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
#'
#' @examples
#' data(pbr28)
#' s1 <- pbr28$tacs[[1]]
#'
#' weights_create(
#'   s1$StartTime/60,
#'   (s1$StartTime + s1$Duration)/60,
#'   radioisotope = "C11",
#'   tac = s1$WB)
#'
#' weights_create(
#'   s1$StartTime/60,
#'   (s1$StartTime + s1$Duration)/60,
#'   radioisotope = "C11", method=9,
#'   tac = s1$WB)
weights_create <- function(t_start, t_end, tac,
                           radioisotope = c("C11", "O15", "F18"),
                           method = 2,
                           minweight = 0.5,
                           minweight_risetopeak = FALSE,
                           weight_checkn = 5,
                           radionuclide_halflife = NULL) {

  radioisotope <- match.arg(radioisotope, c("C11", "O15", "F18"))

  tac <- ifelse(tac < 0, 0, tac)

  tac_uncor <- decay_uncorrect(t_start, t_end, tac, radioisotope, radionuclide_halflife)
  durations <- t_end - t_start
  t_tac <- t_start + durations/2

  corrections <- tac_uncor / tac
  if( is.nan(corrections[1]) ) corrections[1] <- 1

  if (!is.null(radionuclide_halflife)) {
    hl <- radionuclide_halflife
  } else {
    hl <- dplyr::case_when(
      radioisotope == "C11" ~ 20.4,
      radioisotope == "O15" ~ 2.05,
      radioisotope == "F18" ~ 109.8
    )
  }


  # pmaxes here to prevent denominators sending weights to infinity
  calcweights <- dplyr::case_when(
    method == 1 ~ durations / pmax(tac_uncor, 0.01*max(tac_uncor), na.rm=TRUE),
    method == 2 ~ sqrt(durations*tac_uncor),
    method == 3 ~ sqrt(durations) / pmax(tac, 0.01*max(tac), na.rm=TRUE),
    method == 4 ~ sqrt(durations),
    method == 5 ~ durations * exp( (-1*log(2)) / hl ),
    method == 6 ~ durations / pmax(tac, 0.01*max(tac, na.rm=TRUE), na.rm=TRUE),
    method == 7 ~ durations,
    method == 8 ~ durations^2 / pmax(tac_uncor, 0.01*max(tac_uncor, na.rm=TRUE), na.rm=TRUE),
    # method == 9 ~ (durations^2 / pmax(tac*durations, 0.01*max(tac*durations))) * corrections^2
    method == 9 ~ (durations^2 / pmax(tac, 0.01*max(tac), na.rm=TRUE) * durations) * corrections^2
  )

  # Fixing before checking
  calcweights[durations==0] <- 0

  # Check to make sure things don't go horribly wrong
  maxweight  <- max(calcweights, na.rm=TRUE)
  maxweights <- tail(calcweights[order(calcweights)], weight_checkn)

  if( min(maxweights, na.rm=TRUE) < 0.5*maxweight ) {
    maxweight <- median(maxweights)
  }

  calcweights[calcweights > maxweight] <- maxweight
  calcweights <- calcweights/maxweight

  # Apply scaling factor
  minweights <- rep(minweight, length(calcweights))

  if( minweight_risetopeak ) {
    # TAC midtimes
    t_tac <- t_start + 0.5*(t_end - t_start)

    # Find the maximum TAC weight
    t_weightpeak <- t_tac[tail(which(calcweights==max(calcweights, na.rm=TRUE)), 1)]

    # Assign times a fraction of the peaktime
    t_peakfrac <- t_tac / t_weightpeak
    t_peakfrac[t_peakfrac > 1] <- 1

    # Define a rising minimum to the peak, and then minweight
    minweights <- minweights * t_peakfrac
    minweights[durations==0] <- 0
  }

  if( any(calcweights < minweights) ) {

    min_calcweight <- min(calcweights[durations!=0], na.rm=TRUE)

    # scale 0 - 1
    calcweights <- calcweights - min_calcweight / (1- min_calcweight)

    # proportion
    calcweights <- minweights + calcweights * (1-minweights)
  }

  # Any dur=0 shouldn't have the scaling factor
  calcweights[durations==0] <- 0

  return(calcweights)
}

#' Create weights from BIDS data
#'
#' This uses the BIDS PET JSON sidecar information and a TAC to create a
#' weights vector.
#'
#' @param petinfo The parsed information from the PET JSON sidecar. This can be
#' extracted simply using jsonlite::fromJSON() on the JSON file, or else at a
#' study level using bids_parse_study()
#' @param tac The time activity curve: preferably for a region with high signal-to-noise ratio.
#' @param method Which method should be used? 1 represents duration^2 /
#'   (tac_uncorrected). 2 represents sqrt(durations*tac_uncorrected). 3
#'   represents duration / tac. 4 represents sqrt(durations). 5 represents
#'   durations * exp((-ln(2)) / halflife ). 6 represents durations /
#'   tac. 7 represents durations.
#' @param minweight The minimum weight. Weights will be calculated as a fraction
#'   between this value and 1. A zero frame with duration=0 will be set to 0
#'   though.
#' @param minweight_risetopeak Should there be a linear rise of the minimum weight to
#'   the peak weight value? This means that values before the maximum weight can
#'   be below the minweight specified. This is helpful for downweighting noisy
#'   values towards the start of the TAC. Defaults to FALSE.
#' @param weight_checkn The number of values of the weights to check to make
#'   sure that things haven't gone terribly wrong. It will check that this
#'   number of weights are all above half of the maximum weight.
#'
#' @return A vector of weights.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' bids_weights_create(studydata$petinfo[[1]],tacdata$Neo.cx)
#' }
weights_create_bids <- function(petinfo, tac, method=2,
                                minweight = 0.7,
                                minweight_risetopeak = FALSE,
                                weight_checkn=5) {

  weights_create(t_start = petinfo$FrameTimesStart/60,
                 t_end = with(petinfo, FrameTimesStart +
                                FrameDuration)/60,
                 tac = tac,
                 radioisotope = petinfo$TracerRadionuclide,
                 method=method, minweight=minweight,
                 weight_checkn=weight_checkn)

}


weights_to_var <- function(weights) {
  var_est <- rep(NA_real_, length(weights))
  valid   <- is.finite(weights) & weights > 0
  var_est[valid] <- 1 / weights[valid]
  var_est
}


#' Transform weights for Logan plot
#'
#' Transforms original TAC weights to weights appropriate for Logan plot
#' regression using the delta method. The Logan plot Y-axis is
#' integral(TAC)/TAC, so the variance propagation accounts for the cumulative
#' integral in the numerator and the TAC value in the denominator.
#'
#' @param t_tac Numeric vector of times for each frame in minutes.
#' @param dur Numeric vector of the time durations of the frames in minutes.
#' @param tac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame.
#' @param weights_original Numeric vector of the original weights assigned to
#'   each frame.
#'
#' @return A numeric vector of transformed weights suitable for Logan plot
#'   regression. NA values indicate frames that cannot be weighted (e.g., first
#'   frame, zero TAC values).
#'
#' @details This function uses the delta method to propagate variance from the
#'   original TAC measurements to the transformed Y values (integral(TAC)/TAC).
#'   This function can be used for both Logan and refLogan models since the
#'   Y-axis transformation is identical.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
weights_Logan_transform <- function(t_tac, dur, tac, weights_original) {

  stopifnot(length(t_tac) == length(dur),
            length(tac)   == length(weights_original),
            length(tac)   == length(t_tac))

  n_frames <- length(tac)

  # Var(tac) up to common scale factor; invalid NLS weights -> NA
  var_tac <- weights_to_var(weights_original)

  # Integral of tissue TAC
  integ     <- cumsum(tac * dur)          # ∫_0^{t_i} tac(τ) dτ
  integ_lag <- c(0, integ[-n_frames])     # integral up to previous frame

  weights_logan <- rep(NA_real_, n_frames)

  for (i in seq_len(n_frames)) {
    # First point or invalid tac -> no meaningful Logan weight
    if (i == 1L || !is.finite(tac[i]) || tac[i] == 0) {
      weights_logan[i] <- NA_real_
      next
    }

    idx_up_to_i <- seq_len(i)
    usable_idx  <- idx_up_to_i[is.finite(var_tac[idx_up_to_i])]

    if (length(usable_idx) == 0L) {
      weights_logan[i] <- NA_real_
      next
    }

    # Delta-method derivatives:
    # For j < i:  dy_i / d tac_j = dur_j / tac_i
    # For j = i:  dy_i / d tac_i = -integ_lag[i] / tac_i^2
    deriv_y_wrt_tac <- dur[usable_idx] / tac[i]

    if (i %in% usable_idx) {
      idx_current <- which(usable_idx == i)
      deriv_y_wrt_tac[idx_current] <- -integ_lag[i] / (tac[i]^2)
    }

    var_y_i <- sum((deriv_y_wrt_tac^2) * var_tac[usable_idx])

    weights_logan[i] <- if (is.finite(var_y_i) && var_y_i > 0) 1 / var_y_i else NA_real_
  }

  return(weights_logan)
}


#' Transform weights for reference Patlak plot
#'
#' Transforms original TAC weights to weights appropriate for reference Patlak
#' plot regression using the delta method. The refPatlak Y-axis is ROI/REF.
#'
#' @param t_tac Numeric vector of times for each frame in minutes.
#' @param reftac Numeric vector of radioactivity concentrations in the reference
#'   tissue for each frame.
#' @param roitac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame.
#' @param weights_original Numeric vector of the original weights assigned to
#'   each frame.
#'
#' @return A numeric vector of transformed weights suitable for reference Patlak
#'   plot regression. NA values indicate frames that cannot be weighted (e.g.,
#'   zero reference TAC values).
#'
#' @details This function uses the delta method to propagate variance from the
#'   original TAC measurements to the transformed Y values (ROI/REF). The
#'   reference TAC is treated as noise-free.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
weights_refPatlak_transform <- function(t_tac, reftac, roitac,
                                        weights_original) {

  stopifnot(length(t_tac) == length(reftac),
            length(roitac) == length(reftac),
            length(roitac) == length(weights_original))

  n_frames <- length(roitac)

  # Var(roitac) up to a common scale factor; invalid weights -> NA
  var_roitac <- weights_to_var(weights_original)

  weights_refpatlak <- rep(NA_real_, n_frames)

  for (i in seq_len(n_frames)) {

    # Need finite roitac, reftac, and variance
    if (!is.finite(roitac[i]) || !is.finite(reftac[i]) ||
        reftac[i] == 0 || !is.finite(var_roitac[i]) || var_roitac[i] <= 0) {
      weights_refpatlak[i] <- NA_real_
      next
    }

    # y_i = roitac_i / reftac_i
    # dy_i / d roitac_i = 1 / reftac_i   (reftac treated as noise-free)
    dy_droitac <- 1 / reftac[i]

    var_y_i <- (dy_droitac^2) * var_roitac[i]

    weights_refpatlak[i] <- if (is.finite(var_y_i) && var_y_i > 0) 1 / var_y_i else NA_real_
  }

  return(weights_refpatlak)
}


#' Transform weights for Patlak plot
#'
#' Transforms original TAC weights to weights appropriate for Patlak plot
#' regression using the delta method. The Patlak Y-axis is TAC/AIF.
#'
#' @param t_tac Numeric vector of times for each frame in minutes.
#' @param tac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame.
#' @param input Data frame containing the blood, plasma, and parent fraction
#'   concentrations over time. This can be generated using the
#'   \code{blood_interp} function. The AIF column will be interpolated to
#'   match t_tac times.
#' @param weights_original Numeric vector of the original weights assigned to
#'   each frame.
#'
#' @return A numeric vector of transformed weights suitable for Patlak plot
#'   regression. NA values indicate frames that cannot be weighted (e.g.,
#'   zero AIF values).
#'
#' @details This function uses the delta method to propagate variance from the
#'   original TAC measurements to the transformed Y values (TAC/AIF). The
#'   arterial input function is treated as noise-free.
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
weights_Patlak_transform <- function(t_tac, tac, input, weights_original) {

  stopifnot(length(tac) == length(t_tac),
            length(tac) == length(weights_original))

  n_frames <- length(tac)

  # Interpolate AIF to TAC times
  aif <- pracma::interp1(input$Time, input$AIF, t_tac, method = "linear")

  # Var(tac) up to a common scale factor; invalid weights -> NA
  var_tac <- weights_to_var(weights_original)

  weights_patlak <- rep(NA_real_, n_frames)

  for (i in seq_len(n_frames)) {

    # Need finite tac, aif, and variance
    if (!is.finite(tac[i]) || !is.finite(aif[i]) ||
        aif[i] == 0 || !is.finite(var_tac[i]) || var_tac[i] <= 0) {
      weights_patlak[i] <- NA_real_
      next
    }

    # y_i = tac_i / aif_i
    # dy_i / d tac_i = 1 / aif_i   (aif treated as noise-free)
    dy_dtac <- 1 / aif[i]

    var_y_i <- (dy_dtac^2) * var_tac[i]

    weights_patlak[i] <- if (is.finite(var_y_i) && var_y_i > 0) 1 / var_y_i else NA_real_
  }

  return(weights_patlak)
}
