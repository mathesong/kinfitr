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
#' @param method Which method should be used? 1 represents duration /
#'   (tac_uncorrected). 2 represents sqrt(durations*tac_uncorrected). 3
#'   represents duration / tac. 4 represents sqrt(durations). 5 represents
#'   durations * exp((-ln(2)) / halflife ). 6 represents durations /
#'   tac. 7 represents durations. 8 represents duration^2 / tac_uncor.
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
#'   tac = s1$WB, minweight_risetopeak=TRUE)
weights_create <- function(t_start, t_end, tac,
                           radioisotope = c("C11", "O15", "F18"),
                           method = 2,
                           minweight = 0.67,
                           minweight_risetopeak = FALSE,
                           weight_checkn = 5) {

  radioisotope <- match.arg(radioisotope, c("C11", "O15", "F18"))

  tac_uncor <- decay_uncorrect(t_start, t_end, tac, radioisotope)
  durations <- t_end - t_start
  t_tac <- t_start + durations/2

  hl <- dplyr::case_when(
    radioisotope == "C11" ~ 20.4,
    radioisotope == "O15" ~ 2.05,
    radioisotope == "F18" ~ 109.8
  )

  calcweights <- dplyr::case_when(
    method == 1 ~ durations / tac_uncor,
    method == 2 ~ sqrt(durations*tac_uncor),
    method == 3 ~ sqrt(durations) / tac,
    method == 4 ~ sqrt(durations),
    method == 5 ~ durations * exp( (-1*log(2)) / hl ),
    method == 6 ~ durations / tac,
    method == 7 ~ durations,
    method == 8 ~ durations^2 / tac_uncor
  )

  # Fixing before checking
  calcweights[durations==0] <- 0

  # Check to make sure things don't go horribly wrong
  maxweight  <- max(calcweights)
  maxweights <- tail(calcweights[order(calcweights)], weight_checkn)

  if( min(maxweights) < 0.5*maxweight ) {
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
    t_weightpeak <- t_tac[which.max(calcweights)]

    # Assign times a fraction of the peaktime
    t_peakfrac <- t_tac / t_weightpeak
    t_peakfrac[t_peakfrac > 1] <- 1

    # Define a rising minimum to the peak, and then minweight
    minweights <- minweights * t_peakfrac
    minweights[durations==0] <- 0
  }

  if( any(calcweights < minweights) ) {

    min_calcweight <- min(calcweights[durations!=0])

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


#' #' Visualise weights
#' #'
#' #' This allows one to visualise the set of weights described to evaluate
#' #' whether they seem sensible.
#' #'
#' #' @param t_start The starting times of the frames in minutes.
#' #' @param t_end The end times of the frames in minutes.
#' #' @param tac The time activity curve.
#' #' @param radioisotope The radioisotope.
#' #' @param method Which method should be used? 1 represents duration^2 /
#' #'   (tac_uncorrected). 2 represents sqrt(durations*tac_uncorrected). 3
#' #'   represents duration / tac. 4 represents sqrt(durations). 5 represents
#' #'   durations * exp((-ln(2)) / halflife ). 6 represents durations /
#' #'   tac. 7 represents durations. Uncorrected refers to decay correction.
#' #' @param minweight The minimum weight. Weights will be calculated as a fraction
#' #'   between this value and 1. A zero frame with duration=0 will be set to 0
#' #'   though.
#' #' @param minweight_risetopeak Should there be a linear rise of the minimum weight to
#' #'   the peak weight value? This means that values before the maximum weight can
#' #'   be below the minweight specified. This is helpful for downweighting noisy
#' #'   values towards the start of the TAC. Defaults to FALSE.
#' #' @param weight_checkn The number of values of the weights to check to make
#' #'   sure that things haven't gone terribly wrong. It will check that this
#' #'   number of weights are all above half of the maximum weight.
#' #' @param se_meanperk Visualise the standard error under the assumption that the
#' #'   noise in the data
#' #'
#' #' @return The vector of weights.
#' #'
#' #' @author Granville J Matheson, \email{mathesong@@gmail.com}
#' #'
#' #' @import ggplot2
#' #'
#' #' @export
#' #'
#' #' @examples
#' #' data(pbr28)
#' #' s1 <- pbr28$tacs[[1]]
#' #'
#' #' weights_vis(
#' #'   t_start = s1$StartTime/60,
#' #'   t_end = (s1$StartTime + s1$Duration)/60,
#' #'   tac = s1$WB,
#' #'   radioisotope = "C11",
#' #'   minweight_risetopeak=TRUE)
#' weights_vis <- function(t_start, t_end, tac,
#'                            radioisotope = c("C11", "O15", "F18"),
#'                            method = 2,
#'                            minweight = 0.67,
#'                            minweight_risetopeak = FALSE,
#'                            weight_checkn = 5,
#'                            se_meanperc = 5) {
#'
#'   radioisotope <- match.arg(radioisotope, c("C11", "O15", "F18"))
#'
#'   tac_uncor <- decay_uncorrect(t_start, t_end, tac, radioisotope)
#'   durations <- t_end - t_start
#'   t_tac <- t_start + durations/2
#'
#'   hl <- dplyr::case_when(
#'     radioisotope == "C11" ~ 20.4,
#'     radioisotope == "O15" ~ 2.05,
#'     radioisotope == "F18" ~ 109.8
#'   )
#'
#'   calcweights <- dplyr::case_when(
#'     method == 1 ~ durations^2 / (tac_uncor),
#'     method == 2 ~ sqrt(durations*tac_uncor),
#'     method == 3 ~ sqrt(durations) / tac,
#'     method == 4 ~ sqrt(durations),
#'     method == 5 ~ durations * exp( (-1*log(2)) / hl ),
#'     method == 6 ~ durations / tac,
#'     method == 7 ~ durations
#'   )
#'
#'   # Fixing before checking
#'   calcweights[durations==0] <- 0
#'
#'   # Check to make sure things don't go horribly wrong
#'   maxweight  <- max(calcweights)
#'   maxweights <- tail(calcweights[order(calcweights)], weight_checkn)
#'
#'   if( min(maxweights) < 0.5*maxweight ) {
#'     maxweight <- median(maxweights)
#'   }
#'
#'   calcweights[calcweights > maxweight] <- maxweight
#'   calcweights <- calcweights/maxweight
#'
#'   # Apply scaling factor
#'   minweights <- rep(minweight, length(calcweights))
#'
#'   if( minweight_risetopeak ) {
#'     # TAC midtimes
#'     t_tac <- t_start + 0.5*(t_end - t_start)
#'
#'     # Find the maximum TAC weight
#'     t_weightpeak <- t_tac[which.max(calcweights)]
#'
#'     # Assign times a fraction of the peaktime
#'     t_peakfrac <- t_tac / t_weightpeak
#'     t_peakfrac[t_peakfrac > 1] <- 1
#'
#'     # Define a rising minimum to the peak, and then minweight
#'     minweights <- minweights * t_peakfrac
#'     minweights[durations==0] <- 0
#'   }
#'
#'   if( any(calcweights < minweights) ) {
#'
#'     min_calcweight <- min(calcweights[durations!=0])
#'
#'     # scale 0 - 1
#'     calcweights <- calcweights - min_calcweight / (1- min_calcweight)
#'
#'     # proportion
#'     calcweights <- minweights + calcweights * (1-minweights)
#'   }
#'
#'   # Any dur=0 shouldn't have the scaling factor
#'   calcweights[durations==0] <- 0
#'
#'
#'
#'   # Visualise
#'
#'   weightdata <- data.frame(
#'     weights = calcweights,
#'     minweights = minweights,
#'     tac = tac / max(tac),
#'     tac_uncor = tac_uncor / max(tac)
#'   )
#'
#'   ggplot()
#'
#' }
