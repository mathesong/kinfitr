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
#' @param method Which method should be used? 1 represents duration^2 /
#'   tac_uncorrected. 2 represents sqrt(durations) * sqrt(tac_uncorrected). 3
#'   represents duration / tac. 4 represents sqrt(durations). 5 represents
#'   durations * exp((-ln(2)) / halflife ). 6 represents durations /
#'   tac_uncorrected. 7 represents durations.
#' @param minweight The minimum weight. Weights will be calculated as a fraction
#'   between this value and 1. A zero frame with duration=0 will be set to 0
#'   though.
#' @param weight_checkn The number of values of the weights to check to make
#'   sure that things haven't gone terribly wrong. It will check that this
#'   number of weights are all above half of the maximum weight.
#'
#' @return The vector of weights.
#' @export
#'
#' @examples
#' data(pbr28)
#' s1 <- pbr28$tacs[[1]]
#'
#' # Assuming the data were not decay-corrected (they are)
#' weights_create(
#'   s1$StartTime/60,
#'   (s1$StartTime + s1$Duration)/60,
#'   radioisotope = "C11",
#'   tac = s1$WB)
weights_create <- function(t_start, t_end, tac,
                                    radioisotope = c("C11", "O15", "F18"),
                                    method = 2,
                                    minweight = 0.7, weight_checkn = 5) {

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
    method == 1 ~ durations^2 / tac_uncor,
    method == 2 ~ sqrt(durations) * sqrt(tac_uncor),
    method == 3 ~ sqrt(durations) / tac,
    method == 4 ~ sqrt(durations),
    method == 5 ~ durations * exp( (-1*log(2)) / hl ),
    method == 6 ~ durations / tac_uncor,
    method == 7 ~ durations
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
  calcweights <- minweight + calcweights * (1-minweight)

  # This shouldn't have the scaling factor
  calcweights[durations==0] <- 0

  return(calcweights)
}

