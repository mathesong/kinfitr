#' Create weights using durations and counts
#'
#' This is a non-standard weighting scheme, which is designed to approximate the
#' weighting scheme used at Karolinska Institutet, which is still awaiting
#' publication. It is advised to calculate weights once for each PET
#' measurement, using a large region (e.g. whole brain). It uses the following
#' equation: sqrt(durations) / sqrt(counts).
#'
#' @param t_start The starting times of the frames in minutes.
#' @param t_end The end times of the frames in minutes.
#' @param tac The time activity curve.
#' @param radioisotope The radioisotope.
#' @param minweight The minimum weight. Weights will be calculated between this
#'   value and 1.
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
#' weights_sqrtdurcount(
#'   s1$StartTime/60,
#'   (s1$StartTime + s1$Duration)/60,
#'   radioisotope = "C11",
#'   tac = s1$WB)
weights_sqrtdurcount <- function(t_start, t_end, tac,
                                    radioisotope = c("C11", "O15", "F18"),
                                    minweight = 0.7, weight_checkn = 5) {

  radioisotope <- match.arg(radioisotope, c("C11", "O15", "F18"))

  tac_uncor <- decay_uncorrect(t_start, t_end, tac, radioisotope)
  durations <- t_end - t_start

  calcweights <- sqrt(durations) * sqrt(tac_uncor)

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

  calcweights[durations==0] <- 0

  return(calcweights)
}

