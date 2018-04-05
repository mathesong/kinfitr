#' Standardised Uptake Value
#'
#' Function to calculate the standardised uptake value and its integral.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param injRad The injected radioactivity.  If not included, this will be set to 1 in case one is using SUV ratios.
#' @param bodymass The body mass of the participant. If not included, this will be set to 1 in case one is using SUV ratios.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#'
#' @return A list with a data frame of the calculated parameters \code{out$par} and a dataframe containing the TACs both of the
#' original data and the transformed values \code{out$tacs}
#'
#' @examples
#' SUV(t_tac, tac, injRad=150, bodymass=85)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

SUV <- function(t_tac, tac, injRad=1, bodymass=1, frameStartEnd) {

  # Tidying

  tidyinput <- tidyinput_art(t_tac, tac, tac, frameStartEnd) # Don't need weights, thus just set to same as tac

  t_tac <- tidyinput$t_tac
  tac <- tidyinput$tac


  # 'Model'

  denominator <- injRad / bodymass

  intSUV <- pracma::trapz(t_tac, tac) / denominator
  suvtac <- tac / denominator


  # Output

  par <- as.data.frame(list(intSUV = intSUV))

  tacs <- data.frame(Time = t_tac, Target = tac, SUV = suvtac)

  out <- list(par = par, tacs = tacs)

  return(out)
}
