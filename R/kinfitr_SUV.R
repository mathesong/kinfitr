#' Standardised Uptake Value
#'
#' Function to calculate the standardised uptake value and its integral.
#'
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param dur_tac Numeric vector of durations for each frame in minutes. This will be used instead of the middle points of the frames if provided.
#' @param injRad The injected radioactivity.  If not included, this will be set to 1 in case one is using SUV ratios.
#' @param bodymass The body mass of the participant. If not included, this will be set to 1 in case one is using SUV ratios.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This can be used to assess time stability for example.
#' @param timeStartEnd Optional. This allows one to specify the beginning and end time point instead of defining the frame numbers using frameStartEnd. This function will restrict the model to all time frames whose t_tac is between the values, i.e. c(0,5) will select all frames with midtimes during the first 5 minutes. Note that this requires t_tac.
#'
#' @return A list with a data frame of the calculated parameters \code{out$par} and a dataframe containing the TACs both of the
#' original data and the transformed values \code{out$tacs}
#'
#' @examples
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' dur_tac <- pbr28$tacs[[2]]$Duration / 60
#' tac <- pbr28$tacs[[2]]$FC
#'
#'
#' fit1 <- SUV(tac, t_tac, injRad = 150, bodymass = 85)
#' fit2 <- SUV(tac, dur_tac = dur_tac, injRad = 150, bodymass = 85)
#' fit3 <- SUV(tac, t_tac = t_tac, dur_tac = dur_tac, injRad = 150, bodymass = 85)
#' fit4 <- SUV(tac, t_tac = t_tac, dur_tac = dur_tac, injRad = 150, bodymass = 85, frameStartEnd = c(1,5))
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

SUV <- function(tac, t_tac = NULL, dur_tac = NULL, injRad = 1, bodymass = 1,
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

  if (!is.null(dur_tac) & is.null(t_tac)) {

    tidyinput <- tidyinput_art(dur_tac, tac, tac,
                               frameStartEnd) # Don't need weights, thus just set to same as tac
    dur_tac <- tidyinput$t_tac

    tac <- tidyinput$tac
  }

  if (is.null(dur_tac) & !is.null(t_tac)) {
    tidyinput <- tidyinput_art(t_tac, tac, tac,
                               frameStartEnd) # Don't need weights, thus just set to same as tac
    t_tac <- tidyinput$t_tac

    tac <- tidyinput$tac
  }

  if (is.null(dur_tac) & is.null(t_tac)) {
    stop("Either t_tac or dur_tac must be provided")
  }

  if (!is.null(dur_tac) & !is.null(t_tac)) {
    tidyinput1 <- tidyinput_art(dur_tac, tac, tac,
                               frameStartEnd) # Don't need weights, thus just set to same as tac
    dur_tac <- tidyinput1$t_tac

    tidyinput2 <- tidyinput_art(t_tac, tac, tac,
                               frameStartEnd) # Don't need weights, thus just set to same as tac
    t_tac <- tidyinput2$t_tac

    tac <- tidyinput1$tac
  }


  # 'Model'

  denominator <- injRad / bodymass

  if (is.null(dur_tac)) {
    intSUV <- pracma::trapz(t_tac, tac) / denominator
  } else {
    intSUV <- sum((tac * dur_tac) / denominator)
  }

  suvtac <- tac / denominator


  # Output

  par <- as.data.frame(list(intSUV = intSUV))

  tacs <- data.frame(Target = tac, SUV = suvtac)

  if (is.null(dur_tac)) {
    tacs$Duration <- NA
  } else {
    tacs$Duration <- dur_tac
  }


  if (is.null(t_tac)) {
    tacs$Time <- NA
  } else {
    tacs$Time <- t_tac
  }

  out <- list(par = par, tacs = tacs)

  return(out)
}
