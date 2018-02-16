#' Simplified Reference Tissue Model
#'
#' Function to fit the SRTM model of Lammertsma and Hume (1996) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame. We include zero at
#' time zero: if not included, it is added.
#' @param roitac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param frameStartEnd Optional. This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#' @param R1.start Optional. Starting parameter for fitting of R1. Default is 1.
#' @param R1.lower Optional. Lower bound for the fitting of R1. Default is 0.
#' @param R1.upper Optional. Upper bound for the fitting of R1. Default is 10.
#' @param k2.start Optional. Starting parameter for fitting of k2. Default is 0.1.
#' @param k2.lower Optional. Lower bound for the fitting of k2. Default is 0.
#' @param k2.upper Optional. Upper bound for the fitting of k2. Default is 1.
#' @param bp.start Optional. Starting parameter for fitting of bp. Default is 1.5.
#' @param bp.lower Optional. Lower bound for the fitting of bp. Default is -10.
#' @param bp.upper Optional. Upper bound for the fitting of bp. Default is 15.
#' @param multstart_iter Number of iterations for starting parameters. Default is 1.
#'   For more information, see \code{\link[nls.multstart]{nls_multstart}}. If
#'   specified as 1 for any parameters, the original starting value will be
#'   used, and the multstart_lower and multstart_upper values ignored.
#' @param multstart_lower Optional. Lower bounds for starting parameters. Defaults
#'   to halfway between start and lower starting parameters.
#' @param multstart_upper Optional. Upper bounds for starting parameters. Defaults
#'   to halfway between start and upper starting parameters.
#' @param printvals Optional. This displays the parameter values for each iteration of the
#' model. This is useful for debugging and changing starting values and upper and lower
#' bounds for parameters.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit},
#' the model weights \code{out$weights}, and a dataframe containing the TACs both of the data and the fitted
#' values \code{out$tacs}.
#'
#' @examples
#' srtm(t_tac, reftac, roitac)
#' srtm(t_tac, reftac, roitac, weights, frameStartEnd = c(1,11), bp.upper=1)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Lammertsma AA, Hume SP. Simplified reference tissue model for PET receptor studies. Neuroimage. 1996 Dec 31;4(3):153-8.
#'
#' @export

srtm <- function(t_tac, reftac, roitac, weights, frameStartEnd,
                 R1.start = 1 , R1.lower = 0 , R1.upper = 10 ,
                 k2.start=0.1 , k2.lower = 0 , k2.upper=1 ,
                 bp.start=1.5 , bp.lower=-10 , bp.upper=15,
                 multstart_iter=1, multstart_lower, multstart_upper,
                 printvals=F) {


  # Tidying

  tidyinput <- tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd)

  t_tac   <- tidyinput$t_tac
  reftac  <- tidyinput$reftac
  roitac  <- tidyinput$roitac
  weights <- tidyinput$weights


  # Parameters

  start <- c( R1 = R1.start, k2 = k2.start, bp = bp.start)
  lower <- c( R1 = R1.lower, k2 = k2.lower, bp = bp.lower)
  upper <- c( R1 = R1.upper, k2 = k2.upper, bp = bp.upper)

  if( length(multstart_iter) == length(start) || length(multstart_iter) == 1 ) {

    ### Missing multstart boundaries
    if(missing(multstart_lower)) {
      multstart_lower = start - (start-lower)/2
    }
    if(missing(multstart_upper)) {
      multstart_upper = start + (upper-start)/2
    }

    ### No multstart for some variables ###
    if( any(multstart_iter==1) ) {
      non_iterable <- which(multstart_iter==1)
      multstart_lower[non_iterable] = start[non_iterable]
      multstart_upper[non_iterable] = start[non_iterable]
    }

    ### Incorrect multstart boundaries
    if( length(multstart_lower) != length(start) || length(multstart_upper) != length(start) ) {
      stop('multstart_lower and multstart_upper should be the same length as the number of parameters')
    }
  } else {
    stop('multstart_iter should be of length 1 or of the same length as the number of parameters')
  }


  # Solution

  if( prod(multstart_iter) == 1 ) {

    output <- minpack.lm::nlsLM(roitac ~ srtm_model(t_tac, reftac, R1, k2, bp),
                                start =  start, lower = lower, upper = upper,
                                weights=weights,
                                control = minpack.lm::nls.lm.control(maxiter = 200),
                                trace=printvals)

  } else {

    output <- nls.multstart::nls_multstart(roitac ~ srtm_model(t_tac, reftac, R1, k2, bp),
                                          supp_errors = 'Y',
                                          start_lower = multstart_lower,
                                          start_upper = multstart_upper,
                                          iter = multstart_iter, convergence_count = FALSE,
                                          lower = lower, upper=upper, modelweights=weights)

  }

  # Output

  tacs <- data.frame(Time = t_tac, Reference = reftac, Target = roitac, Target_fitted=as.numeric(fitted(output)))

  par = as.data.frame(as.list(coef(output)))

  par.se = as.data.frame(as.list(sqrt(abs(vcov(output)[,1]))))
  names(par.se) = paste0(names(par.se), '.se')

  out <- list(par = par, par.se = par.se,
              fit = output, weights = weights, tacs = tacs,
              model='srtm')

  return(out)

}


#' Model: Simplified Reference Tissue Model
#'
#' This is the SRTM model itself by which predicted values are generated.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a zero.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param R1 Parameter value for R1
#' @param k2 Parameter value for k2
#' @param bp Parameter value for bp
#'
#' @return A numeric vector of the predicted values of the TAC in the target region.
#'
#' @examples
#' srtm_model(t_tac, reftac, R1=0.9, k2=0.1, bp=1.5)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Lammertsma AA, Hume SP. Simplified reference tissue model for PET receptor studies. Neuroimage. 1996 Dec 31;4(3):153-8.
#'
#' @export


srtm_model <- function(t_tac, reftac, R1, k2, bp) {

  interptime <- pracma::linspace( min(t_tac) , max(t_tac) , 1024 )
  step <- interptime[2] - interptime[1]

  iref <- pracma::interp1(t_tac, reftac, interptime, method="linear")

  a <- (k2-(R1*k2/(1+bp)))*iref
  b <- exp((-k2/(1+bp))*interptime)

  ND <- R1*iref
  BOUND <- kinfit_convolve(a,b,step)
  tmp <- ND+BOUND

  outtac <- pracma::interp1(interptime, tmp, t_tac)

  return(outtac)
}


#' Plot: Simplified Reference Tissue Model
#'
#' Function to visualise the fit of the SRTM model to data.
#'
#' @param srtmout The output object of the SRTM fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#' @param refname Optional. The name of the Reference Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' plot_srtmfit(srtmout)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_srtmfit <- function(srtmout, roiname, refname) {

  measured <- data.frame(Time = srtmout$tacs$Time,
                       Reference = srtmout$tacs$Reference,
                       ROI.measured = srtmout$tacs$Target,
                       Weights = weights(srtmout$fit))

  fitted <- data.frame(Time = srtmout$tacs$Time,
                         ROI.fitted = srtmout$tacs$Target_fitted,
                         Weights = weights(srtmout$fit))

  if(missing(roiname)) {roiname = 'ROI'}
  if(missing(refname)) {refname = 'Reference'}

  measured = plyr::rename(measured, c('ROI.measured' = paste0(roiname, '.measured'),
                                      'Reference' = refname))

  fitted = plyr::rename(fitted, c('ROI.fitted' = paste0(roiname, '.fitted')) )

  tidymeasured <- tidyr::gather(measured, key=Region, value=Radioactivity,
                                -Time, -Weights, factor_key = F)

  tidyfitted <- tidyr::gather(fitted, key=Region, value=Radioactivity,
                              -Time, -Weights, factor_key = F)


  Region <- forcats::fct_inorder(factor(c(tidymeasured$Region, tidyfitted$Region)) )

  myColors <- RColorBrewer::brewer.pal(3,"Set1")
  names(myColors) <- levels(Region)
  colScale <- scale_colour_manual(name = "Region",values = myColors)

  outplot <- ggplot(tidymeasured, aes(x=Time, y=Radioactivity, colour=Region)) +
    geom_point(data=tidymeasured, aes(shape='a', size=Weights)) +
    geom_line(data=tidyfitted) +
    guides(shape=FALSE, color=guide_legend(order=1)) + colScale +
    scale_size(range=c(1,3))

  return(outplot)
}
