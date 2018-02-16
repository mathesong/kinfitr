#' Simultaneous Estimation of Non-Displaceable Binding (SIME)
#'
#' Function to fit the SIME Model of Ogden et al (2015) to data to estimate Vnd for a set of TACs.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param tacdf Named dataframe of TACs in wide format, i.e. each TAC should be a column.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param roiweights Optional. Numeric vector of the weights assigned to each ROI in the fitting. If not specified, uniform weights
#' will be used.
#' @param inpshift Optional. The number of minutes by which to shift the timing of the input data frame forwards or backwards.
#' If not specified, this will be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB_fixed Optional. The blood volume fraction.  If not specified, this will be set to 0.05. This can be fitted using 1TCM or 2TCM.
#' @param twotcmstart Optional. The function can fit a 2TCM model to one of the ROIs and use the estimated k2, k3 and k4 as starting
#' parameters for the rest of the fits. If left alone, these parameters will be specified as below. If one wishes to run the 2TCM to
#' start off, use a numeric value to specify which column of \code{tacdf} to use for fitting this: best to use the largest ROI.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#' @param k2.start Optional. Starting parameter for fitting of k2. Default is 0.1.
#' @param k2.lower Optional. Lower bound for the fitting of k2. Default is 0.
#' @param k2.upper Optional. Upper bound for the fitting of k2. Default is 0.5.
#' @param k3.start Optional. Starting parameter for fitting of k3. Default is 0.1.
#' @param k3.lower Optional. Lower bound for the fitting of k3. Default is 0.
#' @param k3.upper Optional. Upper bound for the fitting of k3. Default is 0.5.
#' @param k4.start Optional. Starting parameter for fitting of k4. Default is 0.1.
#' @param k4.lower Optional. Lower bound for the fitting of k4. Default is 0.
#' @param k4.upper Optional. Upper bound for the fitting of k4. Default is 0.5.
#'
#' @return A list with a data frame of the fitted parameter \code{out$par}, the dataframe containing the times and TACs \code{out$tacs},
#' the mean cost values after fitting (after ROI weighting) \code{out$fitvals}, the ROI cost values after fitting (before ROI
#' weighting) \code{out$roifits}, the blood input data frame after time shifting \code{input}, a vector of the weights \code{out$weights},
#' a vector of the ROI weights \code{out$roiweights}, the inpshift value used \code{inpshift} and the vB value used \code{out$vB_fixed}.
#'
#'
#' @examples
#' Vndgrid <- seq(from=0, to=5, by=0.1)
#' SIMEout <-  SIME(t_tac, tacdf, input, Vndgrid, weights = weights,
#'                  inpshift = onetcmout$par$inpshift, vB_fixed = onetcmout$par$vB)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
#'
#' @references Ogden RT, Zanderigo F, Parsey RV. Estimation of in vivo nonspecific binding in positron emission tomography studies without requiring a reference region. NeuroImage. 2015 Mar 31;108:234-42.
#'
#' @importFrom dplyr "%>%"

SIME <- function(t_tac, tacdf, input, Vndgrid, weights, roiweights,
                 inpshift = 0, vB_fixed, twotcmstart, frameStartEnd,
                 k2.start = 0.1 , k2.lower = 0 , k2.upper = 0.5 ,
                 k3.start = 0.1 , k3.lower = 0 , k3.upper = 0.5 ,
                 k4.start = 0.1 , k4.lower = 0 , k4.upper = 0.5 ) {

  # Tidying

  if(missing(weights)) {
    weights = rep(1, nrow(tacdf))
  }

  lengths <- c(length(t_tac), nrow(tacdf), length(weights))
  if(!all(lengths == lengths[1])) {
    stop('The lengths of the times, TACs and/or weights are not equal')
  }

  if(!missing(frameStartEnd)) {
    t_tac <- t_tac[ frameStartEnd[1] : frameStartEnd[2] ]
    tacdf <- tacdf[ frameStartEnd[1] : frameStartEnd[2] , ]
  }

  if(min(t_tac) < 0) {
    stop('There are negative times in the TAC')
  }

  if(missing(roiweights) == T) {
    roiweights = rep(1, ncol(tacdf))
  }
  roiweights = roiweights / max(roiweights)

  if(length(roiweights) != ncol(tacdf)) {
    stop('The number of ROIs and roiweights do not match')
  }

  if(min(t_tac) > 0) {
    t_tac = c(0, t_tac)
    tacdf = rbind(0, tacdf)
    weights = c(0, weights)
  }

  newvals = shift_timings_df(t_tac, tacdf, input, inpshift)

  t_tac <- newvals$t_tac
  tacdf <- newvals$tacdf

  input <- newvals$input

  t_inp <- newvals$input$Time


  # Parameters

  start <- c(k2 = k2.start, k3 = k3.start, k4 = k4.start)
  lower <- c(k2 = k2.lower, k3 = k3.lower, k4 = k4.lower)
  upper <- c(k2 = k2.upper, k3 = k3.upper, k4 = k4.upper)


  # 2tcm Starting Parameters

  if(!missing(twotcmstart)) {
    twotcmout <- twotcm(t_tac, tac = tacdf[,twotcmstart], input,
                        weights=weights, inpshift = inpshift,
                        vB_fixed = vB_fixed)
    start[1] <- twotcmout$par$k2
    start[2] <-  twotcmout$par$k3
    start[3] <-  twotcmout$par$k4
    if(missing(vB_fixed)) { vB_fixed <- twotcmout$par$vB }
  } else { if(missing(vB_fixed)) { vB_fixed <- 0.05 } }

  # Solution

  tacdf$Time <- t_tac
  tacdf$weights <- weights

  tidytacs <- tidyr::gather(tacdf, key=Region, value=Radioactivity, -Time, -weights)
  frames <- nrow(tidytacs)

  tidytacs <- tidytacs[rep(1:nrow(tidytacs),times = length(Vndgrid)),]
  tidytacs$Vnd <- rep(Vndgrid, each=frames)

  Cost <- tidytacs %>%
    dplyr::group_by(Vnd, Region) %>%
    dplyr::do(gridCost = SIMEroi(t_tac=.$Time, tac=.$Radioactivity, input=input, Vnd=.$Vnd[1],
                                 vB_fixed=vB_fixed, weights = .$weights,
              k2.start = start[1] , k2.lower = lower[1] , k2.upper = upper[1] ,
              k3.start = start[2] , k3.lower = lower[2] , k3.upper = upper[2] ,
              k4.start = start[3] , k4.lower = lower[3] , k4.upper = upper[3])) %>%
    dplyr::mutate(gridCost = as.numeric(gridCost))



  # Calculating SSmean

  Cost_wide <- tidyr::spread(Cost, Region, gridCost)
  roiweighting <- t(replicate(length(Vndgrid), roiweights))

  SSmean = rowMeans( Cost_wide[,-1] * roiweighting )


  # Output

  Vnd = Vndgrid[which.min(SSmean)]

  par = as.data.frame(list(Vnd = Vnd))

  tacs = data.frame(Time = t_tac)
  tacs <- cbind(tacs, tacdf)

  fitvals = data.frame(Vnd = Vndgrid, SSmean = SSmean)

  out <- list( par = par, tacs = tacs, fitvals = fitvals, roifits = Cost_wide, input = input,
               weights=weights, roiweights = roiweights, inpshift = inpshift, vB = vB_fixed, model='SIME')

  return(out)


}

#' Cost Function: SIME
#'
#' Function to obtain the cost for a given ROI and Vnd for a TAC.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame including a zero at time
#' zero.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function. It should already be shifted if a shift is desired, using \code{shift_timings}.
#' @param Vnd The specified Vnd value to calculate the cost for for the given ROI.
#' @param vB_fixed Optional. The blood volume fraction.  If not specified, this will be set to 0.05. This can be fitted using 1TCM or 2TCM.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param k2.start Optional. Starting parameter for fitting of k2. Default is 0.1.
#' @param k2.lower Optional. Lower bound for the fitting of k2. Default is 0.
#' @param k2.upper Optional. Upper bound for the fitting of k2. Default is 0.5.
#' @param k3.start Optional. Starting parameter for fitting of k3. Default is 0.1.
#' @param k3.lower Optional. Lower bound for the fitting of k3. Default is 0.
#' @param k3.upper Optional. Upper bound for the fitting of k3. Default is 0.5.
#' @param k4.start Optional. Starting parameter for fitting of k4. Default is 0.1.
#' @param k4.lower Optional. Lower bound for the fitting of k4. Default is 0.
#' @param k4.upper Optional. Upper bound for the fitting of k4. Default is 0.5.
#'
#' @return If the fit converged for the given Vnd value, the function will return the weighted sum of squares of the
#' residuals. If the fit did not converge for the given Vnd value, the function will return NA.
#'
#'
#' @examples
#' SIMEroi(t_tac, tac, input, Vnd=5, vB_fixed=0.05, weights = weights)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Ogden RT, Zanderigo F, Parsey RV. Estimation of in vivo nonspecific binding in positron emission tomography studies without requiring a reference region. NeuroImage. 2015 Mar 31;108:234-42.
#'
#' @export

SIMEroi <- function(t_tac, tac, input, Vnd, vB_fixed, weights,
                    k2.start = 0.1 , k2.lower = 0 , k2.upper = 0.5 ,
                    k3.start = 0.1 , k3.lower = 0 , k3.upper = 0.5 ,
                    k4.start = 0.1 , k4.lower = 0 , k4.upper = 0.5 ) {


  start <- c(k2 = k2.start, k3 = k3.start, k4 = k4.start)
  lower <- c(k2 = k2.lower, k3 = k3.lower, k4 = k4.lower)
  upper <- c(k2 = k2.upper, k3 = k3.upper, k4 = k4.upper)

  pred <-   tryCatch( {
    minpack.lm::nlsLM(tac ~ SIME_model(t_tac, input, Vnd, k2, k3, k4, vB=vB_fixed),
                      start = start, lower = lower, upper = upper,
                      weights=weights, control = minpack.lm::nls.lm.control(maxiter = 200)) },
    error = function(e) NA )

  if(is.na(pred[1])) {SSw = NA} else {
    SSw <-  sum( weights(pred) * ( residuals(pred)^2 ) )
    return(SSw)
  }
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
#' SIME_model(t_tac, input, Vnd=5, k2=0.1, k3=0.05, k4=0.04, vB=0.05)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Ogden RT, Zanderigo F, Parsey RV. Estimation of in vivo nonspecific binding in positron emission tomography studies without requiring a reference region. NeuroImage. 2015 Mar 31;108:234-42.
#'
#' @export

SIME_model <- function(t_tac, input, Vnd, k2, k3, k4, vB) {

  blood <- input$blood
  plasma <- input$plasma
  parentfrac <- input$parentfrac

  i_inp <- plasma*parentfrac

  interptime <- input$Time
  step <- interptime[2] - interptime[1]

  K1 <- k2*Vnd

  delta <- sqrt((k2+k3+k4)^2 - 4*k2*k4)
  th1   <- (k2+k3+k4+delta) / 2
  th2   <- (k2+k3+k4-delta) / 2
  ph1   <- K1*(th1-k3-k4) / delta
  ph2   <- K1*(th2-k3-k4) / (-delta)

  a = ph1*exp(-th1*interptime) + ph2*exp(-th2*interptime)
  b = i_inp

  i_outtac <- kinfit_convolve(a,b,step)

  # Correction for vB
  i_outtac <- i_outtac*(1-vB) + vB*i_inp

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
#' plot_SIMEfit(SIMEout)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_SIMEfit <- function(SIMEout) {

  fitvals <- SIMEout$fitvals

  roifit <- tidyr::gather(SIMEout$roifits, Region, Cost, -Vnd)

  minmax <- list(min = min(c(fitvals$SSmean, roifit$Cost)),
                 max = max(c(fitvals$SSmean, roifit$Cost)) )


  outplot <- ggplot(roifit, aes(x = Vnd, y = Cost)) +
    geom_point(aes(colour = Region)) +
    geom_line(data=fitvals, aes(x = Vnd, y = SSmean), size=1) +
    geom_vline(xintercept = SIMEout$par$Vnd, linetype="dashed") +
    xlab(expression(V[ND])) +
    scale_y_log10("Cost",
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))

  return(outplot)

}
